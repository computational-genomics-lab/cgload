#!/usr/bin/env Rscript
# txdbmaker in the loader comparison: build a TxDb from the GFF3, then write out
# what it holds.
#
#   Rscript benchmark/run_txdb.R GFF OUTDIR [CIRCULAR_SEQIDS]
#
# CIRCULAR_SEQIDS is a comma-separated list of the sequences the file marks
# Is_circular=true (truth.py records them); makeTxDbFromGFF otherwise knows
# only names such as chrM and MT. Writes into OUTDIR:
#
#   versions.txt     R, txdbmaker, GenomicFeatures and rtracklayer versions, and
#                    which package's makeTxDbFromGFF ran
#   types.txt        the GFF3 types it imports and treats as gene, transcript, exon,
#                    CDS and stop codon, from the package itself
#   transcripts.tsv  one row per transcript: name, sequence, strand, span, type, genes
#   exons.tsv        one row per exon of each transcript, with its rank
#   cds.tsv          one row per CDS piece of each transcript
#   genes.tsv        one row per gene ID the TxDb holds
#   summary.json     counts, seconds spent in makeTxDbFromGFF, warnings
#   txdb.sqlite      the TxDb itself (saveDb), for the query timings
#
# Every warning and message makeTxDbFromGFF gives is copied to standard error,
# so the comparison can tell an announced drop from a silent one.

suppressPackageStartupMessages(library(GenomicFeatures))
# makeTxDbFromGFF moved from GenomicFeatures to txdbmaker in Bioconductor 3.19.
# The benchmark's environment installs txdbmaker; an older installation that
# has only GenomicFeatures still runs, and versions.txt says which one did.
maker <- if (requireNamespace("txdbmaker", quietly = TRUE)) "txdbmaker" else "GenomicFeatures"
make_txdb <- getExportedValue(maker, "makeTxDbFromGFF")
version_of <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) as.character(packageVersion(pkg)) else "not installed"
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: run_txdb.R GFF OUTDIR [CIRCULAR_SEQIDS]")
gff <- args[[1]]
out <- args[[2]]
circular <- if (length(args) >= 3 && nzchar(args[[3]])) strsplit(args[[3]], ",")[[1]] else character(0)
dir.create(out, showWarnings = FALSE, recursive = TRUE)

writeLines(c(
  R.version.string,
  paste("txdbmaker", version_of("txdbmaker")),
  paste("GenomicFeatures", version_of("GenomicFeatures")),
  paste("rtracklayer", version_of("rtracklayer")),
  paste("makeTxDbFromGFF from", maker)
), file.path(out, "versions.txt"))

# Which GFF3 types makeTxDbFromGFF imports and treats as gene, transcript, exon,
# CDS and stop codon: the lists inside the package, read from its namespace so
# that the comparison follows the version that ran. Read before the long step,
# so a version that keeps them elsewhere stops at once, saying so.
type_lists <- c(gene = ".GENE_TYPES", transcript = ".TX_TYPES", exon = ".EXON_TYPES",
                cds = ".CDS_TYPES", stop_codon = ".STOP_CODON_TYPES",
                imported = "GFF_FEATURE_TYPES")
type_lines <- character(0)
for (kind in names(type_lists)) {
  ns <- asNamespace(maker)
  if (!exists(type_lists[[kind]], envir = ns, inherits = FALSE)) {
    message("benchmark script error: ", maker, " has no ", type_lists[[kind]],
            "; update run_txdb.R for this version")
    quit(status = 3)
  }
  values <- unique(get(type_lists[[kind]], envir = ns, inherits = FALSE))
  type_lines <- c(type_lines, paste(kind, paste0(maker, "::", type_lists[[kind]]), values,
                                    sep = "\t"))
}
writeLines(type_lines, file.path(out, "types.txt"))

warnings_seen <- 0L
started <- Sys.time()
txdb <- withCallingHandlers(
  make_txdb(gff, format = "gff3", circ_seqs = circular),
  warning = function(w) {
    warnings_seen <<- warnings_seen + 1L
    message("WARNING: ", conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))

columns <- c("tx_id", "tx_name", "gene_id")
if ("tx_type" %in% columns(txdb) || "TXTYPE" %in% columns(txdb)) columns <- c(columns, "tx_type")
tx <- tryCatch(transcripts(txdb, columns = columns),
               error = function(e) transcripts(txdb, columns = c("tx_id", "tx_name", "gene_id")))
genes_of <- vapply(as.list(mcols(tx)$gene_id), function(g) paste(sort(unique(g)), collapse = ","),
                   character(1))
tx_type <- if ("tx_type" %in% names(mcols(tx))) as.character(mcols(tx)$tx_type) else rep("", length(tx))
tx_table <- data.frame(
  tx_id = mcols(tx)$tx_id, tx_name = mcols(tx)$tx_name,
  seqid = as.character(seqnames(tx)), strand = as.character(strand(tx)),
  start = start(tx), end = end(tx), type = tx_type, genes = genes_of,
  stringsAsFactors = FALSE
)
write.table(tx_table, file.path(out, "transcripts.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE, na = "")

parts <- function(grl, name_column) {
  if (length(grl) == 0) {
    return(data.frame(tx_id = integer(0), rank = integer(0), seqid = character(0),
                      strand = character(0), start = integer(0), end = integer(0),
                      name = character(0)))
  }
  flat <- unlist(grl, use.names = FALSE)
  owners <- rep(as.integer(names(grl)), lengths(grl))
  meta <- mcols(flat)
  data.frame(
    tx_id = owners,
    rank = if ("exon_rank" %in% names(meta)) meta$exon_rank else NA_integer_,
    seqid = as.character(seqnames(flat)), strand = as.character(strand(flat)),
    start = start(flat), end = end(flat),
    name = if (name_column %in% names(meta)) as.character(meta[[name_column]]) else NA_character_,
    stringsAsFactors = FALSE
  )
}
exons_by_tx <- exonsBy(txdb, by = "tx", use.names = FALSE)
cds_by_tx <- cdsBy(txdb, by = "tx", use.names = FALSE)
write.table(parts(exons_by_tx, "exon_name"), file.path(out, "exons.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "")
write.table(parts(cds_by_tx, "cds_name"), file.path(out, "cds.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "")

gene_ids <- tryCatch(keys(txdb, keytype = "GENEID"), error = function(e) character(0))
write.table(data.frame(gene_id = gene_ids), file.path(out, "genes.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

saveDb(txdb, file.path(out, "txdb.sqlite"))

summary <- sprintf(
  '{"tool": "txdbmaker", "maker": "%s", "version": "%s", "seconds": %.3f, "transcripts": %d, "exons": %d, "cds": %d, "genes": %d, "warnings": %d, "circular": "%s"}',
  maker, version_of(maker), seconds, length(tx),
  length(unlist(exons_by_tx)), length(unlist(cds_by_tx)), length(gene_ids), warnings_seen,
  paste(circular, collapse = ",")
)
writeLines(summary, file.path(out, "summary.json"))
cat(summary, "\n")
