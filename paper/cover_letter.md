25 November 2026

The Editor, *BMC Bioinformatics*

Dear Editor,

We submit "cgload: genome annotations loaded with a completeness check and
served through a read-only web portal" for consideration as a Research article
in *BMC Bioinformatics*.

Most biologists meet a genome annotation through a web portal built on a
database, and the database is filled by a loader that can drop or alter
features without a message. cgload covers both steps. It loads an assembly and
its annotation into SQLite or PostgreSQL in one transaction, compares what it
stored with an independent count of the input file, and rolls back a load that
falls short. One command, `cgload portal`, then takes the genomes to a
read-only portal, installing BLAST+ if the computer has none: search by gene
name, product or functional term, gene pages with a genome browser, sequences
and proteins, BLAST against the genomes and their proteins, and the record of
how each genome was loaded.

We compared cgload with gffutils, txdbmaker and python-chado on 30 RefSeq
genomes (18.0 million features), under a protocol fixed and tagged before the
first run. cgload kept every feature in every field compared and declared its
12 repairs; each of the other loaders stored some features in another form, or
not at all, without a message. We are explicit about what cgload's own check
cannot see: it compares counts, so it catches a lost or invented feature but
not a renamed one or a changed value. One portal then served all 30 genomes
from a single database, with a median of {{portal:gene_page_median}} ms for a
gene's page, and every BLAST search that finished found its query's own
protein or place among its hits. It was built against the problems genome
portals are known for: it opens the database read-only, binds every value a
visitor sends, starts no program but BLAST+, with arguments it chooses, and
loads no script from another site, and we tested it with injection strings
and thousands of random requests. We also report where
cgload is weaker: it needs far more memory than gffutils, answers single
queries from the command line more slowly, and serves the portal from one
process.

We believe the work suits *BMC Bioinformatics* because it offers an
open-source tool that takes a group from annotation files to a checked
database and a portal its users can work with, and an open benchmark of tools
that many groups already use. The protocol, scripts and results are public, so
others can put the same questions to further loaders and files.

The software is open source under the MIT license, at
https://github.com/computational-genomics-lab/cgload, and reviewers can
install it and run the portal on their own computer anonymously from there.
The source code, with the benchmark protocol, scripts and results, is also
supplied as Additional file 7.

On journal policy: we used large language models, as the Methods section "Use
of large language models" and Additional file 5 describe. Claude (Anthropic)
wrote code, including the portal, benchmark scripts and a draft of the text
under our direction, and checked drafts and the submission files for errors;
ChatGPT (OpenAI) also reviewed drafts for errors. All authors have reviewed
the code, the results and the text, and take full responsibility for them.

The manuscript has not been published and is not under consideration
elsewhere. All authors have approved it and declare no competing interests.

<!-- The submitted letter also gives every author's e-mail address; they are
left out of this public copy. -->

Yours sincerely,

Sucheta Tripathy, on behalf of all authors\
CSIR-Indian Institute of Chemical Biology, Kolkata, India\
tsucheta@gmail.com
