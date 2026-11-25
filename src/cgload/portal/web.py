"""The portal's web application: pages for people, JSON for programs (D-111).

Read-only by construction. Every route answers GET (and HEAD); nothing takes a
form that changes state, so there is no session, no cookie and nothing for a
cross-site request to forge. What a request may say is bounded before it is
used: a search is at most ``MAX_QUERY_LENGTH`` characters, a region at most
``MAX_REGION_TEXT``, a page number at most ``MAX_PAGE``, and a file is served
only if it is one of the four files of an assembly's directory that a genome
browser reads. Values from the database reach a page through Jinja's
autoescaping, and the one drawing made in Python escapes its labels.

Every response carries a Content-Security-Policy that allows scripts from the
portal itself and nowhere else, and no inline script, so text a file put in an
attribute cannot run in a reader's browser even if it got past the escaping.
"""

from __future__ import annotations

import logging
import re
import threading
from urllib.parse import quote

from flask import (
    Flask,
    Response,
    abort,
    jsonify,
    redirect,
    render_template,
    request,
    send_file,
    url_for,
)
from flask.json.provider import DefaultJSONProvider
from markupsafe import Markup
from werkzeug.exceptions import HTTPException

from cgload import __version__
from cgload.db.export_source import StoredFeature, load_context
from cgload.export import export_gff3
from cgload.portal.data import (
    ANNOTATION_FILENAME,
    REGION_LIMIT,
    FamilyView,
    Genome,
    PortalData,
    display_name,
    exceptions,
    genetic_code,
)
from cgload.portal.model import gene_model_svg
from cgload.portal.search import MAX_QUERY_LENGTH, RANK_LABELS
from cgload.portal.sequence import GENETIC_CODES, fasta_text
from cgload.query import QueryError, resolve_region

log = logging.getLogger("cgload.portal")

MAX_REGION_TEXT = 300
MAX_PAGE = 1000
PER_PAGE = 50
#: Bases either side of a feature in its browser view.
FLANK = 0.15

CONTENT_SECURITY_POLICY = "; ".join([
    "default-src 'self'",
    "script-src 'self'",
    # igv.js sets styles from script; no stylesheet comes from anywhere else.
    "style-src 'self' 'unsafe-inline'",
    "img-src 'self' data: blob:",
    "connect-src 'self'",
    "font-src 'self' data:",
    "worker-src 'self' blob:",
    "object-src 'none'",
    "base-uri 'none'",
    "form-action 'self'",
    "frame-ancestors 'none'",
])
SECURITY_HEADERS = {
    "Content-Security-Policy": CONTENT_SECURITY_POLICY,
    "X-Content-Type-Options": "nosniff",
    "X-Frame-Options": "DENY",
    "Referrer-Policy": "same-origin",
    "Cross-Origin-Opener-Policy": "same-origin",
    "Cross-Origin-Resource-Policy": "same-origin",
    "Permissions-Policy": "camera=(), microphone=(), geolocation=(), interest-cohort=()",
}

#: Where a cross-reference prefix leads. Only accessions of the expected form
#: become links, and each is quoted into the URL.
XREF_LINKS: dict[str, tuple[str, str]] = {
    "GeneID": (r"\d+", "https://www.ncbi.nlm.nih.gov/gene/{}"),
    "Genbank": (r"[A-Z]{1,6}_?\d+(?:\.\d+)?", "https://www.ncbi.nlm.nih.gov/nuccore/{}"),
    "GenBank": (r"[A-Z]{1,6}_?\d+(?:\.\d+)?", "https://www.ncbi.nlm.nih.gov/nuccore/{}"),
    "HGNC": (r"HGNC:\d+", "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{}"),
    "MIM": (r"\d+", "https://www.omim.org/entry/{}"),
    "UniProtKB/Swiss-Prot": (r"[A-Z0-9]{6,10}", "https://www.uniprot.org/uniprotkb/{}"),
    "UniProtKB/TrEMBL": (r"[A-Z0-9]{6,10}", "https://www.uniprot.org/uniprotkb/{}"),
    "FLYBASE": (r"FBgn\d+", "https://flybase.org/reports/{}"),
    "FlyBase": (r"FBgn\d+", "https://flybase.org/reports/{}"),
    "SGD": (r"S\d+", "https://www.yeastgenome.org/locus/{}"),
    "MGI": (r"MGI:\d+", "https://www.informatics.jax.org/marker/{}"),
    "ZFIN": (r"ZDB-GENE-[\d-]+", "https://zfin.org/{}"),
    "WormBase": (r"WBGene\d+", "https://wormbase.org/species/c_elegans/gene/{}"),
    "Ensembl": (r"ENS[A-Z]*\d+(?:\.\d+)?", "https://www.ensembl.org/id/{}"),
    "miRBase": (r"MI\d+", "https://www.mirbase.org/hairpin/{}"),
    "taxon": (r"\d+", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={}"),
    "GO": (r"GO:\d{7}", "https://www.ebi.ac.uk/QuickGO/term/{}"),
    "InterPro": (r"IPR\d{6}", "https://www.ebi.ac.uk/interpro/entry/InterPro/{}"),
}
#: Functional annotation: analysis -> (accession form, URL).
TERM_LINKS: dict[str, tuple[str, str]] = {
    "GO": (r"GO:\d{7}", "https://www.ebi.ac.uk/QuickGO/term/{}"),
    "Pfam": (r"PF\d{5}", "https://www.ebi.ac.uk/interpro/entry/pfam/{}"),
    "InterPro": (r"IPR\d{6}", "https://www.ebi.ac.uk/interpro/entry/InterPro/{}"),
    "KEGG_ko": (r"(?:ko:)?(?P<acc>K\d{5})", "https://www.kegg.jp/entry/{}"),
    "KEGG_Pathway": (r"(?:ko|map)\d{5}", "https://www.kegg.jp/pathway/{}"),
    "KEGG_Module": (r"M\d{5}", "https://www.kegg.jp/module/{}"),
    "KEGG_Reaction": (r"R\d{5}", "https://www.kegg.jp/entry/{}"),
    "EC": (r"\d+\.\d+\.\d+\.\d+", "https://enzyme.expasy.org/EC/{}"),
    "PANTHER": (r"PTHR\d+(?::SF\d+)?",
                "https://www.pantherdb.org/panther/family.do?clsAccession={}"),
}


def _link(form: str, template: str, value: str) -> str | None:
    match = re.fullmatch(form, value)
    if match is None:
        return None
    accession = match.groupdict().get("acc") or value
    return template.format(quote(accession, safe=":"))


def xref_link(value: str) -> str | None:
    prefix, _, accession = value.partition(":")
    rule = XREF_LINKS.get(prefix)
    if rule is None:
        return None
    if prefix == "GO":  # GO:0005524 is the accession whole
        accession = value
    return _link(rule[0], rule[1], accession)


def term_link(analysis: str, accession: str) -> str | None:
    rule = TERM_LINKS.get(analysis)
    return None if rule is None else _link(rule[0], rule[1], accession)


feature_name = display_name


def feature_note(item: StoredFeature) -> str:
    for key in ("product", "description", "Note"):
        values = item.attributes.get(key)
        if values:
            return values[0]
    return ""


def feature_json(item: StoredFeature) -> dict:
    return {
        "row": item.row_id,
        "id": item.source_id,
        "name": feature_name(item),
        "type": item.feature_type,
        "source_type": item.source_type,
        "seqid": item.seqid,
        "start": item.start,
        "end": item.end,
        "strand": item.strand,
        "segments": [[segment.start, segment.end, segment.phase] for segment in item.segments],
        "parent": item.parent_source_id,
        "attributes": item.attributes,
    }


def genome_json(genome: Genome) -> dict:
    return {
        "id": genome.assembly_id,
        "organism": genome.organism,
        "strain": genome.strain,
        "assembly": genome.name,
        "accession": genome.accession,
        "sequences": genome.sequences,
        "features": genome.features,
        "rows": genome.rows,
        "loaded_by": genome.cgload_version,
        "loaded_at": genome.loaded_at.isoformat() if genome.loaded_at else None,
        "searchable": genome.indexed is not None,
        "files": {name: url_for("data_file", assembly_id=genome.assembly_id, name=name)
                  for name in ("assembly.fa", "assembly.fa.fai", ANNOTATION_FILENAME,
                               ANNOTATION_FILENAME + ".tbi")
                  if genome.path(name) is not None},
    }


class SafeJSONProvider(DefaultJSONProvider):
    """JSON with <, >, & and ' written as escapes, so a value from a file is
    inert even where a browser is made to read JSON as HTML."""

    sort_keys = False

    def dumps(self, obj, **kwargs) -> str:
        text = super().dumps(obj, **kwargs)
        return (text.replace("<", "\\u003c").replace(">", "\\u003e")
                .replace("&", "\\u0026").replace("'", "\\u0027"))


#: Row and assembly ids are 32-bit integers in the schema; a larger number in a
#: URL names nothing, and must not reach the driver as an overflowing value.
MAX_ID = 2**31 - 1


def create_app(portal: PortalData, *, title: str = "cgload portal",
               count_in_background: bool = True) -> Flask:
    """The Flask application serving ``portal``."""
    app = Flask(__name__)
    app.json_provider_class = SafeJSONProvider
    app.json = SafeJSONProvider(app)
    app.config.update(
        MAX_CONTENT_LENGTH=4096,
        JSON_SORT_KEYS=False,
        SEND_FILE_MAX_AGE_DEFAULT=3600,
        TEMPLATES_AUTO_RELOAD=False,
    )
    app.extensions["cgload_portal"] = portal

    if count_in_background:
        threading.Thread(target=_count, args=(portal,), daemon=True,
                         name="cgload-portal-count").start()

    @app.context_processor
    def _globals():
        return {"portal_title": title, "version": __version__, "feature_name": feature_name,
                "feature_note": feature_note, "xref_link": xref_link, "term_link": term_link,
                "rank_labels": RANK_LABELS, "max_query": MAX_QUERY_LENGTH,
                "genomes_all": portal.genomes}

    @app.after_request
    def _headers(response: Response) -> Response:
        for name, value in SECURITY_HEADERS.items():
            response.headers.setdefault(name, value)
        if request.path.startswith(("/g/", "/api/")) and "/files/" not in request.path:
            response.headers.setdefault("Cache-Control", "no-store")
        return response

    # -- arguments -----------------------------------------------------------

    def text_arg(name: str, limit: int) -> str:
        value = request.args.get(name, "")
        if len(value) > limit:
            abort(400, description=f"'{name}' is at most {limit} characters")
        return value.strip()

    def page_arg() -> int:
        raw = request.args.get("page", "1")
        if not raw.isdigit() or not 1 <= int(raw) <= MAX_PAGE:
            abort(400, description=f"'page' is a number from 1 to {MAX_PAGE}")
        return int(raw)

    def genome_or_404(assembly_id: int) -> Genome:
        try:
            return portal.genome(assembly_id)
        except LookupError as exc:
            abort(404, description=str(exc))

    def family_or_404(genome: Genome, row: int) -> FamilyView:
        if row > MAX_ID:
            abort(404, description=f"no feature {row} in {genome.label}")
        try:
            return portal.family(genome, row)
        except LookupError as exc:
            abort(404, description=str(exc))

    def feature_url(genome: Genome, item: StoredFeature) -> str:
        return url_for("feature_page", assembly_id=genome.assembly_id, row=item.row_id)

    def browser_view(genome: Genome, seqid: str, start: int, end: int) -> dict | None:
        files = {name: genome.path(name) for name in
                 ("assembly.fa", "assembly.fa.fai", ANNOTATION_FILENAME,
                  ANNOTATION_FILENAME + ".tbi")}
        if files["assembly.fa"] is None or files["assembly.fa.fai"] is None:
            return None
        pad = max(int((end - start + 1) * FLANK), 50)
        view = {
            "name": genome.label,
            "fasta": url_for("data_file", assembly_id=genome.assembly_id, name="assembly.fa"),
            "fai": url_for("data_file", assembly_id=genome.assembly_id,
                           name="assembly.fa.fai"),
            "locus": f"{seqid}:{max(start - pad, 1)}-{end + pad}",
            "locate": url_for("locate", assembly_id=genome.assembly_id) + "?q=$FEATURE$",
        }
        if files[ANNOTATION_FILENAME] and files[ANNOTATION_FILENAME + ".tbi"]:
            view["gff"] = url_for("data_file", assembly_id=genome.assembly_id,
                                  name=ANNOTATION_FILENAME)
            view["tbi"] = url_for("data_file", assembly_id=genome.assembly_id,
                                  name=ANNOTATION_FILENAME + ".tbi")
        return view

    # -- pages -----------------------------------------------------------------

    @app.get("/")
    def index():
        return render_template("index.html", genomes=portal.genomes())

    @app.get("/about")
    def about():
        return render_template("about.html")

    @app.get("/search")
    def search_all():
        query = text_arg("q", MAX_QUERY_LENGTH)
        groups = []
        if query:
            for genome in portal.genomes():
                try:
                    result, rows, total = portal.search(genome, query, per_page=5, limit=200)
                except QueryError as exc:
                    abort(400, description=str(exc))
                if rows:
                    groups.append((genome, rows, total, result.truncated))
        return render_template("search_all.html", query=query, groups=groups)

    @app.get("/g/<int:assembly_id>/")
    def genome_page(assembly_id: int):
        genome = genome_or_404(assembly_id)
        return render_template(
            "genome.html", genome=genome, types=portal.counts_by_type(assembly_id),
            sequences=portal.longest_sequences(assembly_id),
            analyses=portal.analyses(assembly_id), files=portal.source_files(assembly_id),
            receipt=genome.receipt())

    @app.get("/g/<int:assembly_id>/search")
    def search_page(assembly_id: int):
        genome = genome_or_404(assembly_id)
        query = text_arg("q", MAX_QUERY_LENGTH)
        page = page_arg()
        if not query:
            return redirect(url_for("genome_page", assembly_id=assembly_id))
        try:
            result, rows, total = portal.search(genome, query, page=page, per_page=PER_PAGE)
        except QueryError as exc:
            abort(400, description=str(exc))
        if total == 1 and page == 1 and rows:
            return redirect(feature_url(genome, rows[0][1]))
        pages = max((total + PER_PAGE - 1) // PER_PAGE, 1)
        return render_template("search.html", genome=genome, query=query, result=result,
                               rows=rows, total=total, page=page, pages=pages)

    @app.get("/g/<int:assembly_id>/id/<path:source_id>")
    def by_source_id(assembly_id: int, source_id: str):
        genome = genome_or_404(assembly_id)
        if len(source_id) > MAX_QUERY_LENGTH * 3:
            abort(400, description="that ID is too long")
        try:
            row = portal.row_of_id(genome, source_id)
        except LookupError as exc:
            abort(404, description=str(exc))
        return redirect(url_for("feature_page", assembly_id=assembly_id, row=row))

    @app.get("/g/<int:assembly_id>/f/<int:row>")
    def feature_page(assembly_id: int, row: int):
        genome = genome_or_404(assembly_id)
        family = family_or_404(genome, row)
        item = family.selected
        try:
            sequences = portal.sequences_of(family, item)
        except (FileNotFoundError, ValueError, KeyError) as exc:
            log.warning("no sequence for %s: %s", item.source_id, exc)
            sequences = {}
        # The drawing escapes every label itself (model.py), and a test feeds it
        # a name that is a script (tests/test_portal_security.py).
        drawing = gene_model_svg(family, lambda member: feature_url(genome, member))
        model = Markup(drawing)  # noqa: S704  # nosec B704
        code = genetic_code(item) if item.feature_type == "CDS" else None
        return render_template(
            "feature.html", genome=genome, family=family, item=item, model=model,
            sequences=sequences, notes=exceptions(item),
            code=(code, GENETIC_CODES[code][2]) if code else None,
            browser=browser_view(genome, family.top.seqid, family.top.start, family.top.end),
            members=family.depth_first(), feature_url=lambda m: feature_url(genome, m))

    @app.get("/g/<int:assembly_id>/region")
    def region_page(assembly_id: int):
        genome = genome_or_404(assembly_id)
        text = text_arg("r", MAX_REGION_TEXT)
        if not text:
            return redirect(url_for("genome_page", assembly_id=assembly_id))
        try:
            region, listed, more = portal.region(genome, text)
        except QueryError as exc:
            abort(400, description=str(exc))
        return render_template(
            "region.html", genome=genome, region=region, listed=listed, more=more,
            too_long=region.end - region.start + 1 > REGION_LIMIT, region_limit=REGION_LIMIT,
            browser=browser_view(genome, region.seqid, region.start, region.end),
            feature_url=lambda m: feature_url(genome, m))

    # -- downloads -------------------------------------------------------------

    @app.get("/g/<int:assembly_id>/f/<int:row>/<kind>.fa")
    def feature_fasta(assembly_id: int, row: int, kind: str):
        genome = genome_or_404(assembly_id)
        family = family_or_404(genome, row)
        if kind not in ("genomic", "spliced", "cds", "protein"):
            abort(404, description=f"no {kind} sequence")
        try:
            sequences = portal.sequences_of(family, family.selected)
        except (FileNotFoundError, ValueError, KeyError) as exc:
            abort(404, description=f"the sequence cannot be read: {exc}")
        if kind not in sequences:
            abort(404, description=f"a {family.selected.feature_type} has no {kind} sequence")
        item = family.selected
        name = item.source_id or f"row{item.row_id}"
        header = (f"{name} {kind} {item.seqid}:{item.start}-{item.end}({item.strand}) "
                  f"{genome.organism}")
        return _download(fasta_text(header, sequences[kind]), f"{_safe(name)}.{kind}.fa")

    @app.get("/g/<int:assembly_id>/f/<int:row>/family.gff3")
    def family_gff3(assembly_id: int, row: int):
        genome = genome_or_404(assembly_id)
        family = family_or_404(genome, row)
        with portal.engine.connect() as connection:
            context = load_context(connection, assembly_id)
        text, _ = export_gff3(
            family.members, version=context.cgload_version,
            vocabulary_hash=context.vocabulary_hash,
            sequence_regions={family.top.seqid: context.sequence_regions[family.top.seqid]},
            comments=[f"cgload portal: the family of {family.top.source_id or 'a feature'}"])
        name = family.top.source_id or f"row{family.top.row_id}"
        return _download(text, f"{_safe(name)}.gff3")

    @app.get("/g/<int:assembly_id>/region.fa")
    def region_fasta(assembly_id: int):
        genome = genome_or_404(assembly_id)
        text = text_arg("r", MAX_REGION_TEXT)
        try:
            with portal.engine.connect() as connection:
                region = resolve_region(connection, assembly_id, text)
        except QueryError as exc:
            abort(400, description=str(exc))
        if region.end - region.start + 1 > REGION_LIMIT:
            abort(400, description=f"a region download is at most {REGION_LIMIT:,} bp; "
                                   "the whole assembly is a file of its own")
        try:
            bases = portal.sequence(genome).fetch(region.seqid, region.start, region.end)
        except (FileNotFoundError, ValueError, KeyError) as exc:
            abort(404, description=f"the sequence cannot be read: {exc}")
        header = f"{region.seqid}:{region.start}-{region.end} {genome.organism}"
        return _download(fasta_text(header, bases),
                         f"{_safe(region.seqid)}_{region.start}-{region.end}.fa")

    @app.get("/g/<int:assembly_id>/files/<name>")
    def data_file(assembly_id: int, name: str):
        genome = genome_or_404(assembly_id)
        path = genome.path(name)
        if path is None:
            abort(404, description="no such file")
        response = send_file(path, conditional=True, download_name=name,
                             mimetype="application/octet-stream")
        response.headers["Accept-Ranges"] = "bytes"
        return response

    # -- JSON ----------------------------------------------------------------

    @app.get("/healthz")
    def health():
        with portal.engine.connect() as connection:
            connection.exec_driver_sql("SELECT 1")
        return jsonify(status="ok", version=__version__, genomes=len(portal.genomes()))

    @app.get("/api/v1/genomes")
    def api_genomes():
        return jsonify(genomes=[genome_json(genome) for genome in portal.genomes()])

    @app.get("/api/v1/genomes/<int:assembly_id>/search")
    def api_search(assembly_id: int):
        genome = genome_or_404(assembly_id)
        query = text_arg("q", MAX_QUERY_LENGTH)
        page = page_arg()
        try:
            result, rows, total = portal.search(genome, query, page=page, per_page=PER_PAGE)
        except QueryError as exc:
            abort(400, description=str(exc))
        return jsonify(query=query, total=total, truncated=result.truncated,
                       searchable=result.indexed, page=page,
                       results=[{"matched_by": RANK_LABELS[match.rank], "field": match.field,
                                 "feature": feature_json(item),
                                 "top": None if top is None else feature_json(top)}
                                for match, item, top in rows])

    @app.get("/api/v1/genomes/<int:assembly_id>/features/<int:row>")
    def api_feature(assembly_id: int, row: int):
        genome = genome_or_404(assembly_id)
        family = family_or_404(genome, row)
        return jsonify(feature=feature_json(family.selected), top=feature_json(family.top),
                       family=[feature_json(item) for item in family.members],
                       annotations={str(key): value for key, value in
                                    family.annotations.items()})

    @app.get("/api/v1/genomes/<int:assembly_id>/region")
    def api_region(assembly_id: int):
        genome = genome_or_404(assembly_id)
        text = text_arg("r", MAX_REGION_TEXT)
        try:
            region, listed, more = portal.region(genome, text)
        except QueryError as exc:
            abort(400, description=str(exc))
        return jsonify(region={"seqid": region.seqid, "start": region.start, "end": region.end},
                       features=[feature_json(item) for item in listed], truncated=more)

    @app.get("/g/<int:assembly_id>/locate")
    def locate(assembly_id: int):
        """For the genome browser's search box: where the best match lies."""
        genome = genome_or_404(assembly_id)
        query = text_arg("q", MAX_QUERY_LENGTH)
        if not query:
            return jsonify([])
        try:
            _, rows, _ = portal.search(genome, query, per_page=1, limit=50)
        except QueryError:
            return jsonify([])
        return jsonify([{"chromosome": item.seqid, "start": item.start, "end": item.end,
                         "name": feature_name(item)} for _, item, _ in rows])

    # -- errors --------------------------------------------------------------

    @app.errorhandler(HTTPException)
    def http_error(exc: HTTPException):
        if request.path.startswith("/api/") or request.path.endswith("/locate"):
            return jsonify(error=exc.description, status=exc.code), exc.code
        return render_template("error.html", code=exc.code, name=exc.name,
                               message=exc.description), exc.code

    @app.errorhandler(Exception)
    def server_error(exc: Exception):
        log.error("error serving %s", request.path, exc_info=exc)
        if request.path.startswith("/api/"):
            return jsonify(error="internal error", status=500), 500
        return render_template("error.html", code=500, name="Internal error",
                               message="Something went wrong; it has been logged."), 500

    return app


def _count(portal: PortalData) -> None:
    try:
        portal.count_features()
    except Exception:  # noqa: BLE001 - counting is a convenience; the pages work without it
        log.exception("counting features failed")


def _safe(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name)[:120] or "feature"


def _download(text: str, filename: str) -> Response:
    response = Response(text, mimetype="text/plain")
    response.headers["Content-Disposition"] = f'attachment; filename="{filename}"'
    return response
