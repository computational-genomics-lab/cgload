# The portal

`cgload serve` turns a cgload database into a website: people search genomes by
name, open a gene with its transcripts drawn to scale, browse any region in a
genome browser and download sequences; programs get the same answers as JSON.
It only reads. Decisions D-111 to D-114 in `docs/decisions.md` say why it is
built the way it is.

## Two commands

```console
$ pip install "cgload[portal]"
$ cgload index                       # once after loading: the name index
$ cgload serve                       # http://127.0.0.1:8000/
```

`cgload index` builds the search index of every loaded assembly (or of those
given with `--assembly-id`) and checks it by counting, as `load` checks a load:
for each attribute it indexes, the features that carry the attribute are counted
in the database and in the index, and an assembly whose counts differ is left
without an index. Run it again after loading more; it rebuilds only what it is
asked to. Without an index, the portal still finds features by exact ID and by
functional term, and says that names are not searchable.

`cgload serve` takes the same `--db-url` (or `CGLOAD_DB_URL`) as every other
command, and `--data-dir` if the data directory has moved since loading: the
portal reads each assembly, and the annotation files a genome browser needs,
from there.

| Option | Default | |
|---|---|---|
| `--host` | `127.0.0.1` | `0.0.0.0` listens on every interface |
| `--port` | `8000` | |
| `--threads` | `8` | requests served at once |
| `--title` | `cgload portal` | shown on every page |
| `--url-prefix` | none | the path behind a web server, e.g. `/genomes` |

## What people can do

- **Search** one genome or all of them, by a gene's name or a synonym, a locus
  tag, a protein or transcript ID, a word of its product or description, a
  cross-reference (`GeneID:672`, or just `672`), or a functional accession that
  eggNOG-mapper or InterProScan assigned (`GO:0005524`, `K00799`, `PF00069`,
  `IPR000719`, `2.7.11.1`). Exact IDs come first, then exact names in any
  case, then names that begin with the search, then features whose text has
  every word, then functional matches; within each, genes before transcripts
  before CDS.
- **Open a gene**: its transcripts drawn to scale (exons, coding parts,
  introns); the genome browser at its locus; its attributes, with
  cross-references linked to NCBI Gene, HGNC, UniProt, FlyBase, SGD, MGI,
  QuickGO, InterPro, Pfam, KEGG and others; its functional annotation; the
  whole family as a table; and its sequences.
- **Sequences**: a feature's genomic span; a transcript spliced; a CDS from its
  first whole codon, and its protein, translated with the NCBI genetic code
  the file names in `transl_table` (D-113). Where the file declares an
  exception a plain translation does not apply, the page says so.
- **Browse a region** (`chr:start-end`): the genome browser (igv.js) reading
  the assembly and the annotation export that `load` wrote, and a table of the
  genes there.
- **Download** a gene's family as GFF3, sequences as FASTA, a region up to
  5 Mb, or the whole assembly and annotation.
- **See how a genome was loaded**: the load's receipt (`manifest.json`): which
  files, their checksums, the checks `load` ran and their outcome, and any
  coordinates cgload repaired.

## For programs

| URL | Returns |
|---|---|
| `/api/v1/genomes` | every assembly, its counts and its file URLs |
| `/api/v1/genomes/<id>/search?q=...&page=n` | matches, 50 a page, each with its family's top feature |
| `/api/v1/genomes/<id>/features/<row>` | a feature, its whole family and its functional annotation |
| `/api/v1/genomes/<id>/region?r=chr:start-end` | the top-level features in a region |
| `/g/<id>/files/<name>` | `assembly.fa`, `assembly.fa.fai`, `annotation.gff3.gz`, `.tbi`, with byte ranges |
| `/healthz` | `status: ok`, the version and the number of genomes, while the database answers |

A feature's `row` is the database row of its first segment, stable for as long
as the assembly is not reloaded. `/g/<id>/id/<ID>` redirects from a file's own
ID to the feature's page.

## Putting it on a network

`cgload serve` runs waitress, a production WSGI server, and is safe to expose
as it is; it does not speak HTTPS, and should sit behind a web server that
does. With nginx:

```nginx
server {
    listen 443 ssl;
    server_name genomes.example.org;
    # ssl_certificate ...; ssl_certificate_key ...;
    client_max_body_size 4k;           # the portal takes no uploads
    location / {
        proxy_pass http://127.0.0.1:8000;
        proxy_set_header Host $host;
        proxy_read_timeout 120s;
    }
}
```

And a systemd unit that starts it with the machine and restarts it if it stops:

```ini
[Unit]
Description=cgload portal
After=network.target

[Service]
User=cgload
Environment=CGLOAD_DB_URL=sqlite:////srv/cgload/cgload.db
ExecStart=/srv/cgload/env/bin/cgload serve --data-dir /srv/cgload/data --threads 8
Restart=on-failure
ProtectSystem=strict
PrivateTmp=true
NoNewPrivileges=true

[Install]
WantedBy=multi-user.target
```

`ProtectSystem=strict` makes the whole file system, the database and data
directory with it, read-only to the process, so even a fault in the portal
could not change them; `PrivateTmp` gives it the scratch space the server
uses to send large files. On PostgreSQL, give the portal a role that can only
`SELECT`.

Access control, where the data are not public, belongs to the web server
(`auth_basic`, or the institution's single sign-on); the portal has no
accounts.

### Serving more people at once

One `cgload serve` is one process. Its threads share Python's interpreter
lock, so it answers about as many requests a second as one core allows
(`benchmark/PROTOCOL.md` section 13 measures it). The portal only reads and
keeps no state, so more processes behind the same web server serve more
people: run one `cgload serve --port 800N` per process, from a systemd
template unit `cgload@.service` whose command ends `--port %i`, and list them
in nginx:

```nginx
upstream cgload {
    server 127.0.0.1:8001;
    server 127.0.0.1:8002;
    server 127.0.0.1:8003;
    server 127.0.0.1:8004;
}
# and in the server block: proxy_pass http://cgload;
```

Readers can share one SQLite file safely; on PostgreSQL each process keeps
its own connections.

## How it is kept safe

- It opens the database read-only: SQLite in `mode=ro` with `query_only`,
  PostgreSQL in read-only transactions. A test reads every page and checks the
  database file is unchanged.
- It answers GET and HEAD only, and has no session or cookie.
- Every value from a request is a bound parameter, and its length is bounded.
- Text from the annotation files is escaped wherever it is shown; a test loads
  a file whose names are scripts and checks every page.
- Every response carries a Content-Security-Policy that allows scripts from the
  portal itself only, and no inline script.
- Of each assembly's directory, four files may be read, by name.
- CI runs the security tests on Linux, macOS and Windows, and bandit on the
  code. The benchmark (PROTOCOL.md section 13) sends injection strings and
  3,000 random requests to a running server that holds all 30 genomes.

## What it does not do

It does not edit, upload, align (no BLAST), or keep accounts, and it fetches
nothing from other sites: igv.js is served by the portal itself (D-114). A
protein is cgload's translation of the CDS rows, not a curated sequence.
