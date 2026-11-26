# The portal

`cgload serve` turns a cgload database into a website: people search genomes by
name, open a gene with its transcripts drawn to scale, browse any region in a
genome browser and download sequences; programs get the same answers as JSON.
It only reads. Decisions D-111 to D-114 in `docs/decisions.md` say why it is
built the way it is.

## Two commands

The portal comes with cgload; nothing else to install. After loading:

```console
$ cgload index                       # once after loading: the name index
$ cgload serve
cgload serve: the portal is running, read-only, with 1 genome.
  Open it on this computer:   http://127.0.0.1:8000/
  or from another computer:   http://10.0.0.239:8000/
  Anyone who can reach this computer can open it. To keep it to this computer: --host 127.0.0.1;
  to serve it beyond your network, put a web server with HTTPS in front (docs/portal.md).
  Stop it with Ctrl-C.
```

Open one of the addresses it prints in a web browser. It listens on every
network interface of the computer, so the other computers of your network can
open it too (D-115); `--host 127.0.0.1` keeps it to this computer alone.

**If something is in the way:**

- *Port 8000 is taken* by another program: `cgload serve` uses the next free
  port and says which. If you asked for a port with `--port`, it tells you the
  port is taken instead, and you choose another.
- *Another computer gets "connection refused" or waits forever:* the server's
  firewall is closing the port. On Ubuntu, for example,
  `sudo ufw allow from 10.0.0.0/24 to any port 8000 proto tcp` opens it to a
  network of that range only. Also check that the portal was not started with
  `--host 127.0.0.1`.
- *Names find nothing:* run `cgload index` (and again after loading more).
- *A genome's sequences and browser are missing:* the data directory has moved
  since loading; give it with `--data-dir`.
- *Unpublished data:* anyone who can reach the computer can read the portal.
  Keep it to this computer with `--host 127.0.0.1`, or put it behind a web
  server that asks for a password (below).

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
| `--host` | `0.0.0.0`, every interface | `127.0.0.1`: this computer only, as behind a web server |
| `--port` | `8000` | the next free port if 8000 is taken, unless `--port` is given |
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
does, which then is the only thing that reaches it: start it with
`--host 127.0.0.1`. With nginx:

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
ExecStart=/srv/cgload/env/bin/cgload serve --host 127.0.0.1 --port 8000 --data-dir /srv/cgload/data
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
lock, so it answers about as many requests a second as one core allows: on the
benchmark's server, with the 30 genomes in one database, about 40 a second,
whether 8 or 32 users asked at once (`benchmark/PROTOCOL.md` section 13).
The portal only reads and keeps no state, so more processes behind the same web server serve more
people: run one `cgload serve --host 127.0.0.1 --port 800N` per process, from
a systemd template unit `cgload@.service` whose command ends `--port %i`, and
list them in nginx:

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
- It listens on every interface unless told `--host 127.0.0.1`, and says so
  when it starts; who may read the data is decided by who can reach the
  computer, or by the web server in front of it.
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
