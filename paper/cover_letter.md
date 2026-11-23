23 November 2026

The Editor, *BMC Bioinformatics*

Dear Editor,

We submit "cgload: loading genome annotations into a database and checking
that every feature arrived" for consideration as a Research article in
*BMC Bioinformatics*.

Genome annotations in GFF3 are routinely loaded into databases, and a loader
can drop or alter features without a message. cgload loads an assembly and its
annotation into SQLite or PostgreSQL in one transaction, compares what it
stored with an independent count of the input file, and rolls back a load that
falls short. We compared it with gffutils, txdbmaker and python-chado on 30
RefSeq genomes (18.0 million features), under a protocol fixed and tagged
before the first run. cgload kept every feature in every field compared and
declared its 12 repairs. We are explicit about what its own check can and
cannot see: it compares counts, so it catches a lost or invented feature but
not a renamed one or a changed value. Each of the other loaders stored some
features in another form, or not at all, without a message: gffutils 1,072,125
with the setting its documentation offers for such files, and python-chado
430,836 in the 21 genomes it finished. We also report where cgload is weaker:
it needs far more memory than gffutils and answers queries more slowly.

We believe the work suits *BMC Bioinformatics* because it offers both an
open-source tool and an open benchmark of tools that many groups already use to
turn annotation files into databases. The protocol, scripts and results are
public, so others can put the same questions to further loaders and files.

The software is open source under the MIT license, at
https://github.com/computational-genomics-lab/cgload, and reviewers can
install and test it anonymously from there. The source code, with the
benchmark protocol, scripts and results, is also supplied as Additional file 6.

On journal policy: we used large language models, as the Methods section "Use
of large language models" and Additional file 5 describe. Claude (Anthropic)
wrote code, benchmark scripts and a draft of the text under our direction, and
checked drafts and the submission files for errors; ChatGPT (OpenAI) also
reviewed drafts for errors. All authors have reviewed the code, the results and
the text, and take full responsibility for them.

The manuscript has not been published and is not under consideration
elsewhere. All authors have approved it and declare no competing interests.

<!-- The submitted letter also gives every author's e-mail address; they are
left out of this public copy. -->

Yours sincerely,

Sucheta Tripathy, on behalf of all authors\
CSIR-Indian Institute of Chemical Biology, Kolkata, India\
tsucheta@gmail.com
