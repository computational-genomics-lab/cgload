[Date]

The Editor, *BMC Bioinformatics*

Dear Editor,

We submit "cgload: loading genome annotations into a database and checking
that every feature arrived" for consideration as a Software article.

Genome annotations in GFF3 are routinely loaded into databases, and a loader
can drop or alter features without a message. cgload loads an assembly and its
annotation into SQLite or PostgreSQL in one transaction, compares what it
stored with an independent count of the input file, and rolls back a load that
falls short. We compared it with gffutils, txdbmaker and python-chado on 30
RefSeq genomes (18.0 million features), under a protocol fixed and tagged
before the first run. cgload kept every feature and declared its 12 repairs.
Each of the other loaders stored some features in another form, or not at all,
without a message: gffutils 1,072,125 with the setting its documentation
offers for such files, and python-chado 430,836 in the 21 genomes it finished.
We also report where cgload is weaker: it needs far more memory than gffutils
and answers queries more slowly.

The software is open source under the MIT license, at
https://github.com/computational-genomics-lab/cgload, and reviewers can
install and test it anonymously from there. The source of the version
described is also supplied as Additional file 6. The protocol, scripts and
results of the benchmark are in the same repository.

A large language model (Claude, Anthropic) was used to write code, benchmark
scripts and a draft of the text, under our direction. This is documented in
the manuscript, in the section "Use of a large language model", and in
Additional file 5. All authors have reviewed the code, the results and the
text, and take full responsibility for them.

The manuscript has not been published and is not under consideration
elsewhere. All authors have approved it and have no competing interests.

Yours sincerely,

Sucheta Tripathy (corresponding author)\
CSIR-Indian Institute of Chemical Biology, Kolkata, India\
[e-mail]
