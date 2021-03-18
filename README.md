## BaCelLo: prediction of protein subcellular localization

#### Publication

Andrea Pierleoni, Pier Luigi Martelli, Piero Fariselli, Rita Casadio, BaCelLo: a balanced subcellular localization predictor, Bioinformatics, Volume 22, Issue 14, 15 July 2006, Pages e408-e416.

#### Requirements
- python=3.8.3
- numpy=1.18.1
- biopython

#### Installation and configuration

Before running bacello you need to set and export a variable named BACELLO_HOME to point to the program installation dir:
```
$ export BACELLO_HOME='/path/to/bacello'
```

#### Usage

To run the program:

```
./bacello.py -f fasta_file -p blast_pssm_file -k [A,P,F] -o out_file

```

where fasta_file contains the input protein sequence in FASTA format, blast_pssm_file is the PSSM file as generated by PSI-BLAST (-out_ascii_pssm) aligning the input sequence against the UniProtKB/SwissProt database (to be run externally) and the last parameter is the kingdom the input protein belongs to [A: animal, P: plant, F: fungi].

Running on a testing animal sequence:

```
./bacello.py -f testdata/T96060011003.fasta -p testdata/T96060011003.pssm -k A -o testdata/T96060011003.gff3
```
