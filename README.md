
About GAZE
==========

GAZE is a tool for creating your own custom gene prediction
system. It assembles arbitrary gene prediction signal and content data
(given in the General Feature Format) into predictions of complete
gene structure. The assembly is performed according to a user-defined
configuration file (in XML), which includes the gene structure model to
specify how the candidate gene structures are validated and
scored. See the "Documentation" section of this file for pointers
towards other useful resources for finding out about GAZE.

Dependencies
============

GAZE has one primary dependency, the expat library for XML parsing:

http://sourceforge.net/projects/expat/


Installing
==========

make

[results in executable "gaze"]


Running
=======

gaze <options> sequence_name1 sequence_name2 ...

Some important options are:

-help                    :   Print a full list of command-line options
-structure_file <fname>  :   The GAZE configuration file
-dna_file <fname>        :   DNA sequence file (in Fasta)
-gff_file                :   Prediction data file (in GFF v2) 

The structure file argument is compulsory. The DNA and GFF file
arguments can be specified multiple times with different files,
allowing GAZE to process several sequences in one run. Appropriate
data for a given sequence name is the first DNA file entry with a
Fasta header matching the name (in the case of the DNA file) and all
GFF lines with a "sequence identifier" field matching the name (in the
case of the GFF files). A full user guide for GAZE can be found in the 
"docs" directory of this distribution.



Documentation
=============

For an introduction to the rationale behind gaze, please refer
to the following publication (which is also the preferred way
of citing GAZE):
Howe K.L., Chothia T. and Durbin R.
GAZE: a generic framework for the integration of gene prediction data
by dynamic programming.
Genome Research 12:1418-1427

More comprehensive documentation and examples for GAZE can be found in 
my PhD thesis:

ftp://ftp.sanger.ac.uk/pub/resources/theses/howe/

Kevin Howe
Wellcome Trust Sanger Institute
February 2003

