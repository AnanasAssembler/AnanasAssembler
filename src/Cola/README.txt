Cola consists of an efficient implementation of a collection of sequence
alignment algorithms, extending the Smith-Waterman and Smith-Waterman-Gap-Affine
methods by the ability to apply a scoring function that is an arbitrary
function of the size of consecutive nucleotide matches.

The main point for a user to run cola is through runCola which accepts the following arguments:
-q  : Query sequence in FASTA format
-t  : Target sequence in FASTA format
-a  : Aligner type - Choose a number from 1 to 4 for the following modes - 1 : NSGA , 2 : NS , 3 : SWGA, 4 : SW 
-o  : Aligner gap open penalty
-m  : Aligner mismatch penalty
-e  : Aligner gap extension penalty
-all: align all to all sequences (default is false)

N.B. Only the first 3 arguments are compulsory and the rest are optional. 
If the aligner penalty parameters are not provided by user, system defaults will be used. 

See following examples for more detail:

To run in one of 4 modes:
1) Nonlinear Scoring Gap-Affine (NSGA): 
./runCola -t samples/human.X.part.superShort.fasta -q samples/dog.X.part.superShort.fasta -a 1 -o 200 -m 8 -e 20

2) Nonlinear Scoring (NS): 
./runCola -t samples/human.X.part.superShort.fasta -q samples/dog.X.part.superShort.fasta -a 2 -o 100 -m 8

3) Smith-Waterman Gap-Affine (SWGA): 
./runCola -t samples/human.X.part.superShort.fasta -q samples/dog.X.part.superShort.fasta -a 3 -o 6 -m 1 -e 1

4) Smith-Waterman (SW): 
./runCola -t samples/human.X.part.superShort.fasta -q samples/dog.X.part.superShort.fasta -a 4 -o 1 -m 1


See Cola manuscript for more detail on the underlying algorithm and methodology.
