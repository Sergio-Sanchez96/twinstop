# Twinstop
Twinstop is a SECIS-independent pipeline useful to identify selenoproteins.
It is based on the evolutionary conservation around the UGA-coding selenocysteine between two homologous 
selenoproteins from two closely related transcriptomes.

Twinstop pipeline can be divided in two phases:
1. Whole CDS transcripts alignments
2. Identification of selenoprotein candidates

First phase starts with a tBLASTx where the two nucleotide transcriptomes are faced to get high scored local
alignments of the homologous regions. This generates milions of results, some with more than one ORFs others
with incomplete ORF. But, more importantly, a lot of redundancy produced by having repeated sequences in 
different ORFs.

So, Twinstop will ensure having the whole sequence of all homologous genes, using different strategies to get 
the best transcript of each gene:
- By keeping only the best scored ORF in each tBLASTx hit
- By keeping only the best scored hit, among those that overlap, per pair of transcripts (Q/S)
- By extending the ORFs both up/downstream until completing the CDS or finishing the transcript

Second phase, to identify bona-fide selenoproteins we use different filters:
- UGA filter, hits with an aligned TGA in the nucleotide sequences
- Conservation filter, hits with an up/downstream score greater than 50/150 (based on controls)
- Mutation filter, hits with more than 5 aa changes both up/downstream

The annotation of candidates is done by blastp between the protein sequences of the queries 
against UniRef50 database.
