# GenGTools.Extractor
Supplementary tool for extracting sequences with high protein-coding potential from genome assembly graph.

Requirements: python2, sklearn, numpy, pandas, matplotlib

Usage: GenGTools-extractor.py input_graph.fastg<STRING> reference.fasta<STRING> CDS.fasta<STRING> RATIO<FLOAT> tracks.fasta<STRING> tracks_len<INT>

input_graph.fastg - .fastg-file, created, for example, with SPAdes assembler.

reference.fasta - linear genome of clothest species

CDS.fasta - CDS-sequences of clothest species with ncbi-specified header

For example:
  >lcl|NT_033779.5_cds_NP_787955.2_1 [gene=CG11023] [locus_tag=Dmel_CG11023] [db_xref=FLYBASE:FBpp0289914,GeneID:33155] [protein=uncharacterized protein, isoform C] [protein_id=NP_787955.2] [location=join(7680..8116,8193..8589,8668..9276)] [gbkey=CDS]

RATIO - parameter for splitting training set into the training dataset and test dataset

tracks.fasta - output fasta file

tracks_len - length of the output sequences
