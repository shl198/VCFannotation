### This project aims to create a tool to annotate short variants in vcf format

#### Prepare database files
First, we need to prepare some files that serve as databases for this tool to retrieve information from. The 
1. Download raw reference and annotation data
position information of each mRNA: The steps to get this file is:
    (a) Get whole genome sequence of human_g1k_v37 from ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz.
    (b) Get the annotation file for this genome from ftp://ftp.ensembl.org/pub/grch37/release-95/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz
2. Generate files that are easier to use
    (a) Here we use tool gff3ToGenePred from https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=2&ved=2ahUKEwiRqPjng_TfAhW0HjQIHXLTCakQFjABegQICBAB&url=http%3A%2F%2Fhgdownload.cse.ucsc.edu%2Fadmin%2Fexe%2Flinux.x86_64%2Fgff3ToGenePred&usg=AOvVaw0O7-M4IOZ1AR9u7OeYANfQ to transfer gff3 to genepred format. This allows more easier way to parse the gene positions. The way to run this is **gff3ToGenePred input.gff3 output.txt **
    (b) Use gffread to extract full mRNA sequence fasta file using gffread from http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.10.6.Linux_x86_64.tar.gz. The command is: **gffread input.gff3 -g human_g1k_v37.fasta -w human_g1k_v37.rna.fa**
    
#### Strategy to annotate the vcf file
1. an interval tree was built using start and end position of all the mRNAs.
2. for each variation, it checks if the mutation locates in the mRNA or not
3. If mutation locates in the mRNA, then decide if it locates in UTR5, CDS, intron or UTR3
4. If it locates in CDS, then further annotation if it's non synonymous or synonymouse SNV, frame shift insertion/deletion or not.
5. Parse each line to extract the coverage support.
6. Use Exac database to retrieve allele frequency information. from here

#### Annotate the vcf file
The code works in python 3. First you need to install some python packages:
* intervaltree: Used to create a tree object to store position information of mRNA sequences
* Biopython: Used to translate dna sequence to protein sequence
* requests: Used to trace the information of Exac

#### How to run this pipeline
1. prepare the following files:
	(a) input raw vcf file <br />
	(b) pref file prepared from gff3 file <br />
	(c) full mRNA sequence fasta file <br />

2. In terminal run the command:
	**python annotate.py -i input.vcf -p human_g1k_v37.pred -f human_g1k_v37.rna.fa -o anno_vcf_anno.vcf**

3. add allele frequency information:
	**python add_allele_frequency.py -i anno_vcf_anno.vcf -o anno_vcf_anno_freq.vcf**

The final annotation file anno_vcf_anno_freq.vcf includes original vcf columns and the following addition columns:
* variat_type: variation type, eg: synonymous mutation
* variant_detail: details of the nucleotide change and peptide change
* normal_cov: read coverage at the mutation position for normal sample
* vaf5_cov: read coverage at the mutation position for patient sample
* normal_mut_cov: reads support mutation in normal sample
* vaf5_mut_cov: reads support mutation in patient sample    
* normal_mut_cov_percent: percentage of reads support mutation in normal sample
* vaf5_mut_cov_percent: percentage of reads support mutation in patient sample
* allele_frequency: allele frequency from ExAC database