# MPRA Oligo-Tag Pipeline

Dependencies
  * BWA MEM (http://bio-bwa.sourceforge.net/)
  * FLASH (https://sourceforge.net/projects/flashpage/)
  * Preseq (http://smithlabresearch.org/software/preseq/)
  
Required Setup
  * BWA, FLASH & Preseq executables need to be in your PATH 
  * Update run.VectorReconstruction_MPRA.sh so INSTALL_PATH points to the scripts folder
  * Fasta file of oligo sequences needs to be indexed using BWA
  
Syntax

run.VectorReconstruction_MPRA.sh [Oligo Sequence Fasta] [ID for Output files] [CPU Threads] [Illumina FASTQ Read 1] [Illumina FASTQ Read 2]


**WARNING**
This pipeline can change without notice. I will try my best to make it backwards compatable.
