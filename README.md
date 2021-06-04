# ENCODE_FADS_MPRA

Download data from ENCODE DCC


##Prereqs
Conda
git

#########
# Create Env
conda create --name FADS_ENCODE python=3.6
conda install --name FADS_ENCODE -c bioconda flash2=2.2.00 minimap2=2.17 preseq perl-text-levenshteinxs bioconductor-deseq2=1.28.0 bioconductor-rtracklayer bioconductor-genomicranges bioconductor-bsgenome.hsapiens.ucsc.hg19 r-ggplot2 r-gridextra r-tidyr r-splitstackshape

conda activate FADS_ENCODE

#########
# Download files
xargs -n1 -a OL13_encode_download.txt -I URL curl -O -L URL
git clone https://github.com/tewhey-lab/ENCODE_FADS_MPRA.git

# Setup
conda activate FADS_ENCODE

INSTALL_PATH="../"
THREADS=20
MAX_MEMORY=30


# Run
mkdir oligo_tag
cd oligo_tag

# Step to create barcode/Oligo connections
sbatch -q batch -N 1 -n 22 --mem 33G -t 24:00:00 -o slurm-%A_%a.out --wrap "../MPRA_Oligo-Tag_pipeline/run.VectorReconstruction_MPRA.sh ../files/ENCFF474GEU.fasta.gz OL13_FADS $THREADS ../files/ENCFF148NVC.fastq.gz ../files/ENCFF103XEY.fastq.gz"

##Tag Seq

mkdir tag_seq
cd tag_seq

touch tmp.out

gzip -dc ../files/ENCFF425XBZ.fastq.gz ../files/ENCFF212MOM.fastq.gz ../files/ENCFF163MTX.fastq.gz ../files/ENCFF315UZK.fastq.gz | perl ${INSTALL_PATH}/MPRA_Tag_Analysis/matchadapter_TagRead.pl -A -H 20 TCTAGAGGTTCGTCG OL13_FADS_K562_rep1
gzip -dc ../files/ENCFF671NYX.fastq.gz ../files/ENCFF734AEJ.fastq.gz ../files/ENCFF952FCD.fastq.gz ../files/ENCFF721KDZ.fastq.gz | perl ${INSTALL_PATH}/MPRA_Tag_Analysis/matchadapter_TagRead.pl -A -H 20 TCTAGAGGTTCGTCG OL13_FADS_K562_rep2
gzip -dc ../files/ENCFF552JYE.fastq.gz | perl ${INSTALL_PATH}/MPRA_Tag_Analysis/matchadapter_TagRead.pl -A -H 20 TCTAGAGGTTCGTCG OL13_FADS_K562_rep3
gzip -dc ../files/ENCFF994NKI.fastq.gz | perl ${INSTALL_PATH}/MPRA_Tag_Analysis/matchadapter_TagRead.pl -A -H 20 TCTAGAGGTTCGTCG OL13_FADS_K562_rep4

gzip -dc ../files/ENCFF379VAH.fastq.gz ../files/ENCFF321WKZ.fastq.gz ../files/ENCFF967KGU.fastq.gz ../files/ENCFF148OPN.fastq.gz | perl ${INSTALL_PATH}/MPRA_Tag_Analysis/matchadapter_TagRead.pl -A -H 20 TCTAGAGGTTCGTCG OL13_FADS_plasmid_rep1
gzip -dc ../files/ENCFF626IQN.fastq.gz ../files/ENCFF727THA.fastq.gz ../files/ENCFF670SCJ.fastq.gz ../files/ENCFF609IJK.fastq.gz | perl ${INSTALL_PATH}/MPRA_Tag_Analysis/matchadapter_TagRead.pl -A -H 20 TCTAGAGGTTCGTCG OL13_FADS_plasmid_rep2
gzip -dc ../files/ENCFF438QVC.fastq.gz ../files/ENCFF227BNM.fastq.gz ../files/ENCFF752NUR.fastq.gz ../files/ENCFF117ETI.fastq.gz | perl ${INSTALL_PATH}/MPRA_Tag_Analysis/matchadapter_TagRead.pl -A -H 20 TCTAGAGGTTCGTCG OL13_FADS_plasmid_rep3
gzip -dc ../files/ENCFF067AUY.fastq.gz ../files/ENCFF598MFS.fastq.gz ../files/ENCFF577HFU.fastq.gz ../files/ENCFF222NLL.fastq.gz | perl ${INSTALL_PATH}/MPRA_Tag_Analysis/matchadapter_TagRead.pl -A -H 20 TCTAGAGGTTCGTCG OL13_FADS_plasmid_rep4

cat OL13_FADS_K562_rep1.match | perl ${INSTALL_PATH}/MPRA_Tag_Analysis/associate_tags.pl stdin ${INSTALL_PATH}/oligo_tag/OL13_FADS.merged.rc.match.enh.mapped.barcode.ct.parsed tmp.out > OL13_FADS_K562_rep1.tag
cat OL13_FADS_K562_rep2.match | perl ${INSTALL_PATH}/MPRA_Tag_Analysis/associate_tags.pl stdin ${INSTALL_PATH}/oligo_tag/OL13_FADS.merged.rc.match.enh.mapped.barcode.ct.parsed tmp.out > OL13_FADS_K562_rep2.tag
cat OL13_FADS_K562_rep3.match | perl ${INSTALL_PATH}/MPRA_Tag_Analysis/associate_tags.pl stdin ${INSTALL_PATH}/oligo_tag/OL13_FADS.merged.rc.match.enh.mapped.barcode.ct.parsed tmp.out > OL13_FADS_K562_rep3.tag
cat OL13_FADS_K562_rep4.match | perl ${INSTALL_PATH}/MPRA_Tag_Analysis/associate_tags.pl stdin ${INSTALL_PATH}/oligo_tag/OL13_FADS.merged.rc.match.enh.mapped.barcode.ct.parsed tmp.out > OL13_FADS_K562_rep4.tag

cat OL13_FADS_plasmid_rep1.match | perl ${INSTALL_PATH}/MPRA_Tag_Analysis/associate_tags.pl stdin ${INSTALL_PATH}/oligo_tag/OL13_FADS.merged.rc.match.enh.mapped.barcode.ct.parsed tmp.out > OL13_FADS_plasmid_rep1.tag
cat OL13_FADS_plasmid_rep2.match | perl ${INSTALL_PATH}/MPRA_Tag_Analysis/associate_tags.pl stdin ${INSTALL_PATH}/oligo_tag/OL13_FADS.merged.rc.match.enh.mapped.barcode.ct.parsed tmp.out > OL13_FADS_plasmid_rep2.tag
cat OL13_FADS_plasmid_rep3.match | perl ${INSTALL_PATH}/MPRA_Tag_Analysis/associate_tags.pl stdin ${INSTALL_PATH}/oligo_tag/OL13_FADS.merged.rc.match.enh.mapped.barcode.ct.parsed tmp.out > OL13_FADS_plasmid_rep3.tag
cat OL13_FADS_plasmid_rep4.match | perl ${INSTALL_PATH}/MPRA_Tag_Analysis/associate_tags.pl stdin ${INSTALL_PATH}/oligo_tag/OL13_FADS.merged.rc.match.enh.mapped.barcode.ct.parsed tmp.out > OL13_FADS_plasmid_rep4.tag

# Create count matrix
perl ${INSTALL_PATH}/MPRA_Tag_Analysis/compile_bc.pl -ECMS -A 0.05 ${INSTALL_PATH}/MPRA_Tag_Analysis/sample_list.txt OL13_FADS_K562_Counts.out >  OL13_FADS_K562_Counts.log

# Create final processed data 
cd ../count_analysis
Rscript --vanilla ../MPRA_Tag_Analysis/FADS_MPRA_Analysis.R ../tag_seq/OL13_FADS_K562_Counts.out OL13_FADS

cut -f1-11 OL13_FADS_K562_20210512.bed |  awk '(NR>1 && $7 !~ "NA"){if($10~"NA"){$10=0;$11=0};$5=0;print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' | sort -k1,1 -k2,2n > OL13_FADS_Tile_K562.hg19.enc.bed


