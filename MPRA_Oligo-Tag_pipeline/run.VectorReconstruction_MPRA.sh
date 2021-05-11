INSTALL_PATH=" /projects/tewher/bin/MPRA_Oligo-Tag_pipeline/scripts"

#ID="<Output Prefix>"
#READ_A="<Read 1 fastq | fastq.gz>"
#READ_B="<Read 2 fastq | fastq.gz>"
#REF="<Reference fasta>"
#THREADS=<number of cpus>

REF=$1
ID=$2
THREADS=$3
READ_A=$4
READ_B=$5

flash2 -r 150 -f 274 -s 20 -o ${ID}.merged -t $THREADS $READ_A $READ_B > ${ID}.merged.log
perl ${INSTALL_PATH}/fq2RCfa.pl ${ID}.merged.extendedFrags.fastq > ${ID}.merged.rc.fasta
perl ${INSTALL_PATH}/matchadapter.pl ${ID}.merged.rc.fasta ${ID}.merged.rc
awk '{print ">"$1"#"$5"\n"$4}' ${ID}.merged.rc.match > ${ID}.merged.rc.match.enh.fa
minimap2 --for-only -Y --secondary=no -t $THREADS -m 10 -n 1 --end-bonus 12 -O 5 -E 1 -k 10 -2K50m --MD --eqx --cs=long -c -a $REF ${ID}.merged.rc.match.enh.fa > ${ID}.merged.rc.match.enh.sam 2> ${ID}.merged.rc.match.enh.log
perl ${INSTALL_PATH}/SAM2MPRA.pl -C -B ${ID}.merged.rc.match.enh.sam ${ID}.merged.rc.match.enh.mapped  > ${ID}.merged.rc.match.enh.mapped.log 2>&1

grep PASS ${ID}.merged.rc.match.enh.mapped  | sort -S16G -k4 > ${ID}.merged.rc.match.enh.mapped.enh.pass.sort
sort -S16G -k2 ${ID}.merged.rc.match.enh.mapped > ${ID}.merged.rc.match.enh.mapped.barcode.sort
perl ${INSTALL_PATH}/Ct_seq.pl ${ID}.merged.rc.match.enh.mapped.barcode.sort 2 4 > ${ID}.merged.rc.match.enh.mapped.barcode.ct
perl ${INSTALL_PATH}/Ct_seq.pl ${ID}.merged.rc.match.enh.mapped.enh.pass.sort 4 2 > ${ID}.merged.rc.match.enh.mapped.enh.pass.ct
awk '{ct[$4]++}END{for (i in ct)print i "\t" ct[i]}' ${ID}.merged.rc.match.enh.mapped.barcode.ct | sort -k1n > ${ID}.merged.rc.match.enh.mapped.barcode.ct.hist
preseq lc_extrap -H ${ID}.merged.rc.match.enh.mapped.barcode.ct.hist -o ${ID}.merged.rc.match.enh.mapped.barcode.ct.hist.preseq -s 25000000 -n 1000 -e 1000000000
awk '($5 == 0)' ${ID}.merged.rc.match.enh.mapped.barcode.ct | awk '{ct[$2]++;cov[$2]+=$4;}END{for(i in ct)print i "\t" ct[i] "\t" cov[i]}' > ${ID}.merged.rc.match.enh.mapped.barcode.ct.plothist
perl ${INSTALL_PATH}/parse_map.pl ${ID}.merged.rc.match.enh.mapped.barcode.ct > ${ID}.merged.rc.match.enh.mapped.barcode.ct.parsed
Rscript --vanilla ${INSTALL_PATH}/plot_barcode_stats.r $ID $REF
