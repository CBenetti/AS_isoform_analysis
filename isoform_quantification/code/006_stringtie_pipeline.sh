cd ../STAR/fastq
dirl=$('ls')
cd ../../CB
for d in ${dirl[@]}
do
echo $d
/mnt/oncog/analysis/cell_bank/dataset/v2/stringtie/stringtie -x chrM -G data/gencode.v44.annotation.gtf -A out/$d/gene_abund.tab -o out/$d/output.gtf ../STAR/fastq/$d/Aligned.sortedByCoord.out.bam
done


