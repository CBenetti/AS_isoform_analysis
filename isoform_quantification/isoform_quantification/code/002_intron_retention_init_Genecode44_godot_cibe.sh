
###Initialization (only once, already run)
#docker pull cloxd/irfinder:2.0

###Genome build (only once, already run)
#docker run -u 1018:1005 -w /data/cell_bank/IRFinder -v /mnt/oncog/analysis/cell_bank/dataset/v2/CB/out/IRFinder/:/data/cell_bank/IRFinder cloxd/irfinder:2.0 \BuildRefDownload -r Refdir_Genecode44_Ens110 ftp://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz  
#docker run -u 1018:1005 -w /data/cell_bank/IRFinder -v /mnt/oncog/analysis/cell_bank/dataset/v2/CB/out/IRFinder/:/data/cell_bank/IRFinder cloxd/irfinder:2.0 BuildRefProcess -r ./Refdir_Genecode44_Ens110


####Go trough every directory
cd ../STAR/fastq
dirl=$('ls')
cd ../../CB
for d in ${dirl[@]}
do

#dirl=("COGA5.STAR/" "COGA5L.STAR/" "HDC114.STAR/" "HDC114S17.STAR/" "HROC24.STAR/" "HROC24S4.STAR/" "HROC277.STAR/" "HROC277MET2.STAR/" "LIM1215.STAR/" "LIM1215S22.STAR/" "LS180.STAR/" "LS180S23.STAR/"  "LS411N.STAR/" "LS411NS23.STAR/" "SNU1040.STAR/" "SNU1040S11.STAR/" "SNU1235.STAR/" "SNU1235S21.STAR/")
#length_cycle=(${#dirl[@]})
#for ((i=1; i<=$length_cycle; i++))
#for ((i=1; i<=2; i++))
#do  
	echo "sample: $d";
	###BAM header adaptation;
	echo "BAM header adaptation"
	samtools view -H ../STAR/fastq/$d/Aligned.sortedByCoord.out.ribo.ex.bam > tmp/header.sam;
	#linux: 
	sed -e 's/chr//g' tmp/header.sam > tmp/headerens.sam;
	#sed 's/chr//g' STAR/fastq/${dirl[$i]}header.sam > STAR/fastq/${dirl[$i]}headerens.sam;
	samtools reheader tmp/headerens.sam ../STAR/fastq/$d/Aligned.sortedByCoord.out.ribo.ex.bam > tmp/Aligned.goodheader.out.bam;
	###IRFinder
	echo "IRFinder";
	docker run -u 1018:1005 -w /data/cell_bank -v /mnt/oncog/analysis/cell_bank/dataset/v2/CB/tmp/:/data/cell_bank/data/ -v /mnt/oncog/analysis/cell_bank/dataset/v2/CB/out/IRFinder/out:/data/cell_bank/out -v /mnt/oncog/analysis/cell_bank/dataset/v2/CB/out/IRFinder/Refdir_Genecode44_Ens110/:/data/cell_bank/genome/ cloxd/irfinder:2.0 BAM -r genome -d out/$d data/Aligned.goodheader.out.bam;
	###IRselection
	#echo "IR selection";
	#cd IRFinder/Results/$dirl[$i];
	#Rscript '../../../code/002_irselection.R';
	#cd ../../../;
done
