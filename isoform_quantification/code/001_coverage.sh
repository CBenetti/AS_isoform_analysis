cd ../STAR/fastq/
dirl=$('ls')
cd ../../CB
#dirl=("C10.STAR/" "C10S1.STAR/" "COGA5.STAR/" "COGA5L.STAR/" "HDC114.STAR/" "HDC114S17.STAR/" "HROC24.STAR/" "HROC24S4.STAR/" "HROC277.STAR/" "HROC277MET2.STAR/" "LIM1215.STAR/" "LIM1215S22.STAR/" "LS180.STAR/" "LS180S23.STAR/"  "LS411N.STAR/" "LS411NS23.STAR/" "SNU1040.STAR/" "SNU1040S11.STAR/" "SNU1235.STAR/" "SNU1235S21.STAR/")
#length=${#dirl[@]}
for d in ${dirl[@]}
do
	counts=$(grep 'Uniquely mapped reads number' ../STAR/fastq/${d}/Log.final.out | cut -d'|' -f2- );
	echo -e  $d"\t"$counts >> out/coverage.txt;
done
