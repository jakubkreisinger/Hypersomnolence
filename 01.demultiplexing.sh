# files with primer info
ADAPT_PATH=$(echo /media/kreising/DATA/data/Radka_Janet/External_data/primers.fasta)
MATRIX_PATH=$(echo /media/kreising/DATA/data/Radka_Janet/External_data/matrix.txt)

cd /media/kreising/DATA/data/Radka_Janet/00.RAW_DATA
mkdir ../demultiplexed
rm ../demultiplexed/*


#skewer - demultiplexing 
SAMPLES=$(ls -1| grep "R[12].fastq.gz"|sed -r 's/R[12].fastq.gz//'|sort| uniq)
for i in $SAMPLES; do skewer -x $ADAPT_PATH -M $MATRIX_PATH -b -m head -k 35 -d 0 -t 8 "$i"R1.fastq.gz "$i"R2.fastq.gz -o ../demultiplexed/"$i"; done >> ./log


cd ../demultiplexed/
gzip *fastq

cd ..
mkdir 02A.DEMULTI.16s

cd ./demultiplexed

#skewer - primer trimming
SAMPLES_16S=$(ls -1| grep "assigned-F_"|sed -r 's/-pair[12].fastq.gz//'|sort| uniq)
for i in $SAMPLES_16S; do skewer -x NNNNCCTACGGGNGGCWGCAG -y GACTACHVGGGTATCTAATCC  -m head -k 35 -d 0 -t 8 "$i"-pair1.fastq.gz "$i"-pair2.fastq.gz -o ../02A.DEMULTI.16s/"$i"; done >> ../02A.DEMULTI.16s/log.trimmed.head


cd ../02A.DEMULTI.16s
gzip *fastq

