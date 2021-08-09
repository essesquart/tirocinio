#!/bin/sh

genome_mit=`cat genome_mit.db`
genome_mit_file=`cat genome_mit.fa`
#genome_nuc=`cat genome_nuc.db`

for infile in *.trim.kraken2.fastq

	do outfile=$infile\_Remove.sam
	file=$infile\_Remove.fastq
	outfile1=$infile\_Remove.bam
        prefile1=$infile\_Pre1_Remove.sam
	prefile2=$infile\_Pre2_Remove.bam
	postfile1=$infile\_Post1_Remove.bam
	postfile2=$infile\_Post2_Remove.fastq
	postfile3=$infile\_Post3_Remove.sam
	postfile3=$infile\_Post.sam

	mito=$infile\_Multi_MitoUnique.bam
	outfile2=$infile\_Multi_MitoNuclear.bam

	bowtie2 -x $genome_mit/$genome_mit -q $infile -S $outfile -N 1 -i C,1 -L 18

	samtools view -b -S $outfile > $outfile1

	samtools flagstat $outfile1 >> Log_Mapped.txt

	samtools view -b -F 4 $outfile1 > $mito

	bedtools bamtofastq -i $outfile2 -fq $file

	bowtie2 -x $genome_nuc/$genome_nuc -q $file -S $prefile1 -i C,1 -L 22

	samtools view -b -S $prefile1 > $prefile2

	samtools flagstat $prefile2 >> Log_Mapped2.txt

	samtools view -b -f 4 $prefile2 > $postfile1

	bedtools bamtofastq -i $postfile1 -fq $postfile2

	bowtie2 -x $genome_mit/$genome_mit -q $postfile2 -S $postfile3 -N 1 -i C,1 -L 18

	samtools view -b -F 4 $postfile3 > $mito

	rm *Remove.sam
	rm *Remove.bam
	rm *Remove.fastq

	find *Unique.bam > Name|grep mapped\ \( Log_Mapped.txt > Percentage		
	find *Nuclear.bam > Name1|grep mapped\ \( Log_Mapped1.txt > Percentage1
	find *Unique.bam > Name2|grep mapped\ \( Log_Mapped2.txt > Percentage2

	paste Name Percentage > ResultsMito.txt
        paste Name1 Percentage1 > ResultsMitoNuclear.txt
	paste Name2 Percentage2 > ResultsMito.txt

done

if [ -d BedFiles ]
	then rm -r BedFiles
	fi
mkdir BedFiles
cp $genome_mit/$genome_mit_file BedFiles
cp *Unique.bam BedFiles

cd BedFiles
for infile in *.bam

	do outfile=$infile\_Map.bed

	outfile1=$infile\_Sort.bed

	outfile2=$infile\_Cluster.csv

	bedtools bamtobed -i $infile > $outfile

	bedtools sort -i $outfile > $outfile1

	bedtools cluster -i $outfile1 > $outfile2

	rm *Sort.bed

done
cd ..
