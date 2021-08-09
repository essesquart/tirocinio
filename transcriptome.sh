#!/bin/sh

#SRA_entries va cambiato in base a quali organismi stiamo analizzando (e anche in base al loro sesso)
#transcriptome e UTR_db vanno anche cambiati in base all'organismo.

#questo script deve girare all'interno della cartella Transcriptome, con all'interno il file 3UTRef.qualcosa.fasta

SRA_entries="SRR5811966"
layout="paired"
transcriptome="Noge_kraken2"
UTR_db="3UTRef.Inv.fasta"
num_threads=12
adapter_file="illumina.adapters.fa"
kraken2_db="kraken_20190620"
minlen=75

fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files $SRA_entries

cp /home/diegocarli/$adapter_file .

if [ layout == "paired" ]

	then for fastq in *_1.fastq
		do java -jar /opt/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads $num_threads $fastq ${fastq%_1.fastq}_2.fastq -baseout ${fastq%_1.fastq}_trimmed ILLUMINACLIP:$adapter_file:2:30:10 AVGQUAL:20 LEADING:3 TRAILING:3 SLIDINGWINDOW:25:33 MINLEN:$minlen
		done

	for trimmed_fastq in *_trimmed_*
		do awk '{if ($1~/@SRR/) print "@"$2"/1"; else print $0}' $trimmed_fastq > ${trimmed_fastq%_*P}_fixed_${trimmed_fastq:(-3)}
		done

	cat *_1P > trimmed_fixed_1P
	cat *_2P > trimmed_fixed_2P

	kraken2 --db /media/storage/dbs/$kraken2_db --threads $num_threads --output kraken2.out --paired --use-names --report kraken2.report --classified-out classified_out# classified_out# --unclassified-out unclassified_out# unclassified_out# trimmed_fixed_*P

	fastqc unclassified_out_*

	Trinity --seqType fq --max_memory 10G --full_cleanup --no_normalize_reads --CPU $num_threads --left unclassified_out_1 --right unclassified_out_2 --output "$transcriptome"_trinity

	else for fastq in *.fastq
		do java -jar /opt/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads $num_threads $fastq ${fastq%.fastq}.trim.fastq ILLUMINACLIP:$adapter_file:2:30:10 AVGQUAL:20 LEADING:3 TRAILING:3 SLIDINGWINDOW:25:33 MINLEN:$minlen
		done

	for fastq in *.trim.fastq
        	do awk '{if ($1~/@SRR/) print "@"$2"/1"; else print $0}' $fastq > ${fastq%.fastq}.fixed.fastq
        	done

	for fastq in *.trim.fixed.fastq
		do kraken2 --db /media/storage/dbs/$kraken2_db --threads $num_threads --output kraken2.out --use-names --report kraken2.report --classified-out ${fastq%.fixed.fastq}.cont.fastq --unclassified-out ${fastq%.fixed.fastq}.kraken2.fastq $fastq
	        done

	fastqc *.trim.kraken2.fastq

	if [ -f file.list ]
		then rm file.list
		fi
	for fastq in *.trim.kraken2.fastq
		do echo ,$fastq >> file.list
		done
	sed -i ':a;N;$!ba;s/\n//g' file.list
	sed -i 's/,//1' file.list
	fastq_file=`cat file.list`

	Trinity --seqType fq --max_memory 10G --full_cleanup --no_normalize_reads --CPU $num_threads --single $fastq_file --output "$transcriptome"_trinity

	rm file.list

	fi

rm $adapter_file

if [ -d 3UTR ]
	then rm -rf 3UTR
	fi
mkdir 3UTR
cd 3UTR
#Ricorda che l'opzione di default è il codice genetico standard. Per cambiarlo, devi modificare il file perl '3UTR_orf_20170816.pl' inserendo, per ogni funzione translate, '-codontable_id => 2', per esempio (codice mitocondriale).
perl /opt/ExUTR-v0.1.0/bin/3UTR_orf_20190703.pl -i ../"$transcriptome"_trinity.Trinity.fasta -d /media/storage/dbs/swissprot/swissprot -a $num_threads -o "$transcriptome"_3UTR -l un
perl /opt/ExUTR-v0.1.0/bin/3UTR_ext_20170816.pl -i1 ../"$transcriptome"_trinity.Trinity.fasta -i2 "$transcriptome"_3UTR_orfs.fa -a $num_threads -o "$transcriptome"_3UTR.fasta -x 2500 -m 1 -d ../$UTR_db
