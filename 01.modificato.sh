#!/bin/sh

output="NoGe_results"

genome_mit="CrGi_mit_doubled"
genome_mit_file="CrassoGigas.mit.fasta"
#genome_nuc="CrGi_nuc
#genome_nuc_file="CrassoGigas.nuc.fasta"

num_threads=12

#prossime righe per riprendere tutto quello che serve dalle varie cartelle, da togliere nello script finale, serve solamente nelle varie prove.
cp /home/diegocarli/scripts/* .

SRA_entries="SRR13528976 SRR13528975 SRR13528974"
adapter_file="illumina.adapters.fa"
kraken2_db="kraken_20210616"

#Dice a R se controllare o no gli allineamenti Bowtie2.
alignment_check="FALSE"

#Stabilisce una dimensione per i cluster di sequenze (individuati da USEARCH) da tenere in quanto molto rappresentati (e quindi molto trascritti!).

cluster_threshold=200

#Stabilisce il nome dei file di configurazione per R. Cambiando questa opzione, va cambiata anche corrispondentemente negli script R.

conf_temp=conf.temp
genomecov_temp=genomecov.temp

#Scrive il nome di base degli output in un file che verrà cancellato poi alla fine.

echo $output > output.basename
echo $genome_mit > genome_mit.db
echo $genome_mit_file > genome_mit.fa
echo $genome_nuc > genome_nuc.db

if [ -d R ]
	then rm -r R
	fi

cp ../Scripts/duplicate_genome.R .

mkdir R
cd R

 Rscript ../duplicate_genome.R
cd ..
rm -r R

#genome_mit è l'indice sulla quale bowtie2 va a scrivere?
bowtie2-build -f $genome_mit_file $genome_mit
bowtie2-build -f $genome_nuc_file $genome_nuc

if [ -d $genome_mit ]
	then rm -rf $genome_mit
	fi
mkdir $genome_mit

if [ -d $genome_nuc ]
	then rm -rf $genome_nuc
	fi
mkdir $genome_nuc

mv $genome_mit.* $genome_mit
if [ -f $genome_mit_file ]
	then mv $genome_mit_file $genome_mit
	fi
mv $genome_nuc.* $genome_nuc
if [ -f $genome_nuc_file ]
	then mv $genome_nuc_file $genome_nuc
	fi

fastq-dump --defline-seq '@$sn[_cds	$rn]/$ri' --split-files $SRA_entries

#commento perchè gli adapter file ce li ho già
cp /home/federicoplazzi/$adapter_file .
#cp /home/diegocarli/$adapter_file .

for fastq in *.fastq
	do java -jar /opt/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads $num_threads $fastq ${fastq%.fastq}.trim.fastq ILLUMINACLIP:$adapter_file:2:30:10 AVGQUAL:20 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20
	done
rm $adapter_file

for fastq in *.trim.fastq
	do kraken2 --db /media/storage/dbs/$kraken2_db --threads $num_threads --use-names --report kraken2.report --classified-out ${fastq%.trim.fastq}.cont.fastq --unclassified-out ${fastq%.trim.fastq}.trim.kraken2.fastq --output kraken2.out $fastq
	done

fastqc *.trim.kraken2.fastq

num_fastq=`ls *.trim.kraken2.fastq | grep -c trim\.kraken2\.fastq`
num_curr_fastq=0

while [ $num_curr_fastq -lt $num_fastq ]
	do num_curr_fastq=$(( num_curr_fastq+1 ))
	curr_fastq=`ls *.trim.kraken2.fastq | sed -n ''$num_curr_fastq' p'`
	sed -i '1~4 s/^@/@'$num_curr_fastq'_/g' $curr_fastq
	done

sh 01.core.modificato.sh

mv genome_mit.fa BedFiles
cp remove_duplicate_region.R BedFiles

cd BedFiles

#Prepara la cartella di R

mkdir R

#Concatena i file .bed e ne costruisce uno unico, che mette in ordine.

cat *.bed > total.bed
bedtools sort -i total.bed > total.sorted.bed
rm total.bed

#Esegue uno script R che va a cercare l'ultima riga in cui ci sono mappature di read che cominciano prima della fine del genoma mitocondriale non duplicato e toglie tutte le righe successive.

cd R
Rscript ../remove_duplicate_region.R
cd ..

rm total.sorted.bed

cd ..

genome_mit="mit_emp"
genome_mit_file=`cat genome_out.out`

bowtie2-build -f $genome_mit_file $genome_mit

if [ -d $genome_mit ]
	then rm -rf $genome_mit
	fi
mkdir $genome_mit
mv $genome_mit.* $genome_mit
mv $genome_mit_file $genome_mit

echo $genome_mit > genome_mit.db
echo $genome_mit_file > genome_mit.fa

sh 01.core.modificato.sh

mv genome_mit.fa BedFiles

rm genome_mit.db
rm genome_nuc.db
rm genome_out.out

cp 02.sh BedFiles
cp 03.sh BedFiles
cp *.R BedFiles

cd BedFiles

mkdir R

#Prepara un po' di nomi.

output_fa=$output.fa
output_conf=$conf_temp
genomecov_conf=$genomecov_temp
output_centroids=$output.centroids.fa
output_selected_centroids=$output.centroids.selected.fa
output_clusters=$output.clusters
clusters_bed=$output.clusters.bedfiles
clusters_fasta=$output.clusters.fasta

#Rimuove, se c'è, il file FASTA con le read allineate sul mitocondrio, ma non sul nucleare.

if [ -f $output_fa ]
	then rm $output_fa
	fi

#Scrive il nome di questo FASTA da fare nel file di configurazione da fare per R.

echo $output_fa > $output_conf

#Rimuove, se c'è, il file FASTA con i centroidi selezionati (quelli dei cluster con almeno 200 sequenze).

if [ -f $output_selected_centroids ]
	then rm $output_selected_centroids
	fi

#Scrive il nome di questo FASTA da fare nel file di configurazione da fare per R.

echo $output_selected_centroids >> $output_conf

#Scrive il nome della cartella in cui finiranno i FASTA dei cluster selezionati e in cui finiranno concatenati in un unico FASTA nel file di configurazione da fare per R..

echo $clusters_fasta >> $output_conf

#Scrive anche il nome del file con il genoma mitocondriale nel file di configurazione per R.

echo $genome_mit_file >> $output_conf

#Cancella, se ci sono, le cartelle di output di R e le rifà.

if [ -d output.R ]
	then rm -rf output.R
	fi
mkdir output.R

#Scrive i nomi di tutti i CSV nel file di configurazione di R e poi scrive (con bamtools) un unico ('>>') FASTA con le sequenze delle read che hanno allineato sul mitocondrio, ma non sul nucleare.

for bam in *.bam
	do echo $bam\_Cluster.csv >> $output_conf
	samtools fasta -n $bam >> $output_fa
	done

#Cancella, se ci sono, le cartelle in cui mettere i cluster e i relativi file BED e le rifà.

if [ -d $output_clusters ]
	then rm -rf $output_clusters
	fi
mkdir $output_clusters
if [ -d $clusters_bed ]
	then rm -rf $clusters_bed
	fi
mkdir $clusters_bed
if [ -d $clusters_fasta ]
	then rm -rf $clusters_fasta
	fi
mkdir $clusters_fasta

#Trova i cluster con CD-HIT. -M 0 indica di usare tutta la memoria disponibile; -c è il livello di similarità GLOBALE (vengono tutti 100%); -r 0 indica il +/+.

sed -i 's/>.\+:\([0-9]\+\):\([0-9]\+\):\([0-9]\+\)$/>\1_\2_\3/g' $output_fa

cd-hit-est -i $output_fa -o $output_centroids -T $num_threads -M 0 -c 0.99 -r 0

#Costruisce un file BED nella cartella apposta solo per i cluster più grandi di '$cluster_threshold' sequenze e mette i relativi FASTA nella cartella apposta.

cp ../clusters1.py .
python clusters1.py $output_fa $output_centroids.clstr $cluster_threshold $clusters_fasta

#Costruisce un file BED nella cartella apposta solo per i cluster più grandi di '$cluster_threshold' sequenze.

cd $clusters_fasta
for cluster in `ls`
	do num_sequences=`grep -c ">" $cluster`
		for s in `grep ">" $cluster | sed 's/>//g' | sed 's/_/:/g'`
			do grep -w $s ../*bed >> ../$clusters_bed/${cluster%"fa"}"bed"
		done
	done
cd ..


#Concatena i FASTA dei cluster selezionati, in modo che R poi possa andare a vedere la distribuzione delle lunghezze delle read.

cd $clusters_fasta
cat *.fa > complete.fa
cd ..

#Fa girare uno script R che produce il file di genoma che serve a bedtools per calcolare la coverage.

cd R

Rscript ../genome.R

cd ..

#Calcola la coverage dei singoli cluster (quelli selezionati). Per sapere il nome del file di genoma, legge la prima parola di quello scritto da R.

genome_file=`awk '{print $1}' *.genome`.genome

cd $clusters_bed
sed -i 's/.\+_Multi_MitoUnique\.bam_Map\.bed:\(.\+\)/\1/g' *.bed
for bed in *.bed
	do if [ $bed != "total.bed" ]
		then bedtools sort -i $bed > ${bed%.bed}.sorted.bed
		bedtools genomecov -d -i ${bed%.bed}.sorted.bed -g ../$genome_file > ${bed%.bed}.genomecov
		bedtools genomecov -5 -d -i ${bed%.bed}.sorted.bed -g ../$genome_file > ${bed%.bed}.genomecov.5
		bedtools genomecov -3 -d -i ${bed%.bed}.sorted.bed -g ../$genome_file > ${bed%.bed}.genomecov.3
		fi
	done
cd ..

#Scrive un FASTA con i centroidi dei cluster selezionati.

for s in `grep ">" $output_centroids | sed 's/>//g' | sed 's/_/:/g'`
	do grepfile=`echo $s | sed 's/:/_/g' | sed 's/\//_/g'`
	grep -w $s $clusters_bed/*.sorted.bed > $grepfile
	if [ -s $grepfile ]
		then grep -w -A 1 ">"$grepfile $output_centroids >> $output_selected_centroids
		fi
	rm $grepfile
	done

for s in `grep ">" $output_centroids`
	do control=`grep $s $clusters_fasta/*fa`
	if [ -n "$control" ]
		then grep -w -A 1 $s $output_centroids >> $output_selected_centroids
	fi
	done

#Prepara un secondo file di configurazione per R, questa volta per lo script che fa i grafici delle coverage dei cluster selezionati.

echo $clusters_bed > $genomecov_conf

cd $clusters_bed
ls *.sorted.bed >> ../$genomecov_conf
cd ..

sed -i 's/\.sorted\.bed//g' $genomecov_conf

#Entra nella cartella di R e fa girare gli script di R.

cd R

echo $alignment_check > alignment_check.in

Rscript ../length.distribution.R
mv Rplots.pdf ../output.R/length.hist.pdf

rm alignment_check.in

#Lo script R genomecov.R, oltre a produrre i grafici della coverage, segnala anche i doppioni nei nomi delle read. A quel punto, bisogna modificare a mano i file .bed eliminando il doppione (o i doppioni...) fuori luogo.

Rscript ../genomecov.R
mv Rplots.pdf ../output.R/genomecov.pdf

cd ..

#A questo punto, bisogna controllare a occhio il file genomecov.pdf nella cartella degli output di R e scegliere i cluster "interessanti", scrivendone i numeri nel file 'selected.list'.
#Segue lo script 02.sh...
