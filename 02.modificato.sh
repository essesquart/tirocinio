#!/bin/sh

#Stabilisce il nome dei file di configurazione per R. Cambiando questa opzione, va cambiata anche corrispondentemente negli script R.

conf_temp=conf.temp
genomecov_temp=genomecov.temp

#Prepara un po' di nomi.
#Si fa partire questo script all'interno di BedFiles(?)

output=`cat ../output.basename`
output_selected_centroids=$output.centroids.selected.fa
clusters_bed=$output.clusters.bedfiles
clusters_fasta=$output.clusters.fasta

#Segue dallo script 01.sh...
#A questo punto, bisogna aver controllato a occhio il file genomecov.pdf nella cartella degli output di R e scegliere i cluster "interessanti", scrivendone i numeri nel file 'selected.list'.
#Costruisce il file 'smithRNAs.fa', dove mette i centroidi scelti a mano dal file 'selected.list' e li annota con il cluster e le posizioni sul genoma, che prende dai bedfile.

if [ -f smithRNAs.fa ]
	then rm smithRNAs.fa
	fi
for c in `cat selected.list`
	do grep ">" $clusters_fasta/$c.fa > seqs.list
	grep -A 1 -f seqs.list $output_selected_centroids > smithRNA.fa
	sed -i 's/>//g' smithRNA.fa
	sed -i 's/_/:/g' smithRNA.fa
	seq=`head -1 smithRNA.fa`
	pos_start=`grep -w -m 1 $seq $clusters_bed/$c.sorted.bed | awk '{print $2}'`
	pos_start=$(( pos_start+1 ))
	pos_end=`grep -w -m 1 $seq $clusters_bed/$c.sorted.bed | awk '{print $3}'`
	strand=`grep -w -m 1 $seq $clusters_bed/$c.sorted.bed | awk '{print $6}'`
	if [ $strand == "+" ]
		then sed -i '1 s/.\+/>smithRNA_C'$c'_'$pos_start'_'$pos_end'/g' smithRNA.fa
		else sed -i '1 s/.\+/>smithRNA_C'$c'_'$pos_end'_'$pos_start'/g' smithRNA.fa
		fi
	cat smithRNA.fa >> smithRNAs.fa
	done
rm seqs.list
rm smithRNA.fa

#A questo punto, a mano, bisogna andare a pescare le regioni da cui trascrivono gli smithRNA da usare come pre-smithRNA (l'intero tRNA o l'intera UR, per esempio). Si prosegue poi con lo script '03.sh'.

#Blastenna gli smithRNA contro il trascrittoma assemblato con Trinity e calcola le energie libere dell'ibrido RNA-RNA. Le corrispondenze esatte vanno greppate, perché per le sequenze corte blastn funziona male.

#Per prima cosa, usiamo revseq (di EMBOSS) per reversare e complementare i seed ('-sbegin1 4 -send1 10') e gli smithRNA interi.

revseq smithRNAs.fa -tag FALSE -sbegin1 4 -send1 10 -outseq RC_seed_smithRNAs.fa
revseq smithRNAs.fa -tag FALSE -outseq RC_smithRNAs.fa

#Contiamo quanti smith abbiamo.

num_smith=`wc -l smithRNAs.fa | awk '{print $1}'`
num_smith=$(( num_smith/2 ))
for s in `seq 1 $num_smith`

	#Estrae il nome di uno smithRNA e su quello poi lavora.

	do curr_smith=`head -n $(( s*2-1 )) smithRNAs.fa | tail -n 1`
	curr_smith=`echo $curr_smith | sed 's/>//g'`
	if [ -f $curr_smith.BLAST.target ]
		then rm $curr_smith.BLAST.target
		fi
	if [ -f $curr_smith.selected.UTR.fa ]
		then rm $curr_smith.selected.UTR.fa
		fi
	if [ -f $curr_smith.targets ]
		then rm $curr_smith.targets
		fi

	#Estrae seed e RNA intero reversato e complementato.

	curr_seed=`head -n $(( s*2 )) RC_seed_smithRNAs.fa | tail -n 1`
	curr_RNA=`head -n $(( s*2 )) RC_smithRNAs.fa | tail -n 1`

	#Cerca il seed nelle 3' UTR: questo ci deve essere esattamente, per cui, per evitare storie con blastn, usiamo il grep. Costruisce il file .seed, il FASTA delle 3' UTR che soddisfano questo primo requisito.

	grep -B 1 $curr_seed ../../Transcriptome/3UTR/*_3UTR.fasta > $curr_smith.seed
	sed -i 's/--//g' $curr_smith.seed
	sed -i 's/\([^\t]*\).*/\1/g' $curr_smith.seed
	sed -i 's/^>.\+/#####&#####/g' $curr_smith.seed
	sed -i ':a;N;$!ba;s/\n//g' $curr_smith.seed
	sed -i 's/#####/\n/g' $curr_smith.seed
	sed -i '1 d' $curr_smith.seed

	#Per sicurezza, greppa nel file .seed anche lo smithRNA intero (reversato e complementato) e lo scrive nel file .target: in realtà, di solito, tutti i trascritti che finiscono qui li trova anche BLAST.

	grep -B 1 $curr_RNA $curr_smith.seed > $curr_smith.target
	sed -i 's/--//g' $curr_smith.target
	sed -i 's/^>.\+/#####&#####/g' $curr_smith.target
	sed -i ':a;N;$!ba;s/\n//g' $curr_smith.target
	sed -i 's/#####/\n/g' $curr_smith.target
	sed -i '1 d' $curr_smith.target

	#Costruisce un database BLAST sul file .seed e ci cerca lo smithRNA (solo in Plus/Minus - '-strand minus') con l'opzione per le sequenze corte ('-task blastn-short').
	#I risultati con almeno 11 match vengono salvati nel file .BLAST.target.

	makeblastdb -in $curr_smith.seed -dbtype nucl -title "$curr_smith"_3UTR -out "$curr_smith"_3UTR
	head -n $(( s*2 )) smithRNAs.fa | tail -n 2 > $curr_smith.query.fa
	blastn -task blastn-short -db "$curr_smith"_3UTR -query $curr_smith.query.fa -strand minus -out $curr_smith.temp.BLAST.target -outfmt '6 qseqid sseqid pident nident length mismatch gapopen qstart qend sstart send sstrand evalue'
	num_lines=`wc -l $curr_smith.temp.BLAST.target | awk '{print $1}'`
	for h in `seq 1 $num_lines`
		do if [ `head -n $h $curr_smith.temp.BLAST.target | tail -n 1 | awk '{print $4}'` -ge 11 ]
			then head -n $h $curr_smith.temp.BLAST.target | tail -n 1 >> $curr_smith.BLAST.target
			fi
		done

	#Recupera le sequenze delle 3' UTR che hanno superato entrambi i test (match perfetto del seed; almeno 11 match in generale) e le scrive nel file .selected.UTR.fa.

	for U in `awk '{print $2}' $curr_smith.BLAST.target`
		do grep -A 1 $U ../../Transcriptome/3UTR/*_3UTR.fasta >> $curr_smith.selected.UTR.fa
		sed -i 's/\([^\t]*\).*/\1/g' $curr_smith.selected.UTR.fa
		done

	#Fa andare PITA e RNAhybrid (se nel file .selected.UTR.fa c'è scritto qualcosa!).

	perl /usr/local/PITA-v6/pita_prediction.pl -utr $curr_smith.selected.UTR.fa -mir $curr_smith.query.fa -flank_up 3 -flank_down 15 -prefix $curr_smith
	if [ -s $curr_smith.selected.UTR.fa ]
		then RNAhybrid -f 3,10 -e -20 -p 0.05 -s 3utr_fly -t $curr_smith.selected.UTR.fa -q $curr_smith.query.fa > $curr_smith.RNAhybrid_results.out
		fi

	#Mette insieme i risultati. Per prima cosa, selezioniamo le hit con meno di -9 kJ di ddG di PITA usando un miniscript R perché sono offeso con bash.

	echo $curr_smith > curr_smith.in
	cd R
	Rscript ../pita_results.R
	cd ..

	#Chi è lo re?

	num_lines=`wc -l "$curr_smith"_pita_results.tab | awk '{print $1}'`
	if [ $num_lines -ge 2 ]
		then for t in `seq 2 $num_lines`
			do curr_PITA_target=`head -n $t "$curr_smith"_pita_results.tab | tail -n 1 | awk -F"\t" '{print $1}'`
			too_long_grep=`grep -c "target too long: $curr_PITA_target" $curr_smith.RNAhybrid_results.out`
			curr_grep=`grep -c $curr_PITA_target $curr_smith.RNAhybrid_results.out`
			if [ $too_long_grep -eq 0 -a $curr_grep -ne 0 ]
				then grep -A 1 $curr_PITA_target ../../Transcriptome/3UTR/*_3UTR_orfs.fa >> $curr_smith.targets

				#Toglie i duplicati nei file dei target, perché ad ARGOT non piacciono.

				cat $curr_smith.targets | paste - - > tab_$curr_smith.targets
				sort -u tab_$curr_smith.targets > sorted_tab_$curr_smith.targets
				awk '{print $1"\n"$2}' sorted_tab_$curr_smith.targets > $curr_smith.target.fa
				fi
			done
		fi
	done

#Scrive i nomi delle colonne nel file dei risultati di BLAST, per rendere più semplice la lettura e il controllo.

sed -i '1 i qseqid\tsseqid\tpident\tnident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tsstrand\tevalue' *.BLAST.target

#Scrive il file definitivo degli smithRNA, lasciando solo quelli che hanno passato tutti i filtri.

if [ -f validated.smithRNAs.fa ]
	then rm validated.smithRNAs.fa
	fi
for validated in `ls *.target.fa | sed 's/\.target\.fa//g'`
	do grep -A 1 $validated smithRNAs.fa >> validated.smithRNAs.fa
	done
mv smithRNAs.fa putative.smithRNAs.fa
mv validated.smithRNAs.fa smithRNAs.fa

#Cancella i file temporanei.

rm *.targets
rm curr_smith.in
rm *.temp.BLAST.target
rm *.nhr
rm *.nin
rm *.nsq
rm RC*
rm *.seed
rm *.selected.UTR.fa
rm *.query.fa
rm *_pita_results_targets.tab

#Cancella i file temporanei: il file di configurazione di R, il file con il genoma mitocondriale con un nome standard, il file di genoma e la cartella con i FASTA dei cluster selezionati.

#rm -r $clusters_fasta
