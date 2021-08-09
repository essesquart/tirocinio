#!/bin/sh

#Stabilisce la temperatura di folding del RNA

RNAfoldT=20

#Prepara un po' di nomi: tutti i file si chiameranno 'output' e la sigla della specie.

output=`cat ../output.basename`
output_mit=`cat genome_mit.fa`
pre_dG=$output.RNAfold

#Segue dallo script 03.sh...
#A questo punto, bisogna aver preparato il file 'pre_smithRNA_regions.list' con i nomi giusti dei file dei pre_smithRNA, sulla falsariga dei nomi degli smith veri e propri.
#Bisogna che nel file ci sia una tabella in tre colonne: nome dei pre_smith (come quello degli smith, ma con 'pre_' davanti), inizio sul genoma mitocondriale e fine.
#Va a pescare nel genoma le parti richieste e le ripiega con RNAfold - occhio alla temperatura di folding!

if [ -f pre_smithRNAs.fa ]
	then rm pre_smithRNAs.fa
	fi
num_pre_smith=`wc -l pre_smithRNA_regions.list | awk '{print $1}'`
pre_smith=1
while  [ $pre_smith -le $num_pre_smith ]
	do boundaries=`sed -n ''$pre_smith' p' pre_smithRNA_regions.list`
	name=`echo $boundaries | awk '{print $1}'`
	start=`echo $boundaries | awk '{print $2}'`
	end=`echo $boundaries | awk '{print $3}'`
	if [ $start -lt $end ]
		then extractalign $output_mit -outseq C$pre_smith.pre.fa -regions "$start $end"
		else extractalign $output_mit -outseq C$pre_smith.pre.fa -regions "$end $start"
		revseq C$pre_smith.pre.fa -tag FALSE -outseq RC_C$pre_smith.pre.fa
		rm C$pre_smith.pre.fa
		mv RC_C$pre_smith.pre.fa C$pre_smith.pre.fa
		fi
	sed -i '1 d' C$pre_smith.pre.fa
	sed -i '1 i >'$name'_'$start'_'$end'' C$pre_smith.pre.fa
	sed -i 's/ /_/g' C$pre_smith.pre.fa
	cat C$pre_smith.pre.fa >> pre_smithRNAs.fa
	pre_smith=$(( pre_smith+1 ))
	done
rm *.pre.fa

if [ -d $pre_dG ]
	then rm -r $pre_dG
	fi
mkdir $pre_dG
RNAfold -T $RNAfoldT pre_smithRNAs.fa > pre_smithRNAs.dG
mv *.ps $pre_dG

#Cancella i file temporanei: il file di configurazione di R, il file con il genoma mitocondriale con un nome standard, il file di genoma e la cartella con i FASTA dei cluster selezionati.

rm ../output.basename
rm genome_mit.fa
