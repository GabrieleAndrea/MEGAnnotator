#!/bin/bash

# MEGAnnotator: a Multi-threaded Enhanced prokaryotic Genome Annotator
# Copyright (C) 2015 Lugli Gabriele Andrea

# MEGAnnotator is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# MEGAnnotator is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#Technology selection
TECHNOLOGY=$(zenity --list --radiolist --title "ASAP" --text "What kind of technology do you have used?" --hide-header --column "Select" --column "Technology" FALSE "454" FALSE "IonTorrent" FALSE "Illumina" FALSE "Other")
if [ $TECHNOLOGY = '454' ] ; then
	echo $TECHNOLOGY technology selected
	if INPUTFASTQ=$(zenity --file-selection --title="Select input fastq file" --file-filter='*.fastq *.fq')
		then echo $INPUTFASTQ input selected
		else zenity --warning --title="ASAP stopped" --text="ASAP has been interrupted.\nPlease restart the script and follow the instructions."; exit
	fi
elif [ $TECHNOLOGY = 'IonTorrent' ] ; then
	echo $TECHNOLOGY technology selected
	if INPUTFASTQ=$(zenity --file-selection --title="Select input fastq file" --file-filter='*.fastq *.fq')
		then echo $INPUTFASTQ input selected
		else zenity --warning --title="ASAP stopped" --text="ASAP has been interrupted.\nPlease restart the script and follow the instructions."; exit
	fi
elif [ $TECHNOLOGY = 'Illumina' ] ; then
	echo $TECHNOLOGY technology selected
	READS=$(zenity --list --radiolist --title "ASAP" --text "What kind of data do you have?" --hide-header --column "Select" --column "Analyses" FALSE "Single-Read" FALSE "Paired-End")
	if [ $READS = 'Single-Read' ] ; then
		if INPUTFASTQ=$(zenity --file-selection --title="Select input fastq file" --file-filter='*.fastq *.fq')
			then echo $INPUTFASTQ input selected
			else zenity --warning --title="ASAP stopped" --text="ASAP has been interrupted.\nPlease restart the script and follow the instructions."; exit
		fi
	elif [ $READS = 'Paired-End' ] ; then
		if INPUTFASTQ=$(zenity --file-selection --title="Select input fastq files" --multiple --file-filter='*.fastq *.fq')
			then echo $INPUTFASTQ input selected
			else zenity --warning --title="ASAP stopped" --text="ASAP has been interrupted.\nPlease restart the script and follow the instructions."; exit
		fi
	else
		zenity --error --title="ASAP Error" --text="No option were chosed.\nPlease restart the script and follow the instructions."; exit
	fi
elif [ $TECHNOLOGY = 'Other' ] ; then
	if INPUTFASTQ=$(zenity --file-selection --title="Select input fastq file" --file-filter='*.fastq *.fq')
		then echo $INPUTFASTQ input selected
		else zenity --warning --title="ASAP stopped" --text="ASAP has been interrupted.\nPlease restart the script and follow the instructions."; exit
	fi
else
	zenity --error --title="ASAP Error" --text="No option were chosed.\nPlease restart the script and follow the instructions."; exit
fi

#Project name
if PROJECT=$(zenity --entry --title="Project Name" --text="Enter the project name.\nPlease, enter alphanumeric characters only\n(e.g. AH17 - clone05 - Coli - 1349).")
	then echo $PROJECT project selected
	else zenity --warning --title="ASAP stopped" --text="ASAP has been interrupted.\nPlease restart the script and follow the instructions."; exit
fi

#Reads length
if RLENGTH=$(zenity --entry --title="Reads Length" --text="Enter the length of the sequenced reads.\nPlease, enter numeric characters only\n(e.g. 100 - 150 - 250 - 400).")
	then echo $RLENGTH reads length selected
	else zenity --warning --title="ASAP stopped" --text="ASAP has been interrupted.\nPlease restart the script and follow the instructions."; exit
fi

#Number of threads
if NTHREADS=$(zenity --entry --title="Threads" --text="Enter the number of simultaneous threads\n(based on your machine specification).\nPlease, enter numeric characters only\n(e.g. 2 - 8 - 24 - 64).")
	then echo $NTHREADS treads selected
	else zenity --warning --title="ASAP stopped" --text="ASAP has been interrupted.\nPlease restart the script and follow the instructions."; exit
fi

#RAM threshold
if MAXMEM=$(zenity --entry --title="RAM" --text="Enter the GigaByte threshold of the RAM\n(based on your machine specification).\nPlease, enter numeric characters only\n(e.g. 30 - 60 - 120 - 250).")
	then echo $MAXMEM GigaByte of RAM selected
	else zenity --warning --title="ASAP stopped" --text="ASAP has been interrupted.\nPlease restart the script and follow the instructions."; exit
fi

#Databases selection
if RAPSEARCHDB=$(zenity --file-selection --title="Select RAPSearch database file")
	then echo $RAPSEARCHDB database selected
	else zenity --warning --title="ASAP stopped" --text="ASAP has been interrupted.\nPlease restart the script and follow the instructions."; exit
fi

lib/./metagenome_counter.sh &

#Files selection
mkdir temp ;
echo $PROJECT > temp/project_name ;
echo $INPUTFASTQ > temp/fastq_data ;
if [ ! -z "$REFERENCE" ] ; then
	echo $REFERENCE > temp/reference_genome ;
fi
echo $RLENGTH > temp/reads_length ;
echo $NTHREADS > temp/n_threads ;
echo $MAXMEM > temp/ram_threshold ;
if [ ! -z "$READS" ] ; then
	echo $READS > temp/reads_type ;
else
	echo Single-Read > temp/reads_type ;
	READS=$(cat temp/reads_type) ;
fi
echo $TECHNOLOGY > temp/thecnology_type ;
echo $RAPSEARCHDB > temp/database_rapsearch ;

#Organizing input
mkdir input ;

if [ $RLENGTH -ge 240 ] ; then
	echo 21,33,55,77,99,127 > temp/khmer_parameters ;
elif [ $RLENGTH -ge 190 -a $RLENGTH -le 239 ] ; then
	echo 21,33,55,77,99 > temp/khmer_parameters ;
elif [ $RLENGTH -ge 140 -a $RLENGTH -le 189 ] ; then
	echo 21,33,55,77 > temp/khmer_parameters ;
elif [ $RLENGTH -ge 100 -a $RLENGTH -le 139 ] ; then
	echo 21,33,55 > temp/khmer_parameters ;
elif [ $RLENGTH -le 99 ] ; then
	echo 21,33 > temp/khmer_parameters ;
fi
KMERPAR=$(cat temp/khmer_parameters) ;

#=====================================================================================================================================================================================#
echo "#§#§#	Phase 1: Metagenome Assembly	#§#§#" ;

if [ $READS = 'Single-Read' ] ; then
	cp $INPUTFASTQ input/input.fastq ;
	bin/SPAdes/bin/spades.py -t $NTHREADS --disable-gzip-output --cov-cutoff 5 -m $MAXMEM -k $KMERPAR --careful -s input/input.fastq -o ${PROJECT}_assembly ;
else
	cp temp/fastq_data temp/paired_data ;
	sed -i 's/|/\n/' temp/paired_data ;
	split -l 1 temp/paired_data temp/paired ;
	PAIRED1=$(cat temp/pairedaa) ;
	PAIRED2=$(cat temp/pairedab) ;
	cp $PAIRED1 input/input1.fastq ;
	cp $PAIRED2 input/input2.fastq ;
	bin/SPAdes/bin/spades.py -t $NTHREADS --disable-gzip-output --cov-cutoff 5 -m $MAXMEM -k $KMERPAR --careful -1 input/input1.fastq -2 input/input2.fastq -o ${PROJECT}_assembly ;
fi

mkdir results ;
cp ${PROJECT}_assembly/contigs.fasta results/${PROJECT}_multifasta.fasta ;

#=====================================================================================================================================================================================#
echo "#§#§#     Phase 2: ORFs prediction     #§#§#" ;
bin/./prodigal.linux -f gff -a temp/aaORFs.fasta -p meta -i results/${PROJECT}_multifasta.fasta -o temp/ORFs.gff ;

#=====================================================================================================================================================================================#
echo "#§#§#     Phase 3: RapSearch2 vs nrDB     #§#§#" ;
bin/./rapsearch -q temp/aaORFs.fasta -d $RAPSEARCHDB -o temp/file.rapsearch -s f -e 0.0001 -l 20 -z $NTHREADS -b 10 -v 10 ;

#locus_tag definition
grep -v "^#" temp/file.rapsearch.m8 > temp/file.rapsearch ;
gawk '{print $1}' temp/file.rapsearch > temp/first_column ;
gawk '!x[$0]++' temp/first_column > temp/locus_tag ;

#output modding
sed 's/\ \[/\n/g' temp/file.rapsearch > temp/rapsearch_tbl ;
grep -f temp/locus_tag temp/rapsearch_tbl > temp/rapsearch_tbl0 ;
grep -i -v "MULTISPECIES:" temp/rapsearch_tbl0 > temp/rapsearch_tbl1 ;
grep -i -v "conserved protein" temp/rapsearch_tbl1 > temp/rapsearch_tbl2 ;
grep -i -v "putative" temp/rapsearch_tbl2 > temp/rapsearch_tbl3 ;
grep -i -v "predicted" temp/rapsearch_tbl3 > temp/rapsearch_tbl4 ;
grep -i -v "uncharacterized" temp/rapsearch_tbl4 > temp/rapsearch_tbl5 ;
grep -i -v "unknown" temp/rapsearch_tbl5 > temp/rapsearch_tbl6 ;
grep -i -v "hypothetical" temp/rapsearch_tbl6 > temp/rapsearch_tbl7 ;
grep -i -v "unnamed" temp/rapsearch_tbl7 > temp/rapsearch_tbl8 ;
grep -i -v "membrane protein" temp/rapsearch_tbl8 > temp/rapsearch_tbl9 ;
grep -i -v "secreted protein" temp/rapsearch_tbl9 > temp/rapsearch_tbl10 ;
grep -i -v "RecName" temp/rapsearch_tbl10 > temp/rapsearch_tbl11 ;
grep -i -v "hipothetical" temp/rapsearch_tbl11 > temp/rapsearch_tbl12 ;
grep -i -v "hipothetycal" temp/rapsearch_tbl12 > temp/rapsearch_tbl13 ;
grep -i -v "hypothetycal" temp/rapsearch_tbl13 > temp/rapsearch_tbl14 ;
grep -i -v "hypotetical" temp/rapsearch_tbl14 > temp/rapsearch_tbl15 ;
grep -i -v "membrane spanning" temp/rapsearch_tbl15 > temp/rapsearch_tbl16 ;
grep -i -v "probable" temp/rapsearch_tbl16 > temp/rapsearch_tbl17 ;
grep -i -v "possible" temp/rapsearch_tbl17 > temp/rapsearch_tbl18 ;
grep -i -v "domain protein" temp/rapsearch_tbl18 > temp/rapsearch_tbl19 ;
grep -v "COG" temp/rapsearch_tbl19 > temp/rapsearch_tbl20 ;
grep -i -v "| protein$" temp/rapsearch_tbl20 > temp/rapsearch_tbl21 ;
grep -i -v "Crystal Structure" temp/rapsearch_tbl21 > temp/rapsearch_tbl22 ;
grep -i -v "PF..... family protein" temp/rapsearch_tbl22 > temp/rapsearch_tbl23 ;
grep -i -v "UPF.... protein" temp/rapsearch_tbl23 > temp/rapsearch_tbl24 ;
grep -i -v "| ORF" temp/rapsearch_tbl24 > temp/rapsearch_tbl25 ;
grep -i -v "LOW QUALITY PROTEIN: " temp/rapsearch_tbl25 > temp/rapsearch_tbl26 ;
grep -i -v "uncharacterised" temp/rapsearch_tbl26 > temp/rapsearch_tbl27 ;
grep -i -v "unknwon" temp/rapsearch_tbl27 > temp/rapsearch_tbl28 ;
grep -i -v "unusual" temp/rapsearch_tbl28 > temp/rapsearch_tbl29 ;
grep -i -v "possibl" temp/rapsearch_tbl29 > temp/rapsearch_tbl30 ;
grep -i -v "low molecular" temp/rapsearch_tbl30 > temp/rapsearch_tbl31 ;
grep -i -v "macro domain" temp/rapsearch_tbl31 > temp/rapsearch_tbl32 ;
grep -i -v "small protein" temp/rapsearch_tbl32 > temp/rapsearch_tbl33 ;
grep -i -v "cupin" temp/rapsearch_tbl33 > temp/rapsearch_tbl34 ;
grep -i -v "repeat, PF....." temp/rapsearch_tbl34 > temp/rapsearch_tbl35 ;
grep -i -v "repeat, UPF...." temp/rapsearch_tbl35 > temp/rapsearch_tbl36 ;
grep -i -v "repeat, TIGR....." temp/rapsearch_tbl36 > temp/rapsearch_tbl37 ;
grep -i -v "PF..... repeat" temp/rapsearch_tbl37 > temp/rapsearch_tbl38 ;
grep -i -v "UPF.... repeat" temp/rapsearch_tbl38 > temp/rapsearch_tbl39 ;
grep -i -v "TIGR..... repeat" temp/rapsearch_tbl39 > temp/rapsearch_tbl40 ;
grep -i -v "conserved repeat" temp/rapsearch_tbl40 > temp/rapsearch_tbl41 ;
grep -i -v "hupothetical" temp/rapsearch_tbl41 > temp/rapsearch_tbl42 ;
grep -i -v "DUF... family protein" temp/rapsearch_tbl42 > temp/rapsearch_tbl43 ;
grep -i -v "DUF.... family protein" temp/rapsearch_tbl43 > temp/rapsearch_tbl44 ;
grep -i -v "TIGR..... family protein" temp/rapsearch_tbl44 > temp/file_greppato_rapsearch_temp ;

#takes the first hit of rapsearch
grep "| " temp/file_greppato_rapsearch_temp > temp/file_greppato_rapsearch ;
cat temp/locus_tag | xargs -I{} grep -m 1 -w {} temp/file_greppato_rapsearch > temp/first_hit_rapsearch_temp.txt ;
sed -i 's/\t/\ \t/g' temp/first_hit_rapsearch_temp.txt ;
sed -i 's/|\ /\t/' temp/first_hit_rapsearch_temp.txt ;
gawk -F"\t" '{print $1,"\t",$3}' temp/first_hit_rapsearch_temp.txt > temp/first_hit_rapsearch.txt ;
sed -i 's/\ \t\ /\t/g' temp/first_hit_rapsearch.txt ;

#all the locus tag generation
grep "^>" temp/aaORFs.fasta > temp/all_the_locus_tag ;
sed -i 's/\#/\t/' temp/all_the_locus_tag ;
gawk -F"\t" '{print $1}' temp/all_the_locus_tag > temp/locus_tag_TOT ;
sed -i 's/>//' temp/locus_tag_TOT ;

#=====================================================================================================================================================================================#
echo "#§#§#     Phase 4: merging annotation     #§#§#" ;
while read LTAG
do
let count=$count+1
grep -w $LTAG temp/first_hit_rapsearch.txt > temp/file${count}_NCBI.single
if [ -s temp/file${count}_NCBI.single ] ; then
	grep -m 1 "^" temp/file${count}_NCBI.single >> temp/merged_annotation
else
	echo "$LTAG 	hypothetical protein" >> temp/merged_annotation ;
fi
done < temp/locus_tag_TOT

#gff editing
cat temp/ORFs.gff | gawk '{print $9}' > temp/ORFs2.gff ;
sed -i 's/\_/\t/' temp/ORFs2.gff ;
sed -i 's/;partial=/\t/' temp/ORFs2.gff ;
cat temp/ORFs2.gff | gawk '{print $2}' > temp/ORFs3.gff ;
cat temp/ORFs.gff | gawk '{print $1}' > temp/ORFs4.gff ;
pr -mts_ temp/ORFs4.gff temp/ORFs3.gff > temp/gff_locus_tag ;
sed -i 's/\#_//g' temp/gff_locus_tag ;
sed -i 's/gff-version//g' temp/gff_locus_tag ;
sed -i 's/$/\ /' temp/gff_locus_tag ;
sed -i 's/^\ /NULLNULLNULL/' temp/gff_locus_tag ;

while read LTAG2
do
let count2=$count2+1
grep -w $LTAG2 temp/merged_annotation > temp/file${count2}_gff.single
if [ -s temp/file${count2}_gff.single ] ; then
	grep -m 1 "^" temp/file${count2}_gff.single >> temp/merged_gff
else
	echo "aaaaa" >> temp/merged_gff
fi
done < temp/gff_locus_tag

#gff generation
cat temp/merged_gff | gawk -F"\t" '{print $2}' > temp/merged_gff_C2 ;
pr -mts@ temp/ORFs.gff temp/merged_gff_C2 > temp/annotated_ORFs.gff ;
sed -i 's/\@$//' temp/annotated_ORFs.gff ;
sed 's/\@/product=/' temp/annotated_ORFs.gff > results/${PROJECT}_annotated_ORFs.gff ;
cp temp/merged_annotation results/${PROJECT}_annotation_list ;

#=====================================================================================================================================================================================#
echo "#§#§#     Phase 5: gbk generation and rRNA prediction     #§#§#" ;
mkdir temp/single_gffs ;
mkdir temp/single_contigs ;
mkdir temp/single_gbks ;
mkdir temp/single_contigs_length ;
mkdir temp/single_rRNA ;

#genbank creation with rRNA prediction
grep "^>" results/${PROJECT}_multifasta.fasta > temp/contigs_name ;
sed -i 's/>//' temp/contigs_name ;
sed "s/\tID=/\tlocus_tag=${PROJECT}_/" results/${PROJECT}_annotated_ORFs.gff > temp/annotated_ORFs_locus_tag.gff ;
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < results/${PROJECT}_multifasta.fasta > temp/multifasta_1line.fasta ;

while read LTAG3
do
let count3=$count3+1
	grep -P "$LTAG3\t" temp/annotated_ORFs_locus_tag.gff > temp/single_gffs/contig_${count3}.gff
	grep -A 1 "$LTAG3$" temp/multifasta_1line.fasta > temp/single_contigs/contig_${count3}.fasta
	seqret -sequence temp/single_contigs/contig_${count3}.fasta -feature -fformat gff -fopenfile temp/single_gffs/contig_${count3}.gff -osformat genbank -auto -outseq temp/single_gbks/contig_${count3}_tmp.gbk
	grep -v "/note" temp/single_gbks/contig_${count3}_tmp.gbk > temp/single_gbks/contig_${count3}.gbk

	rnammer -S bac temp/single_contigs/contig_${count3}.fasta -m tsu,ssu,lsu -gff temp/single_rRNA/rRNA_${count3}.gff
	grep -v "#" temp/single_rRNA/rRNA_${count3}.gff > temp/single_rRNA/rRNA_contig_${count3} ;
	if [ -s temp/single_rRNA/rRNA_contig_${count3} ] ; then
		sed -i 's/5s_rRNA/locus_tag=XXX_XrRNA5S;product=5s ribosomal RNA/g' temp/single_rRNA/rRNA_contig_${count3}
		sed -i 's/16s_rRNA/locus_tag=XXX_XrRNA16S;product=16s ribosomal RNA/g' temp/single_rRNA/rRNA_contig_${count3}
		sed -i 's/23s_rRNA/locus_tag=XXX_XrRNA23S;product=23s ribosomal RNA/g' temp/single_rRNA/rRNA_contig_${count3}
		cat temp/single_rRNA/rRNA_contig_${count3} | gawk '{print $3,"\t",$4,"\t",$5,"\t",$9}' > temp/single_rRNA/rRNA_contig_${count3}_awkapped
		cat temp/single_rRNA/rRNA_contig_${count3} | gawk '{print $7}' > temp/single_rRNA/rRNA_contig_${count3}_strend
		count6=0

		while read LTAG6
		do
		let count6=$count6+1
		awk "NR==$count6" temp/single_rRNA/rRNA_contig_${count3}_awkapped > temp/single_rRNA/rRNAline_${count6}_of_contig_${count3}.rrna
			if [ "$LTAG6" = "+" ] ; then
				sed -i 's/$/\"/' temp/single_rRNA/rRNAline_${count6}_of_contig_${count3}.rrna
				sed -i 's|rRNA\ \t\ |\ \ \ \ \ rRNA\ \ \ \ \ \ \ \ \ \ \ \ |' temp/single_rRNA/rRNAline_${count6}_of_contig_${count3}.rrna
				sed -i "s|\ \t\ locus_tag=XXX_X|\n\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ /locus_tag=\"${PROJECT}_|" temp/single_rRNA/rRNAline_${count6}_of_contig_${count3}.rrna
				sed -i 's|;|\"\n\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ /|' temp/single_rRNA/rRNAline_${count6}_of_contig_${count3}.rrna
				sed -i 's|product=|product="|' temp/single_rRNA/rRNAline_${count6}_of_contig_${count3}.rrna
				sed -i 's|\ \t\ |..|' temp/single_rRNA/rRNAline_${count6}_of_contig_${count3}.rrna
			elif [ "$LTAG6" = "-" ] ; then
				sed -i 's/$/\"/' temp/single_rRNA/rRNAline_${count6}_of_contig_${count3}.rrna
				sed -i 's|rRNA\ \t\ |\ \ \ \ \ rRNA\ \ \ \ \ \ \ \ \ \ \ \ complement(|' temp/single_rRNA/rRNAline_${count6}_of_contig_${count3}.rrna
				sed -i "s|\ \t\ locus_tag=XXX_X|)\n\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ /locus_tag=\"${PROJECT}_|" temp/single_rRNA/rRNAline_${count6}_of_contig_${count3}.rrna
				sed -i 's|;|\"\n\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ /|' temp/single_rRNA/rRNAline_${count6}_of_contig_${count3}.rrna
				sed -i 's|product=|product="|' temp/single_rRNA/rRNAline_${count6}_of_contig_${count3}.rrna
				sed -i 's|\ \t\ |..|' temp/single_rRNA/rRNAline_${count6}_of_contig_${count3}.rrna
			fi
		done < temp/single_rRNA/rRNA_contig_${count3}_strend

		cat temp/single_rRNA/*.rrna > temp/premature_contig_${count3}_rRNA.gff
		sed -i 's/5s/5s ribosomal RNA/' temp/premature_contig_${count3}_rRNA.gff
		sed -i 's/16s/16s ribosomal RNA/' temp/premature_contig_${count3}_rRNA.gff
		sed -i 's/23s/23s ribosomal RNA/' temp/premature_contig_${count3}_rRNA.gff
		rm temp/single_rRNA/*.rrna

		if [ -s temp/single_gffs/contig_${count3}.gff ] ; then
			cat temp/single_gbks/contig_${count3}.gbk | head -n 2 > temp/head_${count3}.gbk
			cat temp/single_gbks/contig_${count3}.gbk | tail -n +3 > temp/tail_${count3}.gbk
			cat temp/head_${count3}.gbk temp/premature_contig_${count3}_rRNA.gff temp/tail_${count3}.gbk > temp/single_gbks/contig_${count3}_with_rRNA.gbk
		else
			cat temp/single_gbks/contig_${count3}.gbk | head -n 1 > temp/head_${count3}.gbk
			echo "FEATURES             Location/Qualifiers" >> temp/head_${count3}.gbk
			cat temp/single_gbks/contig_${count3}.gbk | tail -n +2 > temp/tail_${count3}.gbk
			cat temp/head_${count3}.gbk temp/premature_contig_${count3}_rRNA.gff temp/tail_${count3}.gbk > temp/single_gbks/contig_${count3}_with_rRNA.gbk
		fi
	fi
	
	if [ $count3 -eq 1 ] ; then
		if [ -e temp/single_gbks/contig_${count3}_with_rRNA.gbk ] ; then
			mv temp/single_gbks/contig_${count3}_with_rRNA.gbk temp/single_gbks/catted_union.gbk
		else
			mv temp/single_gbks/contig_${count3}.gbk temp/single_gbks/catted_union.gbk
		fi
	else
		if [ -e temp/single_gbks/contig_${count3}_with_rRNA.gbk ] ; then
			cat temp/single_gbks/catted_union.gbk temp/single_gbks/contig_${count3}_with_rRNA.gbk > temp/single_gbks/catted_gbk
			union -sequence temp/single_gbks/catted_gbk -sformat genbank -outseq temp/single_gbks/catted_union.gbk -osformat genbank -feature Y -auto
			rm temp/single_gbks/catted_gbk
		else
			cat temp/single_gbks/catted_union.gbk temp/single_gbks/contig_${count3}.gbk > temp/single_gbks/catted_gbk
			union -sequence temp/single_gbks/catted_gbk -sformat genbank -outseq temp/single_gbks/catted_union.gbk -osformat genbank -feature Y -auto
			rm temp/single_gbks/catted_gbk
		fi
	fi
done < temp/contigs_name

grep -v "/partial" temp/single_gbks/catted_union.gbk > temp/union.gbk ;
sed -i "s/LOCUS       Contig_01/LOCUS       ${PROJECT}/" temp/union.gbk ;

#fasta features editing
readseq -a -f=2 results/${PROJECT}_multifasta.fasta > temp/multigbk_fromfasta.gbk ;
grep "^LOCUS" temp/multigbk_fromfasta.gbk > temp/contigs_locus ;
sed -i 's/\ \ \ \ \ \ \ /\t/g' temp/contigs_locus ;
cat temp/contigs_locus | gawk '{print $3}' > temp/contigs_length ;
echo 1 > temp/contigs_starts ;

while read LTAG4
do
let count4=$count4+1
	head -$count4 temp/contigs_length > temp/single_contigs_length/${LTAG4}_length
	cat temp/single_contigs_length/${LTAG4}_length | awk '{sum+=$1} END {print sum}' > temp/single_contigs_length/${LTAG4}_end
	TMPVAR=$(cat temp/single_contigs_length/${LTAG4}_end)
	echo $TMPVAR >> temp/contigs_ends
	let TMPVAR=$TMPVAR+1
	echo $TMPVAR >> temp/contigs_starts
done < temp/contigs_name

sed -i '$ d' temp/contigs_starts ;
echo "#${PROJECT}_fasta_features_gff" > temp/fasta_features_tmp.gff ;

while read LTAG5
do
let count5=$count5+1
	evenodd=$(($count5%2))
	if [ $evenodd -eq 0 ] ; then
		echo "gff_seqname	artemis	fasta_record	0000	9999	.	+	.	label=$LTAG5;colour=11" >> temp/fasta_features_tmp.gff
	else
		echo "gff_seqname	artemis	fasta_record	0000	9999	.	+	.	label=$LTAG5;colour=10" >> temp/fasta_features_tmp.gff
	fi
done < temp/contigs_name

#fasta features file generation
sed -i '1s/^/\n/' temp/contigs_starts ;
sed -i '1s/^/\n/' temp/contigs_ends ;
cat temp/fasta_features_tmp.gff | gawk '{print $1,"\t",$2,"\t",$3}' > temp/fasta_features_1-3.gff ;
cat temp/fasta_features_tmp.gff | gawk '{print $6,"\t",$7,"\t",$8,"\t",$9}' > temp/fasta_features_6-9.gff ;
pr -mts temp/fasta_features_1-3.gff temp/contigs_starts temp/contigs_ends temp/fasta_features_6-9.gff > temp/fasta_features.gff ;

cat temp/fasta_features.gff | tail -n +2 > temp/fasta_features_edit.gff ;
sed -i 's|gff_seqname\ \t\ artemis\ \t\ |\ \ \ \ \ |' temp/fasta_features_edit.gff ;
sed -i 's|fasta_record\t|fasta_record\ \ \ \ |' temp/fasta_features_edit.gff ;
sed -i 's|\t.\ \t\ +\ \t\ .\ \t\ |\n\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ /|' temp/fasta_features_edit.gff ;
sed -i 's|;|\n\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ /|' temp/fasta_features_edit.gff ;
sed -i 's|\t|..|' temp/fasta_features_edit.gff ;

#=====================================================================================================================================================================================#
echo "#§#§#     Phase 6: genbank finalization and tRNA prediction     #§#§#" ;
union -sequence results/${PROJECT}_multifasta.fasta -sformat fasta -snucleotide -outseq temp/allnt.fasta -osformat fasta -auto ;
tRNAscan-SE -B -o results/${PROJECT}_tRNA temp/allnt.fasta ;

#tRNa genes editing
cat results/${PROJECT}_tRNA | gawk '{print $3,"\t",$4,"\t",$5}' > temp/tRNA_awkapped ;
sed '1,3d' temp/tRNA_awkapped > temp/tRNA_sedded ;
sort -k1 -n temp/tRNA_sedded > temp/tRNA_sorted ;
cat temp/tRNA_sorted | gawk '{print $1}' > temp/tRNA_column1 ;
cat temp/tRNA_sorted | gawk '{print $2}' > temp/tRNA_column2 ;
cat temp/tRNA_sorted | gawk '{print $3}' > temp/tRNA_column3 ;
sed -i 's/^/tRNA-/' temp/tRNA_column3 ;

while read LTAG7
do
let count7=$count7+1
awk "NR==$count7" temp/tRNA_column1 > temp/tRNA_column1_${count7}
awk "NR==$count7" temp/tRNA_column2 > temp/tRNA_column2_${count7}
NUMONE=$(cat temp/tRNA_column1_${count7})
NUMTWO=$(cat temp/tRNA_column2_${count7})
	if [ "$NUMONE" -gt "$NUMTWO" ] ; then
		echo -e "     tRNA            complement(${NUMTWO}..${NUMONE})\n                     /locus_tag=\"${PROJECT}_tRNA${count7}\"\n                     /product=\"${LTAG7}\"" > temp/tRNA_${count7}.trna
	elif [ "$NUMTWO" -gt "$NUMONE" ] ; then
		echo -e "     tRNA            ${NUMONE}..${NUMTWO}\n                     /locus_tag=\"${PROJECT}_tRNA${count7}\"\n                     /product=\"${LTAG7}\"" > temp/tRNA_${count7}.trna
	fi
done < temp/tRNA_column3

cat temp/*.trna > temp/premature_tRNA.gff ;

#genbank file generation
cat temp/union.gbk | head -n 2 > temp/head_union.gbk ;
cat temp/union.gbk | tail -n +3 > temp/tail_union.gbk ;
cat temp/head_union.gbk temp/fasta_features_edit.gff temp/premature_tRNA.gff temp/tail_union.gbk > results/${PROJECT}.gbk ; sleep 2s ;

#additional file formats
seqret results/${PROJECT}.gbk results/${PROJECT}.gff3 --feature -osf gff3 ;
seqret results/${PROJECT}.gbk results/${PROJECT}.embl --feature -osf embl ;
java -cp bin/readseq.jar run -feat=CDS,rRNA,tRNA -format=19 results/${PROJECT}.gbk -o results/${PROJECT}.xml ;

#cleaning
rm -R input ;
rm -R temp ;
rm results/${PROJECT}_annotated_ORFs.gff ;
rm results/${PROJECT}_tRNA ;
mv results ${PROJECT}_results ;

echo "#§#§#     Prediction Complete     #§#§#" ;

