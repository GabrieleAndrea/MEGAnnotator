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

#Info
FILE=lib/INFO
if zenity --text-info --title="MEGAnnotator" --filename=$FILE --checkbox="I read and accept the terms." --width=550 --height=650
	then echo Starting MEGAnnotator settings selection
	else zenity --warning --title="MEGAnnotator stopped" --text="MEGAnnotator has been interrupted.\nPlease restart the script and follow the instructions."; exit
fi

#Analysis selection
ANALYSIS=$(zenity --list --radiolist --title "MEGAnnotator" --text "What kind of analysis do you want to perform?" --hide-header --column "Select" --column "Analyses" FALSE "Genomic_Assembly" FALSE "Metagenomic_Assembly" FALSE "Prediction_and_Annotation_only")
if [ $ANALYSIS = 'Genomic_Assembly' ] ; then
	echo $ANALYSIS selected
elif [ $ANALYSIS = 'Metagenomic_Assembly' ] ; then
	echo $ANALYSIS selected ;
	./metagenomic.sh ;
	exit
elif [ $ANALYSIS = 'Prediction_and_Annotation_only' ] ; then
	echo $ANALYSIS selected ;
	./annotation.sh ;
	exit
else
	zenity --error --title="MEGAnnotator Error" --text="No option were chosed.\nPlease restart the script and follow the instructions."; exit
fi

#Technology selection
TECHNOLOGY=$(zenity --list --radiolist --title "MEGAnnotator" --text "What kind of technology do you have used?" --hide-header --column "Select" --column "Technology" FALSE "454" FALSE "IonTorrent" FALSE "Illumina" FALSE "Other")
if [ $TECHNOLOGY = '454' ] ; then
	echo $TECHNOLOGY technology selected
	if INPUTFASTQ=$(zenity --file-selection --title="Select input fastq file" --file-filter='*.fastq *.fq')
		then echo $INPUTFASTQ input selected
		else zenity --warning --title="MEGAnnotator stopped" --text="MEGAnnotator has been interrupted.\nPlease restart the script and follow the instructions."; exit
	fi
elif [ $TECHNOLOGY = 'IonTorrent' ] ; then
	echo $TECHNOLOGY technology selected
	if INPUTFASTQ=$(zenity --file-selection --title="Select input fastq file" --file-filter='*.fastq *.fq')
		then echo $INPUTFASTQ input selected
		else zenity --warning --title="MEGAnnotator stopped" --text="MEGAnnotator has been interrupted.\nPlease restart the script and follow the instructions."; exit
	fi
elif [ $TECHNOLOGY = 'Illumina' ] ; then
	echo $TECHNOLOGY technology selected
	READS=$(zenity --list --radiolist --title "MEGAnnotator" --text "What kind of data do you have?" --hide-header --column "Select" --column "Analyses" FALSE "Single-Read" FALSE "Paired-End")
	if [ $READS = 'Single-Read' ] ; then
		if INPUTFASTQ=$(zenity --file-selection --title="Select input fastq file" --file-filter='*.fastq *.fq')
			then echo $INPUTFASTQ input selected
			else zenity --warning --title="MEGAnnotator stopped" --text="MEGAnnotator has been interrupted.\nPlease restart the script and follow the instructions."; exit
		fi
	elif [ $READS = 'Paired-End' ] ; then
		if INPUTFASTQ=$(zenity --file-selection --title="Select input fastq files" --multiple --file-filter='*.fastq *.fq')
			then echo $INPUTFASTQ input selected
			else zenity --warning --title="MEGAnnotator stopped" --text="MEGAnnotator has been interrupted.\nPlease restart the script and follow the instructions."; exit
		fi
	else
		zenity --error --title="MEGAnnotator Error" --text="No option were chosed.\nPlease restart the script and follow the instructions."; exit
	fi
elif [ $TECHNOLOGY = 'Other' ] ; then
	if INPUTFASTQ=$(zenity --file-selection --title="Select input fastq file" --file-filter='*.fastq *.fq')
		then echo $INPUTFASTQ input selected
		else zenity --warning --title="MEGAnnotator stopped" --text="MEGAnnotator has been interrupted.\nPlease restart the script and follow the instructions."; exit
	fi
else
	zenity --error --title="MEGAnnotator Error" --text="No option were chosed.\nPlease restart the script and follow the instructions."; exit
fi

#Project name
if PROJECT=$(zenity --entry --title="Project Name" --text="Enter the project name.\nPlease, enter alphanumeric characters only\n(e.g. AH17 - clone05 - Coli - 1349).")
	then echo $PROJECT project selected
	else zenity --warning --title="MEGAnnotator stopped" --text="MEGAnnotator has been interrupted.\nPlease restart the script and follow the instructions."; exit
fi

#Number of threads
if NTHREADS=$(zenity --entry --title="Threads" --text="Enter the number of simultaneous threads\n(based on your machine specification).\nPlease, enter numeric characters only\n(e.g. 2 - 8 - 24 - 64).")
	then echo $NTHREADS treads selected
	else zenity --warning --title="MEGAnnotator stopped" --text="MEGAnnotator has been interrupted.\nPlease restart the script and follow the instructions."; exit
fi

#Reference genome
zenity --question --title="Reference" --text="Did you have a reference genome?"
if [ $? = 0 ] ; then
	if REFERENCE=$(zenity --file-selection --title="Select reference fasta files" --file-filter='*.fasta *.fsa *.fna')
		then echo $REFERENCE reference selected
		else zenity --warning --title="MEGAnnotator stopped" --text="MEGAnnotator has been interrupted.\nPlease restart the script and follow the instructions."; exit
	fi
else
	echo No reference has been selected
fi

#Contigs selection parameters
if [ $ANALYSIS = 'Genomic_Assembly' ] ; then
	zenity --question --title="Contigs selection parameters" --ok-label="Parameters are fine" --cancel-label="Edit Parameters" --text="Contigs selection parameters:\n- Minimum contig length = 1000\n- Minimum reads per contig 100\n- Minimum contig coverage = 33% of the final assemble average\nEdit values if necessary."
	if [ $? = 0 ] ; then
		CLENGTH=1000 ;
		RCONTIG=100 ;
		echo Minimum contigs length $CLENGTH ;
		echo Minimum reads per contig $RCONTIG ;
	else 
		if CLENGTH=$(zenity --entry --title="Minimum contigs length" --text="Enter the minimum contigs length.\nPlease, enter numeric characters only\n(e.g. 2000 - 1000 - 200).")
			then echo Minimum contigs length $CLENGTH
			else zenity --warning --title="MEGAnnotator stopped" --text="MEGAnnotator has been interrupted.\nPlease restart the script and follow the instructions."; exit
		fi
		if RCONTIG=$(zenity --entry --title="Minimum reads per contig" --text="Enter the minimum reads per contig.\nPlease, enter numeric characters only\n(e.g. 500 - 100 - 10).")
			then echo Minimum contigs length $RCONTIG
			else zenity --warning --title="MEGAnnotator stopped" --text="MEGAnnotator has been interrupted.\nPlease restart the script and follow the instructions."; exit
		fi
		if zenity --info --title="Minimum contig coverage" --text="The minimum contigs coverage is calculated at the end of the assembly.\nThe variable cannot be changed."
			then echo Minimum contig coverage = 33% of the final assemble average
			else zenity --warning --title="MEGAnnotator stopped" --text="MEGAnnotator has been interrupted.\nPlease restart the script and follow the instructions."; exit
		fi
	fi
fi

#Databases selection
if RAPSEARCHDB=$(zenity --file-selection --title="Select RAPSearch database file")
	then echo $RAPSEARCHDB database selected
	else zenity --warning --title="MEGAnnotator stopped" --text="MEGAnnotator has been interrupted.\nPlease restart the script and follow the instructions."; exit
fi
if PFAMDB=$(zenity --file-selection --title="Select Pfam-A.hmm database file")
	then echo $PFAMDB database selected
	else zenity --warning --title="MEGAnnotator stopped" --text="MEGAnnotator has been interrupted.\nPlease restart the script and follow the instructions."; exit
fi

lib/./genome_counter.sh &

#Files selection
mkdir temp ;
echo $PROJECT > temp/project_name ;
echo $INPUTFASTQ > temp/fastq_data ;
if [ ! -z "$REFERENCE" ] ; then
	echo $REFERENCE > temp/reference_genome ;
fi
echo $NTHREADS > temp/n_threads ;
echo $ANALYSIS > temp/assembly_type ;
echo $READS > temp/reads_type ;
echo $TECHNOLOGY > temp/thecnology_type ;
echo $RAPSEARCHDB > temp/database_rapsearch ;
echo $PFAMDB > temp/database_pfam ;

#Organizing input
mkdir input ;
if [ -e temp/reference_genome ] ; then 
	cp $REFERENCE input/reference.fasta
fi
cp bin/mira4.0.2 mira ;
cp temp/project_name temp/project_name_edited ;
sed -i '1iproject =' temp/project_name_edited ;
sed -i ':a;N;$!ba;s/\n/ /g' temp/project_name_edited ;

if [ $TECHNOLOGY = '454' ] ; then
	cp lib/454_manifest.conf manifest.conf ;
	sed -i '5r temp/project_name_edited' manifest.conf ;
	cp $INPUTFASTQ input/input.fastq ;
elif [ $TECHNOLOGY = 'IonTorrent' ] ; then
	cp lib/iontorrent_manifest.conf manifest.conf ;
	sed -i '5r temp/project_name_edited' manifest.conf ;
	cp $INPUTFASTQ input/input.fastq ;
elif [ $TECHNOLOGY = 'Illumina' ] ; then
	if [ $READS = 'Single-Read' ] ; then
		cp lib/illumina_single_manifest.conf manifest.conf ;
		sed -i '5r temp/project_name_edited' manifest.conf ;
		cp $INPUTFASTQ input/input.fastq ;
	else
		cp lib/illumina_paired_manifest.conf manifest.conf ;
		sed -i '5r temp/project_name_edited' manifest.conf ;
		cp temp/fastq_data temp/paired_data ;
		sed -i 's/|/\n/' temp/paired_data ;
		split -l 1 temp/paired_data temp/paired ;
		PAIRED1=$(cat temp/pairedaa) ;
		PAIRED2=$(cat temp/pairedab) ;
		cp $PAIRED1 input/input1.fastq ;
		cp $PAIRED2 input/input2.fastq ;
		zenity --question --title="Reads Length" --ok-label="Length is fine" --cancel-label="Edit Length" --text="Illumina reads are set as 250 bases.\nEdit value if necessary."
		if [ $? = 0 ] ; then
			echo manifest.conf ready
		else 
			gedit manifest.conf ;
			if zenity --info --title="MEGAnnotator" --ok-label="Start MEGAnnotator" --text="Start MEGAnnotator when you are ready."
				then echo Starting MEGAnnotator
				else zenity --warning --title="MEGAnnotator stopped" --text="MEGAnnotator has been interrupted.\nPlease restart the script and follow the instructions."; 
					rm -R input ; rm -R temp ; rm mira ; rm manifest.conf ; exit
			fi
		fi
	fi
else
	cp lib/other_manifest.conf manifest.conf ;
	sed -i '5r temp/project_name_edited' manifest.conf ;
	cp $INPUTFASTQ input/input.fastq ;
	if zenity --info --title="MEGAnnotator" --ok-label="Open manifest.conf" --text="Please edit the manifest.conf file.\nThen come back to MEGAnnotator."
		then echo Editing manifest.conf
		else zenity --warning --title="MEGAnnotator stopped" --text="MEGAnnotator has been interrupted.\nPlease restart the script and follow the instructions."; 
			rm -R input ; rm -R temp ; rm mira ; rm manifest.conf ; exit
	fi
	gedit manifest.conf ;
	if zenity --info --title="MEGAnnotator" --ok-label="Start MEGAnnotator" --text="Start MEGAnnotator when you are ready."
		then echo Starting MEGAnnotator
		else zenity --warning --title="MEGAnnotator stopped" --text="MEGAnnotator has been interrupted.\nPlease restart the script and follow the instructions."; 
			rm -R input ; rm -R temp ; rm mira ; rm manifest.conf ; exit
	fi
fi

#=====================================================================================================================================================================================#
echo "#§#§#	Phase 1: Genome Assembly	#§#§#" ;
./mira -t $NTHREADS manifest.conf > mira.log ;

cd ${PROJECT}_assembly ;
rm -R ${PROJECT}_d_tmp ; rm -R ${PROJECT}_d_chkpt ;
cd .. ;
rm manifest.conf ; rm mira ;

#=====================================================================================================================================================================================#
echo "#§#§#     Phase 2: Contigs selection     #§#§#" ;
mkdir results ;
grep -v "IUPAC" ${PROJECT}_assembly/${PROJECT}_d_info/${PROJECT}_info_contigstats.txt > temp/contigs.txt ;
cut -f1 temp/contigs.txt > temp/contigs2.txt ;

#average calculation
grep "total coverage" ${PROJECT}_assembly/${PROJECT}_d_info/${PROJECT}_info_assembly.txt | sed 's/\ \ Avg\.\ total coverage\:\ //' > temp/average.txt ;
AVERAGE=$(cat temp/average.txt) ;
sed -i 's/\./\t/' temp/average.txt ;
cut -f1 temp/average.txt > temp/average2.txt ;
AVERAGE2=$(cat temp/average2.txt) ;
AVE=$((AVERAGE2/3)) ;

#miraconvert
bin/./miraconvert -t fasta -n temp/contigs2.txt -x $CLENGTH -y $AVE -z $RCONTIG ${PROJECT}_assembly/${PROJECT}_d_results/${PROJECT}_out.maf temp/${PROJECT}_LargeContigs_out > temp/ec.log ;

#=====================================================================================================================================================================================#
echo "#§#§#     Phase 3: Alignment vs. reference genome     #§#§#" ;
cp -a bin/mauve_2.3.1/. temp/ ;
if [ -e input/reference.fasta ] ; then
	cd temp ;
	java -Xmx2g -cp Mauve.jar org.gel.mauve.contigs.ContigOrderer -output mauve -ref ../input/reference.fasta -draft ${PROJECT}_LargeContigs_out_AllStrains.unpadded.fasta ;
	cd mauve ;
	MAUVEFOLDER=$(find . -exec stat \{} --printf="%y\t%n\n" \; | sort -n -r | head -1 | awk '{print $4}' | sed 's/\.\///') ;
	cd $MAUVEFOLDER ;
	cp *.tab ../../unpadded_contigs.tab ;
	mkdir ../../../results/mauve_aligment ;
	cp * ../../../results/mauve_aligment ;
	cd .. ; cd .. ; cd .. ;
else
	echo No reference identified - Proceeding with improvement of quality output
fi

#reordering contigs
if [ -e input/reference.fasta ] ; then
	grep -A 1000 "Ordered Contigs" temp/unpadded_contigs.tab > temp/ordered_tab ;
	sed -i '/^\s*$/d' temp/ordered_tab ;
	grep -B 1000 "Contigs with conflicting" temp/ordered_tab > temp/ordered_tab2 ;

	if [ -s temp/ordered_tab2 ] ; then
		grep -v "Contigs with conflicting" temp/ordered_tab2 > temp/ordered_tab_ok
	else
		cp temp/ordered_tab temp/ordered_tab_ok
	fi

	grep -v "Ordered Contigs" temp/ordered_tab_ok > temp/ordered_tab_ok2 ;
	grep -v "type" temp/ordered_tab_ok2 > temp/ordered_tab_ok3 ;
	cat temp/ordered_tab_ok3 | gawk '{print $2,"\t",$4}' > temp/ordered_tab_ok4 ;
	cat temp/ordered_tab_ok3 | gawk '{print $2}' > temp/ordered_tab_names ;
	grep "complement" temp/ordered_tab_ok4 > temp/ordered_complement ;
	grep "forward" temp/ordered_tab_ok4 > temp/ordered_forward ;
	cat temp/ordered_complement | gawk '{print $1}' > temp/ordered_complement_names ;
	cat temp/ordered_forward | gawk '{print $1}' > temp/ordered_forward_names ;
	sed -i 's/$/\ /' temp/ordered_complement_names ;
	sed -i 's/$/\ /' temp/ordered_forward_names ;
	sed -i 's/$/\ /' temp/ordered_tab_names ;

	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < temp/${PROJECT}_LargeContigs_out_AllStrains.unpadded.fasta > temp/${PROJECT}_contigs_1line.fasta ;
	sed -i 's/$/\ /' temp/${PROJECT}_contigs_1line.fasta ;

	while read COMPLEMENT
	do
	grep -w -A 1 $COMPLEMENT temp/${PROJECT}_contigs_1line.fasta >> temp/complement_contigs.fasta
	done < temp/ordered_complement_names

	while read FORWARD
	do
	grep -w -A 1 $FORWARD temp/${PROJECT}_contigs_1line.fasta >> temp/forward_contigs.fasta
	done < temp/ordered_forward_names

	sed -i 's/\ $//' temp/complement_contigs.fasta ;

	readseq -r -a -format=8 temp/complement_contigs.fasta > temp/complement_contigs_rotated.fasta ;
	sed -i 's/,/\t/' temp/complement_contigs_rotated.fasta ;
	cat temp/complement_contigs_rotated.fasta | gawk '{print $1}' > temp/complement_contigs_rotated2.fasta ;
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < temp/complement_contigs_rotated2.fasta > temp/complement_contigs_rotated3.fasta ;
	sed -i 's/$/\ /' temp/complement_contigs_rotated3.fasta ;
	cat temp/forward_contigs.fasta temp/complement_contigs_rotated3.fasta > temp/tot_contigs.fasta ;

	while read TOTAL
	do
	grep -w -A 1 $TOTAL temp/tot_contigs.fasta >> temp/${PROJECT}_multifasta.fasta
	done < temp/ordered_tab_names

	sed -i 's/\ $//' temp/${PROJECT}_multifasta.fasta ;
fi

if [ -f temp/unpadded_contigs.tab ] ; then
	cp temp/${PROJECT}_multifasta.fasta results/${PROJECT}_multifasta.fasta
else
	cp temp/${PROJECT}_LargeContigs_out_AllStrains.unpadded.fasta results/${PROJECT}_multifasta.fasta
fi

#=====================================================================================================================================================================================#
echo "#§#§#     Phase 4: Improvement of quality output     #§#§#" ;
mkdir temp/snp ;
mkdir results/improvement_quality_results ;

#Reads alignments
cp results/${PROJECT}_multifasta.fasta temp/snp/templato.fasta ;
bwa index temp/snp/templato.fasta ;

if [ -e input/input1.fastq ] ; then
	bwa mem -t $NTHREADS temp/snp/templato.fasta input/input1.fastq input/input2.fastq > temp/snp/${PROJECT}_aligned.sam ; sleep 1s ;
	samtools view -@ $NTHREADS -S -h -b temp/snp/${PROJECT}_aligned.sam > temp/snp/${PROJECT}.bam ; sleep 1s ;
	samtools sort -@ $NTHREADS temp/snp/${PROJECT}.bam  temp/snp/sorted${PROJECT} ; sleep 1s ;
	samtools index temp/snp/sorted${PROJECT}.bam ; sleep 1s ;
	samtools mpileup -f temp/snp/templato.fasta temp/snp/sorted${PROJECT}.bam > temp/snp/sorted${PROJECT}mpileup.bam ;
else
	bwa mem -t $NTHREADS temp/snp/templato.fasta input/input.fastq > temp/snp/${PROJECT}_aligned.sam ; sleep 1s ;
	samtools view -@ $NTHREADS -S -h -b temp/snp/${PROJECT}_aligned.sam > temp/snp/${PROJECT}.bam ; sleep 1s ;
	samtools sort -@ $NTHREADS temp/snp/${PROJECT}.bam  temp/snp/sorted${PROJECT} ; sleep 1s ;
	samtools index temp/snp/sorted${PROJECT}.bam ; sleep 1s ;
	samtools mpileup -f temp/snp/templato.fasta temp/snp/sorted${PROJECT}.bam > temp/snp/sorted${PROJECT}mpileup.bam ;
fi

#SNPs & INDELs with VarScan
cp -a bin/VarScan.v2.3.6.jar temp/ ;
cd temp ;
java -jar VarScan.v2.3.6.jar mpileup2cns snp/sorted${PROJECT}mpileup.bam --min-coverage 8 --min-reads2 2 --min-avg-qual 15 --min-var-freq 0.89 -p-value 99e-02 --variants > snp/${PROJECT}_0.89_validation.tab ;
java -jar VarScan.v2.3.6.jar mpileup2cns snp/sorted${PROJECT}mpileup.bam --min-coverage 8 --min-reads2 2 --min-avg-qual 15 --min-var-freq 0.89 -p-value 99e-02 --output-vcf 1 --variants > snp/${PROJECT}_0.89_validation.vcf ;
cd .. ;
mv temp/snp/${PROJECT}_0.89_validation.tab results/improvement_quality_results/${PROJECT}_0.89_validation.tab ;
mv temp/snp/${PROJECT}_0.89_validation.vcf results/improvement_quality_results/${PROJECT}_0.89_validation.vcf ;
bgzip results/improvement_quality_results/${PROJECT}_0.89_validation.vcf -c > results/improvement_quality_results/${PROJECT}_0.89_validation.vcf.gz ;
tabix -p vcf results/improvement_quality_results/${PROJECT}_0.89_validation.vcf.gz ;

#vcf editing ;
cp results/improvement_quality_results/${PROJECT}_0.89_validation.vcf results/improvement_quality_results/${PROJECT}_0.89_validation_noIUPC.vcf ;
sed -i 's/\tR\t/\tN\t/' results/improvement_quality_results/${PROJECT}_0.89_validation_noIUPC.vcf ;
sed -i 's/\tY\t/\tN\t/' results/improvement_quality_results/${PROJECT}_0.89_validation_noIUPC.vcf ;
sed -i 's/\tS\t/\tN\t/' results/improvement_quality_results/${PROJECT}_0.89_validation_noIUPC.vcf ;
sed -i 's/\tW\t/\tN\t/' results/improvement_quality_results/${PROJECT}_0.89_validation_noIUPC.vcf ;
sed -i 's/\tK\t/\tN\t/' results/improvement_quality_results/${PROJECT}_0.89_validation_noIUPC.vcf ;
sed -i 's/\tM\t/\tN\t/' results/improvement_quality_results/${PROJECT}_0.89_validation_noIUPC.vcf ;
sed -i 's/\tB\t/\tN\t/' results/improvement_quality_results/${PROJECT}_0.89_validation_noIUPC.vcf ;
sed -i 's/\tD\t/\tN\t/' results/improvement_quality_results/${PROJECT}_0.89_validation_noIUPC.vcf ;
sed -i 's/\tH\t/\tN\t/' results/improvement_quality_results/${PROJECT}_0.89_validation_noIUPC.vcf ;
sed -i 's/\tV\t/\tN\t/' results/improvement_quality_results/${PROJECT}_0.89_validation_noIUPC.vcf ;

#GATK improved generation
cp -a bin/CreateSequenceDictionary.jar temp/ ;
cp -a bin/GenomeAnalysisTK.jar temp/ ;
cd temp ;
java -jar CreateSequenceDictionary.jar R=snp/templato.fasta O=snp/templato.dict ;
java -Xmx2g -jar GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R snp/templato.fasta -o ../results/improvement_quality_results/${PROJECT}_improved.fasta --variant ../results/improvement_quality_results/${PROJECT}_0.89_validation_noIUPC.vcf ;
cd .. ;
sed -i 's/>/>Contig_/' results/improvement_quality_results/${PROJECT}_improved.fasta ;
sed -i 's/>Contig_1$/>Contig_01/' results/improvement_quality_results/${PROJECT}_improved.fasta ;
sed -i 's/>Contig_2$/>Contig_02/' results/improvement_quality_results/${PROJECT}_improved.fasta ;
sed -i 's/>Contig_3$/>Contig_03/' results/improvement_quality_results/${PROJECT}_improved.fasta ;
sed -i 's/>Contig_4$/>Contig_04/' results/improvement_quality_results/${PROJECT}_improved.fasta ;
sed -i 's/>Contig_5$/>Contig_05/' results/improvement_quality_results/${PROJECT}_improved.fasta ;
sed -i 's/>Contig_6$/>Contig_06/' results/improvement_quality_results/${PROJECT}_improved.fasta ;
sed -i 's/>Contig_7$/>Contig_07/' results/improvement_quality_results/${PROJECT}_improved.fasta ;
sed -i 's/>Contig_8$/>Contig_08/' results/improvement_quality_results/${PROJECT}_improved.fasta ;
sed -i 's/>Contig_9$/>Contig_09/' results/improvement_quality_results/${PROJECT}_improved.fasta ;
cp results/improvement_quality_results/${PROJECT}_improved.fasta results/${PROJECT}_improved.fasta ;

#contigs reference file
grep ">" results/${PROJECT}_multifasta.fasta > temp/original_contigs_grep ;
sed -i 's/>//' temp/original_contigs_grep ;
grep ">" results/${PROJECT}_improved.fasta > temp/new_contigs_grep ;
sed -i 's/>//' temp/new_contigs_grep ;
pr -mts temp/new_contigs_grep temp/original_contigs_grep > results/contigs_reference_tab ;

#=====================================================================================================================================================================================#
echo "#§#§#     Phase 5: ORFs prediction     #§#§#" ;
bin/./prodigal.linux -f gff -a temp/aaORFs.fasta -i results/${PROJECT}_improved.fasta -o temp/ORFs.gff ;

#=====================================================================================================================================================================================#
echo "#§#§#     Phase 6: RapSearch2 vs nrDB     #§#§#" ;
bin/./rapsearch -q temp/aaORFs.fasta -d $RAPSEARCHDB -o temp/file.rapsearch -s f -e 0.0001 -l 20 -z $NTHREADS -b 10 -v 10  ;

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
grep -i -v "similar" temp/rapsearch_tbl44 > temp/rapsearch_tbl45 ;
grep -i -v "fragment" temp/rapsearch_tbl45 > temp/rapsearch_tbl46 ;
grep -i -v "TIGR..... family protein" temp/rapsearch_tbl46 > temp/file_greppato_rapsearch_temp ;

#takes the first hit of rapsearch
grep "| " temp/file_greppato_rapsearch_temp > temp/file_greppato_rapsearch ;
cat temp/locus_tag | xargs -I{} grep -m 1 -w {} temp/file_greppato_rapsearch > temp/first_hit_rapsearch_temp.txt ;
sed -i 's/\t/\ \t/g' temp/first_hit_rapsearch_temp.txt ;
sed -i 's/|\ /\t/' temp/first_hit_rapsearch_temp.txt ;
gawk -F"\t" '{print $1,"\t",$3}' temp/first_hit_rapsearch_temp.txt > temp/first_hit_rapsearch.txt ;
sed -i 's/\ \t\ /\t/g' temp/first_hit_rapsearch.txt ;

#=====================================================================================================================================================================================#
echo "#§#§#     Phase 7: pfam prediction     #§#§#" ;
hmmscan --cpu $NTHREADS -E 1e-10 --tblout temp/pfam_annotation_tbl -o temp/pfam_annotation $PFAMDB temp/aaORFs.fasta ;

#output modding
grep -v "^#" temp/pfam_annotation_tbl > temp/pfam_annotation_tbl2 ;
sed -i 's/\ \ \ /\ /g' temp/pfam_annotation_tbl2 ;
sed -i 's/\ \ /\ /g' temp/pfam_annotation_tbl2 ;
sed -i 's/\ \ \ /\ /g' temp/pfam_annotation_tbl2 ;
sed -i 's/\ \ /\ /g' temp/pfam_annotation_tbl2 ;
sed -i 's/\ \ \ /\ /g' temp/pfam_annotation_tbl2 ;
sed -i 's/\ \ /\ /g' temp/pfam_annotation_tbl2 ;
sed -i 's/\ \ \ /\ /g' temp/pfam_annotation_tbl2 ;
sed -i 's/\ \ /\ /g' temp/pfam_annotation_tbl2 ;
sed -i 's/\ \ \ /\ /g' temp/pfam_annotation_tbl2 ;
sed -i 's/\ \ /\ /g' temp/pfam_annotation_tbl2 ;
sed -i 's/\ /\t/g' temp/pfam_annotation_tbl2 ;
cut -f 3,19-50 temp/pfam_annotation_tbl2 > temp/pfam_annotation_tbl3 ;
sed -i 's/\t/\ /g' temp/pfam_annotation_tbl3 ;
sed -i 's/\ /\t/' temp/pfam_annotation_tbl3 ;
sed -i 's/domain$/domain-containing protein/' temp/pfam_annotation_tbl3 ;
grep -v "DUF" temp/pfam_annotation_tbl3 > temp/pfam_annotation_tbl4 ;
grep -v "UPF" temp/pfam_annotation_tbl4 > temp/pfam_annotation_tbl5 ;
grep -i -v "protein of unknown function" temp/pfam_annotation_tbl5 > temp/pfam_annotation_tbl6 ;
grep -i -v "domain of unknown function" temp/pfam_annotation_tbl6 > temp/pfam_annotation_tbl7 ;
grep -i -v "small basic protein" temp/pfam_annotation_tbl7 > temp/pfam_annotation_tbl8 ;
grep -i -v "conserved hypothetical protein" temp/pfam_annotation_tbl8 > temp/pfam_annotation_tbl9 ;
grep -i -v "repeated domains containing protein" temp/pfam_annotation_tbl9 > temp/pfam_annotation_tbl10 ;
grep -i -v "repeat domain-containing protein" temp/pfam_annotation_tbl10 > temp/pfam_annotation_tbl11 ;
grep -i -v "Uncharacterised" temp/pfam_annotation_tbl11 > temp/pfam_annotation_tbl12 ;
grep -i -v "Uncharacterized" temp/pfam_annotation_tbl12 > temp/file_greppato_pfam ;

#locus_tag definition
gawk '{print $1}' temp/file_greppato_pfam > temp/first_column_pfam ;
gawk '!x[$0]++' temp/first_column_pfam > temp/locus_tag_pfam ;

#takes the first hit of pfam
cat temp/locus_tag_pfam | xargs -I{} grep -m 1 -w {} temp/file_greppato_pfam > temp/first_hit_pfam.txt ;
sed -i 's/\t/\ \t/g' temp/first_hit_pfam.txt ;

#all the locus tag generation
grep "^>" temp/aaORFs.fasta > temp/all_the_locus_tag ;
sed -i 's/\#/\t/' temp/all_the_locus_tag ;
gawk -F"\t" '{print $1}' temp/all_the_locus_tag > temp/locus_tag_TOT ;
sed -i 's/>//' temp/locus_tag_TOT ;

#=====================================================================================================================================================================================#
echo "#§#§#     Phase 8: merging annotation     #§#§#" ;
while read LTAG
do
let count=$count+1
grep -w $LTAG temp/first_hit_rapsearch.txt > temp/file${count}_NCBI.single
grep -w $LTAG temp/first_hit_pfam.txt  > temp/file${count}_PFAM.single
if [ -s temp/file${count}_NCBI.single ] ; then
	grep -m 1 "^" temp/file${count}_NCBI.single >> temp/merged_annotation
elif [ -s temp/file${count}_PFAM.single ] ; then
	grep -m 1 "^" temp/file${count}_PFAM.single >> temp/merged_annotation
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

#=====================================================================================================================================================================================#
echo "#§#§#     Phase 9: gbk generation and rRNA prediction     #§#§#" ;
mkdir temp/single_gffs ;
mkdir temp/single_contigs ;
mkdir temp/single_gbks ;
mkdir temp/single_contigs_length ;
mkdir temp/single_rRNA ;

#genbank creation with rRNA prediction
grep "^>" results/${PROJECT}_improved.fasta > temp/contigs_name ;
sed -i 's/>//' temp/contigs_name ;
sed "s/\tID=/\tlocus_tag=${PROJECT}_/" results/${PROJECT}_annotated_ORFs.gff > temp/annotated_ORFs_locus_tag.gff ;
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < results/${PROJECT}_improved.fasta > temp/improved_1line.fasta ;

while read LTAG3
do
let count3=$count3+1
	grep -P "$LTAG3\t" temp/annotated_ORFs_locus_tag.gff > temp/single_gffs/contig_${count3}.gff
	grep -A 1 "$LTAG3$" temp/improved_1line.fasta > temp/single_contigs/contig_${count3}.fasta
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
readseq -a -f=2 results/${PROJECT}_improved.fasta > temp/multigbk_fromfasta.gbk ;
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
echo "#§#§#     Phase 10: genbank finalization and tRNA prediction     #§#§#" ;
union -sequence results/${PROJECT}_improved.fasta -sformat fasta -snucleotide -outseq temp/allnt.fasta -osformat fasta -auto ;
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

#cleaning
rm -R input ;
rm -R temp ;
rm results/${PROJECT}_annotated_ORFs.gff ;
rm results/${PROJECT}_improved.fasta ;
rm results/${PROJECT}_multifasta.fasta ;
rm results/${PROJECT}_tRNA ;
rm results/improvement_quality_results/${PROJECT}_0.89_validation_noIUPC.vcf ;
rm results/improvement_quality_results/${PROJECT}_0.89_validation_noIUPC.vcf.idx ;
mv results/contigs_reference_tab results/improvement_quality_results/. ;
mv mira.log ${PROJECT}_assembly/. ;
mv results ${PROJECT}_results ;

echo "#§#§#     Prediction Complete     #§#§#" ;

