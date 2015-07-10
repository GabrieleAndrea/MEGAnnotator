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

(
echo "5"
echo "# Phase0: Allocating Variables"
#========================================================================================================================================================
while [ ! -f "input/input.fasta" ]
do
	sleep 1
done
echo "10"
echo "# Phase1: Genes Prediction"
#========================================================================================================================================================
while [ ! -f "temp/file.rapsearch.tmp0" ]
do
	sleep 1
done
echo "20"
echo "# Phase2: Genes Annotation"
#========================================================================================================================================================
while [ ! -f "temp/first_hit_rapsearch.txt" ]
do
	sleep 1
done
echo "50"
echo "# Phase3: Motives Prediction"
#========================================================================================================================================================
while [ ! -f "temp/locus_tag_TOT" ]
do
	sleep 1
done
echo "75"
echo "# Phase4: Merging Annotation"
#========================================================================================================================================================
while [ ! -f "temp/contigs_name" ]
do
	sleep 1
done
echo "80"
echo "# Phase5: Genbank generation and rRNA prediction"
#========================================================================================================================================================
while [ ! -f "temp/allnt.fasta" ]
do
	sleep 1
done
echo "90"
echo "# Phase6: Genbank finalization and tRNA prediction"
#========================================================================================================================================================
while [ ! -f "temp/tail_union.gbk" ]
do
	sleep 1
done
echo "99"
echo "# MEGAnnotator ends" ; sleep 1
echo "100"
) | zenity --progress --auto-close --no-cancel --title="MEGAnnotator" --text="Setting tings up..." --percentage=0
if [ "$?" = -1 ] ; then
	zenity --error --text="MEGAnnotator ends improperly.\nPlease verify the data quality."
fi
zenity --info --title="MEGAnnotator" --text="Annotation completed"
