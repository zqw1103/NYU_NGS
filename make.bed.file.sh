grep "chr1\s" genome.file.txt > chr1.genome.file.txt

bedtools makewindows -g chr1.genome.file.txt -w 500 > chr1_500bp.bed
