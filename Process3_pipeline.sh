function f_compare_two_gene_list {
	############################################
	# Analysis of two gene list
	# Nan
	# @PGRP
	############################################

		str1=`echo $1 | awk -F. '{print $1}'`;
		str2=`echo $2 | awk -F. '{print $1}'`;
		str="$str1-$str2"

		cat $1 | sort -u > file1_sorted_tmp;
		cat $2 | sort -u > file2_sorted_tmp;
     # # of genes in list1
		awk 'END{print NR}' file1_sorted_tmp > "Num-$str1.txt";

     # # of genes in list2 
     	awk 'END{print NR}' file2_sorted_tmp > "Num-$str2.txt";

     # genes Overlapping
     	grep -f file1_sorted_tmp file2_sorted_tmp | sort -u > "Genes-Overlapping-$str.txt";
		awk 'END{print NR}' "Genes-Overlapping-$str.txt" > "Num-Genes-Overlapping-$str.txt";

     # in list1 not list 2
     	# diff -f "Genes-Overlapping-$str.txt" file1_sorted_tmp > tmp
     	# awk '{if ($2 ~/AT/) print $2}' tmp |sort -u > "Genes-$str1-Not-$str2.txt"

     	grep -v -f "Genes-Overlapping-$str.txt" file1_sorted_tmp > "Genes-$str1-Not-$str2.txt"
		#rm tmp
		awk 'END{print NR}' "Genes-$str1-Not-$str2.txt" > "Num-Genes-$str1-Not-$str2.txt"


     # in list2 not list1
      #  diff -f "Genes-Overlapping-$str.txt" $2 > tmp
     #	awk '{if ($2 ~/AT/) print $2}' tmp |sort -u > "Genes-$str2-Not-$str1.txt"
	#	rm tmp
	#	awk 'END{print NR}' "Genes-$str2-Not-$str1.txt" > "Num-Genes-$str2-Not-$str1.txt"
     	grep -v -f "Genes-Overlapping-$str.txt" file2_sorted_tmp > "Genes-$str2-Not-$str1.txt"
		#rm tmp
		awk 'END{print NR}' "Genes-$str2-Not-$str1.txt" > "Num-Genes-$str2-Not-$str1.txt"

}

# Analysis on derived DEGs v.s. EIN3-R
awk -F  "." '{print $1}' ANan-Gene-list-DEG.txt \
	| sort -u > ANan-DEGs.txt
f_compare_two_gene_list "ANan-DEGs.txt" "EIN3-R.txt"

# Analysis on derived DEGs v.s. DEGs from Chang's paper
cat all.down all.up | sort -u> kat.txt
f_compare_two_gene_list "ANan-Gene-list-DEG.txt" "kat.txt"

# Analysis on derived up-regulated DEGs v.s. DEGs up-regulated from Chang's paper
cat all.up > all-up.txt
f_compare_two_gene_list "ANan-Gene-list-DEG-up-regulated.txt" "all-up.txt"
rm all-up.txt

# Analysis on derived down-regulated DEGs v.s. DEGs down-regulated from Chang's paper
cat all.down > all-down.txt
f_compare_two_gene_list "ANan-Gene-list-DEG-down-regulated.txt" "all-down.txt"
rm all-down.txt

# extract profiles of DEGs.
grep -f ANan-Gene-list-DEG.txt kat-rpkm-expression.csv > Profiles-ANan-DEGs.csv

# extract profiles of up and down regulated genes.
grep -f ANan-Gene-list-DEG-up-regulated.txt kat-rpkm-expression.csv > Profiles-ANan-up-regulated.csv
grep -f ANan-Gene-list-DEG-down-regulated.txt kat-rpkm-expression.csv > Profiles-ANan-down-regulated.csv

# analysis on genes shown in ANan up regulated but not all-up.txt.
# analysis on genes shown in ANan down regulated but not all-down.txt.
grep -f Genes-ANan-Gene-list-DEG-down-regulated-Not-all-down.txt kat-rpkm-expression.csv > Profies-Genes-ANan-Gene-list-DEG-down-regulated-Not-all-down.csv
grep -f Genes-ANan-Gene-list-DEG-up-regulated-Not-all-up.txt kat-rpkm-expression.csv > Profiles-Genes-ANan-Gene-list-DEG-up-regulated-Not-all-up.csv




