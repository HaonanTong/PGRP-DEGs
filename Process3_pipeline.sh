function f_compare_two_gene_list {
	############################################
	# Analysis of two gene list
	# Nan
	# @PGRP
	############################################

		str1=`echo $1 | awk -F. '{print $1}'`;
		str2=`echo $2 | awk -F. '{print $1}'`;
		str="$str1-$str2"

		mkdir "Result_of_Process3/$str"

		cat $1 | sort -u > Result_of_Process3/$str/file1_sorted_tmp;
		cat $2 | sort -u > Result_of_Process3/$str/file2_sorted_tmp;
     # # of genes in list1
		awk 'END{print NR}' Result_of_Process3/$str/file1_sorted_tmp > "Result_of_Process3/$str/Num-$str1.txt";

     # # of genes in list2 
     	awk 'END{print NR}' Result_of_Process3/$str/file2_sorted_tmp > "Result_of_Process3/$str/Num-$str2.txt";

     # genes Overlapping
     	grep -f Result_of_Process3/$str/file1_sorted_tmp Result_of_Process3/$str/file2_sorted_tmp | sort -u > "Result_of_Process3/$str/Genes-Overlapping-$str.txt";
		awk 'END{print NR}' "Result_of_Process3/$str/Genes-Overlapping-$str.txt" > "Result_of_Process3/$str/Num-Genes-Overlapping-$str.txt";

     # in list1 not list 2
     	# diff -f "Genes-Overlapping-$str.txt" file1_sorted_tmp > tmp
     	# awk '{if ($2 ~/AT/) print $2}' tmp |sort -u > "Genes-$str1-Not-$str2.txt"

     	grep -v -f "Result_of_Process3/$str/Genes-Overlapping-$str.txt" Result_of_Process3/$str/file1_sorted_tmp > "Result_of_Process3/$str/Genes-$str1-Not-$str2.txt"
		#rm tmp
		awk 'END{print NR}' "Result_of_Process3/$str/Genes-$str1-Not-$str2.txt" > "Result_of_Process3/$str/Num-Genes-$str1-Not-$str2.txt"


     # in list2 not list1
      #  diff -f "Genes-Overlapping-$str.txt" $2 > tmp
     #	awk '{if ($2 ~/AT/) print $2}' tmp |sort -u > "Genes-$str2-Not-$str1.txt"
	#	rm tmp
	#	awk 'END{print NR}' "Genes-$str2-Not-$str1.txt" > "Num-Genes-$str2-Not-$str1.txt"
     	grep -v -f "Result_of_Process3/$str/Genes-Overlapping-$str.txt" Result_of_Process3/$str/file2_sorted_tmp > "Result_of_Process3/$str/Genes-$str2-Not-$str1.txt"
		#rm tmp
		awk 'END{print NR}' "Result_of_Process3/$str/Genes-$str2-Not-$str1.txt" > "Result_of_Process3/$str/Num-Genes-$str2-Not-$str1.txt"

		rm Result_of_Process3/$str/file1_sorted_tmp
		rm Result_of_Process3/$str/file2_sorted_tmp
}

mkdir 'Result_of_Process3'

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
head -1 kat-rpkm-expression.csv > Result_of_Process3/Profiles-ANan-DEGs.csv
grep -f ANan-Gene-list-DEG.txt kat-rpkm-expression.csv >> Result_of_Process3/Profiles-ANan-DEGs.csv

# extract profiles of up and down regulated genes.
head -1 kat-rpkm-expression.csv > Result_of_Process3/Profiles-ANan-up-regulated.csv
grep -f ANan-Gene-list-DEG-up-regulated.txt kat-rpkm-expression.csv >> Result_of_Process3/Profiles-ANan-up-regulated.csv
head -1 kat-rpkm-expression.csv > Profiles-ANan-down-regulated.csv
grep -f ANan-Gene-list-DEG-down-regulated.txt kat-rpkm-expression.csv >> Result_of_Process3/Profiles-ANan-down-regulated.csv

# analysis on genes shown in ANan up regulated but not all-up.txt.
# analysis on genes shown in ANan down regulated but not all-down.txt.
head -1 kat-rpkm-expression.csv > Result_of_Process3/Profies-Genes-ANan-Gene-list-DEG-down-regulated-Not-all-down.csv
grep -f Result_of_Process3/ANan-Gene-list-DEG-down-regulated-all-down/Genes-ANan-Gene-list-DEG-down-regulated-Not-all-down.txt kat-rpkm-expression.csv >> Result_of_Process3/Profiles-Genes-ANan-Gene-list-DEG-down-regulated-Not-all-down.csv
head -1 kat-rpkm-expression.csv > Result_of_Process3/Profiles-Genes-ANan-Gene-list-DEG-up-regulated-Not-all-up.csv
grep -f Result_of_Process3/ANan-Gene-list-DEG-up-regulated-all-up/Genes-ANan-Gene-list-DEG-up-regulated-Not-all-up.txt kat-rpkm-expression.csv >> Result_of_Process3/Profiles-Genes-ANan-Gene-list-DEG-up-regulated-Not-all-up.csv




