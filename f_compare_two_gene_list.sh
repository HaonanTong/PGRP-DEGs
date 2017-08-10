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

# cd ~/Google\ Drive/PGRP/Nan/Data/20170526/
# f_compare_two_gene_list "Gene-list-fdr2.txt" "EIN3-R.txt"
