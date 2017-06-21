# Analysis on derived DEGs v.s. EIN3-R
awk -F  "." '{print $1}' ANan-Gene-list-DEG.txt \
	| sort -u > ANan-DEGs.txt
f_compare_two_gene_list "ANan-DEGs.txt" "EIN3-R.txt"

# Analysis on derived DEGs v.s. DEGs from Chang's paper
cat all.down all.up | sort -u> kat.txt
f_compare_two_gene_list "ANan-Gene-list-DEG.txt" "kat.txt"