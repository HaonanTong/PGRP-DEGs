mkdir Figures
# convert white space into comma
sed -E 's/[[:space:]]+/,/g' kat.all.rpkm > kat-rpkm-expression.csv 
