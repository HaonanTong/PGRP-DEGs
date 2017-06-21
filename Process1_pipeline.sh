cd ~/Google\ Drive/PGRP/Nan/Data/20170620g
# convert white space into comma
# cat kat.all.rpkm > kat-rpkm-expression
sed -E 's/[[:space:]]+/,/g' kat.all.rpkm > kat-rpkm-expression.csv 
