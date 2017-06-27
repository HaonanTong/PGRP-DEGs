echo '=============================================================='
echo 'Process1: Generating .csv file from txt file in white space...'
bash Process1_pipeline.sh 
echo 'End Process1'
echo '==============================================================/n/n'

echo '=============================================================='
echo 'Process2: deriving DEGs and analyzing...'
export PATH=$PATH:/Applications/MATLAB_R2017a.app/bin/
matlab  -r -nodesktop Process2_main
echo 'End Process2'
echo '==============================================================/n/n'

echo '=============================================================='
echo 'Process3: Comparing DEGs deriving from Nan v.s. Chang's
bash Process3_pipeline.sh 
echo 'End Process3'
echo '==============================================================/n/n'

echo '=============================================================='
echo 'Process4: Ploting and Concluding...'
matlab -r -nodesktop Process4_main
echo 'End Process4'
echo '==============================================================/n/n'