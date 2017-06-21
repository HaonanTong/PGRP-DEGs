%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
% @CopyRight Haonan Tong
% PGRP
% Verify Gene Set EIN3-R
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
T = readtable('kat-rpkm-expression.csv',...
     'ReadVariableNames',true);
% summary(T);
Data = table2array(T(:,2:end));
name = table2array(T(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T-test Table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T0 = Data(:,1:3)';
T1 = Data(:,4:6)';
T2 = Data(:,7:9)';
T3 = Data(:,10:12)';
T4 = Data(:,13:15)';
T5 = Data(:,16:18)';
T6 = Data(:,19:21)';

sig = 0.05;

[h1,p1] = ttest2(T0,T1,sig);
[h2,p2] = ttest2(T0,T2,sig);
[h3,p3] = ttest2(T0,T3,sig);
[h4,p4] = ttest2(T0,T4,sig);
[h5,p5] = ttest2(T0,T5,sig);
[h6,p6] = ttest2(T0,T6,sig);

htable = [h1 ; h2 ; h3 ; h4 ; h5 ; h6 ]';
htable(isnan(htable))=0;
ptable = [p1 ; p2 ; p3 ; p4 ; p5 ; p6 ]';

% T.Properties.VariableNames
p_table = array2table(ptable,'VariableNames',{'p1','p2','p3','p4','p5','p6'},'RowNames',...
    name);
writetable( p_table, 'TABLE_ptable.csv','WriteRowNames',true);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derive mean of Expression data set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
newData = [];
for i = 1:3:21%7 time points; 3 replicates;
   newData = [newData sum(Data(:,i:i+2),2)];
end
newData = 1/3*newData;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rpkm > 1 as a condition of
% significantly expressed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SEMatrix = newData > 1;
[m,n] = size(SEMatrix);
tmp = [];
for i = 2 : n
    tmp = [ tmp SEMatrix(:,i) | SEMatrix(:,1) ];
end
SEMatrix = tmp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the 50% difference as NC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DEMatrix = log2(newData./repmat(newData(:,1),1,size(newData,2)));
DEMatrix_abs = abs(log2(newData./repmat(newData(:,1),1,size(newData,2))));
DEMatrix_abs = DEMatrix_abs(:,2:end);
DEtable = DEMatrix_abs >= 0.58;

DE_table = array2table(DEMatrix_abs,'VariableNames',{'r1','r2','r3','r4','r5','r6'},'RowNames',...
    name);
writetable( DE_table, 'TABLE_DE_table.csv','WriteRowNames',true);

%go up and go down
DE_up_table = DEMatrix(:,2:end) > 0;
DE_down_table = DEMatrix(:,2:end) < 0;

GroundtruthTable = DEtable .* htable .* SEMatrix;
GroundtruthTable_up = GroundtruthTable .* DE_up_table;
GroundtruthTable_down = GroundtruthTable .* DE_down_table;

indx = sum(GroundtruthTable , 2) > 0;% log2(1.5) = 0.58
indx_up = sum(GroundtruthTable_up,2)>0;
indx_down =sum(GroundtruthTable_down,2)>0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output gene list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEGs
fp = fopen('ANan-Gene-list-DEG.txt','wt');
fprintf(fp, '%s\n', name{indx,:});
fclose(fp);

% up-regulated
fp = fopen('ANan-Gene-list-DEG-up-regulated.txt','wt');
fprintf(fp, '%s\n', name{indx_up,:});
fclose(fp);

% down-regulated
fp = fopen('ANan-Gene-list-DEG-down-regulated.txt','wt');
fprintf(fp, '%s\n', name{indx_down,:});
fclose(fp);







