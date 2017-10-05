T = readtable('ANan-DEGs.txt','ReadRowNames',0,'ReadVariableNames',0);
T_ethylene = readtable('GO0009723.txt','ReadRowNames',0,'ReadVariableNames',1);
T_auxin = readtable('GO0009733.txt','ReadRowNames',0,'ReadVariableNames',1);
T_cellwall = readtable('GO0071555.txt','ReadRowNames',0,'ReadVariableNames',1);

ID = table2cell(T);
ng = length(ID);
GO0009723 = zeros(ng,1);%ethylene
GO0009733 = zeros(ng,1);%auxin
GO0071555 = zeros(ng,1);%cellwall
T_GO_ASSOC = table(GO0009723,GO0009733,GO0071555);
T_GO_ASSOC.Properties.RowNames=ID;

T_GO_ASSOC.GO0009723(table2array(T_ethylene),:)=1;%ethylene
T_GO_ASSOC.GO0009733(table2array(T_auxin),:)=1;%auxin
T_GO_ASSOC.GO0071555(table2array(T_cellwall),:)=1;%cellwall

writetable(T_GO_ASSOC,'TABLE_GO_ASSOC.csv','WriteRowNames',1,'WriteVariableNames',1);