function [ fig, ngene, expr, agis ] = f_plotTable( csv, isVariableNames )
% [ fig, ngene, expr, agis ] = plotTable( csv, isVariableNames )
% fig - figure of profiles in csv file
% ngene - # of genes in the file
% expr - profiles of genes in RPKM at each time point
% agis - gene list
% csv - exel file with 3 replicates for each time point, totally 6 time
% point
%       Example of .csv file opened in MS Excell might look as follows:
%       |    AIG    |  R1T1  |  R2T1  |  R3T1  | ... 
%       |AT1G01010.1| 1.2691 | 1.7789 | 3.0794 | ...
%       |    ...    |  ...   |  ...   |   ...  | ...
%       |AT5G67360.1| 150.93 | 243.19 |   ...  | ...
%       |    ...    |  ...   |  ...   |   ...  | ... 
% in this case isVariable set to be 1;
% isVariableName - 1 if ".csv" file contains Variable names. 0 otherwise.
%
if nargin ~= 2
    fprintf('\ntwo parameter required for the function\n');
    fprintf(' path of .csv file as the first parameter\n');
    fprintf(' if the file contains VariableNames as second parameter\n');
    fprintf(' example:  [ fig, ngene, expr, agis ] = plotTable(''kat-rpkm-expression.csv'', 1)\n\n');
    return;
end

%% Read File
if isVariableNames == 1
    T = readtable(csv,...
     'ReadVariableNames',true);
else
    T = readtable(csv,...
     'ReadVariableNames',false);
end
 % summary(T);
Data = table2array(T(:,2:end));
agis = table2array(T(:,1));
[ngene,~] = size(Data);
%% Plot Data Generation
expr = [];
for i = 1:3:21%7 time points; 3 replicates;
   expr = [expr sum(Data(:,i:i+2),2)];
end
expr = 1/3*expr;
tmp = [];
for i = 2:7
   tmp = [tmp log2( expr(:,i)./expr(:,1) )];
end
plotData = tmp;
plotData = [ zeros(size(plotData,1),1) plotData ];

%% Plot
fig = figure;x = 0 : 1 : 6;
hold on;axis([0 6 -2 5])
plot(x, plotData,'Color','[.4,.4,.4]');
plot(x,mean( plotData),'Color','r','LineWidth',4);
xticks(0:6)
xticklabels({'0','0.25','0.5','1','4','12','24'})
title(sprintf( 'plot of expression file\n "%s"', csv))
xlabel('Ethylene treatment(hrs)');
ylabel('Expression-log2ratio(reference at 0 hrs)');
set(gca,'fontsize',14);
end

