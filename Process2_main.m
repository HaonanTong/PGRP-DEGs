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

%%%%%%%%%%%%%%%%%%%%%%%
% T-test Table
%%%%%%%%%%%%%%%%%%%%%%%
% e.g.
% [h,p] = ttest([.001,.001,.001],[.019,.023,.0189]);
% h = 1 different mean value;
% p < 0.05;

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



%%%%%%%%%%%%%%%%%%%%%%%
% ratio-test Table
%%%%%%%%%%%%%%%%%%%%%%%
newData = [];
for i = 1:3:21%7 time points; 3 replicates;
   newData = [newData sum(Data(:,i:i+2),2)];
end
newData = 1/3*newData;
% Test using log2
% newData = log2(newData);
% CW
%newData = ( newData(:,2:end)+1 )./( repmat(newData(:,1),1,6) + 1 );
%newData = log2(newData);

%%%%%%%%%%%%%%%%%%%%%%%
% rpkm > 1 as a condition of
% significantly expressed
%%%%%%%%%%%%%%%%%%%%%%%
SEMatrix = newData > 1;
[m,n] = size(SEMatrix);
tmp = [];
for i = 2 : n
    tmp = [ tmp SEMatrix(:,i) | SEMatrix(:,1) ];
end
SEMatrix = tmp;

% the 50% difference
DEMatrix = abs(log2(newData./repmat(newData(:,1),1,size(newData,2))));
DEMatrix = DEMatrix(:,2:end);
DEtable = DEMatrix >= 0.58;

DE_table = array2table(DEMatrix,'VariableNames',{'r1','r2','r3','r4','r5','r6'},'RowNames',...
    name);
writetable( DE_table, 'TABLE_DE_table.csv','WriteRowNames',true);

GroundtruthTable = DEtable .* htable .* SEMatrix;


% [~, ntime] = size(newData2);
% DEMatrix2 = [];
% for i = 2 : ntime
%     DEMatrix2 = [DEMatrix2 abs(newData2(:,i)-newData2(:,i-1))./newData2(:,i-1) ];
% end
% indx2 = sum(DEMatrix2 > 0.5 , 2) > 0;
% writetable(cell2table(name(indx2)),'Genes-mRNA-ethylene-regulated2.txt','WriteVariableNames',false);
%%
%[~, ntime] = size(newData);
%DEMatrix = [];
%for i = 2: ntime
%    DEMatrix = [DEMatrix abs(newData(:,i)-newData(:,i-1))./newData(:,i-1) ];
%end
indx = sum(GroundtruthTable , 2) > 0;% log2(1.5) = 0.58

% Output names
%writetable(cell2table(name(indx)),'Genes-final-kat.txt','WriteVariableNames',false);


%#########
fp = fopen('ANan-Gene-list-DEG.txt','wt');
fprintf(fp, '%s\n', name{indx,:});
fclose(fp);

%% PLOT Chang
% plot profiles after t-test
[fig1, ngenes1] = plotTable('Profiles-Gene-list-t-test.csv');
title(sprintf('profiles after t-test, \n with %d genes',ngenes1),...
    'FontSize',14)
print(fig1,'./Figures/fig1','-dpng');

% plot profiles of DEGs after Chang's procedures
[fig2, ngenes2] = plotTable('Profiles-Genes-mRNA-ethylene-regulated.csv');
title(sprintf('profiles of DEGs after Chang procedures, \n with %d genes',ngenes2),...
    'FontSize',14)
print(fig2,'./Figures/fig2','-dpng');

% plot profiles of overlapping genes of DEGs after Chang's procedures v.s. ENI3-R
[fig3, ngenes3] = plotTable('Profiles-Genes-overlapping-EIN3-R-mRNA-ethylene-regulated.csv');
title(sprintf('profiles of overlapping genes of DEGs after Chang procedures v.s. ENI3-R \n with %d genes'...
    ,ngenes3),'FontSize',14)
print(fig3,'./Figures/fig3','-dpng');

% plot profiles of genes in EIN3-R but not in DEGs of Chang's procedures
[fig4, ngenes4] = plotTable('Profiles-Genes-EIN3-R-Not-mRNA-ethylene-regulated.csv');
title(sprintf('profiles of genes in EIN3-R but not in DEGs of Chang procedures, \n with %d genes',ngenes4),...
    'FontSize',14)
print(fig4,'./Figures/fig4','-dpng');

%% PLOT Nan
% plot profiles after fdr
[fig5, ngenes] = plotTable('Profiles-Gene-list-fdr.csv');
title(sprintf('profiles after fdr, \n with %d genes',ngenes),...
    'FontSize',14)
print(fig5,'./Figures/fig5','-dpng');

% plot profiles of overlapping genes of EIN3-R vs FDR DEs
[fig6, ngenes] = plotTable('Profiles-Genes-overlapping-EIN3-R-FDR.csv');
title(sprintf('profiles of overlapping genes of EIN3-R vs FDR DEs, \n with %d genes',ngenes),...
    'FontSize',14)
print(fig6,'./Figures/fig6','-dpng');

% profiles of genes in EIN3-R but not in FDR procedures
[fig7, ngenes] = plotTable('Profiles-Genes-EIN3-R-Not-FDR.csv');
title(sprintf('profiles of genes in EIN3-R but not in FDR procedures \n with %d genes'...
    ,ngenes),'FontSize',14)
print(fig7,'./Figures/fig7','-dpng');





