%% Grep TF and nTF from DEGs at different time point
cd Process7
% unix grep TFs and ~TFs

myDir = './'; %gets directory
myFiles = dir(fullfile(myDir,'DEGs-time-period*.csv')); %gets all wav files in struct
for i = 1 : length(myFiles)
    % TFs that showns DE patterns
    unix(sprintf('head -1 %s > TFs-%s',myFiles(i).name,myFiles(i).name));
    unix(sprintf('fgrep -f TFs.txt %s >> TFs-%s',myFiles(i).name,myFiles(i).name));

    % genes not TFs but show DE patterns
    unix(sprintf('fgrep -v -f TFs.txt %s >> nTFs-%s',myFiles(i).name,myFiles(i).name));
end
cd ..

%% Plot and anlyze different period
ngene = zeros(6,6);
% ngene(i,1) - total go down - Down
% ngene(i,2) - total go up - Up
% ngene(i,3) - nTFs go down - Down TFs
% ngene(i,4) - nTFs go up - Up TFs
% ngene(i,5) - TFs go down - Down nTFs
% ngene(i,1) - TFs go up - Up nTFs
% at time phase i.
for i = 1 : 6 % time period
    myDir = sprintf('./Process7/period%d',i); %gets directory
    myFiles = dir(fullfile(myDir,'*.csv')); %gets all wav files in struct
    for j = 1 : length(myFiles)
        [ ~, ngene(i,j)] = f_plotTable2(sprintf('%s/%s',myDir,myFiles(j).name),[],'Mean Plot');
    end
end

T_cluster_summary = array2table(ngene,'RowNames',{'T0.25','T0.5','T1','T4','T12','T24'},'VariableNames',...
    {'Down','Up','Down_TFs', 'Up_TFs', 'Down_nTFs', 'Up_nTFs'});
writetable(T_cluster_summary,'./Process7/TimePhaseAnalysis.csv','WriteRowNames',true,'WriteVariableNames',true);

T_up = T_cluster_summary(:,{'Up','Up_nTFs','Up_TFs'});
writetable(T_up,'./Process7/TimePhaseAnalysis_Up.csv','WriteRowNames',true,'WriteVariableNames',true);

T_down = T_cluster_summary(:,{'Down','Down_nTFs','Down_TFs'});
writetable(T_down,'./Process7/TimePhaseAnalysis_Down.csv','WriteRowNames',true,'WriteVariableNames',true);

%% Barplot
% up
fig1 = figure;
b_up = bar(table2array(T_up(:,2:3)),'stacked');
b_up(2).FaceColor = 'red';%TFs
b_up(1).FaceColor = 'blue';%nTFs
leg={'nTFs','TFs'};
legend(leg,'FontSize',14,'Location','best');
str = {'.25','.5','1','4','12','24'};
set(gca, 'XTickLabel',str, 'XTick',1:numel(str))
xlabel('time(hrs)');
title({'Bar plot for number of genes',...
    'that activated in up-regulated tendency at different time point'},'FontSize',14);
print(fig1,'./Process7/up','-dpng');

% down
fig2 = figure;
b_down = bar(table2array(T_down(:,2:3)),'stacked');
b_down(2).FaceColor = 'red';%TFs
b_down(1).FaceColor = 'blue';%nTFs
leg={'nTFs','TFs'};
legend(leg,'FontSize',14,'Location','northwest');
str = {'.25','.5','1','4','12','24'};
set(gca, 'XTickLabel',str, 'XTick',1:numel(str))
xlabel('time(hrs)');
title({'Bar plot for number of genes','that activated in down-regulated tendency'},'FontSize',14);
print(fig2,'./Process7/down','-dpng');



