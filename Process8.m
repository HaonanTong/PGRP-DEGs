%% Gene Activition at different time point.
% Genes activation time point.
% File that contains infor of genes that show DE pattern at specifice time point
myDir = './'; %gets directory
myFiles = dir(fullfile(myDir,'DEGs-time-period*.csv')); %gets all wav files in struct

T = cell(6,1);
for i = 1 : 6
    % read Tablei that contents gene set show DE pattern at i time point.
    T{i} = readtable(myFiles(i).name,...
     'ReadVariableNames',true,'ReadRowNames',false);
end

for i = 1 : 5
    % compare Tablei v.s. afterwards tablei+1 i+2 ... 6.
    % Remove overlapping genes in afterwards tablei+1 i+2 ... 6
    for j = i+1 : 6
        [C,ia,ib] = intersect(T{i},T{j});
        T{j}(ib,:)=[];
    end
end

for i = 1 : 6
     % Output tables that contains gene set that activated at the same time
     writetable(T{i},sprintf('./Process8/DEGs-time-Activation%d.csv',i),'WriteRowNames',false,'WriteVariableNames',true);
end

% Plot
myDir = './Process8/'; %gets directory
myFiles = dir(fullfile(myDir,'DEGs-time-Activation*.csv')); %gets all wav files in struct
ngene = zeros(6,1);
for i = 1 : 6
    % [ fig, ngene(i), expr,plotData, agis, agis_new ] = f_plotTable2(sprintf('%s%s',myDir,myFiles(i).name),[],'Normalized');
    % f_plotTable2(sprintf('%s%s',myDir,myFiles(i).name),[],'Mean Plot');
end
%
% T_Summary = array2table(ngene,'RowNames',{'.25','.5','1','4','12','24'},'VariableNames',{'ngene'});
% writetable(T_Summary,'./Process8/DEGs-Activation-Summary.csv','WriteRowNames',true,'WriteVariableNames',true);

%% Patterns of Activation(Up/Down regulated compared to original time point.
for Tp = 2 : 7
    Timepoint = Tp;% 2-7
    T = readtable(sprintf('./Process8/DEGs-time-Activation%d.csv',Timepoint-1),...
         'ReadVariableNames',true,'ReadRowNames',true);
    % summary(T);
    Data = table2array(T);
    name = T.Properties.RowNames;

    T0 = Data(:,1:3)';
    T1 = Data(:,4:6)';
    T2 = Data(:,7:9)';
    T3 = Data(:,10:12)';
    T4 = Data(:,13:15)';
    T5 = Data(:,16:18)';
    T6 = Data(:,19:21)';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Derive mean of Expression data set
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mData = [];
    for i = 1:3:21%7 time points; 3 replicates;
       mData = [mData sum(Data(:,i:i+2),2)];
    end
    mData = 1/3*mData;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % the 50% difference as NC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DEi = log2(mData(:,Timepoint)./mData(:,1));
    DEi_up = DEi >= 0.58;
    DEi_down = DEi <= -0.58;

    T_up = T(DEi_up,:);
    T_down = T(DEi_down,:);

    writetable(T_up,sprintf('./Process8/DEGs-time-Activation%d-up.csv',Tp-1),...
            'WriteRowNames',true);

    writetable(T_down ,sprintf('./Process8/DEGs-time-Activation%d-down.csv',Tp-1),...
            'WriteRowNames',true);

    f_plotTable2( sprintf('./Process8/DEGs-time-Activation%d-up.csv',Tp-1), [], 'Mean Plot' );
    f_plotTable2( sprintf('./Process8/DEGs-time-Activation%d-down.csv',Tp-1), [], 'Mean Plot' );
    f_plotTable2( sprintf('./Process8/DEGs-time-Activation%d-up.csv',Tp-1), [], 'Normalized' );
    f_plotTable2( sprintf('./Process8/DEGs-time-Activation%d-down.csv',Tp-1), [], 'Normalized' );
end

%% TFs and nTFs analysis
% Grep TF and nTF from DEGs at different time point
cd Process8
% unix grep TFs and ~TFs

myDir = './'; %gets directory
myFiles = dir(fullfile(myDir,'DEGs-time-Activation*.csv')); %gets all wav files in struct
for i = 1 : length(myFiles)
    % TFs that showns DE patterns
    unix(sprintf('head -1 %s > TFs-%s',myFiles(i).name,myFiles(i).name));
    unix(sprintf('fgrep -f TFs.txt %s >> TFs-%s',myFiles(i).name,myFiles(i).name));

    % genes not TFs but show DE patterns
    unix(sprintf('fgrep -v -f TFs.txt %s >> nTFs-%s',myFiles(i).name,myFiles(i).name));
end

for i = 1 : 6
    mkdir(sprintf('period%d',i));
    unix(sprintf('mv *%d* period%d',i,i));
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
    myDir = sprintf('./Process8/period%d',i); %gets directory
    myFiles = dir(fullfile(myDir,'*.csv')); %gets all wav files in struct
    for j = 1 : length(myFiles)
        [ ~, ngene(i,j)] = f_plotTable2(sprintf('%s/%s',myDir,myFiles(j).name),[],'Mean Plot');
    end
end

T_cluster_summary = array2table(ngene,'RowNames',{'T0.25','T0.5','T1','T4','T12','T24'},'VariableNames',...
    {'Down','Up','Down_TFs', 'Up_TFs', 'Down_nTFs', 'Up_nTFs'});
writetable(T_cluster_summary,'./Process8/TimeActivationAnalysis.csv','WriteRowNames',true,'WriteVariableNames',true);

T_up = T_cluster_summary(:,{'Up','Up_nTFs','Up_TFs'});
writetable(T_up,'./Process8/TimeActivationAnalysis_up.csv','WriteRowNames',true,'WriteVariableNames',true);

T_down = T_cluster_summary(:,{'Down','Down_nTFs','Down_TFs'});
writetable(T_down,'./Process8/TimeActivationAnalysis_down.csv','WriteRowNames',true,'WriteVariableNames',true);

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
