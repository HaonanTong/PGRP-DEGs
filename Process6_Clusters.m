% Nan_Cluster(i,j): i = 1 : 6, indicates gene set that activates at different time
% point. j = 1 : 2, indicates up-regulated or down-regulated respectively
for Tp = 2 : 7
    Timepoint = Tp;% 2-7
    T = readtable(sprintf('DEGs-time-period%d.csv',Timepoint-1),...
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

    writetable(T_up,sprintf('DEGs-time-period%d-up.csv',Tp-1),...
            'WriteRowNames',true);

    writetable(T_down ,sprintf('DEGs-time-period%d-down.csv',Tp-1),...
            'WriteRowNames',true);

    f_plotTable2( sprintf('DEGs-time-period%d-up.csv',Tp-1), [], 'Mean Plot' );
    f_plotTable2( sprintf('DEGs-time-period%d-down.csv',Tp-1), [], 'Mean Plot' );
end
