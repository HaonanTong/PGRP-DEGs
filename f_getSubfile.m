function  [sub_Name, sub_Table, ngene] = f_getSubfile( origin_csv, isVariableNames, indx, file_expr_idx, file_agis_idx, file_Num_agis_idx )
%   [sub_Name, sub_Table, ngene] = f_getSubfile( origin_csv, isVariableNames, indx, file_expr_idx, file_agis_idx, isWriteFile )
% sub_Name: gene list of selected gene from original file based on index
% information
% sub_Table: expresion profile of sub_Name gene list
% ngene: # of gene seleted
% orgin_csv: original expression profile data
% isVariableNames: 1 if the table of origin_csv contains Varible names
% indx: array of logical or positive integers
% file_expr_idx: expression profile of sub gene list
% file_agis_idx: gene list that is selected
% file_Num_agis_idx: # of genes that is selected
% example: [sub_Name, sub_Table, ngene] = f_getSubfile( 'kat-rpkm-expression.csv', 1, [1 2 3], 'profiles_subgl.csv', 'subgl.txt' )
% example: [sub_Name, sub_Table, ngene] = f_getSubfile( 'Profiles-ANan-down-regulated.csv', 0, [1 2 3] )
sub_Name = NaN;
sub_Table = NaN;
ngene = NaN;
if nargin ~= 6 && nargin ~= 3
    help f_getSubfile
    return;
end


if isVariableNames == 1
    T = readtable(origin_csv,...
     'ReadVariableNames',true);
else
    T = readtable(origin_csv,...
     'ReadVariableNames',false);
end

ngene = length(indx);

name = table2array(T(:,1));

sub_Name = name(indx);
sub_Table = T(indx,:);

if nargin == 6
    fp = fopen(file_agis_idx,'wt');
    fprintf(fp, '%s\n', sub_Name{:,:});
    fclose(fp);
    if isVariableNames == 1
        writetable( sub_Table, file_expr_idx,'WriteRowNames',...
            true,'WriteVariableNames',true);
    else
        writetable( sub_Table, file_expr_idx,'WriteRowNames',...
            true,'WriteVariableNames',false);
    end
    
    fp = fopen(file_Num_agis_idx,'wt');
    fprintf(fp, '%d\n', ngene);
    fclose(fp);    
end

end

