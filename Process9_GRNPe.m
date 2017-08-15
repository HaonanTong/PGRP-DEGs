%% interpolation
myDir = './Process9/'; %gets directory
myFiles = dir(fullfile(myDir,'*.csv')); %gets all wav files in struct
Dir_Interp = './Process9/Interp/';
mkdir(Dir_Interp);
for i = 1 : length(myFiles)
    f_interp(sprintf('%s%s',myDir,myFiles(i).name));
end 
%% discritization
mkdir ./Process9/Interp/Dscrtz
dscritz_Dir = './Process9/Interp/Dscrtz';
myFiles = dir(fullfile(Dir_Interp,'*.csv')); %gets all wav files in struct
 
n_levels = 3 ;
for i = 1 : length(myFiles)
    f_discritize(sprintf('%s%s',Dir_Interp,myFiles(i).name),n_levels);
end