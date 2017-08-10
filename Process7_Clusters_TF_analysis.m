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