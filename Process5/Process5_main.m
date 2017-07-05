f_plotTable( 'Profiles-ANan-DEGs.csv', [], 'Mean Plot' );

for i =1 : 6
    f_plotTable( sprintf('DEGs-time-period%d.csv',i),...
        [], 'Mean Plot' );
end 
