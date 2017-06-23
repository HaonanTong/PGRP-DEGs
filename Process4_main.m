%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the profiles of up and down regulated genes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure1
[ fig, ngenes ] = plotTable( ...
    'Profies-Genes-ANan-Gene-list-DEG-down-regulated-Not-all-down.csv', 0 );
title(sprintf('Profies of Genes shown in Nan down regulated gene list but Not Chang, \n with %d genes',ngenes),...
    'FontSize',14)
print(fig,'./Figures/fig1','-dpng');

% Figure2
[ fig, ngenes ] = plotTable( ...
    'Profiles-ANan-up-regulated.csv', 0 );
title(sprintf('Profies of Genes shown in Nan up regulated gene Analyisis, \n with %d genes',ngenes),...
    'FontSize',14)
print(fig,'./Figures/fig2','-dpng');

% Figure3
[ fig, ngenes ] = plotTable( ...
    'Profiles-ANan-down-regulated.csv', 0 );
title(sprintf('Profies of Genes shown in Nan down regulated gene Analysis, \n with %d genes',ngenes),...
    'FontSize',14)
print(fig,'./Figures/fig3','-dpng');

% Figure4
[ fig, ngenes ] = plotTable( ...
    'Profiles-Genes-ANan-Gene-list-DEG-up-regulated-Not-all-up.csv', 0 );
title(sprintf('Profies of Genes shown in Nan up regulated gene list but Not Chang, \n with %d genes',ngenes),...
    'FontSize',14)
print(fig,'./Figures/fig4','-dpng');