[ fig, ngene, expr,plotData, agis, agis_new ] = ...
    f_plotTable2( 'kat-rpkm-expression.csv', [], 'Mean Plot' );
title(sprintf('Expresson profiles of ALL %d genes in Arabidopsis thaliana',...
    ngene),'FontSize',15)
ylim([-8,+8])
print(fig,'./Figures/ALLGENEPLOTS','-dpng');
