function [ fig, G ] = Analyze_GroundTruthTable(GroundTruthTable)
%Analyze_GroundTruthTable plots a histogram of cardinalities at different
%time periods
%   GroundTruthTable       - array of binary-value for a certain condition 
%   fig - Hist plot of cardinalities of groundtruth table at different time
%   period
%   G - cardinalities for each time point;
    % Creating a new figure
    fig = figure;
    clf;
    % Array of cardinalities
    G = sum(GroundTruthTable,1);
    % Number of bars
    lg = length(G);
    % Histogram plot
    bar(1:lg,G,'BarWidth',1);
    % Labeling bars
    for r = 1:lg
        text(r-0.15, G(r)+0.01*sum(G), ['$\mathcal{G}_' int2str(r) '$'],...
            'Interpreter', 'Latex');
    end
    % Setting labels
    xlim(0.5 + [0 lg])
    %set(gca, 'XTick', 0.5:(lg + 0.5));
    %set(gca, 'XTickLabel', lbls);
    grid on
    xlabel('t');
    ylabel('Cardinality');
end