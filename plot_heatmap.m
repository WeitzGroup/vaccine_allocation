function plot_heatmap(xvals, yvals, cdata, lim_lower, lim_upper)
    figure;
    h = heatmap(xvals, yvals, cdata','CellLabelColor','none','ColorLimits',[lim_lower, lim_upper]);
    h.Colormap = parula;
    %h.Title = 'Reduction in A';
    xlabel('Percentage vaccinated daily ($\%$)')
    ylabel('$\kappa$')
    set(gca,'FontSize',20);
    h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
    h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
    h.NodeChildren(3).Title.Interpreter = 'latex';
    h.NodeChildren(3).YDir='normal'; 
end