%% Produce a Bland-Altman plot for two measures
%  Also known as mean-difference plots. The title, axis labels, and legend
%  can be produced by uncommenting the appropriate lines.
%
%  6 December, 2021

function f = blandAltmanPlot(x,y,xString,yString,titleString)

average = mean([x y],2);
difference = x-y;

f = figure;
hold on;
scatter(average,difference,100,'filled','HandleVisibility','off');
yline(mean(difference),'-k','LineWidth',4,'DisplayName','Mean Difference');
yline(mean(difference)+(2*std(difference)),'--k','LineWidth',2,'DisplayName','Mean Difference +- 2 SD');
yline(mean(difference)-(2*std(difference)),'--k','LineWidth',2,'HandleVisibility','off');
%lgd = legend('location','Best');
%xlabel(xString);
%ylabel(yString);
%title(titleString);

yRange = max(abs(ylim));
ylim([-yRange yRange]);

ax = gca;
ax.FontSize = 48;
ax.FontWeight = 'bold';
ax.LineWidth = 10;
%lgd.FontSize = 24;
f.WindowState = 'maximized';
hold off;

end