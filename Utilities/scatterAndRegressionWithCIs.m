% Creates a scatter plot between x and y, performs a linear regression and
% then displays all the results.

function [f,mdl] = scatterAndRegressionWithCIs(x,y,n,xString,yString,titleString,taskName)
mdl = fitlm(x,y,'linear');
xVals = linspace(min(x),max(x));
[mdlVals,mdlCIs] = predict(mdl,xVals','Alpha',0.05);
titleString = strcat([titleString,' ',taskName]);

f = figure;
hold on;

scatter(x,y,'filled');
plot(xVals,mdlVals,'Color',[0.8500 0.3250 0.0980],'LineWidth',2);
plot(xVals,mdlCIs(:,1),'--','Color',[0.8500 0.3250 0.0980]);
plot(xVals,mdlCIs(:,2),'--','Color',[0.8500 0.3250 0.0980]);

dataString = sprintf('Data (n = %i)',n);
lrString = sprintf('Linear Regression (R^{2} = %1.2f)',mdl.Rsquared.Adjusted);
lgd = legend(dataString,lrString,'95% CI','location','Best');

title(titleString);
xlabel(xString);
ylabel(yString);
ax = gca;
ax.FontSize = 32;
ax.FontWeight = 'bold';
ax.LineWidth = 2.5;
lgd.FontSize = 24;
lgd.FontWeight = 'normal';
f.WindowState = 'maximized';
hold off;
end