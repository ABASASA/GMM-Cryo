%% Info
% This is a box/candle plot of the data)


function [fig] = BoxPlotAsaf(fig, xData, yData1, plotProp)
%% Create Labels
labels = cell(length(xData),1);
for i = 1 : length(xData)
    labels{i} = num2str(round(xData(i),4));
end

%% plot
hold on;
boxplot(yData1', 'Labels', labels, 'Whisker',0,'OutlierSize',0.1);
hold on;
hAx=gca;                                   % retrieve the axes handle
xtk=hAx.XTick;                             % and the xtick values to plot() at...
hold on
plot(xtk, mean(yData1,2), plotProp)
end