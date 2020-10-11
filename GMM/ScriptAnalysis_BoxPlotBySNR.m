% Params
SNRs= [1];% 0.5 0.6];
ParametrString = 'Radial sigma';
folder = '061020 - BallsNew_L15_Beta1_TruncatedVol_NoReg_IdealGMM_AnalyticBias_Singels/Radial_sigma_1-10Exp';
numberOfExp = 10;
% Init
proportionsRelError = zeros(length(SNRs), numberOfExp);
proportionsFSC = zeros(length(SNRs), numberOfExp);
errorAndFscLS = zeros(length(SNRs), numberOfExp, 2);
errorAndFscGMM = zeros(length(SNRs), numberOfExp, 2);

for iSNR = 1 : length(SNRs) 
    %% Load
    for  t1 = 1 : numberOfExp 
     t = t1 + (iSNR-1) * numberOfExp;
    load(['E:\Workspace\לימודים\תזה\CryoEM with GMM\cryo_code\ResultsGMM\' folder '\Exp' num2str(t) '\data.mat'])
    proportionsRelError(iSNR,t1) =  rel_errLS ./ rel_errGMM;
    proportionsFSC(iSNR,t1) = resAEstFSCLS ./ resAEstFSCGMM;
    
    errorAndFscLS(iSNR, t1, :) = [rel_errLS, resAEstFSCLS];
    errorAndFscGMM(iSNR, t1, :) = [rel_errGMM, resAEstFSCGMM];
    close all;
    t
    end

end
%% Display
figureToSavePath = ['E:\Workspace\לימודים\תזה\CryoEM with GMM\cryo_code\ResultsGMM\' folder  '\figures\'];
mkdir(figureToSavePath);
for iSNR = 1 : length(SNRs)
       
    %% Estimator
    ratioRelError = squeeze(errorAndFscLS(iSNR,:,1) ./ errorAndFscGMM(iSNR,:,1));
    ratioFSC = squeeze(errorAndFscLS(iSNR,:,2) ./ errorAndFscGMM(iSNR,:,2));

    
    meanRatioRelError = mean(ratioRelError);
    medianRatioRelError = median(ratioRelError);

    meanRatioFSC = mean(ratioFSC);
    medianRatioFSC = median(ratioFSC);
    %% Display
    fig1 = figure;
    subplot(2,1,1);
    histogram(ratioRelError, 20);
    title([ ParametrString ': ' num2str(SNRs(iSNR)) ' Rellative Error of: LS/GMM, Mean: ' num2str(round(meanRatioRelError,2))...
       '  Median: ' num2str(round(medianRatioRelError,2))])

    subplot(2,1,2);
    histogram(ratioFSC, 20);
    title([ParametrString ': ' num2str(SNRs(iSNR)) ' FSC_{50} of: LS/GMM, Mean: ' num2str(round(meanRatioFSC,2))...
       '  Median: ' num2str(round(medianRatioFSC,2))])
    saveas(fig1, [figureToSavePath 'Hist_' ParametrString '_' num2str(SNRs(iSNR)) '.jpg']);
    saveas(fig1, [figureToSavePath 'Hist_' ParametrString '_' num2str(SNRs(iSNR)) '.fig']);

end

for iGraph = 1 : 2
    fig = figure;

    xlabel(ParametrString);
    ylabel('Ratio');
    if (iGraph == 1)
        [fig] = BoxPlotAsaf(fig, SNRs, proportionsRelError, 'b*-');
        title('Ratio Relative Error - LS / GMM');
    elseif (iGraph == 2)
        [fig] = BoxPlotAsaf(fig, SNRs, proportionsFSC, 'b*-');
        title('Ratio FSC_{50}  - LS / GMM');
    end
    hold on;
    plot(1: size(proportionsFSC,1), ones(size(proportionsFSC,1),1),'k--');
    xlabel(ParametrString);
    ylabel('Ratio');
    
    if (iGraph == 1)
        saveas(fig, [figureToSavePath 'Ratio_Rel_Error.jpg']);
        saveas(fig, [figureToSavePath 'Ratio_Rel_Error.fig']);
    elseif ( iGraph == 2)
        saveas(fig, [figureToSavePath 'Ratio_FSC50.jpg']);
        saveas(fig, [figureToSavePath 'Ratio_FSC50.fig']);
    end
end
fig = figure;
hold on;
subplot(2,1,1);
plot(SNRs, mean(squeeze(errorAndFscLS(:,:,1)),2), 'r*-',...
     SNRs, mean(squeeze(errorAndFscGMM(:,:,1)),2), 'b*-');
title('Mean Rel. Error')
legend('LS', 'GMM');
xlabel(ParametrString);
subplot(2,1,2);
plot(SNRs, std(squeeze(errorAndFscLS(:,:,1)),[],2), 'r*-',...
     SNRs, std(squeeze(errorAndFscGMM(:,:,1)),[],2), 'b*-');
title('std Rel. Error')
legend('LS', 'GMM');
xlabel(ParametrString);
saveas(fig, [figureToSavePath 'Rel_Error.jpg']);
saveas(fig, [figureToSavePath 'Rel_Error.fig']);

figure;
hold on;
subplot(2,1,1);
plot(SNRs, mean(squeeze(errorAndFscLS(:,:,2)),2), 'r*-',...
     SNRs, mean(squeeze(errorAndFscGMM(:,:,2)),2), 'b*-');
hold on;
title('Mean FSC_{50}')
legend('LS', 'GMM');
xlabel(ParametrString);
subplot(2,1,2);
hold on;

plot(SNRs, std(squeeze(errorAndFscLS(:,:,2)),[],2), 'r*-',...
     SNRs, std(squeeze(errorAndFscGMM(:,:,2)),[],2), 'b*-');
title('std FSC_{50}')
legend('LS', 'GMM');
xlabel(ParametrString);
saveas(fig, [figureToSavePath 'FSC_50.jpg']);
saveas(fig, [figureToSavePath 'FSC_50.fig']);



%%
if (length(SNRs)  == 1)
   disp(['Mean LS: Rel. Error: ', num2str(mean(squeeze(errorAndFscLS(:,:,1)),2))]); 
   disp(['STD LS: Rel. Error: ', num2str(std(squeeze(errorAndFscLS(:,:,1)),[],2))]); 
    
   
   disp(['Mean GMM: Rel. Error: ', num2str(mean(squeeze(errorAndFscGMM(:,:,1)),2))]); 
   disp(['STD GMM: Rel. Error: ', num2str(std(squeeze(errorAndFscGMM(:,:,1)),[],2))]);
   
   
   disp(['Mean LS: Rel. Error: ', num2str(mean(squeeze(errorAndFscLS(:,:,2)),2))]); 
   disp(['STD LS: Rel. Error: ', num2str(std(squeeze(errorAndFscLS(:,:,2)),[],2))]); 
    
   
   disp(['Mean GMM: Rel. Error: ', num2str(mean(squeeze(errorAndFscGMM(:,:,2)),2))]); 
   disp(['STD GMM: Rel. Error: ', num2str(std(squeeze(errorAndFscGMM(:,:,2)),[],2))]); 
   
   
end