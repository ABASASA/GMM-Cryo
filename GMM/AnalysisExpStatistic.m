foler = '270920 - BallsNew_L15_Beta1_TruncatedVol_NoReg_SNR_01_05_2_IdealGMM-small_Stat';
numberOfExp = 50;
ratioFSC = zeros(numberOfExp, 1);
ratioRelError = zeros(numberOfExp, 1);
%% Load
for  t1 = 1 : numberOfExp 
 t = t1 + 2 * numberOfExp;
load(['E:\Workspace\לימודים\תזה\CryoEM with GMM\cryo_code\ResultsGMM\' foler '\Exp' num2str(t) '\data.mat'])
ratioRelError(t1) =  rel_errLS ./ rel_errGMM;
ratioFSC(t1) = resAEstFSCLS ./ resAEstFSCGMM;
close all;
t
end
%% Estimator
meanRatioRelError = mean(ratioRelError);
medianRatioRelError = median(ratioRelError);

meanRatioFSC = mean(ratioFSC);
medianRatioFSC = median(ratioFSC);
%% Display
fig = figure;
subplot(2,1,1);
histogram(ratioRelError, 20);
title(['Rellative Error of: LS/GMM, Mean: ' num2str(round(meanRatioRelError,2))...
   '  Median: ' num2str(round(medianRatioRelError,2))])

subplot(2,1,2);
histogram(ratioFSC, 20);
title(['FSC_{50} of: LS/GMM, Mean: ' num2str(round(meanRatioFSC,2))...
   '  Median: ' num2str(round(medianRatioFSC,2))])
