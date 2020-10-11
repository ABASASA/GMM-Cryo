function [rel_err, resAEstFSC, figFSC] = AnalytsisResultsFCS_Alignment_relError(AOriginal, AEst, filepath,...
                estimationMethodNameString, isPadding, AFull, gridSize, beta, delta, nameit,c, figCurrent)

if ~exist('figCurrent', 'var')
    figCurrent = figure;
end

if isPadding
    A_Est_padded = AFull;
    A_padded     = AFull;
    for j=1:length(AFull{1})
        if j<=length(AEst{1})
            A_padded{1}{j}     = AOriginal{1}{j};
            A_Est_padded{1}{j} = AEst{1}{j};
        else
            A_padded{1}{j}     = zeros(size(AFull{1}{j}));
            A_Est_padded{1}{j} = zeros(size(AFull{1}{j}));

        end
    end
    AOriginal = A_padded;
    AEst = A_Est_padded;
end

%inverse prolates stransform  TO DO AGAIN WITH OLD PACKAGE
vol_Est  = pswf_t_b_3d(AEst, gridSize, c, delta);
vol      = pswf_t_b_3d(AOriginal,     gridSize, c, delta);
%     A = pswf_t_f_3d(vol, beta, delta);
% print out volumes
VisualVol(vol_Est,[filepath, estimationMethodNameString,'_',nameit]);
VisualVol(vol, [filepath,'vol_',nameit]);
[~, ~, volREst] = cryo_align_densities(vol, vol_Est,1 ,1);
VisualVol(volREst, [filepath,'rot_', estimationMethodNameString, '_',nameit]);
rel_err       = norm(volREst(:)-vol(:))/norm(vol(:));

[resAEstFSC,figFSC] = plotFSC_v2(vol,volREst, 0.5,1, estimationMethodNameString, figCurrent);

saveas(figFSC, [filepath,'_FSC', estimationMethodNameString, '_',nameit], 'fig')
saveas(figFSC, [filepath,'_FSC', estimationMethodNameString, '_',nameit], 'jpg')

end