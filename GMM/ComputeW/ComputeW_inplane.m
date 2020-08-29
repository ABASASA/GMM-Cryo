function [W, Omega] = ComputeW_inplane( A, B, proj_PSWF, weight, total_N,...
                            gamma, Gamma_mat, sign_mat, C_tensor, vecBoolMoments)

%% Compute analytic moments
m1     = FirstMoment_PSWF_v2(A, B, Gamma_mat, sign_mat);
m2     = SecondMoment_PSWF_v2(A, B, gamma, C_tensor, Gamma_mat);

m1Size = size(m1);
m2Size = size(m2);

msAnalyticVec = [m1(:); m2(:)];
msAnalyticVec = msAnalyticVec(vecBoolMoments); % Trimed

%% Compute the diff between aanlytic to each projection
reminders = zeros(length(msAnalyticVec), total_N);
counter = 0;
for i = 1 : total_N
    if weight(i) == 0
        continue;
    end
    counter = counter + 1;
    [m1_hatTmp, m2_hatTmp] = EstimateMoments(proj_PSWF(:,:,i), weight(i), 1, gamma);
    % Trim
    vecMsEmpiricTmp = [m1_hatTmp(:); m2_hatTmp(:)];
    vecMsEmpiricTmp = vecMsEmpiricTmp(vecBoolMoments);
    if(sum(isnan(vecMsEmpiricTmp(:))))
        disp('a');
        EstimateMoments(proj_PSWF(:,:,i), weight(i), 1, gamma);
    end
    reminders(:,i) = msAnalyticVec -  vecMsEmpiricTmp;
end
reminders = reminders(:,1:counter);
%% compute cov
Omega=  cov(reminders.');
AA = svd(Omega);
% semilogy(AA);
Omega = Omega + AA(800) * eye(size(Omega));
W = pinv(Omega);

end