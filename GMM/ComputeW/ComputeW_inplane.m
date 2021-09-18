% Compute the weights.
% Only work for inplane
function [W, Omega, OmegaCut] = ComputeW_inplane( A, B, proj_PSWF, weight, total_N,...
                            gamma, Gamma_mat, sign_mat, C_tensor, vecBoolMoments, svPercent, Bias)

%% Compute analytic moments
m1     = FirstMoment_PSWF_v2(A, B, Gamma_mat, sign_mat);
m2     = SecondMoment_PSWF_v2(A, B, gamma, C_tensor, Gamma_mat);

m1Size = size(m1);
m2Size = size(m2);

msAnalyticVec = [m1(:); m2(:)];
msAnalyticVec = msAnalyticVec(vecBoolMoments); % Trimed

%% Compute the diff between aanlytic to each projection
reminders = zeros(length(msAnalyticVec), total_N);
reminders2 = zeros(length(msAnalyticVec), total_N);

counter = 0;
for i = 1 : total_N
    if weight(i) == 0
        continue;
    end
    counter = counter + 1;
    [m1_hatTmp, m2_hatTmp] = EstimateMoments(proj_PSWF(:,:,i), weight(i), 1, gamma, Bias);
    % Trim
    vecMsEmpiricTmp = [m1_hatTmp(:); m2_hatTmp(:)];
    vecMsEmpiricTmp = vecMsEmpiricTmp(vecBoolMoments);
    if(sum(isnan(vecMsEmpiricTmp(:))))
        disp('a');
        EstimateMoments(proj_PSWF(:,:,i), weight(i), 1, gamma);
    end
%     reminders(:,i) = msAnalyticVec -  vecMsEmpiricTmp;
    reminders(:,i) = vecMsEmpiricTmp;

end
reminders = reminders(:,1:counter);
%% compute cov
Omega=  cov(reminders.');

OmegaCut = Omega;
% AA = svd(Omega);
% semilogy(AA);
% Omega = Omega + AA(800) * eye(size(Omega));
% tic;
W = inv(Omega);
% W = pinv(Omega);
% toc1  = toc;

% tic
% [U,S,V] = svds(Omega, floor(svPercent * size(Omega,1)));
% W = V * diag(1./ diag(S))  * U';

% [U,S,V] = svd(Omega);
% Sd = diag(S);
% Sinv = Sd;
% Sinv(Sd <= Sd(1) * svPercent) = 0;
% OmegaCut = U * diag(Sinv) * V';
% Sinv(Sd > Sd(1) * svPercent) = 1./Sinv(Sd > Sd(1) * svPercent);
% W = V * diag(Sinv)  * U';

% OmegaCut = Omega + eye(size(Omega)) *  Sd(1) * svPercent;
% W = inv(OmegaCut);

% toc2 = toc;

% err = norm(W-W1,'fro')

% disp(['Time inv: ' num2str(toc1)])
% disp(['Time svds: ' num2str(toc2)])

end