function [m1_hat, m2_hat] = EstimateMoments(proj_PSWF, weight, total_N, gamma)
%% average to attain the empirical moments
% fprintf('Averaging...');

% initialize  variables
perm_list = randperm(total_N);
m1_hat = 0;
m2_hat = 0;
for i=1:total_N
    % averaging
    m1_hat = m1_hat + proj_PSWF(:,:,perm_list(i))*weight(perm_list(i));
    current_vec_proj = proj_PSWF(:,:,perm_list(i));
    m2_hat = m2_hat + current_vec_proj(:)*current_vec_proj(:)'*weight(perm_list(i));
end

% second moment initialization
current_m2_check = zeros(max(gamma.ang_idx_2d)+1, nnz(gamma.ang_idx_2d==0),...
                            max(gamma.ang_idx_2d)+1, nnz(gamma.ang_idx_2d==0));

%reshape the second moment and summarize
for m=0:max(gamma.ang_idx_2d)
    for k=0:nnz(gamma.ang_idx_2d==0)-1
        current_m2_check(m+1,k+1,:,:) = reshape(m2_hat(m+1 + (k)*(max(gamma.ang_idx_2d)+1),:),...
            max(gamma.ang_idx_2d)+1, nnz(gamma.ang_idx_2d==0) );
    end
end
m2_hat = current_m2_check/sum(weight(perm_list(1:end)));
% fprintf('...DONE \n');

end