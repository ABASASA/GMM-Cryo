function [m2Bias] =ComputeBiasMatrixEmpricly(B, sigma, gridSize,  g, x_2d, y_2d, gamma, numberOfProjToCalBias, repNum)
beta  = 1;       % Bandlimit ratio
eps_p = 1e-3;    % Prescribed accuracy

if ~exist('numberOfProjToCalBias', 'var')
    numberOfProjToCalBias = 150000;
end
if ~exist('repNum', 'var')
    repNum = 1;
end

%% Compute Projecitons
paramsName = [];%'C1_params' ;


total_N = numberOfProjToCalBias;
m2Bias = 0;
for iRep = 1 : repNum

    [proj ,weight] = ComputeProjection (total_N, gridSize, paramsName, B, g, x_2d, y_2d);

    noise = randn(size(proj)) * sigma;
    [proj_PSWFNoiseOnly] = Projection2PSWF(noise, beta, eps_p, gamma);
    %%
    m2_hat = 0;
    for i=1:total_N
        % averaging
        current_vec_proj = proj_PSWFNoiseOnly(:,:,i);    
        m2_hat = m2_hat + current_vec_proj(:) * current_vec_proj(:)' * weight(i);
    end

    m2Emp = zeros(max(gamma.ang_idx_2d)+1, nnz(gamma.ang_idx_2d==0),...
                                max(gamma.ang_idx_2d)+1, nnz(gamma.ang_idx_2d==0));
    for m=0:max(gamma.ang_idx_2d)
        for k=0:nnz(gamma.ang_idx_2d==0)-1
            m2Emp(m+1,k+1,:,:) = reshape(m2_hat(m+1 + (k)*(max(gamma.ang_idx_2d)+1),:),...
                 max(gamma.ang_idx_2d)+1, nnz(gamma.ang_idx_2d==0) );
        end
    end
    BiasTmp = m2Emp / sum(weight);
    m2Bias = m2Bias + BiasTmp ./ repNum;
end

end