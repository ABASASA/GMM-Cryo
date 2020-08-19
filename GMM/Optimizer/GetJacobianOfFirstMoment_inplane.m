function [m1Grad] = GetJacobianOfFirstMoment_inplane(x, gamma, Gamma_mat, sign_mat,B)


% back to cell arrays
[A, ~] = VecAB_inplane_to_A_B(x, gamma);

%% Create the general strutre of A and B
zerosA{1} = A{1};
for l = 1 : length(A{1})
    zerosA{1}{l} = zeros(size(A{1}{l}));
end

zerosB = B;
for l = 1 : length(B{1})
    zerosB{l} = zeros(size(B{l}));
end

%% init -  max degree, total size
max_deg = min(length(A{1})-1, size(B,1)-1);  
sizeM   = max(gamma.ang_idx_2d)+1 ;
L       = max(gamma.band_idx_3d);
P  = size(B,1);
m1GradA = cell(max(gamma.ang_idx_2d) + 1, nnz(gamma.ang_idx_2d==0));
m1GradB = cell(max(gamma.ang_idx_2d) + 1, nnz(gamma.ang_idx_2d==0));

%% main loop
for m=0:max(gamma.ang_idx_2d)
    for k=0:nnz(gamma.ang_idx_2d==m)-1
        current_val = 0;
        
        m1GradA{m+1, k+1} = zerosA;
        m1GradB{m+1, k+1} = {};

        tmpGradA = zerosA;
        tmpGradB = zerosB;
        
%         counterIndexesRows = 1;
        for l = abs(m) : L
            if l>P-1
%                 m1GradA{m+1, k+1} = tmpGradA;
%                 m1GradB{m+1, k+1} = tmpGradB;
                continue
            end

            for n = -l : l    % volume is assume to be real, otherwise n=-l:l
                valueGradB = 0;
                for s=0:nnz(gamma.band_idx_3d==l)-1
                  %% compute corrent gamma coeff
                    [coeff_lsmk] = get_Gamma_Coeff(l,s,m,k,gamma);
%                     
%                     tmpGradA{1}{l+1}(s+1,n+1) = ((-1)^(m-n)) *...
%                                         coeff_lsmk * B{l+1}(l+1-n,l+1-m);
%                     valueGradB = valueGradB + ...
%                           ((-1)^(m-n)) * coeff_lsmk * A{1}{l+1}(s+1,n+1);
                    
                    if n>=0
                        tmpGradA{1}{l+1}(s+1,n+1) = tmpGradA{1}{l+1}(s+1,n+1) + (((-1)^(m-n)) *...
                                            coeff_lsmk * (B{l+1}(l+1-n,l+1-m)));
                        valueGradB = valueGradB + ...
                         ((-1)^(m-n)) * coeff_lsmk * A{1}{l+1}(s+1,n+1);
%                         mid_val = mid_val + ...
%                         ((-1)^(m-n))*coeff_lsmk*B{l+1}(l+1-n,l+1-m)*A{1}{l+1}(s+1,n+1);
                    else
                        tmpGradA{1}{l+1}(s+1,abs(n)+1) = tmpGradA{1}{l+1}(s+1,abs(n)+1) + conj( ((-1)^(m-n)) *...
                                            coeff_lsmk * (B{l+1}(l+1-n,l+1-m)) * (-1)^abs(n));%+ conj(A{1}{l+1}(s+1,abs(n)+1));
                        valueGradB = valueGradB + ...
                            ((-1)^(m-n)) * coeff_lsmk * conj(A{1}{l+1}(s+1,abs(n)+1)) * (-1)^abs(n);

%                         mid_val = mid_val + ...
%                         ((-1)^(m-n))*coeff_lsmk*B{l+1}(l+1-n,l+1-m)*(-1)^abs(n)*conj(A{1}{l+1}(s+1,abs(n)+1));
                    end  
                end
                tmpGradB{l+1}(l+1-n,l+1-m) = valueGradB;
%                 tmpGradA{1}{l+1}(s+1,abs(n)+1) = tmpGradA{1}{l+1}(s+1,abs(n)+1) ./ (2*l+1);
            end
        end
        m1GradA{m+1, k+1} = tmpGradA;
        m1GradB{m+1, k+1} = tmpGradB;
    end
end

%% Change the order to Jacobian matrix
m1Grad = zeros(length(x), numel(m1GradA));
% % m1Grad = zeros(1, numel(m1GradA));
% r = 6;
for i = 1 :  numel(m1GradA)
    if isempty(m1GradA{i})
       continue; 
    end
    tmpAGrad = cell(1,1);

    tmpAGrad{1} = m1GradA{i}{1}; 
    % For A and B
    m1Grad(:,i) = A_B_to_vec_inplane(tmpAGrad, m1GradB{i}, sizeM-1);

    % For A
%     cc = A_B_to_vec_inplane(tmpAGrad,[], sizeM-1);
%     m1Grad(:,i) = cc;
end
end