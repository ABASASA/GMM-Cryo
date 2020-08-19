function [val, grad] = GMM_cost_in_plane_uniform(vec_AB, gamma, C_tensor, Gamma_mat , ...
                                        sign_mat, mu1_hat, mu2_hat, W, vecBoolMoments)

% A cost function for (complex) least squares optimization, with gradient!


% balancer between the two terms of mu1 and mu2
% if isempty(weight)
%     weight = 1;
% end

% back to cell-array form
[A,B] = VecAB_inplane_to_A_B(vec_AB, gamma);
TmpABVec = vec_AB;

%% For only A
% [A,~] = VecAB_inplane_to_A_B(vec_AB, gamma);
% TmpABVec = A_B_to_vec_inplane(A, B, max(gamma.ang_idx_2d));
%% For a single element uncomment
% i = 7;
% 
% vec_ABtmp = Avec;
% vec_ABtmp(i) = vec_AB;
% vec_AB = vec_ABtmp;
% 
% [A,~] = VecAB_inplane_to_A_B(vec_AB, gamma);
% TmpABVec = A_B_to_vec_inplane(A, B, max(gamma.ang_idx_2d));
%% 
% first moment
m1     = FirstMoment_PSWF_v2(A, B, Gamma_mat, sign_mat);
inner1 = m1 - mu1_hat;
term1  = norm(inner1(:))^2;   % "vector" norm

% % second moment
m2     = SecondMoment_PSWF_v2(A, B, gamma, C_tensor, Gamma_mat);
inner2 = m2 - mu2_hat;
term2  = norm(inner2(:))^2;   % "vector" norm

g = [inner1(:); inner2(:)]; % moment function

g = g(vecBoolMoments); % Trim

% function evaluation
val = g' * W * g;% + weight*term2;  % just second moment at this point

% gradient evaluation, if needed
if nargout>1
    %% Asaf implentation
%     J = GetJacobianOfFirstMoment_inplane(vec_AB, gamma, Gamma_mat, sign_mat, B);
% %     
% %     J = J(i, :); % UnComment for single element only
%     J = J(1:142, :); % UnComment for A only
% %     
%     grad = 2 * J * g;                   
%     grad = transpose(grad);
%     
    %% Nir's implentation

    innersW = W * g;
    
    expandedInnerW = zeros(size(vecBoolMoments));
    expandedInnerW(vecBoolMoments) = innersW;
    
    innerTmp1 = reshape(expandedInnerW(1:numel(inner1)), size(inner1));
    innerTmp2 = reshape(expandedInnerW(numel(inner1) + 1 : end), size(inner2));

%     innerTmp1 = reshape(W *  inner1(:), size(inner1));
    grad1 = get_grad_inplane_PSWF_GMM( TmpABVec, gamma, C_tensor, Gamma_mat, sign_mat,...
                              innerTmp1, innerTmp2);
% %     grad1 = grad1(i); % uncomment for single element
% %     grad1 = grad1(1:142); % uncomment for only A
    grad = grad1;


    % filtering the (gradient of the) zero angular coefficient to be REAL
    % IF RUNNING GRAD CHECK MUST COMMENT IT!
    start_ind = 1;                        
   j=1; % for j=1:length(A{1})
        s = size(A{1}{j},1);
        l_ind = start_ind:(start_ind+s-1);
        grad(l_ind) = real(grad(l_ind));
       start_ind = start_ind + numel(A{1}{j});
   % end
end

end


