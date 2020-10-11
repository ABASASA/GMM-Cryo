function [t_opt, xEst, xcost, info] =  Iterative_GMM_Inplane(N, P, gamma, C_tensor,...
                    Gamma_mat, sign_mat, m1_hat, m2_hat, vecBoolMoments, proj_PSWF, weight,...
                    svPercent, initial_guess, numMaxIter, numOfStepsGMM, WFirst, Bias)
   t_opt = 0; 
   
   %% First step with ginput weights matrix as weights
   
   [t_optCurrent, xEst, xcost, info] = Optimization_GMM_inplane(N, P, gamma, C_tensor,...
                    Gamma_mat, sign_mat, m1_hat, m2_hat, WFirst, vecBoolMoments, initial_guess, floor(numMaxIter / numOfStepsGMM));
    t_opt = t_opt + t_optCurrent;
    
    
    %% Number of GMM Step
    for i = 2 : numOfStepsGMM
        
        [A, B] = VecAB_inplane_to_A_B(xEst, gamma);
        
        [W, ~] = ComputeW_inplane(A, B, proj_PSWF, weight, size(proj_PSWF,3),...
                                gamma, Gamma_mat, sign_mat, C_tensor, vecBoolMoments, svPercent, Bias);
%       [W, ~,OmegaCut] = ComputeW_inplane(A, B, proj_PSWF, weight, total_N,...
%                                 gamma, Gamma_mat, sign_mat, C_tensor, vecBoolMoments, svPercent, Bias);
        [t_optCurrent, xEst, xcost, info] = Optimization_GMM_inplane(N, P, gamma, C_tensor,...
                    Gamma_mat, sign_mat, m1_hat, m2_hat, W, vecBoolMoments, xEst, 2 * floor(numMaxIter / numOfStepsGMM)); 
        t_opt = t_opt + t_optCurrent;

    end
    
end