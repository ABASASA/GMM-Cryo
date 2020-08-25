function [t_opt, xEst, xcost, info] = Optimization_GMM_inplane(vec_AB_GT, P, gamma, C_tensor,...
                    Gamma_mat, sign_mat, m1_hat, m2_hat, W, vecBoolMoments, initial_guess )

tic


N = length(vec_AB_GT);

%% Define manifold
manifold  = euclideancomplexfactory(N, 1);
problem.M = manifold;
%% Choose cost & grad function
if P==1
%     problem.costgrad = @(x)  LS_cost_w_grad_JUST_VOL_grad_fix(x, B, gamma, C_tensor, Gamma_mat , ...
%         sign_mat, m1_hat, m2_hat, weight);
else
    problem.costgrad = @(x) GMM_cost_in_plane_uniform(x, gamma, C_tensor, ...            
        Gamma_mat, sign_mat, m1_hat, m2_hat, W, vecBoolMoments);
end
%     checkgradient(problem);

options.maxinner = 35;
options.tolgradnorm  = 1e-20;
options.maxiter  =  60;
%% Run optimization
[xEst, xcost, info, ~] = trustregions(problem, initial_guess, options);
t_opt = toc;
% save(['working_space_',nameit])
end