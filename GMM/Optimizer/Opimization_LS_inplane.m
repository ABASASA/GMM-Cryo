function [t_opt, x, xcost, info ] = Opimization_LS_inplane(N, weight, P, gamma, C_tensor, ...            
        Gamma_mat, sign_mat, m1, m2, initial_guess )
% optimizing option 2
manifold  = euclideancomplexfactory(N, 1);
problem.M = manifold;
if P==1
%     problem.costgrad = @(x)  LS_cost_w_grad_JUST_VOL_grad_fix(x, B, gamma, C_tensor, Gamma_mat , ...
%         sign_mat, m1, m2, weight);
else
    problem.costgrad = @(x) LS_cost_in_plane_uniform(x, gamma, C_tensor, ...            
        Gamma_mat, sign_mat, m1, m2, weight);
end

options.maxinner = 35;
options.tolgradnorm  = 1e-20;
options.maxiter  = 60;
[x, xcost, info, ~] = trustregions(problem, initial_guess, options);
t_opt = toc;

end