%% running optimization
clc
option = 2;
tic
if option ==1
    % optimizing option 1
    %     V_ab = A_B_to_VecAB(A, B, gamma.band_idx_3d);
    %     % WITH GT distribution
    %     initial_guess((end-length(B_part)+1):end) = V_ab((end-length(B_part)+1):end);
    %
    %     [A_est, B_est] = LS_matlab_PSWF(initial_guess, m1, m2, gamma, ...
    %         C_tensor, Gamma_mat , sign_mat);
    %     t_opt = toc;
else
    % optimizing option 2

    N = length(vec_AB_GT);
%     N = length(vec_A_GT);
%     N = 1;

    manifold  = euclideancomplexfactory(N, 1);
    problem.M = manifold;
    weight    = numel(m1_hat)/numel(m2_hat);
    weight    = norm(m1_hat(:))^2/norm(m2_hat(:))^2;
    numElinMopment =  numel(m1_hat) + numel(m2_hat);
    weight = rand(numElinMopment, numElinMopment);
    weight = triu(weight,1) + triu(weight).'; 
    weight(1) =100;
    if P==1
        problem.costgrad = @(x)  LS_cost_w_grad_JUST_VOL_grad_fix(x, B, gamma, C_tensor, Gamma_mat , ...
            sign_mat, m1_hat, m2_hat, weight);
    else
       % problem.costgrad = @(x) LS_cost_w_grad(x, gamma, C_tensor, ...            
       %     Gamma_mat, sign_mat, m1, m2, weight);
       Avec = A_B_to_VecAB(A,[], M-1);

        problem.costgrad = @(x) GMM_cost_in_plane_uniform(x, gamma, C_tensor, ...            
            Gamma_mat, sign_mat, m1_hat, m2_hat, weight, B, Avec);
    end
    checkgradient(problem);
    
    options.maxinner = 35;
    options.tolgradnorm  = 1e-20;
    options.maxiter  = 60;
    [x1, xcost1, info1, ~] = trustregions(problem, initial_guess, options);
    x = x1;
    %     options.maxiter  = 200;
    %     options.minstepsize = 1e-20;
    %     % [x, xcost, info, options] = steepestdescent(problem, x, options);
    %     [x, xcost, info, options] = conjugategradient(problem, x1, options);
    
    [A_est, ~] = VecAB_inplane_to_A_B(x, gamma);
end
t_opt = toc;
save(['working_space_',nameit]);