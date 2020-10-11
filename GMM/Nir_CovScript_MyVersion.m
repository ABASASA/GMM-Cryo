% script name: "test_mu2_convergence_v5"
%
% We test the moments code using simulated volume in k-space
%
% Note: the most expensive part is calculating the second moment, with a new code
% we should repeat the test for larger P, for example.
%
% NS, 2019

clear variables; 
% close all; 
clc;

%% Create some probability (in alpha, beta, gamma). This probability will be replaced by a probability with finite expansion in Wigner-D
P  = 4;     % aribitrary, chosen to be small
beta  = 1;       % Bandlimit ratio
delta = 0.99;    % Truncation parameter
eps_p = 1e-3;    % Prescribed accuracy
SNR = 0;
total_N      = 20000;
R1 = randn(3); [R1,~] = qr(R1,0);
R2 = randn(3); [R2,~] = qr(R2,0);
R3 = randn(3); [R3,~] = qr(R3,0);


if P==1
    B{1} = 1;
else

    rho = @(a) 1*exp(-5*norm(eul2rotm(a)-R1,'fro')^2)+ 1*exp(-5*norm(eul2rotm(a)-R2,'fro')^2) + 1*exp(-5*norm(eul2rotm(a)-R3,'fro')^2);
    B   = WignerD_transform(rho, P);  % the coefficients
    %"Project" to positive
    load('SO3_fifteen.mat');
    [AI,AE,Acon] = linear_cons_B2(numel(B)-1,SO3);
    [B,Breal]    = project_B_to_positive2(AI,AE,Acon,B);
    scl = 1/B{1};
    for i=1:numel(B)
        B{i} =  scl*B{i};
    end

end

%% Simulate kspace function
% gridSize = 11;  % number of voxels in each dimenison. take odd.
gridSize = 15;  % number of voxels in each dimenison. take odd.

% Sum of Gaussian
sigma = 200;
% sigma = 100;
T     = (gridSize/5)*[0 0 0; 0.1 0.1 0; -0.2 0.1 0.2; -0.2 0 -0.2]';% ; 0.13 -0.2 0.1;0.1 0 -0.15]';

g     =  @(k) exp(1i*k'*T(:,1)).*exp(-pi^2/sigma*sum(k.*k)).'  + ...
    exp(1i*k'*T(:,2)).*exp(-pi^2/sigma*sum(k.*k)).'  + ...
    exp(1i*k'*T(:,3)).*exp(-pi^2/sigma*sum(k.*k)).' + ...
    exp(1i*k'*T(:,4)).*exp(-pi^2/sigma*sum(k.*k)).'; % + ...
%          exp(1i*k'*T(:,5)).*exp(-pi^2/sigma*4*sum(k.*k)).';
T     = (gridSize/5)*[0 0 0; 0.08 0.1 0; -0.1 0 0.1]';% ; 0.13 -0.2 0.1;0.1 0 -0.15]';
g     =  @(k) exp(1i*k'*T(:,1)).*exp(-pi^2/sigma*sum(k.*k)).'  + ...
    exp(1i*k'*T(:,2)).*exp(-pi^2/sigma*sum(k.*k)).'  + ...
    exp(1i*k'*T(:,3)).*exp(-pi^2/sigma*sum(k.*k)).' ;%+ ...
%% Vol on kspace over cartesian grid
radius   = floor(gridSize/2);
if mod(gridSize,2)==0
    x_1d = (-radius:1:radius-1);%/radius;   % - Even number of points
else
    x_1d = (-radius:1:radius);%/radius;   % - Odd number of points
end

% the 3d grid
[x_3d,y_3d,z_3d] = meshgrid(x_1d,x_1d,x_1d);
vec_3d_grid      = [x_3d(:),y_3d(:),z_3d(:)].'*pi/2;

% parameterization of the 2d slices
x_2d = x_3d(:,:,1);
y_2d = y_3d(:,:,1);

% evaluate the function
volf = g(vec_3d_grid);

% back to real domain
volk = reshape(volf, gridSize, gridSize, gridSize);
vol  = real(fftshift(ifftn(ifftshift(volk)))); % in real domain
vol  = vol*floor(size(x_2d,1)/2);   % Nir Sunday ?
% VisualVol(vol,'take a look')


%% calculate 3D expansion and gamma coefficients


fprintf('Calculating 3D coefficients...');
A = pswf_t_f_3d(vol, beta, delta);
fprintf('DONE \n');

fprintf('Calculating gamma coefficients...');
c     = beta*pi*radius;              % Nyquist bandlimit
gamma = PSWF_2D_3D_T_mat(c, delta, eps_p);
fprintf('DONE \n');

%% checking reconstruct volume from coefficients
vol_hat = pswf_t_b_3d(A, gridSize, c, delta);
% norm(vol(:)-vol_hat(:)) % norm(vol(:)-vol_hat(:))/norm(vol(:))

r    = sqrt(x_3d.^2 + y_3d.^2 + z_3d.^2);
ball = (r <= max(x_3d(:)));
% vol_supp = zeros(size(vol)); vol_supp(ball) = vol(ball); % figure; vol3d('cdata',real(vol_supp))

% figure; vol3d('cdata',vol);     title('Original Volume')
% figure; vol3d('cdata',vol_hat); title('Reconstructed Volume')
v1 = vol(ball);
v2 = vol_hat(ball);
disp(['3D vol reconstruction - relative squared error (inside the ball): ',num2str( (norm(v1(:)-v2(:)).^2)/norm(v1(:)).^2 ) ])


%% calculating moments
fprintf('Calculating moments: \n');


% preprocessing
L = max(gamma.band_idx_3d);  % the overall degree. Common to both 2D and 3D
P = size(B,1);               % expansion length of the distribution
M = max(gamma.ang_idx_2d)+1; % angular size of the moment
fprintf('Run preprocessing...');
tic
fprintf('Run preprocessing...');
tic
[sign_mat, Gamma_mat] = preprocessing_mu1_coef(P, L, gamma);
[C_tensor]            = preprocessing_mu2_coef_V2(gamma, P, M);
tt = toc;
fprintf('DONE in about %d seconds \n', round(tt));

fprintf('Moments calculation...');
% first moment
[ m1_true ] = FirstMoment_PSWF_v1_fixed(A, B, gamma);%Gamma_mat, sign_mat);
% [ m1_true ] = FirstMoment_PSWF_v2(A, B, Gamma_mat, sign_mat);

fprintf('First moment is DONE. \n');

% second moment
fprintf('Mu2 calculation...');
tic
[m2_true] = SecondMoment_PSWF_v2(A, B, gamma, C_tensor, Gamma_mat);
tt = toc;
fprintf('DONE in about %d seconds \n', round(tt));
%% prepare the comparison data

repeat_trial = 5;

fprintf('Preparing projections: \n');

% sampling SO(3))
fprintf('Uniform sampling...')
R = generate_SO3_uniform_array_v2(total_N);
fprintf('DONE \n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projection simulating
fprintf('projection simulating...');
sigma = 0;
[projs ,weight] = ComputeProjection(total_N, gridSize, [], B, g, x_2d, y_2d);
% if SNR ~= 0
%     [projs, ~, ~, sigma] = cryo_addnoise(projs, SNR,'gaussian');
% end

fprintf('DONE \n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% overall energy
norm_m1_square = norm(m1_true(:))^2;
norm_m2_square = norm(m2_true(:))^2;

% comparison points
number_of_check_points = 50;
start_N                = 5;
% check_points           = floor(linspace(start_N, total_N, number_of_check_points));
check_points           = floor(logspace(log10(start_N), log10(total_N), number_of_check_points));

%% moving to prolates coordinates -- define expansion parameters
fprintf('Moving generated projections to PSWF coordinates...');

[proj_PSWF] = Projection2PSWF(projs, beta, eps_p, gamma);

fprintf('...DONE \n')

%% average to attain the empirical moments
% weight = ones(total_N,1);
fprintf('Averaging...');
% initiliaze
relative_m1 = zeros(size(check_points));
relative_m2 = zeros(size(check_points));
absolute_m1 = zeros(size(check_points));
absolute_m2 = zeros(size(check_points));
m1_hat = 0;
m2_hat = 0;
% m1_true = m1_true(2:end,:);
for rt=1:repeat_trial
    % initialize inner variables
    check_counter = 1;
    inner_relative_m1 = zeros(size(check_points));
    inner_relative_m2 = zeros(size(check_points));
    inner_absolute_m1 = zeros(size(check_points));
    inner_absolute_m2 = zeros(size(check_points));
    perm_list = randperm(size(R,3));
    m1_hat = 0;
    m2_hat = 0;
    m1_hat1 = 0;
    m2_hat1 = 0;
    for i=1:total_N
%         % averaging
%         m1_hat1 = m1_hat1 + proj_PSWF(:,:,perm_list(i))*weight(perm_list(i));
%          current_vec_proj = proj_PSWF(:,:,perm_list(i));
%          m2_hat1 = m2_hat1 + current_vec_proj(:)*current_vec_proj(:)'*weight(perm_list(i));
%         [m1_hatnew, m2_hatnew] = EstimateMoments(proj_PSWF(:,:,perm_list(i)),...
%                                 weight(perm_list(i)), 1, gamma, sigma);
                            
%         m1_hat = m1_hat + proj_PSWF(:,:,perm_list(i))*weight(perm_list(i));
%         current_vec_proj = proj_PSWF(:,:,perm_list(i));
%         m2_hat = m2_hat + current_vec_proj(:)*current_vec_proj(:)'*weight(perm_list(i));                    
%         
%         m1_hat = m1_hat + m1_hatnew;
%         m2_hat = m2_hat + m2_hatnew;
        % check point control
        if ismember(i,check_points)
            [current_m1_check, current_m2_check] = EstimateMoments(proj_PSWF(:,:,perm_list(1:i)),...
                                weight(perm_list(1:i)), i, gamma, sigma);
                            
            % first moment
%             current_m1_check     = m1_hat/sum(weight(perm_list(1:i)));
%             current_m1_check = current_m1_check(2:end,:);
%             inner_absolute_m1(check_counter) = norm(m1_true(:)-current_m1_check(:))^2;
            inner_absolute_m1(check_counter) = norm(m1_true(:)-current_m1_check(:))^2;
            inner_relative_m1(check_counter) = inner_absolute_m1(check_counter)/norm_m1_square;
            
%             % second moment initialization
%             current_m2_check = zeros(max(gamma.ang_idx_2d)+1, nnz(gamma.ang_idx_2d==0),max(gamma.ang_idx_2d)+1, nnz(gamma.ang_idx_2d==0));
%             
%             %reshape the second moment and summarize
%             for m=0:max(gamma.ang_idx_2d)
%                 for k=0:nnz(gamma.ang_idx_2d==0)-1
%                     current_m2_check(m+1,k+1,:,:) = reshape(m2_hat(m+1 + (k)*(max(gamma.ang_idx_2d)+1),:),...
%                         max(gamma.ang_idx_2d)+1, nnz(gamma.ang_idx_2d==0) );
%                     current_m2_check(m+1,k+1,m+1,k+1) = current_m2_check(m+1,k+1,m+1,k+1) ...
%                                                         - sigma^2;
%                 end
%             end
%             current_m2_check = m2_hat;
%             current_m2_check = current_m2_check/sum(weight(perm_list(1:i)));% - sigma^2;
            inner_absolute_m2(check_counter) = norm(m2_true(:)-current_m2_check(:))^2;
            inner_relative_m2(check_counter) = inner_absolute_m2(check_counter)/norm_m2_square;
            check_counter = check_counter + 1;        
        end 
    end
    % adding current trial results
    relative_m1 = relative_m1 + (1/repeat_trial)*inner_relative_m1;
    absolute_m1 = absolute_m1 + (1/repeat_trial)*inner_absolute_m1;
    relative_m2 = relative_m2 + (1/repeat_trial)*inner_relative_m2;
    absolute_m2 = absolute_m2 + (1/repeat_trial)*inner_absolute_m2;    
    fprintf('*');
end

fprintf('...DONE \n');


%% loglog
figure; loglog(check_points/1000,relative_m1,'LineWidth',2.5); %legend('\mu_1')
xlabel('Number of samples (in thousands)'); ylabel('Square relative error'); %set(gca,'FontSize',22)
title('First moment')

grid on

figure; loglog(check_points/1000,relative_m2,'r','LineWidth',3); %legend('\mu_2');
xlabel('Number of samples (in thousands)'); ylabel('Square relative error'); %set(gca,'FontSize',22)
title('Second moment')

grid on
