function [projs ,weight] = ComputeProjectionFromAspire (vol, total_N, gridSize, paramsName, B)
%% prepare the comparison data

fprintf('Preparing projections: \n');

% sampling SO(3))
fprintf('Uniform sampling...')
R = generate_SO3_uniform_array_v2(total_N);
fprintf('DONE \n');

% projection simulating
fprintf('projection simulating...');

% projs = cryo_project_gaussian(paramsName, gridSize, 1, R); % uniform at the moment

fprintf('projection simulating...');
% vol    = cryo_gaussian_phantom_3d(paramsName,gridSize,1);

projs = cryo_project(vol, R); % uniform at the moment

weight = zeros(total_N,1);
% rho    = @(x) wignerD_expansion2(B,x);
rho    = @(x) B_to_Density_func(B,x);

for i=1:total_N
    currEul = rotm2eul(R(:,:,i));
    weight(i) = rho(currEul);

%     weight(i) = rho(rotm2eul(R(:,:,i)));
end

weight = real(weight);
weight(weight<0) = 0;
fprintf('DONE \n');

end
