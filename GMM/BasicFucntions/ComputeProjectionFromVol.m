function [projs ,weight] = ComputeProjectionFromVol (total_N, gridSize, vol, B)
%% prepare the comparison data

fprintf('Preparing projections: \n');

% sampling SO(3))
fprintf('Uniform sampling...')
R = generate_SO3_uniform_array_v2(total_N);
fprintf('DONE \n');

% projection simulating
fprintf('projection simulating...');


projs = cryo_project(vol, R); % uniform at the moment
% 
% %   Prepare weights for integration according to distribution
weight = zeros(total_N,1);
rho    = @(x) wignerD_expansion2(B,x);
for i=1:total_N
    weight(i) = rho(rotm2eul(R(:,:,i)));
end
% 
weight = real(weight);
weight(weight<0) = 0;
fprintf('DONE \n');

end
