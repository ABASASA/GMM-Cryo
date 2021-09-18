function [projs ,weight] = ComputeProjection (total_N, B, g, x_2d, y_2d)
%% prepare the comparison data

fprintf('Preparing projections: \n');

% sampling SO(3))
fprintf('Uniform sampling...')
R = generate_SO3_uniform_array_v2(total_N);
fprintf('DONE \n');

% projection simulating
fprintf('projection simulating...');
[projs, weight] = make_projs_v2(g, R, B, x_2d, y_2d);

weight = real(weight);
weight(weight<0) = 0;
fprintf('DONE \n');

end
