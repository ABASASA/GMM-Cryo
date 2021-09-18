function [covaraince, cov_c] = Noise3on3GassianVec(gridSize,L, sigmaScalar)

%% Generate regular grid and basis functions (if needed)
x_1d_grid = -L:1:L;   % - Odd number of points
[x_2d_grid,y_2d_grid] = meshgrid(x_1d_grid,x_1d_grid);
r_2d_grid = sqrt(x_2d_grid.^2 + y_2d_grid.^2);
points_inside_the_circle = (r_2d_grid <= L);
%% Find middle
if mod(gridSize, 2) == 0
    middle = gridSize / 2;
else
    middle = (gridSize - 1) / 2;
end

% quad1 = floor(grisSize/4);
% quad3 = floor(3 * gridsize/4);
%% compute 
step = 0;

% indexes in circle
indexes = 1 : gridSize.^2;
indexes_circ = indexes(points_inside_the_circle);
% init
covaraince = zeros(gridSize.^2);
cov_c = zeros(length(indexes_circ));

for i = 1 : length(indexes_circ)
    i_image = indexes_circ(i);
    
    [row, col] = ind2sub([gridSize,gridSize], i_image);
    for irow = max(row-step,1) : min(row+step,gridSize)
        diff1 = irow - row;
       for icol = max(col-step,1) : min(col+step,gridSize) 
            diff2 = icol - col;
            if (diff1 == 0 & diff2 == 0)
                amplitude =  1;(1/sqrt(2) - sqrt(((row-middle) / gridSize)^2 + ((col-middle) / gridSize)^2 ));
            else
                amplitude = (0.2) *... %/ sqrt(diff1^2 + diff2^2)) * ...
                        (1/sqrt(2) - sqrt(((row-middle) / gridSize)^2 + ((col-middle) / gridSize)^2 ));    
            end
            amplitude = sigmaScalar^2 * amplitude;


            curr_ind_image = sub2ind([gridSize,gridSize], irow, icol);

            if ~points_inside_the_circle(curr_ind_image)
                continue;
            end
            covaraince(curr_ind_image, i_image) = amplitude;
            covaraince(i_image, curr_ind_image) = amplitude;
            
            curr_ind_circ = find(~(indexes_circ - curr_ind_image),1);
            cov_c(i, curr_ind_circ) = amplitude;
            cov_c(curr_ind_circ,i) = amplitude;
       end
    end
    
end
%cov_c = 0.35 * sigmaScalar * cov_c + 0.9* eye(size(cov_c)) * sigmaScalar;
%covaraince = 0.35 * sigmaScalar * covaraince + 0.9 * eye(size(covaraince)) * sigmaScalar;


cov_c = cov_c+  2*abs(min(eig(cov_c))) * eye(size(cov_c));
covaraince = covaraince + 2*abs(min(eig(covaraince))) * eye(size(covaraince));

%% Cut
% 
% covaraince = covaraince(:);
% cov_c = diag(covaraince(points_inside_the_circle));
% covaraince(~points_inside_the_circle) = 0;
% covaraince = diag(covaraince);

end