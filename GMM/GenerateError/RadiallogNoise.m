function [sigma] = RadiallogNoise(gridSize, sigmaScalar)


%% Find middle
if mod(gridSize, 2) == 0
    middle = gridSize / 2;
else
    middle = (gridSize - 1) / 2;
end

%% compute 

sigma = zeros(gridSize, gridSize);

for i = 1 : gridSize
    for j = 1 : gridSize
        
        sigma(i,j) = sqrt( ((i-middle) / gridSize)^2 + ((j-middle) / gridSize)^2 );
    end
end

sigma = sigma * sigmaScalar;

end