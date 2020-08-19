function [vecBoolMoments, boolM1, boolM2] = TrimMoments_inplane(m1Size,...
                        m2Size, gamma, L, P)


boolM1 = true(m1Size);
boolM2 = true(m2Size);

%% Trim M1 zeros

for i = 1 : size(boolM1, 1) - 1
    boolM1(i + 1, nnz(gamma.ang_idx_2d==i) + 1 : end) = 0;
end

%% Trim M2
M = length(gamma.coeff);
for q1 = 1 : m2Size(1)
    for t1 = 1 : m2Size(2)
        if boolM1(q1,t1)
            boolM2(q1, t1, :, :) = boolM1;
        else
            boolM2(q1, t1, :, :) = 0;
        end
    end
end

%% Create index map
vecBoolMoments = [boolM1(:); boolM2(:)];
end