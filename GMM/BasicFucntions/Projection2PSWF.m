function[proj_PSWF] = Projection2PSWF(projs, beta, eps_p, gamma)
%% moving to prolates coordinates -- define expansion parameters
fprintf('Moving generated projections to PSWF coordinates...');

[im_coef, p_struct] = pswf_t_f_2d(projs, floor(size(projs,1)/2), beta, eps_p, 1, []);

% accuracy check for the 2d
imagesP = pswf_t_b_2d( im_coef, p_struct, 1 );
%norm(images(:,:,1)-projs(:,:,1))/norm(projs(:,:,1))
norm(imagesP(:)-projs(:))/norm(projs(:))

proj_PSWF = zeros(max(gamma.ang_idx_2d)+1, nnz(gamma.ang_idx_2d==0), size(im_coef,2));
for m=0:max(gamma.ang_idx_2d)
    for k=0:nnz(gamma.ang_idx_2d==m)-1
        coeff2Dcurr_m = im_coef(p_struct.ang_freq==m,:);
        proj_PSWF(m+1,k+1,:) = coeff2Dcurr_m(k+1,:);
    end
end
fprintf('...DONE (imagesc to PSWF) \n')
end