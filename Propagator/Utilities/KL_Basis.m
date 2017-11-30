function [ KL, sing_vals,cov_mat ] = KL_Basis( Noll_list, r_0, R,dx,N )
% [ KL ] = KL_Basis( Noll_list, r_0, R, N )
%  
% Construct a matrix KL whose columns are sampled Karhunen-Loeve modes for
% Kolmogoroff turbulence, as given by Noll in JOSA, 3 October 1975
%
% Outputs:
% KL - Karhunen-Loeve basis matrix
% sing_vals - singular values associated with modes
%
% Inputs:
% Noll_list - Integer vector of Noll indices to be used
% r_0 - Fried correlation length
% R - Radius of circular aperture
% dx - pixel spacing
% N - number of desired pixels per side of returned KL modes [N^2,1]
%     vectors
%
% TODO:
% Fix normalization?


% Turn on/off debugging features
DEBUG = false;

% Make a Zernike Basis
[Z,~] = Zernike_Basis(Noll_list,ones(length(Noll_list),1),(2*R/dx),N);
Z = double(Z.');

% Compute the analytical Zernike Covariance Matrix from [Noll, josa, 1975]
% cov_mat = NollZernikeCovMat(Noll_list,r_0,R);
cov_mat = Wang_Markey_CovMat(Noll_list,2*R,r_0);

% Eig Approach
eigen_vals = eig(cov_mat);
sing_vals = eigen_vals;
[V,D] = eig(cov_mat);
KL = Z*V*abs(D)^-0.5;

rms_vec = zeros(size(KL,2),1);
% Renormalize to unity rms
for ii = 1:size(KL,2)
    tmp = KL(:,ii);
    tmp2 = tmp(tmp~=0);
    rms_vec(ii) = rms(tmp2);
    KL(:,ii) = KL(:,ii) / rms_vec(ii);
end

if DEBUG
    for ii = 1:size(KL,2)
        tmp = KL(:,ii);
        tmp2 = tmp(tmp~=0);
        rms_vec(ii) = rms(tmp2);
    end
end

if DEBUG
    for ii = 1:size(KL,2)
        subplot(ceil(sqrt(length(Noll_list))),ceil(sqrt(length(Noll_list))),ii)
        imagesc(real(reshape(KL(:,ii),N,N)));
        axis square;
        title(sprintf('KL Mode %d',ii));
        colorbar;
%         drawnow;
    end
end


% % SVD Approach
% sing_vals = svd(cov_mat);
% [U,S,V] = svd(cov_mat);
% KL = Z*U*abs(S)^-0.5;
% 
% if DEBUG
%     for ii = 1:size(V,2)
%         subplot(ceil(sqrt(length(Noll_list))),ceil(sqrt(length(Noll_list))),ii)
%         imagesc(real(reshape(Z(:,ii),N,N)));
%         axis square;
%         title(sprintf('KL Mode %d',ii));
%         colorbar;
% %         drawnow;
%     end
% end






end

