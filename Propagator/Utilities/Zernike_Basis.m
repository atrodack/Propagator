
function [Z,rms_] = Zernike_Basis(Noll_list,amplitude_list,pix_D,N)
% Z = Zernike_Basis(Noll_list,amplitude_list,pix_D,N)
% Make a Matrix whose columns are the Zernike basis functions given by
% Noll_list, with rms given in amplitude_list.
%
% pix_D is the diameter of your pupil in pixels
% N is the grid size of the matrix


if (length(Noll_list) ~= length(amplitude_list))
    error('Must provide amplitude for each basis function!');
end

K = length(Noll_list);

dx = 2/pix_D;
x = ((1:N) - (N/2))*dx;
[X,Y] = meshgrid(x);
[PHI,RHO] = cart2pol(X,Y);
unit_circ = single(RHO<=1);


Z = init_variable(K,N*N,1,'single',0);
rms_ = init_variable(K,1,1,'single',0);

for ii = 1:K
    [n,m] = Noll(Noll_list(ii));
    tmp = Zernike2D(n,m,unit_circ.*RHO,unit_circ.*PHI);
    rms_orig = rms(tmp(tmp~=0));
    tmp = tmp * (1 / rms_orig); 
    vec = amplitude_list(ii) * tmp(:);
    
    Z(ii,:) = vec;
    rms_(ii) = rms(vec(vec~=0));
    
    clear tmp;
    
    % DEBUG
%     tmp = reshape(Z(ii,:),[N,N]);
%     imagesc(tmp);
%     axis square; axis xy;
%     colorbar;
%     title(sprintf('k = %d',ii));
%     drawnow;
%     pause(0.3);
    
end



end