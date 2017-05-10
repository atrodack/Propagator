function [ alpha_n ] = compute_Zernikecoefs( obj, Noll_list, pix_D )
% [ alpha_n ] = compute_Zernikecoefs( obj, Noll_list, pix_D )
%   Detailed explanation goes here


N = size(obj,1);
K = length(Noll_list);
alpha_n = init_variable(K,1,1,'single',0);


dx = 2/pix_D;
x = ((1:N) - (N/2))*dx;
[X,Y] = meshgrid(x);
[PHI,RHO] = cart2pol(X,Y);
unit_circ = single(RHO<=1);


counter = 0;
for ii = 1:K
    [n,m] = Noll(Noll_list(ii));
    tmp = Zernike2D(n,m,RHO,PHI);
    tmp = tmp.*unit_circ;
    rms_orig = rms(tmp(tmp~=0));
    basis = tmp * (1 / rms_orig); 
    tmp = conj(basis) .* obj;
    alpha_n(counter+1) = sum(tmp(:)) / sum(unit_circ(:));
    counter = counter+1;
    
    clear tmp;
    
end


end

