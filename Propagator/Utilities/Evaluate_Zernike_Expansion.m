function [ obj ] = Evaluate_Zernike_Expansion( alphas,Noll_list,pix_D, N )
%[ obj ] = Evaluate_Zernike_Expansion( alphas,Noll_list,pix_D ,N)
%
%   Function for computing the result of a zernike expansion. Takes a list
%   of Zernike coefficients, a list of Noll indicies, the diameter of the
%   unit circle in pixels, and the sampling lattice size.


K = length(alphas);
if ~isequal(K,length(Noll_list))
    error('Need a coefficient for each Zernike');
end


obj = init_variable(N,N,1,'single',0);

dx = 2/pix_D;
x = ((1:N) - (N/2))*dx;
[X,Y] = meshgrid(x);
[PHI,RHO] = cart2pol(X,Y);
unit_circ = single(RHO<=1);


counter = 1;
for ii = 1:K
    [n,m] = Noll(Noll_list(ii));
    tmp = Zernike2D(n,m,RHO,PHI);
    tmp = tmp.*unit_circ;
    rms_orig = rms(tmp(tmp~=0));
    basis = tmp * (1 / rms_orig); 
    tmp = (basis.*alphas(counter));
    obj = obj+tmp';
    counter = counter+1;
    
    clear tmp;
    
end

% obj = rot90(obj,2);
obj = obj';

end

