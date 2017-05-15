function [ alpha_n ] = compute_Zernikecoefs( obj, Noll_list, pix_D, basis_set )
% [ alpha_n ] = compute_Zernikecoefs( obj, Noll_list, pix_D, basis_set )
%   
% Obj: The object to decompose into the Zernike Expansion
% Noll_list: A list of Zernike expansion functions to use given by Noll
%            index
% pix_D: Diameter of unit circle in pixels
% basis_set: Matrix of Zernike modes created by "Zernike_Basis.m"

if nargin < 4
    use_mat = false;
    
else
    use_mat = true;
    
    
end


if ~use_mat
    
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
        tmp = (basis)' .* obj;
        alpha_n(counter+1) = sum(tmp(:)) / sum(unit_circ(:));
        counter = counter+1;
        
%         clear tmp;
        
    end

else
    [Nzerns, ~] = size(basis_set);
    col = basis_set(1,:);
    mask = col(col~=0);
    mask = mask./mask;
    
    alpha_n = init_variable(Nzerns,1,1,'single',0);
    counter = 0;
    
    for ii = 1:Nzerns
        tmp = (basis_set(ii,:))' .* obj(:);
        alpha_n(counter+1) = sum(tmp) / sum(mask);
        counter = counter+1;
        
%         clear tmp;
    end
    
end

