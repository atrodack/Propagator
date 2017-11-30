
function [Z,rms_] = Zernike_Basis(Noll_list,rms_list,pix_D,N)
% Z = Zernike_Basis(Noll_list,rms_list,pix_D,N)
% Make a Matrix whose rows are the Zernike basis functions given by
% Noll_list, with rms given in rms_list.
%
% pix_D is the diameter of your pupil in pixels
% N is the grid size of the matrix

DEBUG = false;


if (length(Noll_list) ~= length(rms_list))
    error('Must provide rms for each basis function!');
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
    tmp = Zernike2D(n,m,RHO,PHI);
    tmp = tmp .* unit_circ;
    
    vec = rms_list(ii) * tmp(:);
    
    Z(ii,:) = vec;
    rms_(ii) = rms(vec(vec~=0));
    %     rms_(ii) = rms(vec);
    
    clear tmp;
    
    if DEBUG
        tmp = reshape(Z(ii,:),[N,N]);
        figure(10000)
        imagesc(tmp);
        axis square; axis xy;
        colorbar;
        title(sprintf('k = %d',ii));
        drawnow;
        pause(0.3);
    end
end



end