function [ f ] = FastInverseFourierT(x,dx, mrows, ncols)
%FastInverseFourierT Two-dimensional discrete Fourier Transform with built in fftshifts.
%   FastFourierT(X,dx) returns the two-dimensional Fourier transform of matrix X.
%   If X is a vector, the result will have the same orientation.
%
%   FastInverseFourierT(X,dx,MROWS,NCOLS) pads matrix X with zeros to size MROWS-by-NCOLS
%   before transforming.
%
%   Class support for input X: 
%      float: double, single
[N,M] = size(x);
if nargin < 3
    if length(dx) == 2
        dy = dx(1);
        dx = dx(2);
    else
        dy = dx;
    end
    f = ifftshift(ifft2(ifftshift(x)))*dx*N*M*dy;
    
elseif nargin == 4
    if length(dx)==2
        dy = dx(1);
        dx = dx(2);
    else
        dy = dx;
    end
    f = ifftshift(ifft2(ifftshift(x),mrows,ncols))*dx*N*M*dy;


end