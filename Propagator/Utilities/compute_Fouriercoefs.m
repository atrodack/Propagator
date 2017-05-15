function [ alpha_n ] = compute_Fouriercoefs( obj, L, Nbasis)
%[ alpha_n ] = compute_Fouriercoefs( obj, L, Nbasis)
%   perform inner products to compute Fourier coefficients
%
%   obj: the object (matrix) to find the coefficients of
%   L: square side length
%   Nbasis: Number of Basis functions to use. Nearest square number to what
%   is asked for is used.

FoVr = L/2;
FoVsize = size(obj,1);


endpoint = round((round(sqrt(Nbasis))-1)/2);
K = length(-endpoint:endpoint)^2;

fprintf('Using %d Fourier Modes\n',K);

alpha_n = zeros(K,1);

% Coordinates and scales
x = linspace(-FoVr,FoVr,FoVsize);
y = x;
[X,Y] = meshgrid(x,y);



counter = 0;
for k = -endpoint:endpoint
    for l = -endpoint:endpoint
        basis = (1/L) * exp(2*pi*1i*(k*X + l*Y)/L);
        tmp = (basis)'.*obj;
        alpha_n(counter+1) = sum(tmp(:))/FoVsize;
        counter = counter+1;
    end
end

end

