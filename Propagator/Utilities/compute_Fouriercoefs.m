function [ alpha_n ] = compute_Fouriercoefs( obj, L, endpoint,FoVr)
%[ alpha_n ] = compute_Fouriercoefs( obj, L, endpoint,FoVr)
%   perform inner products to compute Fourier coefficients
%   obj is the object (matrix)
%   L = square side length
%   endpoint = integer that determines how many modes are created
%   (ie 16 for 1024 modes, running from -15 to 16 for double index)
% FoVr = coordinate radius (probably always just L/2)


FoVsize = size(obj,1);
val = length(-endpoint+1:endpoint)^2;
alpha_n = zeros(val,1);

% Coordinates and scales
x = linspace(-FoVr,FoVr,FoVsize);
y = x;
[X,Y] = meshgrid(x,y);



counter = 0;
for k = -endpoint+1:endpoint
    for l = -endpoint+1:endpoint
        basis = (1/L) * exp(2*pi*1i*(k*X + l*Y)/L);
        tmp = conj(basis)'.*obj;
        alpha_n(counter+1) = sum(tmp(:))/FoVsize;
        counter = counter+1;
    end
end

end

