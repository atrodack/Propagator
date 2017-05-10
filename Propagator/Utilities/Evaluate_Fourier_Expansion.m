function [ obj ] = Evaluate_Fourier_Expansion( alphas, L, endpoint, N)
% [ obj ] = Evaluate_Fourier_Expansion( alphas, L, endpoint, FoVsize)
%   L = square side length
%   endpoint = integer that determines how many modes are created
%   (ie 16 for 1024 modes, running from -15 to 16 for double index)

FoVr = L/2;
obj = zeros(N,N);

% Coordinates and scales
x = linspace(-FoVr,FoVr,N);
y = x;
[X,Y] = meshgrid(x,y);


counter = 1;
for k = -endpoint+1:endpoint
    for l = -endpoint+1:endpoint
        basis = (1/L) * exp((2*pi*1i*(k*X + l*Y)/L));
        tmp = (((basis)).*alphas(counter)) / N;
        obj = obj + tmp';
        counter = counter+1;
    end
end


end

