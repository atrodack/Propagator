function [ obj ] = Evaluate_Fourier_Expansion( alphas, L, Nbasis, N)
% [ obj ] = Evaluate_Fourier_Expansion( alphas, L, endpoint, FoVsize)
% alphas: Fourier Coefficients for basis functions  
% L: square side length
% Nbasis: Number of basis functions to use
% N: Sampling lattice size


FoVr = L/2;
obj = zeros(N,N);

% Coordinates and scales
x = linspace(-FoVr,FoVr,N);
y = x;
[X,Y] = meshgrid(x,y);

endpoint = round((round(sqrt(Nbasis))-1)/2);
K = length(-endpoint:endpoint)^2;

% fprintf('Using %d Fourier Modes\n',K);

counter = 1;
for k = -endpoint:endpoint
    for l = -endpoint:endpoint
        basis = (1/L) * exp((2*pi*1i*(k*X + l*Y)/L));
        tmp = (((basis)).*alphas(counter)) / N;
        obj = obj + tmp';
        counter = counter+1;
    end
end


end

