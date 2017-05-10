function [ F ] = Fourier_Basis(L, endpoint,N )
%[ F ] = Fourier_Basis(  L, endpoint, N )
% Makes a matrix whose rows are normalized Fourier basis functions
% 
% L: side of the square expansion functions are inscribed in
%
% endpoint: integer value index determining highest frequency basis
%           function to make. See code below to see how it relates to the
%           creation of the basis functions
% 
% N: Sampling lattice size in pixels

DEBUG = false;

FoVr = L/2;
K = length(-endpoint+1:endpoint)^2;

% Coordinates and scales
x = linspace(-FoVr,FoVr,N);
y = x;
[X,Y] = meshgrid(x,y);

F = init_variable(K,N*N,1,'single',0);

counter = 1;
for k = -endpoint+1:endpoint
    for l = -endpoint+1:endpoint
        basis = (1/L) * exp(2*pi*1i*(k*X + l*Y)/L);
        F(counter, :) = basis(:);
        
        if DEBUG
            tmp = reshape(F(counter,:),[N,N]);
            figure(10000)
            subplot(1,2,1)
            imagesc(real(tmp));
            axis square; axis xy;
            colorbar;
            title(sprintf('k = %d real part',counter));
            subplot(1,2,2)
            imagesc(imag(tmp));
            axis square; axis xy;
            colorbar;
            title(sprintf('k = %d imaginary part',counter));
            drawnow;
            pause(0.3);
        end
        
        counter = counter+1;
                
    end
end


end

