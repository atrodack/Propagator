function [ F ] = Fourier_Basis(L, Nbasis,N )
%[ F ] = Fourier_Basis(  L, Nbasis, N )
% Makes a matrix whose rows are normalized Fourier basis functions
% 
% L: side of the square expansion functions are created in (width of the
%    matrix in length)
%
% Nbasis: number of basis functions to create. Input should be a square
%         integer. If non-square integer is chosen, the closest square integer is
%         used.
%
% N: Sampling lattice size in pixels

DEBUG = false;

FoVr = L/2;
endpoint = round((round(sqrt(Nbasis))-1)/2);
K = length(-endpoint:endpoint)^2;

fprintf('Using %d Fourier Modes\n',K);

% Coordinates and scales
x = linspace(-FoVr,FoVr,N);
y = x;
[X,Y] = meshgrid(x,y);

F = init_variable(K,N*N,1,'single',0);

counter = 1;
for k = -endpoint:endpoint
    for l = -endpoint:endpoint
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

