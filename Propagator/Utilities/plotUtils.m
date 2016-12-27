function plotUtils(titl,xlab, ylab )
%PLOTUTILS Function to set figure parameters
%   Detailed explanation goes here

if nargin == 0
    axis xy;
    axis square;
    colormap(gray(256));
    colorbar;
    
elseif nargin == 1
    title(titl);
    axis xy;
    axis square;
    colormap(gray(256));
    colorbar;
    
elseif nargin == 3
    title(titl);
    xlabel(xlab);
    ylabel(ylab);
    axis xy;
    axis square;
    colormap(gray(256));
    colorbar;
    
end

end

