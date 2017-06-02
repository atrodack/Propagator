function plotUtils(titl,xlab, ylab )
%PLOTUTILS Function to set figure parameters
%   Detailed explanation goes here

if nargin == 0
    axis xy;
    axis square;
    colormap(gray(256));
    colorbar;
    
elseif nargin == 1
    if length(titl) == 1
        title(titl);
    elseif length(titl) == 3
        title(titl(1),titl(2),title(3));
    end
    axis xy;
    axis square;
    colormap(gray(256));
    colorbar;
    
elseif nargin == 3
    if iscell(titl)
        if length(titl) == 1
            title(titl);
        elseif length(titl) == 3
            title(titl{1},titl{2},titl{3});
        end
    else
        title(titl)
    end
    
    if iscell(xlab)
        if length(xlab) == 1
            xlabel(xlab);
        elseif length(xlab) == 3
            xlabel(xlab{1},xlab{2},xlab{3});
        end
    else
        xlabel(xlab)
    end
    
   if iscell(ylab)
        if length(ylab) == 1
            ylabel(ylab);
        elseif length(ylab) == 3
            ylabel(ylab{1},ylab{2},ylab{3});
        end
    else
        ylabel(ylab)
   end
   
   axis xy;
    axis square;
    colormap(gray(256));
    colorbar;
    
end

end

