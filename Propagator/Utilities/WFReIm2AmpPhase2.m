function [WFamp,WFphase] = WFReIm2AmpPhase2(OS,WFreal,WFimag)
% [WFamp,WFphase] = WFReIm2AmpPhase2(OS,WFreal,WFimag)
% Converts real and imaginary parts of WF to amplitude and
% phase.


if nargin == 2 % Assume user isn't using an OptSys object, save input OS as WFamp, and WFamp as WFphase
    WFimag = WFreal;
    WFreal = OS;    
    numLambdas = size(WFreal,3);
    
    WFamp = zeros(size(WFreal,1),size(WFreal,2),numLambdas);
    WFphase = WFamp;

    for jj = 1:numLambdas
        WFamp(:,:,jj) = sqrt(WFreal(:,:,jj).*WFreal(:,:,jj) + WFimag(:,:,jj).*WFimag(:,:,jj));
        WFphase(:,:,jj) = atan2(WFimag(:,:,jj),WFreal(:,:,jj));
    end
    
elseif nargin == 3
    numLambdas = size(WFreal,3);
    
    WFamp = zeros(size(WFreal,1),size(WFreal,2),numLambdas);
    WFphase = WFamp;

    for jj = 1:numLambdas
        WFamp(:,:,jj) = sqrt(WFreal(:,:,jj).*WFreal(:,:,jj) + WFimag(:,:,jj).*WFimag(:,:,jj));
        WFphase(:,:,jj) = atan2(WFimag(:,:,jj),WFreal(:,:,jj));
    end
    
    %Store result in OS
    OS.WFamp = WFamp;
    OS.WFphase = WFphase;
end

end % of WFReIm2AmpPhase2

