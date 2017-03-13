function [WFamp,WFphase] = WFReIm2AmpPhase(OS,WFreal,WFimag)
% [WFamp,WFphase] = WFReIm2AmpPhase(OS,WFreal,WFimag)
% Converts real and imaginary parts of WF to amplitude and
% phase.


if nargin == 2 % Assume user isn't using an OptSys object, save input OS as WFamp, and WFamp as WFphase
    WFimag = WFreal;
    WFreal = OS;
    sz = size(WFreal);
    WFamp = zeros(sz(1),sz(2));
    WFphase = zeros(sz(1),sz(2));
    for ii = 1:sz(1)
        for jj = 1:sz(2)
            re = WFreal(ii,jj);
            im = WFimag(ii,jj);
            amp = sqrt(re*re + im*im);
            phase = atan2(im,re);
            WFamp(ii,jj) = amp;
            WFphase(ii,jj) = phase;
        end
    end
elseif nargin == 3
    sz = size(WFreal);
    WFamp = zeros(sz(1),sz(2));
    WFphase = zeros(sz(1),sz(2));
    for ii = 1:sz(1)
        for jj = 1:sz(2)
            re = WFreal(ii,jj);
            im = WFimag(ii,jj);
            amp = sqrt(re*re + im*im);
            phase = atan2(im,re);
            WFamp(ii,jj) = amp;
            WFphase(ii,jj) = phase;
        end
    end
    
    % set OS
    OS.WFamp = WFamp;
    OS.WFphase = WFphase;
end

end % of WFReIm2AmpPhase
