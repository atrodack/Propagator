function [WFamp,WFphase] = WFReIm2AmpPhase(OS,WFreal,WFimag)
% [WFamp,WFphase] = WFReIm2AmpPhase(OS,WFreal,WFimag)
% Converts real and imaginary parts of WF to amplitude and
% phase.

if nargin < 2
    if isa(OS,'OptSys') == 0
        error('Argument must be an OptSys object');
    end
    sz = size(OS.WF);
    for ii = 1:sz(1)
        for jj = 1:sz(2)
            re = OS.WFreal(ii,jj);
            im = OS.WFimag(ii,jj);
            amp = sqrt(re*re + im*im);
            phase = atan2(im,re);
            OS.WFamp(ii,jj) = amp;
            OS.WFphase(ii,jj) = phase;
        end
    end
elseif nargin == 2 % Assume user isn't using an OptSys object, save input OS as WFamp, and WFamp as WFphase
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
end

end % of WFReIm2AmpPhase
