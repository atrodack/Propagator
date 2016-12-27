function [WFreal,WFimag] = WFAmpPhase2ReIm(OS,WFamp,WFphase)
% [WFreal,WFimag] = WFAmpPhase2ReIm(OS,WFamp,WFphase)
% sets the real and imaginary parts of the wavefront from
% amplitude and phase components. If no amplitude and phase
% components are given, act on stored ones, and store resutls

if nargin < 2
    if isa(OS,'OptSys') == 0
        error('Argument must be an OptSys object');
    end
    sz = size(OS.WF);
    for ii = 1:sz(1)
        for jj = 1:sz(2)
            amp = OS.WFamp(ii,jj);
            phase = OS.WFphase(ii,jj);
            OS.WFreal(ii,jj) = amp*cos(phase);
            OS.WFimag(ii,jj) = amp*sin(phase);
        end
    end
elseif nargin == 2 % Assume user isn't using an OptSys object, save input OS as WFamp, and WFamp as WFphase
    WFphase = WFamp;
    WFamp = OS;
    sz = size(WFamp);
    WFreal = zeros(sz(1),sz(2));
    WFimag = zeros(sz(1),sz(2));
    for ii = 1:sz(1)
        for jj = 1:sz(2)
            amp = WFamp(ii,jj);
            phase = WFphase(ii,jj);
            WFreal(ii,jj) = amp*cos(phase);
            WFimag(ii,jj) = amp*sin(phase);
        end
    end
elseif nargin == 3
    sz = size(WFamp);
    WFreal = zeros(sz(1),sz(2));
    WFimag = zeros(sz(1),sz(2));
    for ii = 1:sz(1)
        for jj = 1:sz(2)
            amp = WFamp(ii,jj);
            phase = WFphase(ii,jj);
            WFreal(ii,jj) = amp*cos(phase);
            WFimag(ii,jj) = amp*sin(phase);
        end
    end
end
end % of WF2ReIm