function [WFamp,WFphase] = WFReIm2AmpPhase2(OS,WFreal,WFimag)
% [WFamp,WFphase] = WFReIm2AmpPhase2(OS,WFreal,WFimag)
% Converts real and imaginary parts of WF to amplitude and
% phase.

if nargin < 2
    if isa(OS,'OptSys') == 0
        error('Argument must be an OptSys object');
    end

    amp = sqrt(OS.WFreal.*OS.WFreal + OS.WFimag.*OS.WFimag);
    phase = atan2(OS.WFimag,OS.WFreal);
    OS.WFamp = amp;
    OS.WFphase = phase;
    
elseif nargin == 2 % Assume user isn't using an OptSys object, save input OS as WFamp, and WFamp as WFphase
    WFimag = WFreal;
    WFreal = OS;
    
    WFamp = sqrt(WFreal.*WFreal + WFimag.*WFimag);
    WFphase = atan2(WFimag,WFreal);
    
elseif nargin == 3
    WFamp = sqrt(WFreal.*WFreal + WFimag.*WFimag);
    WFphase = atan2(WFimag,WFreal);
end

end % of WFReIm2AmpPhase2

