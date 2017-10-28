function field_ = FraunhoferPropWF2(field, pscale, direction,tiltphase)
% F = FraunhoferProp(F, propdist, pscale, lambda)
% Monochromatic Calculation - doesn't not scale psf size with
% wavelength

if nargin < 3
    direction = 1;
    tiltphase = 1;
end

sz = size(field);
if length(sz) == 2
    if sz(2) == 1
        sz(1) = sqrt(sz(1));
        sz(2) = sz(1);
        vec = true;
    else
        vec = false;
    end
    sz(3) = 1;
end

field_ = init_variable(sz(1),sz(2),sz(3),'single',0);

if vec == 0
    for ii = 1:sz(3)
        if direction == -1
            field_ = gather(fft2_back(field,(206.29508972167969*pscale .* pscale)^-1));
%             [centerPupil,~] = makeShiftPhase(sz(1),sz(3),'single');
            field_ = field_.*tiltphase;
        else
%             [~,centerFocal] = makeShiftPhase(sz(1),sz(3),'single');
            field_ = gather(fft2_fwd(field.*tiltphase,(206.29508972167969*pscale .* pscale)));
        end
    end

end


[amp,pha] = WFReIm2AmpPhase2(real(field_),imag(field_));
field_(:,:,ii) = amp(:,:,ii) .* exp(1i * pha(:,:,ii));

end % FraunhoferProp