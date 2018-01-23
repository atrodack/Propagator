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
%             pscale = params.pscale;
%             lambda = params.lambda;
%             propdist = params.propdist;
%             coeff = exp(1i*propdist*2*pi/lambda) / (1i*lambda*propdist);
            field_ = gather(fft2_back(field,(pscale .* pscale)^-1));
            % [centerPupil,~] = makeShiftPhase(sz(1),sz(3),'single');
            field_ = sz(2)^2.*field_.*tiltphase;
        else
%             pscale = params.pscale;
%             lambda = params.lambda;
%             propdist = params.propdist;
%             coeff = exp(1i*propdist*2*pi/lambda) / (1i*lambda*propdist);
            % [~,centerFocal] = makeShiftPhase(sz(1),sz(3),'single');
            field_ = gather(fft2_fwd(field.*tiltphase,(pscale .* pscale)));
        end
    end
elseif vec == true
%     if direction == -1
%         field_ = gather(fft2_back(field,(pscale .* pscale)^-1));
%         field_ = sz(2)^2.*field_.*tiltphase;
%     else
%         field_ = gather(fft2_fwd(field.*tiltphase,(pscale .* pscale)));
%     end
end



[amp,pha] = WFReIm2AmpPhase2(real(field_),imag(field_));
field_(:,:,ii) = amp(:,:,ii) .* exp(1i * pha(:,:,ii));

end % FraunhoferProp
