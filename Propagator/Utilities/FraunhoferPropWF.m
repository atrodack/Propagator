function field_ = FraunhoferPropWF(field, propdist, pscale, lambda,direction,photonz)
% F = FraunhoferProp(F, propdist, pscale, lambda)
% Monochromatic Calculation - doesn't not scale psf size with
% wavelength

if nargin < 5
    direction = 1;
    photonz = false;
elseif nargin < 6
    photonz = false;
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

Nphotons_in_input_field = zeros(1,sz(3));
if photonz
    for jj = 1:sz(3)
        amp(:,:,jj) = sqrt(real(field(:,:,jj)).*real(field(:,:,jj)) + imag(field(:,:,jj)).*imag(field(:,:,jj)));
        Nphotons_in_input_field(1,jj) = sum(sum(amp(:,:,jj)));
    end
end

k = (2*pi)./lambda;

field_ = init_variable(sz(1),sz(2),sz(3),'single',0);

if vec == 0
    for ii = 1:length(lambda)
        coeff = ((exp(1i.*k(ii)*propdist))/(1i.*(lambda(ii))*propdist));
        if direction == -1
            field_(:,:,ii) = gather(coeff * ifftshift(ifft2(ifftshift(OptSys.halfpixelShift( field(:,:,ii),1 ) ))) .* (pscale .* pscale));
        else
            field_(:,:,ii) = gather(coeff * fftshift(fft2(fftshift(OptSys.halfpixelShift( field(:,:,ii),1 ) ))) .* (sz(1).* sz(2) .* pscale .* pscale));
%             field_(:,:,ii) = gather(coeff * fftshift(fft2(fftshift(( field(:,:,ii) ) ))) .* (sz(1).* sz(2) .* pscale .* pscale));

        end
    end

else
    for ii = 1:length(lambda)
        coeff = ((exp(1i.*k(ii)*propdist))/(1i.*(lambda(ii))*propdist));
        if direction == -1
            field_(:,:,ii) = gather((coeff * ifftshift(ifft2(ifftshift(OptSys.halfpixelShift( reshape(field(:,:,ii),sz(1),sz(2)),1 ) ))) .* (pscale .* pscale)));
        else
            field_(:,:,ii) = gather(coeff * fftshift(fft2(fftshift(OptSys.halfpixelShift( reshape(field(:,:,ii),sz(1),sz(2)),1 ) ))) .* (sz(1).* sz(2) .* pscale .* pscale));
%             field_(:,:,ii) = gather(coeff * fftshift(fft2(fftshift(( reshape(field(:,:,ii),sz(1),sz(2)) ) ))) .* (sz(1).* sz(2) .* pscale .* pscale));
                        
        end
    end
end


[amp,pha] = WFReIm2AmpPhase2(real(field_),imag(field_));
for ii = 1:size(amp,3)
    if photonz
        amp(:,:,ii) = (amp(:,:,ii) / sum(sum(amp(:,:,ii)))) * Nphotons_in_input_field(1,ii);
    end
    field_(:,:,ii) = amp(:,:,ii) .* exp(1i * pha(:,:,ii));
end

end % FraunhoferProp