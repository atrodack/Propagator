function [ psf ] = add_noise( psf,ReadNoise,useShotNoise, useReadNoise )
%ADD_NOISE
%   [ psf ] = add_noise( useShotNoise, useReadNoise )

if useShotNoise
    psf = double(psf);
    psf = double(psf)*(1e-12);
    psf = 1e12 * imnoise(psf,'poisson');
end
   
if useReadNoise
    psf = psf + (randn(size(psf)))*ReadNoise;
end
psf = single(abs(psf));

end

