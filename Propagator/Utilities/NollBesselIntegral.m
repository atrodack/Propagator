function [I_n_n_prime] = NollBesselIntegral(n,n_prime)
% [I_n_n_prime] = NollBesselIntegral(n,n_prime)
%
% Function for analytical evaluation of the integral term in Eq. 25 of
% Noll, Robert J. "Zernike Polynomials and atmospheric turbulence", JOSA, 3
% October, 1975 given in the Appendix of said paper
%
% Output:
% I_n_n_prime - evaluation of the integral
%
% Inputs:
% n - n index of Zernike polynomial i
% n_prime - n index of Zernike polynomial j



if n ~= 0 || n_prime ~= 0
    numerator = gamma(14/3) .* gamma((n+n_prime-(14/3)+3)/2);
    denominator1 = 2^(14/3) * gamma((-n+n_prime+(14/3)+1)/2);
    denominator2 = gamma((n-n_prime+(14/3)+1)/2);
    denominator3 = gamma((n+n_prime+(14/3)+3)/2);
    I_n_n_prime = numerator ./ (denominator1.*denominator2.*denominator3);
else
    I_n_n_prime = 0;

end
    
    
end