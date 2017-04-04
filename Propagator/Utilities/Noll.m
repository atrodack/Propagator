function [n, m] = Noll(j)
% [n,m] = Noll(j)
%
% Author: Jared Males / Rus Belikov
%
% this function returns the Zernike coefficients corresponding to Noll
% index j (http://en.wikipedia.org/wiki/Zernike_polynomials)

n = ceil(-1.5 + sqrt(0.25 + 2*j) - 1e-10);
% 1e-10 is to avoid wrong rounding due to numerical precision

jrem = j - (n.*(n+1)/2+1);
m = (jrem + mod(jrem,2).*abs(mod(n,2)-1) + abs(mod(jrem,2)-1).*mod(n,2)).*(-sign(mod(j,2)-0.5));
end
