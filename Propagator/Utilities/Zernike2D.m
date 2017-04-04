function Z = Zernike2D(n,m,rho,phi)
% Z = Zernike2D(n,m,rho,phi)
%
% Author: Jared Males / Rus Belikov
%
% Return the Zernike mode Z_n^m(rho, phi)
% normalized so that R(1) = 1
% rho and phi can be 2D arrays
%
if m >= 0
    Z = Zernike(n,m,rho).*cos(m*phi);
else
    Z = Zernike(n,-m,rho).*sin(-m*phi);
end

end