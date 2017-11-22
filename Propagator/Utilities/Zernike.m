function R = Zernike(n,m,x)
% R = Zernike(n,m,x)
%
% Author: Jared Males / Rus Belikov
%
% Return the Zernike Polynomial R_n^m of x
% normalized so that R(1) = 1
%
R=0;
if (n-m)/2==round((n-m)/2),
    for l=0:(n-m)/2,
        R=R+x.^(n-2*l)*(-1)^l*gamma(n-l+1)/(gamma(l+1)*gamma((n+m)/2-l+1)*gamma((n-m)/2-l+1));
    end
end
R(isnan(R)) = 0;

% At this point R(1) = 1, which is the standard zernike normalization. The
% next few lines normalize it to unity rms on the unit circle

R_norm=R*sqrt(2*(n+1));
if m == 0,
   R_norm=R_norm/sqrt(2);
end
R = R_norm*2.8884;

end