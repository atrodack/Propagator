function [ CovMat ] = NollZernikeCovMat( Noll_list, r_0, R)
%[ CovMat ] = NollZernikeCovMat( Noll_list, r_0, R)
% 
% Analytical computation of a Zernike matrix representation of the
% Kolmogoroff phase spectrum via calculation of the covariance matrix in
% Fourier Space.
%
% Reference: Noll, Robert J. "Zernike polynomials and atmospheric
% turbulence", JOSA, 3 October, 1975
%
% Output:
% CovMat - The covariance matrix [length(Noll_list), length(Noll_list)]
%
% Inputs:
% Noll_list - Integer vector of Noll indices to be used
% r_0 - Fried correlation length
% R - Radius of circular aperture



    CovMat = zeros(length(Noll_list),length( Noll_list));
    n_ = zeros(length(Noll_list),1);
    m_ = n_;
    
    for jj = 1:length(Noll_list)
        [n_(jj),m_(jj)] = Noll(Noll_list(jj));
    end
    
    for ii = 1:length(Noll_list)
        for jj = 1:length(Noll_list)'
            n = n_(ii);
            n_prime = n_(jj);
            m = m_(ii);
            m_prime = m_(jj);
            
            % Kronecker Delta
            if m == m_prime
                delta_mm = 1;
            else
                delta_mm = 0;
            end
            
            t1 = (0.046/pi).*(R/r_0).^(5/3) .* sqrt((n+1).*(n_prime+1));
            t2 = (-1).^((n+n_prime-2*n)/2) .* delta_mm;
            Inn = NollBesselIntegral(n,n_prime);
            
            CovMat(ii,jj) = t1 .* t2 .* Inn;
            
        end
    end
    



end

