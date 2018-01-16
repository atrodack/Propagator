function [ CovMat ] = Wang_Markey_CovMat( Noll_list, D, r0 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

 CovMat = zeros(length(Noll_list),length( Noll_list));
    n_ = zeros(length(Noll_list),1);
    m_ = n_;
    
    for jj = 1:length(Noll_list)
        [n_(jj),m_(jj)] = Noll(Noll_list(jj));
    end
    
    
    for i = 1:length(Noll_list)
        for i_p = 1:length(Noll_list)'
            n = n_(i);
            n_prime = n_(i_p);
            m = abs(m_(i));
            m_prime = abs(m_(i_p));
            
            % Kronecker Delta
            if m == m_prime
                K_delta = 1;
            else
                K_delta = 0;
            end
            
            expon = (n+n_prime-2*m)/2;
            
            % i-i' even
            if mod(i-i_p,2) == 0
                t1 = 0.0072 * (D/r0)^(5/3);
                t2 = (-1)^expon;
                t3 = sqrt((n+1)*(n_prime+1));
                t4 = pi^(8/3);
                delta_mm = K_delta;
                
                t5_num = gamma(14/3) * gamma((n+n_prime-(5/3))/2);
                t5_den1 = gamma((n-n_prime+(17/3))/2);
                t5_den2 = gamma((n_prime - n + (17/3))/2);
                t5_den3 = gamma((n+n_prime + (23/3))/2);
                
                t5 = t5_num / (t5_den1 * t5_den2 * t5_den3);
                
                CovMat(i,i_p) = t1 * t2 * t3 * t4 * delta_mm * t5;
                
            % i-i' odd    
            elseif mod(i-i_p,2) == 1
                CovMat(i,i_p) = 0;
            end
        end
    end





end

