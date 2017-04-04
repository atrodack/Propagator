function [ pad_im ] = zeropad( input, n )
%pad_im = zeropad(input,n)
%
% Author: Jared Males/Rus Belikov
%
% Pads an image "input" with zeros to increase the size of the matrix by a
% factor of n

N = length(input)/2;
z = zeros(N*(n-1));
zp = zeros(N*(n-1),N*2);

pad_im = [z zp z; zp' input zp'; z zp z];


end

