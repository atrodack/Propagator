function [ f_out ] = fft2_back( f_in, scale)
% [ f_out ] = fft2_fwd( f_in, scale )
%   Detailed explanation goes here




sz = size(f_in);
if length(sz) == 2
    sz(3) = 1;
end


if isa(f_in,'gpuArray')
    flag = true;
    tmp = gather(f_in);
    datatype = class(tmp);
    clear tmp;
else
    datatype = class(f_in);
    flag = false;
end

f_out = init_variable(sz(1),sz(2),sz(3),datatype,0);
if flag
    f_out = gpuArray(f_out);
end

for ii = 1:size(f_in,3)
        f_out(:,:,ii) = ifftshift(ifft2(ifftshift(f_in(:,:,ii))))*scale;
end


end

