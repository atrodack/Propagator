function [ img_out,normalize_val ] = normalize_im( img_in )
%NORMALIZE_IM [ img_out ] = normalize_im( img_in )
%   Detailed explanation goes here

img_out = img_in;
normalize_val = max(max(abs(img_in)));
img_out = img_out ./ normalize_val;

end

