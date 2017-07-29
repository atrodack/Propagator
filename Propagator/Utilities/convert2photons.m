function [ img_out ] = convert2photons( img_in, Nphotons )
%CONVERT2PHOTONS [ img_out ] = convert2photons( img_in, Nphotons )
%   Detailed explanation goes here

Nphotons = double(Nphotons);
img_in = double(img_in);
img_in = img_in ./ max(max(abs(img_in)));


total = sum(img_in(:));
tmp = img_in * ((1e-12) * (Nphotons/total));
img_out = (1e12) * imnoise(tmp,'poisson');


end

