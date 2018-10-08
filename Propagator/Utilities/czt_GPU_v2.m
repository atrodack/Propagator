function g = czt_GPU_v2(x,k,info)
% g = czt_GPU_v2(x,k,info)
%
% see czt.m or czt_v2.m for detailed information. This is an equivalent
% function, but assumes GPU use to handle multiple matrix pages
% simultaneously. Supports multiple wavelengths through use
% of info struct.
%
% x: paged array where the third dimmension is used to store multiple
% fields to perform CZT on at once.
%
% info: struct containing x0, y0, dx0, dy0, dx1, dy1, f, and lambdaList
% x1 = out.x(1)
% y1 = out.y(1)
% dx0 = in.dx
% dy0 = in.dy
% dx1 = out.dx
% dy1 = out.dy
% f = focal length
% lambdaList = vector of wavelengths. If only one wavelength is needed, must be
%              a vector of length size(x,3) with each element equal to that wvl

%% Setup

% Compute CZT parameters
ax = exp(2*pi*1i * info.dx0 * info.x1 ./(info.lambdaList * info.f));
wx = exp(-2*pi*1i * info.dx0 * info.dx1 ./(info.lambdaList * info.f));
ay = exp(2*pi*1i * info.dy0 * info.y1 ./(info.lambdaList * info.f));
wy = exp(-2*pi*1i * info.dy0 * info.dy1 ./(info.lambdaList * info.f));

[m, n, z] = size(x); oldm = m;
if m == 1, x = x(:); [m, n] = size(x); end


%% ------- Length for power-of-two fft.

nfft = 2^nextpow2(m+k-1);

%% ------- Premultiply data.

kk = ( (-m+1):max(k-1,m-1) ).';
kk2 = (kk .^ 2) ./ 2;
nn = (0:(m-1))';
aax = zeros(m,1,z);
aay = aax;

for thisLambda = 1:z
    wwx(:,:,thisLambda) = wx(thisLambda) .^ (kk2);   % <----- Chirp filter is 1./ww
    wwy(:,:,thisLambda) = wy(thisLambda) .^ (kk2);


    aax_ = ax(1,thisLambda) .^ ( -nn );
    aay_ = ay(1,thisLambda) .^ ( -nn );
    aax(:,:,thisLambda) = aax_.*wwx(m+nn,thisLambda);
    aay(:,:,thisLambda) = aay_.*wwy(m+nn,thisLambda);
end

%% Construct DFT Matrices and send to GPU (leave double precision)
X = (1:nfft);
fX = X * (1/nfft);
dft_fwd_y = rot90(exp(1i*-1*2*pi* X.' * fX),2);
dft_fwd_y = gpuArray(dft_fwd_y);
dft_rev_y = conj(dft_fwd_y)/nfft;


%% Do the CZT in the y-direction
y = bsxfun(@times,x,aay);
% fft(y,nfft) pads with trailing zeros - replicate behavior and pad all
% pages at once
y = [y; zeros(nfft-size(y,1),size(y,2),size(y,3))];
v = ([1./wwy(1:(k-1+m),:,:)]);
v = gpuArray([v;zeros(nfft-size(v,1),1,size(v,3))]);
fy = pagefun(@mtimes,dft_fwd_y,y);
fv = pagefun(@mtimes,dft_fwd_y,v);
fy = bsxfun(@times,fy,fv);
g = pagefun(@mtimes,dft_rev_y,fy);
g = (permute(bsxfun(@times, g( m:(m+k-1), :, : ) , wwy( m:(m+k-1),ones(1, n),:)),[2 1 3]));


%% Do the CZt in the x-direction

% Redefine n to size things correctly
[m, n, z] = size(g);

y = bsxfun(@times,g,aax);
% fft(y,nfft) pads with trailing zeros
y = [y; zeros(nfft-size(y,1),size(y,2),size(y,3))];
v = ([1./wwy(1:(k-1+m),:,:)]);
v = gpuArray([v;zeros(nfft-size(v,1),1,size(v,3))]);
fy = pagefun(@mtimes,dft_fwd_y,y);
fv = pagefun(@mtimes,dft_fwd_y,v);
fy = bsxfun(@times,fy,fv);
g = pagefun(@mtimes,dft_rev_y,fy);
g = (permute(bsxfun(@times, g( m:(m+k-1), :, : ) , wwy( m:(m+k-1),ones(1, n),:)),[2 1 3]));

% Clean up
clear dft_fwd_y dft_rev_y y v fy fv;

end
