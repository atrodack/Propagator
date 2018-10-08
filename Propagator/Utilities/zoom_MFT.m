function [MFToutput, MFTdisp] = zoom_MFT(inputValues, inputMask, outputMask, zoom_factor, direct)
% Author:   Justin Knight
% Contact:  jknight@optics.arizona.edu
%
% Matrix Fourier transform to compute a custom sampled version of
% the observation plane specified by outputMask from a rectangle enclosing
% the set of input pixels (inputMask). This method is much faster than
% zoom_DFT.
%
% inputs:       inputValues - input matrix of the values we wish to apply
%               the DFT operation to
%
%               inputMask - input matrix of the binary mask over which we wish to
%               apply the DFT operation to
%
%               outputMask - mask describing the pixels over which we wish
%               to compute a higher sampled version of - we're zero-padding
%
%               zoom_factor - multiplier related to the increased number of
%               samples across the output coordinate plane
%
%               direct - direction of propagation: -1 for forward, +1 for
%               reverse
%
% outputs:      MFToutput - MFT of inputMask consisting of inputValues across the smallest rectangle
%               surrounding ROI from outputMask
%
%               MFTdisp - full scale display of MFT of inputValues over
%               outputMask

[inArrayrow, inArraycol] = find(inputMask);
x_min = min(inArraycol); x_max = max(inArraycol);
y_min = min(inArrayrow); y_max = max(inArrayrow);
[ysize, xsize] = size(inputMask);

[outArrayrow, outArraycol] = find(outputMask);
fx_min = min(outArraycol); fx_max = max(outArraycol);
fy_min = min(outArrayrow); fy_max = max(outArrayrow);
% numptsout = length(outArraycol);

% input coordinates
X = ((1/xsize) * ((x_min:x_max) - 1)) - 0.5;
Y = ((1/ysize) * ((y_min:y_max) - 1)) - 0.5;

% output coordinates
fX = (1/zoom_factor) * (((1/xsize) * ((fx_min:fx_max) - 1)) - 0.5) * xsize;
fY = (1/zoom_factor) * (((1/ysize) * ((fy_min:fy_max) - 1)) - 0.5) * ysize;

% [phasor(coord outer product) * smallest rectangle containing input mask * phasor(coord outer product)]
% = DFT across smallest rectangle containing the output mask.
MFToutput = exp(1i*direct*2*pi*fY.'*Y)*inputValues(y_min:y_max, x_min:x_max)*exp(1i*direct*2*pi*X.'*fX);

MFToutput = MFToutput/zoom_factor;

MFTdisp = zeros(ysize, xsize);

MFTdisp(fy_min:fy_max, fx_min:fx_max) = MFToutput;


end