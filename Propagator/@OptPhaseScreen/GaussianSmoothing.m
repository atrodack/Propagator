function [ screen ] = GaussianSmoothing( PS, ind, SmoothSize, WindLag )
% screen = GaussianSmoothing(PS,ind,scale,transition)
%   Function for applying a Gaussian smoothing to the screen
%
% Adapted from code written by Johanan L. Codona


dx = PS.dx_;
Kwidth = 5*SmoothSize;
SMOOTH = fspecial('gaussian',ceil(Kwidth/dx),SmoothSize/dx);
screen = PS.field_(:,:,ind);
screen = conv2(screen,SMOOTH,'same');

if nargin>2
    screen = circshift(screen,WindLag); % Leaves edge artifacts :(
end

end

