function [ screen ] = GaussianSmoothing( PS, ind, SmoothSize, WindLag, Device, verbose )
% screen = GaussianSmoothing(PS,ind,scale,transition)
%   Function for applying a Gaussian smoothing to the screen
%
% Adapted from code written by Johanan L. Codona

if nargin < 6
    verbose = true;
elseif nargin < 5
    verbose = true;
    Device = 1;
end
nGPUs = gpuDeviceCount;

dx = PS.dx_;
Kwidth = 5*SmoothSize;
% SMOOTH = fspecial('gaussian',ceil(Kwidth/dx),SmoothSize/dx);
SMOOTH = fspecial('gaussian',Kwidth,SmoothSize);
SMOOTH = single(SMOOTH);
screen = PS.field_(:,:,ind);

if nGPUs > 0
    device = gpuDevice(Device);
    SMOOTH = gpuArray(SMOOTH);
    screen = gpuArray(screen);
    
    if verbose
    %     fprintf('***************************************************\n');
    %     fprintf('*         Now Using GPU %s        *\n',device.Name);
    %     fprintf('***************************************************\n\n');
    cprintf('comment','***************************************************\n')
    cprintf('comment','* ');
    cprintf('text','         Now Using GPU ');
    cprintf('-err','%s ',device.Name);
    cprintf('comment','      *\n')
    cprintf('comment','***************************************************\n\n')
    end
end

screen = gather(conv2(screen,SMOOTH,'same'));

if nargin>2
    screen = circshift(screen,WindLag); % Leaves edge artifacts :(
end

end

