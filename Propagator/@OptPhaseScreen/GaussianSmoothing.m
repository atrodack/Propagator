function [ screen ] = GaussianSmoothing( PS, ind, SmoothSize, WindLag, Device )
% screen = GaussianSmoothing(PS,ind,scale,transition)
%   Function for applying a Gaussian smoothing to the screen
%
% Adapted from code written by Johanan L. Codona

if nargin < 5
    Device = 1;
end
nGPUs = gpuDeviceCount;

dx = PS.dx_;
Kwidth = 5*SmoothSize;
SMOOTH = fspecial('gaussian',ceil(Kwidth/dx),SmoothSize/dx);
screen = PS.field_(:,:,ind);

if nGPUs > 0
    device = gpuDevice(Device);
    SMOOTH = gpuArray(SMOOTH);
    screen = gpuArray(screen);
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

screen = gather(conv2(screen,SMOOTH,'same'));

if nargin>2
    screen = circshift(screen,WindLag); % Leaves edge artifacts :(
end

end

