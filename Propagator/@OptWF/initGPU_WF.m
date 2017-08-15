function [ WF ] = initGPU_WF( WF,gpu2use )
%initGPU_WF Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    gpu2use = 1;
end

WF.nGPUs = gpuDeviceCount;

WF.DEVICES = cell(WF.nGPUs,1);

% Initialize Detected GPUs
if WF.nGPUs == 1
    WF.DEVICES{1} = gpuDevice(1);
    WF.useGPU = 1;  % initialize to use GPU if there is one
    
elseif WF.nGPUs == 2
    if gpu2use == 2
        WF.DEVICES{2} = gpuDevice(1);
        WF.DEVICES{1} = gpuDevice(2);
    else
        WF.DEVICES{2} = gpuDevice(2);
        WF.DEVICES{1} = gpuDevice(1);
    end
    
    WF.useGPU = 1;  % initialize to use GPU if there is one
    
elseif WF.nGPUS > 2
    fprintf('Only supports up to 2 GPUs\n');
    n = input('Please Select a GPU to use:  ');
    m = input('Please Select a second GPU to use:  ');
    WF.DEVICES{2} = gpuDevice(m);
    WF.DEVICES{1} = gpuDevice(n);
    WF.nGPUS = 2;
    WF.useGPU = 1;  % initialize to use GPU if there is one
    
else
    warning('GPU:noGPU','No GPU to use!');
    WF.useGPU = 0;
    %                 c = parcluster('local'); % build the 'local' cluster object
    %                 OS.nWorkers = c.NumWorkers; % get the number of CPU cores from it
end


    

end

