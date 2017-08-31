function [ WF ] = GPUify_WF( WF, verbose )
%GPUIFY_WF Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    verbose = true;
end

% Check if GPU use is allowed
if(WF.useGPU ~=1)
    warning('GPU:CPUselected', 'GPU is not being used. One may not exist. If you know it does, run initGPU()');
    return;
end

% Make sure the field is of type float
WF.set_datatype('single');

if WF.nGPUs	 == 1
    % Send field to gpu
    WF.set_field(gpuArray(WF.field_));
    
    
    if verbose
        %     fprintf('***************************************************\n');
        %     fprintf('*         Now Using GPU %s        *\n',WF.DEVICES{1}.Name);
        %     fprintf('***************************************************\n\n');
        cprintf('comment','***************************************************\n')
        cprintf('comment','* ');
        cprintf('text','         Now Using GPU ');
        cprintf('-err','%s ',WF.DEVICES{1}.Name);
        cprintf('comment','      *\n')
        cprintf('comment','***************************************************\n\n')
    end
    %             warning('GPU:PROPNS','Propagation is currently pixel by pixel, and not a matrix multiply. This will be incredibly slow on GPU. Consider using CPU');
elseif WF.nGPUs == 2
    
    % Send field to gpu
    WF.set_field(gpuArray(WF.field_));
    
    if verbose
        %     fprintf('***************************************************\n');
        %     fprintf('*         Now Using GPU %s        *\n',WF.DEVICES{1}.Name);
        %     fprintf('***************************************************\n\n');
        cprintf('comment','***************************************************\n')
        cprintf('comment','* ');
        cprintf('text','         Now Using GPU ');
        cprintf('-err','%s ',WF.DEVICES{1}.Name);
        cprintf('comment','      *\n')
        cprintf('comment','***************************************************\n\n')
    end
    
end

end

