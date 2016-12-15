classdef OptSys < matlab.mixin.Copyable
    %OptSys Class for storing optical system parameters
    %   Detailed explanation goes here
    
    
    %% Properties
    properties(Constant=true, GetAccess='protected')
        
        
    end % of constant properties
    
    properties(GetAccess='public',SetAccess='public')
        name; % give a name to the structure/system
        
        
    end % of public properties
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        numElements_ = 1; % number of elements in system
        f_number_; % F/# of system
        ELEMENTS_; % cell array of element cards
        
        
        % GPU properties
        useGPU; % use a GPU flag
        nGPUs;
        DEVICES;
        
        % Pool Characteristics
        nWorkers; % number of workers for parpool
        
        
        % Misc
        verbose; % print extra info
        interpolate_method;
        
    end % of protected properties
    
    %% Methods
    methods
        
        %% Constructor
        
        function OS = OptSys(numElements)
            % OS = OptSys(numElements)
            
            OS.numElements_ = numElements;
            OS.initOptSys;
            
            
        end % of contructor
        
        
        %% Set Properties
        
        function OS = initOptSys(OS)
            %OS = initializeOptSys(OS)
            
            OS.ELEMENTS_ = cell(OS.numElements_,1);
            OS.initGPU;
            
        end % of initOptSys
            
        function OS = initGPU(OS)
            
            % Count number of GPUs
            OS.nGPUs = gpuDeviceCount;
            OS.DEVICES = cell(OS.nGPUs,1);
            
            % Initialize Detected GPUs
            if OS.nGPUs == 1
                OS.DEVICES{1} = gpuDevice(1);
                OS.useGPU = 1;  % initialize to use GPU if there is one

            elseif OS.nGPUs == 2
                OS.DEVICES{1} = gpuDevice(1);
                OS.DEVICES{2} = gpuDevice(2);
                OS.nWorkers = 2;
                OS.useGPU = 1;  % initialize to use GPU if there is one
                
            elseif OS.nGPUS > 2
                fprintf('Only supports up to 2 GPUs\n');
                n = input('Please Select a GPU to use:  ');
                m = input('Please Select a second GPU to use:  ');
                OS.DEVICES{1} = gpuDevice(n);
                OS.DEVICES{2} = gpuDevice(m);
                OS.nGPUS = 2;
                OS.nWorkers = 2;
                OS.useGPU = 1;  % initialize to use GPU if there is one
                
            else
                warning('GPU:noGPU','No GPU to use!\n');
                OS.useGPU = 0;
                OS
            end
            

          
        end % of initGPU

        
    end % of methods
    
    %% Static Methods
    methods(Static = true)
        
        
    end % of static methods
end


%% IDEAS
% build a cell array of type Element objects
% use info from Element structures to simulate the optical system

% This class should do all the math, unless it is just a simple thing that
% an element should know

% Propagation type to be used currently stored in element
% This might change. All the actual propagation should occur here

% Element type and lens/mirror stored in element. This should help making
% separate propagation methods for certain items doable (if necessary)