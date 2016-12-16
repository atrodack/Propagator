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
        numElements_ = 0; % number of elements in system
        f_number_; % F/# of system
        ELEMENTS_; % cell array of element cards
        
        conjugations_ = zeros(1,1);
        
        
        % GPU properties
        useGPU; % use a GPU flag
        nGPUs;
        DEVICES;
        
        % Pool Characteristics
        nWorkers; % number of workers for parpool
        
        
        % Misc
        verbose = 1; % print extra info
        interpolate_method = [];
        
    end % of protected properties
    
    %% Methods
    methods
        
        %% Constructor
        
        function OS = OptSys()
            % OS = OptSys()
            
            OS.initOptSys;
            
            
        end % of contructor
        
        
        %% Set Properties
        
        function OS = initOptSys(OS)
            %OS = initializeOptSys(OS)
            
            OS.ELEMENTS_ = cell(OS.numElements_,1);
            OS.nGPUs = gpuDeviceCount;
            if OS.nGPUs > 0
                OS.initGPU;
            end
            
        end % of initOptSys
        
        function OS = addElement(OS,elem)
            % OS = addElement(OS,elem)
            % adds and Element object to the optical system
            
            if isa(elem,'OptElement')
                OS.numElements_ = OS.numElements_ + 1;
                OS.ELEMENTS_{OS.numElements_} = elem;
                OS.setConjugations;
                
            else
                error('I do not understand anything but OptElement objects');
            end
        end % of addElement
                
        
        function conjugations = setConjugations(OS)
            % OS = setConjugations(OS)
            % Adds the zpos_ of an added element to the conjugations list
            
            conjugations = zeros(OS.numElements_,1);
            for z = 1:OS.numElements_
                conjugations(z,1) = OS.ELEMENTS_{z}.z_position_;
            end
            OS.conjugations_ = conjugations;
            
        end % of setConjugations
        
        function OS = getElementFNum(OS,element_num)
            % OS = getElementFNum(OS,element_num)
            % sets the f/# of the system to that of an element
            
            OS.f_number_ = OS.ELEMENTS_{element_num}.getFNumber;
        end % of getElementFNum
            
        %% Propagation Methods
        
        
        
        
        %% GPU Methods
        
        function OS = initGPU(OS)
            
            OS.DEVICES = cell(OS.nGPUs,1);
            
            % Initialize Detected GPUs
            if OS.nGPUs == 1
                OS.DEVICES{1} = gpuDevice(1);
                OS.useGPU = 1;  % initialize to use GPU if there is one

            elseif OS.nGPUs == 2
                OS.DEVICES{2} = gpuDevice(2);
                OS.DEVICES{1} = gpuDevice(1);
                OS.nWorkers = 2;
                OS.useGPU = 1;  % initialize to use GPU if there is one
                
            elseif OS.nGPUS > 2
                fprintf('Only supports up to 2 GPUs\n');
                n = input('Please Select a GPU to use:  ');
                m = input('Please Select a second GPU to use:  ');
                OS.DEVICES{2} = gpuDevice(m);
                OS.DEVICES{1} = gpuDevice(n);
                OS.nGPUS = 2;
                OS.nWorkers = 2;
                OS.useGPU = 1;  % initialize to use GPU if there is one
                
            else
                warning('GPU:noGPU','No GPU to use!\n');
                OS.useGPU = 0;
                c = parcluster('local'); % build the 'local' cluster object
                OS.nWorkers = c.NumWorkers % get the number of CPU cores from it
            end
            
        end % of initGPU

        
        function mem = GPUavailableMem(OS)
            % mem = GPUavailableMem(OS)
            % Returns the available memory on the active GPU
            
            mem = OS.DEVICES{1}.AvailableMemory;
                
            
        end % of GPUavailableMem
        
        function OS = useCPU(OS)
            % useCPU(OS)
            % Tell code to not use GPUs. Initialize CPU workers
            
            OS.useGPU = 0;
            c = parcluster('local'); % build the 'local' cluster object
            OS.nWorkers = c.NumWorkers % get the number of CPU cores from it
        end % useCPU    
        
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
% interpolation method is also set in element. this could change as well

% Element type and lens/mirror stored in element. This should help making
% separate propagation methods for certain items doable (if necessary)