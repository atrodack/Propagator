classdef OptWFS < matlab.mixin.Copyable
    %OPTWFS Base Wavefront Sensor Class
    %   Framework properties and methods for Wavefront Sensing:
    %   Children of this class will be actual devices because they are all
    %   so different from eachother. Most properties and methods therefore
    %   will be defined in the device classes.
    
    
    %% Public Properties
    properties(GetAccess='public',SetAccess='public')
        name; % give a name to the structure/system
        verbose = 1; % print extra info
        interpolate_method = 'cubic'; % selects a method of interpolation
    end % of public properties
    
    %% Protected Properties
    properties(GetAccess='public',SetAccess='private')
        
        % Data Type
        default_data_type = 'single'; % data_type of object properties
        
        % WFS Propterties
        estimatedWF_;
        
        % Grid Properties
        N_; %square grid size
        
        
        
    end % of protected properties
    
    
    
    %% Methods
    methods
        
        %% Set Property Methods
        function WFS = setGridSize(WFS,N)
            % WFS = setGridSize(WFS,N)
            WFS.N_ = N;
        end % of setGridSize
        
        function WFS = setEstimatedWF(WFS,EWF)
            %WFS = setEstimatedWF(WFS,EWF)
            WFS.estimatedWF_ = EWF;
        end % of setEstimatedWF
        
        %% Utility Methods
        function WFS = setdatatype(WFS,default_data_type)
            % WFS = datatype(default_data_type)
            % Sets the datatype to use. Do not do anything if already of
            % the correct data type.
            % Currently supported:
            % single
            % double
            % uint8
            
            if nargin < 2
                default_data_type = WFS.default_data_type;
            else
                WFS.default_data_type = default_data_type;
            end
            
            switch default_data_type
                case 'single'
                    if ~isa(WFS.estimatedWF_,'single')
                        WFS.estimatedWF_ = single(WFS.estimatedWF_);
                        if WFS.verbose == 1
                            fprintf('Data Type set to single\n');
                        end
                    end
                    
                case 'double'
                    if ~isa(WFS.estimatedWF_,'double')
                        WFS.estimatedWF_ = double(WFS.estimatedWF_);
                        if WFS.verbose == 1
                            fprintf('Data Type set to double\n');
                        end
                    end
                    
                case 'uint8'
                    if ~isa(WFS.estimatedWF_,'uint8')
                        WFS.estimatedWF_ = uint8(WFS.estimatedWF_);
                        if WFS.verbose == 1
                            fprintf('Data Type set to uint8\n');
                        end
                    end
                    
                otherwise
                    error('I do not understand that data type (yet)!');
            end
        end % of setdatatype
        
        
        %% Plotting Methods
        
        
    end % of methods
    
end

