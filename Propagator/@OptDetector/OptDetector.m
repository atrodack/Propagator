classdef OptDetector < OptElement
    %OptDetector Class for Detector objects
    %   Placeholder for the end of a system no longer with chirp-z
    %   transform!
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        useNoise_;
        exposure_time_;
        
        pixel_size_; % [lambda / D]
        FPregion_; % [lambda / D], -FPregion<theta<FPregion
        
    end
    
    methods
         %% Constructor
        function elem = OptDetector(PROPERTIES)
            % elem = OptDetector(PROPERTIES)
            % PROPERTIES is a 7x1 cell array containing element properties
            % If something is unknown or to be set later, leave that
            % position in A empty:
            %
%           PROPERTIES{1,1} = name                  [string]
%           PROPERTIES{2,1} = diameter              [float in m]
%           PROPERTIES{3,1} = zsag                  [float in m]
%           PROPERTIES{4,1} = FPregion              [float in l/d]
%           PROPERTIES{5,1} = pixel_size            [float in l/d]
%           PROPERTIES{6,1} = exposure_time         [float in seconds]
%           PROPERTIES{7,1} = useNoise              [Boolean]
            
            if size(PROPERTIES) == [5,1]
                if iscell(PROPERTIES) == 1
                    A = PROPERTIES;
                else
                    error('PROPERTIES must be 7x1 cell array!');
                end
            else
                error('PROPERTIES must be a 7x1 cell array!');
            end
            
            elem.set_name(A{1,1});
            elem.set_diameter(A{2,1});
            elem.set_zsag(A{3,1});
            elem.set_FPregion(A{4,1});
            elem.set_pixel_size(A{5,1});
            elem.set_z_position(0); %isn't used, but needed to make code happy
            elem.setdatatype();
            elem.set_propagation_method(1);
            elem.addnewline(2);

        end % of contructor
        
        function elem = set_FPregion(elem, nld)
            % elem = set_FPregion(elem, nld)
            elem.FPregion_ = nld;
            
        end % of set_FPregion
        
        function elem = set_pixel_size(elem,M)
            % elem = set_pixel_size(elem,M)
            
            elem.pixel_size_ = M;
        end % of set_pixel_size
        
        function elem = set_noiseFlag(elem,flag)
            % elem = set_noiseFlag(elem,flag)
            elem.useNoise_ = flag;
            
        end % of set_noiseFlag
        
        function elem = set_exposureTime(elem,dt)
            % elem = set_exposureTime(elem,dt)
            elem.exposure_time_ = dt;
            
        end % of set_exposureTime
            
        function descr = describe(elem)
            % elem = describe(elem)

            elemtype = 'OptDetector';
            
            descr = sprintf('%s:\nElement is an %s',elem.name,elemtype);

        end
            
        
    end
    
end

