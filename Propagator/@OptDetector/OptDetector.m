classdef OptDetector < OptElement
    %OptDetector Class for Detector objects
    %   Placeholder for the end of a system
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        
        
    end
    
    methods
         %% Constructor
        function elem = OptDetector(PROPERTIES)
            % elem = OptDetector()
            % PROPERTIES is a 3x1 cell array containing element properties
            % If something is unknown or to be set later, leave that
            % position in A empty:
            %
%           PROPERTIES{1,1} = name                  [string]
%           PROPERTIES{2,1} = diameter              [float in m]
%           PROPERTIES{3,1} = zsag                  [float in m]
            
            if size(PROPERTIES) == [3,1]
                if iscell(PROPERTIES) == 1
                    A = PROPERTIES;
                else
                    error('PROPERTIES must be 3x1 cell array!');
                end
            else
                error('PROPERTIES must be a 3x1 cell array!');
            end
            
            elem.set_name(A{1,1});
            elem.set_diameter(A{2,1});
            elem.set_zsag(A{3,1});
            elem.set_z_position(0); %isn't used, but needed to make code happy
            elem.setdatatype();
            elem.set_propagation_method(1);
            elem.addnewline(2);

        end % of contructor
        
        function elem = set_isFocal(elem,code)
            % elem = set_isFocal
            % sets the propagation type for focusing:
            % code = 0 --> Fresnel
            % code = 1 --> Fourier Transform
            % code = 2 --> Zoom-FFT
            % code = 3 --> Convolution
            
            elem.isFocal_ = code;
        end % of set_isFocal
            
        function descr = describe(elem)
            % elem = describe(elem)

            elemtype = 'OptDetector';
            
            descr = sprintf('%s:\nElement is an %s',elem.name,elemtype);

        end
            
        
    end
    
end

