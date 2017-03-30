classdef OptDetector < OptElement
    %OptDetector Class for Detector objects
    %   Placeholder for the end of a system
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        noisetype_;
        FPregion_; % [lambda / D], -FPregion<theta<FPregion
        M_;
        
    end
    
    methods
         %% Constructor
        function elem = OptDetector(PROPERTIES)
            % elem = OptDetector()
            % PROPERTIES is a 4x1 cell array containing element properties
            % If something is unknown or to be set later, leave that
            % position in A empty:
            %
%           PROPERTIES{1,1} = name                  [string]
%           PROPERTIES{2,1} = diameter              [float in m]
%           PROPERTIES{3,1} = zsag                  [float in m]
%           PROPERTIES{4,1} = FPregion              [float in l/d]
            
            if size(PROPERTIES) == [5,1]
                if iscell(PROPERTIES) == 1
                    A = PROPERTIES;
                else
                    error('PROPERTIES must be 5x1 cell array!');
                end
            else
                error('PROPERTIES must be a 5x1 cell array!');
            end
            
            elem.set_name(A{1,1});
            elem.set_diameter(A{2,1});
            elem.set_zsag(A{3,1});
            elem.set_FPregion(A{4,1});
            elem.set_M(A{5,1});
            elem.set_z_position(0); %isn't used, but needed to make code happy
            elem.setdatatype();
            elem.set_propagation_method(1);
            elem.addnewline(2);

        end % of contructor
        
        function elem = set_FPregion(elem, nld)
            % elem = set_FPregion(elem, nld)
            elem.FPregion_ = nld;
            
        end % of set_FPregion
        
        function elem = set_M(elem,M)
            % elem = set_M(elem,M)
            
            elem.M_ = M;
        end % of set_M
            
        function descr = describe(elem)
            % elem = describe(elem)

            elemtype = 'OptDetector';
            
            descr = sprintf('%s:\nElement is an %s',elem.name,elemtype);

        end
            
        
    end
    
end

