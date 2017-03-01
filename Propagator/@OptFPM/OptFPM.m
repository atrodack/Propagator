classdef OptFPM < OptElement
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        
    end
    
    methods
                %% Constructor
        function elem = OptFPM(PROPERTIES)
            % elem = OptFPM(PROPERTIES)
            % PROPERTIES is a 8x1 cell array containing element properties
            % If something is unknown or to be set later, leave that
            % position in A empty:
            %
%           PROPERTIES{1,1} = name
%           PROPERTIES{2,1} = material (0-10)
%           PROPERTIES{3,1} = amponly
%           PROPERTIES{4,1} = z_position (m)
%           PROPERTIES{5,1} = diameter (m)
%           PROPERTIES{6,1} = zsag (m) [file path or matrix]
            
            if size(PROPERTIES) == [6,1]
                if iscell(PROPERTIES) == 1
                    A = PROPERTIES;
                else
                    error('PROPERTIES must be 6x1 cell array!');
                end
            else
                error('PROPERTIES must be a 6x1 cell array!');
            end
            
            elem.set_name(A{1,1});
            elem.set_material(A{2,1});
            elem.set_amponly(A{3,1});
            elem.set_z_position(A{4,1});
            elem.set_diameter(A{5,1});
            elem.set_zsag(A{6,1});

        end % of contructor
        
        
        function elem = set_amponly(elem,flag)
            % elem = set_amponly(elem,flag)
            % set if the mask is amplitude or phase mask
%            (1 for amp, 0 for phase).
            
            elem.amponly = flag;
            
        end % of set_focal_length
    end
    
end

