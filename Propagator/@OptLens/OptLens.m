classdef OptLens < OptElement
    %OptLens Class for lenses
    %   Detailed explanation goes here
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        focalLength;
        isFocal;
        
        
    end
    
    methods
            %% Constructor
        function elem = OptLens(PROPERTIES)
            % elem = OptLens(PROPERTIES)
            % PROPERTIES is a 8x1 cell array containing element properties
            % If something is unknown or to be set later, leave that
            % position in A empty:
            %
%           PROPERTIES{1,1} = name
%           PROPERTIES{2,1} = material (0-10)
%           PROPERTIES{3,1} = focalLength
%           PROPERTIES{4,1} = z_position (m)
%           PROPERTIES{5,1} = diameter (m)
%           PROPERTIES{6,1} = zsag (m) [file path or matrix]
%           PROPERITES{7,1} = isFocal
            
            if size(PROPERTIES) == [7,1]
                if iscell(PROPERTIES) == 1
                    A = PROPERTIES;
                else
                    error('PROPERTIES must be 7x1 cell array!');
                end
            else
                error('PROPERTIES must be a 7x1 cell array!');
            end
            
            elem.set_name(A{1,1});
            elem.set_material(A{2,1});
            elem.set_focal_length(A{3,1});
            elem.set_z_position(A{4,1});
            elem.set_diameter(A{5,1});
            elem.set_zsag(A{6,1});
            elem.set_isFocal(A{7,1});

        end % of contructor
        
        
        
        
        function elem = set_focal_length(elem,f)
            % elem = set_focal_length(elem,f)
            % sets the focal length of the element to f (in meters)
            % If the element is a mask, use this as a flag to set if the mask
            % is amplitude or phase mask (1 for amp, 0 for phase).
            
            elem.focalLength = f;
            
        end % of set_focal_length
        
        function elem = set_isFocal(elem,val)
            % elem = set_isFocal
            % sets the propagation type:
            % val = 0 --> Fresnel
            % val = 1 --> Fourier Transform
            % val = 2 --> Zoom-FFT
            elem.isFocal = val;
        end % of set_isFocal
        
        function fl = getFocalLength(elem)
            % fl = getFocalLength(elem)
            % returns the value stored in the focal_length_
            
            fl = elem.focal_length_;
        end % of getFocalLength
        
        function isFocal = getPropagationMethod(elem)
            % propagation_method = getPropagationMethod(elem)
            % returns the method to be used for propagation
            
            isFocal = elem.isFocal;
            if isFocal == 0
                descr = 'Fresnel Propagation';
            elseif isFocal == 1
                descr = 'Fourier Transform for focusing';
            elseif isFocal == 2
                descr = 'Subsampled Fourier Transform for focusing to FPM';
            end
            if elem.verbose == 1
                fprintf('The Propagation Method is %s\n',descr);
            end
        end % of getPropagationMethod
        
        %% Utilities
        
        function fnum = getFNumber(elem)
            % fnum = getFNumber(elem)
            % returns the f/# of the element
            
            fnum = elem.getFocalLength / elem.getDiameter;
        end % of getFNumber
        
    end
    
end

