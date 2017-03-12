classdef OptMirror < OptElement
    %OptMirror Class for Mirrors
    %   Detailed explanation goes here
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        focalLength;
        
    end
    
    methods
            %% Constructor
        function elem = OptMirror(PROPERTIES)
            % elem = OptMirror(PROPERTIES)
            % PROPERTIES is a 6x1 cell array containing element properties
            % If something is unknown or to be set later, leave that
            % position in A empty:
            %
%           PROPERTIES{1,1} = name                  [string]
%           PROPERTIES{2,1} = focalLength           [float in m]       
%           PROPERTIES{3,1} = isFocal_              [code {see below}]
%           PROPERTIES{4,1} = z_position            [float in m]
%           PROPERTIES{5,1} = diameter              [float in m]
%           PROPERTIES{6,1} = zsag                  [file path or matrix]
            
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
            elem.set_material(2); %Mirror should only be allowed to be a mirror
            elem.set_focal_length(A{2,1});
            elem.set_isFocal(A{3,1});
            elem.set_z_position(A{4,1});
            elem.set_diameter(A{5,1});
            elem.set_zsag(A{6,1});
            elem.setdatatype();
            elem.addnewline(2);

        end % of contructor
        
        
        %% Methods for setting properties
        function elem = set_focal_length(elem,f)
            % elem = set_focal_length(elem,f)
            % sets the focal length of the element to f (in meters)
            % If the element is a mask, use this as a flag to set if the mask
            % is amplitude or phase mask (1 for amp, 0 for phase).
            
            elem.focalLength = f;
            
        end % of set_focal_length
        
        function elem = set_isFocal(elem,code)
            % elem = set_isFocal
            % sets the propagation type for focusing:
            % code = 0 --> Fresnel
            % code = 1 --> Fourier Transform
            % code = 2 --> Zoom-FFT
            % code = 3 --> Convolution
            
            elem.isFocal_ = code;
        end % of set_isFocal
        
        
        %% Methods for probing properties
        function fl = getFocalLength(elem)
            % fl = getFocalLength(elem)
            % returns the value stored in the focal_length_
            
            fl = elem.focal_length_;
        end % of getFocalLength
        
        function isFocal = getPropagationMethod(elem)
            % propagation_method = getPropagationMethod(elem)
            % returns the method to be used for propagation
            
            isFocal = elem.isFocal_;
            if isFocal == 0
                descr = 'Fresnel Propagation';
            elseif isFocal == 1
                descr = 'Fourier Transform for focusing';
            elseif isFocal == 2
                descr = 'Zoom Fourier Transform for focusing to FPM';
            elseif isFocal == 3
                descr = 'Using Convolution to Apply FPM';
            end
            if elem.verbose == 1
                fprintf('The Propagation Method is %s\n',descr);
            end
        end % of getPropagationMethod
        
        
         %% Propagation Methods
        % Should only be used by the propagation calls from OptSys
        
        function WFout = ApplyElement(elem,WFin,lambda,n0)
            % OS = ApplyElement(elem,WFin,lambda)
            % Applys element to current wavefront
            
            % Get number of wavelengths
            numLambdas = length(lambda);
            % Get wavelength dependent phase factors
            elem.set_phasefactor(lambda,n0);
            
            % Initialize WFout
            WFout = WFin; % done this way to preserve data type...
                          % actual values will be overwritten
            
            % Edit zsag_ to be same dimension as WFin
            if numLambdas ~= 1
                if size(elem.zsag_,3) ~= numLambdas % only do this if necessary
                    elem.Cubify(numLambdas);
                end
            end
            
            for ii = 1:numLambdas
                WFout(:,:,ii) = WFin(:,:,ii) .* exp(-1i * elem.phasefac_(ii) .* elem.zsag_(:,:,ii));
            end
            
            
        end % ApplyElement
                %% Utility methods
        
        function fnum = getFNumber(elem)
            % fnum = getFNumber(elem)
            % returns the f/# of the element
            
            fnum = elem.getFocalLength / elem.getDiameter;
        end % of getFNumber
        
        
        
        function descr = describe(elem)
            % elem = describe(elem)

            elemtype = 'OptMirror';
            zpos = abs(elem.z_position_);
            sz = size(elem.zsag_);
            
            descr = sprintf('%s:\nElement is an %s at %0.3f meters downstream\nGrid Size: [%d %d]',elem.name,elemtype,zpos, sz(1), sz(2));
        end % of describe
        
    end
    
end

