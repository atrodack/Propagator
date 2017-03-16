classdef OptFPM < OptElement
    %OptFPM Class for Focal Plane Masks
    %   Detailed explanation goes here
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        amponly = false;
        
    end
    
    methods
                %% Constructor
        function elem = OptFPM(PROPERTIES)
            % elem = OptFPM(PROPERTIES)
            % PROPERTIES is a 7x1 cell array containing element properties
            % If something is unknown or to be set later, leave that
            % position in A empty:
            %
%           PROPERTIES{1,1} = name                  [string]
%           PROPERTIES{2,1} = material              [0-10]       
%           PROPERTIES{3,1} = amponly               [bool]
%           PROPERTIES{4,1} = isFocal_              [0-3]
%           PROPERTIES{5,1} = z_position            [float in m]
%           PROPERTIES{6,1} = diameter              [float in m]
%           PROPERTIES{7,1} = mask                  [file path or matrix]
            
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
            elem.set_amponly(A{3,1});
            elem.set_propagation_method(A{4,1});
            elem.set_z_position(A{5,1});
            elem.set_diameter(A{6,1});
            elem.set_zsag(A{7,1});
            elem.setdatatype();
            elem.addnewline(2);
            
        end % of contructor
        
        %% Methods for setting properties
        
        
        function elem = set_amponly(elem,flag)
            % elem = set_amponly(elem,flag)
            % set if the mask is amplitude or phase mask
%            (1 for amp, 0 for phase).
            
            elem.amponly = flag;
            
        end % of set_focal_length
        
        
        %% Propagation Methods
        % Should only be used by the propagation calls from OptSys
        
        
        function WFout = ApplyElement(elem,WFin,lambda,n0,pscale,propdist)
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
            
            if elem.propagation_method_ == 0
                error('Not in a focal plane');
                
            elseif elem.propagation_method_ == 1
                
                if elem.amponly ~= false % amplitude mask
                    WFout = WFin .* elem.zsag_;
                else % phase mask
                    for ii = 1:numLambdas
                        WFout(:,:,ii) = WFin(:,:,ii) .* exp(-1i * elem.phasefac_(ii) .* elem.zsag_(:,:,ii));
                    end
                end
                
                fprintf('Moving back to pupil plane after FPM application\n');
                WFout = OptSys.move2PP(WFout,propdist,pscale,lambda);
                
            elseif elem.propagation_method_ == 2
                
            elseif elem.propagation_method_ == 3
                
                
            end
            
        end % ApplyElement
        
        %% Utility Methods
        
        function descr = describe(elem)
            % elem = describe(elem)

            elemtype = 'OptFPM';
            zpos = abs(elem.z_position_);
            sz = size(elem.zsag_);
            
            descr = sprintf('%s:\nElement is an %s\nGrid Size: [%d %d]',elem.name,elemtype, sz(1), sz(2));

        end
    end
    
end

