classdef OptPupil < OptElement
    %OptPupil Class for System Pupils and masks
    %   This is intended to be a binary mask class. I will leave a hook
    %   that allows it to be a phase mask, but leave that to the user to
    %   add.
    
    
    % These properties are within the class OptPupil Scope.
    % This means that the value/s stored in them are accessible to the
    % outside world, but only the class object itself can modify them
    properties(GetAccess = 'public', SetAccess = 'protected')
        amponly = true; % flag for determining if Mask is phase or amplitude
                        % Defaults to amplitude
        
    end
    
    methods
        %% Constructor
        function elem = OptPupil(PROPERTIES)
            % elem = OptPupil(PROPERTIES)
            % PROPERTIES is a 6x1 cell array containing element properties
            % If something is unknown or to be set later, leave that
            % position in A empty:
            %
%           PROPERTIES{1,1} = name                  [string]
%           PROPERTIES{2,1} = material              [0-10]       
%           PROPERTIES{3,1} = amponly               [bool]
%           PROPERTIES{4,1} = z_position            [float in m]
%           PROPERTIES{5,1} = diameter              [float in m]
%           PROPERTIES{6,1} = mask                  [file path or matrix]
            
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
            elem.setdatatype();

        end % of contructor
        
        %% Set Pupil Specific Properties
        
        function elem = set_amponly(elem,flag)
            % elem = set_amponly(elem,flag)
            % set if the mask is amplitude or phase mask
%            (1 for amp, 0 for phase).
            
            elem.amponly = flag;
            
        end % of set_amponly
        
        %% Propagation Methods
        % Should only be used by the propagation calls from OptSys
        
        
        function WFout = ApplyElement(elem,WFin,lambda)
            % OS = ApplyElement(elem,WFin,lambda)
            % Applys element to current wavefront
            
            % Get number of wavelengths
            numLambdas = length(lambda);
            elem.set_phasefactor(lambda);
            
            % Edit zsag_ to be same dimension as WFin
            if numLambdas ~= 1
                if size(elem.zsag_,3) ~= numLambdas
                    if strcmp( elem.default_data_type, 'single')
                        tmp = single(zeros(elem.gridsize_(1),elem.gridsize_(2),numLambdas));
                        WFout = tmp;
                    elseif strcmp( elem.default_data_type, 'double')
                        tmp = double(zeros(elem.gridsize_(1),elem.gridsize_(2),numLambdas));
                        WFout = tmp;
                    elseif strcmp( elem.default_data_type, 'uint8')
                        tmp = uint8(zeros(elem.gridsize_(1),elem.gridsize_(2),numLambdas));
                        WFout = tmp;
                    end
                        
                    for ii = 1:numLambdas
                        tmp(:,:,ii) = elem.zsag_;
                    end
                    elem.set_zsag(tmp);
                    clear tmp;
                end
            end
            
            if elem.amponly ~= false % amplitude mask
                WFout = WFin .* elem.zsag_;
            else % phase mask
                elem.set_phasefactor(lambda);
                for ii = 1:numLambdas
                    WFout(:,:,ii) = WFin(:,:,ii) .* exp(-1i * elem.phasefac_(ii) .* elem.zsag_(:,:,ii));
                end
            end
            
        end % ApplyElement
        
        
        %% Overloaded Operators
        
        function descr = describe(elem)
            % elem = describe(elem)
            if elem.verbose == 1
                flag = true;
                elem.verbose = 0;
            end
            elemtype = 'OptPupil';
            zpos = abs(elem.z_position_);
            sz = size(elem.zsag_);
            
            descr = sprintf('%s:\nElement is an %s at %0.3f meters downstream\nGrid Size: [%d %d]',elem.name,elemtype,zpos, sz(1), sz(2));
            if flag == true
                elem.verbose = 1;
            end
        end
        
    end
    
end

