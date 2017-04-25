classdef OptWF < matlab.mixin.Copyable
    %OPTWF Class for Optical Wavefronts / Fields
    %   Detailed explanation goes here
    
    properties(Constant=true, GetAccess='protected')
        
        
    end % of constant properties
    
    properties(GetAccess='public',SetAccess='public')
        name; % give a name to the structure/system
        savefile = 0;% save wavefronts as FITS
        verbose = 0; % print extra info
        seed; % RNG seed value
        
        
    end % of public properties
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        
        gridsize_; %
        lambda0_; % central wavelength
        lambda_array_; % vector of lambdas to use
        
        default_data_type = 'single';
        
        dx_;
        
        
        field_;         % Property for storing the field
        
    end
    
    
    
    
    
    
    methods
        
        function F = OptWF(PROPERTIES)
            % elem = OptWF(PROPERTIES)
            % PROPERTIES is a 5x1 cell array containing element properties
            % If something is unknown or to be set later, leave that
            % position empty.
            %
%           PROPERTIES{1,1} = name                  [string]
%           PROPERTIES{2,1} = gridsize              [float in m]
%           PROPERTIES{3,1} = dx                    [float in m]
%           PROPERTIES{4,1} = lambda_0              [float in m]
%           PROPERTIES{5,1} = field                 [float in ]
            

            
            if size(PROPERTIES) >= [5,1]
                if iscell(PROPERTIES) == 1
                    A = PROPERTIES;
                else
                    error('PROPERTIES must be 5x1 cell array!');
                end
            else
                error('PROPERTIES must be a 5x1 cell array!');
            end
            
            if size(PROPERTIES) == [5,1]
                % Constructing a Field
                F.set_name(A{1,1});
                F.set_gridsize(A{2,1});
                F.set_dx(A{3,1});
                F.set_Central_wavelength(A{4,1});
                F.set_field(A{5,1});
            elseif size(PROPERTIES) == [10,1]
                % Constructing a Phase Screen, let it's constructor handle it
            end
                
            
        end % of constructor
        
        %% Property Setting Utilities
        
        function F = set_name(F,name)
            % elem = set_name(elem,name)
            F.name = name;
        end % of set_name
        
        function F = set_Central_wavelength(F,lambda)
            % elem = set_Central_wavelength(elem,lambda)
            
            F.lambda0_ = lambda;
        end % of set_refLambda
        
        function F = set_field(F,field)
            % elem = set_field(elem,field)
            
            F.field_ = field;
            F.set_datatype;
            
        end % of set_field
        
        function F = set_seed(F,seed)
            % elem = set_seed(elem,seed)
            
            F.seed = seed;
        end % of set_seed
        
        function F = set_gridsize(F,N)
            % elem = set_N(elem,N)
            F.gridsize_ = N;
        end % of set_N
        
        function F = set_dx(F,dx)
            % elem = set_dx(elem,dx)
            F.dx_ = dx;
        end % of set_dx
        
        %% Utilities
        
        function F = set_datatype(F,default_data_type)
            % elem = datatype(default_data_type)
            % Sets the datatype to use. Do not do anything if already of
            % the correct data type.
            % Currently supported:
            % single
            % double
            % uint8
            
            if nargin < 2
                default_data_type = F.default_data_type;
            else
                F.default_data_type = default_data_type;
            end
            
            switch default_data_type
                case 'single'
                    if ~isa(F.field_,'single')
                        F.field_ = single(F.field_);
                        if F.verbose == 1
                            fprintf('Data Type set to single\n');
                        end
                    end
                    
                case 'double'
                    if ~isa(F.field_,'double')
                        F.OS.field_ = double(F.field_);
                        if F.verbose == 1
                            fprintf('Data Type set to double\n');
                        end
                    end
                    
                case 'uint8'
                    if ~isa(F.field_,'uint8')
                        F.field_ = uint8(F.field_);
                        if F.verbose == 1
                            fprintf('Data Type set to uint8\n');
                        end
                    end
                    
                otherwise
                    error('I do not understand that data type (yet)!');
            end
        end % of setdatatype

        
        function Phi = phase(elem,lambda,ind)
            % Phi = phase(elem,lambda,ind)
            
            if nargin < 2
                lambda = elem.lambda0_;
            end
            if nargin < 3
                ind = 1;
            end
            
            k = (2*pi) / lambda(ind);
            
            
            Phi = k .* elem.field_(:,:,ind);
        end % of phase
        
        function Psi = phasor(elem,lambda,ind)
            % Psi = phasor(elem, lambda,ind)
            
            if nargin < 2
                lambda = elem.lambda0_;
            end
            if nargin < 3
                ind = 1;
            end
            
            Psi = exp(-1i * elem.phase(lambda,ind));
        end % of phasor
        
        function rv = complexNormalRV(elem,sigma, gridsize)
            % rv = complexGaussianRV(sigma)
            % Returns a complex Normal random variable with zero-mean,
            % sigma standard deviation
            
            sz = [gridsize, gridsize];
            if(isempty(elem.seed))
                rv = (sigma/sqrt(2)) * (randn(sz) + 1i * randn(sz));
            else
                rng(elem.seed)
                rv = (sigma/sqrt(2)) * (randn(sz) + 1i * randn(sz));
                rng('default');
            end
        end % of complexNormalRV
        
        function show(F,fignum)
            % show(elem,fignum)
            % Plots matrix stored in field_
            if nargin < 2
                figure;
            else
                figure(fignum);
            end
            
            sz = size(F.field_);
            if length(sz) == 2
                sz(3) = 1;
            end
            
            for ii = 1:sz(3)
                imagesc(F.field_(:,:,ii))
                plotUtils(sprintf('Displacement for %s, wavelength %d',F.name,ii));
                drawnow;
            end
            
        end % of show
        
        function [KX,KY,KR,kx] = Kcoords2D(F)
            % [KX,KY,KR,kx] = Kcoords(OS)
            
            N = F.N_;
            dk = (2*pi) ./ (N*F.dx_);
            kx = ((1:N) - (N/2))*dk;
            [KX,KY] = meshgrid(kx);
            KR = sqrt(KX.^2 + KY.^2);
            
        end % of Kcoords2D
        
    end
    
end

