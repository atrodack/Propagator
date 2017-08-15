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
        amp_;
        pha_;
        
        % GPU properties
        useGPU = 0;         % use a GPU flag
        nGPUs;              % Number of GPUs MATLAB finds
        DEVICES;            % A list of those GPUs
        
    end
    
    
    
    
    
    
    methods
        
        %% Constructor
        
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
                F.set_Central_Wavelength(A{4,1});
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
        
        function F = set_lambda_array(F,lambdalist)
            % OS = set_lambda_array(OS,lambdalist)
            
            if nargin < 2
                if isempty(F.lambda0_) == 0
                    F.lambda_array_ = F.lambda0_;
                else
                    error('Need to include lambdalist or set central wavelength');
                end
            else
                F.lambda_array_ = lambdalist;
                if isempty(F.lambda0_)
                    F.set_Central_Wavelength();
                end
            end
        end % of set_lambda_array
        
        
        function F = set_Central_Wavelength(F,lambda)
            % OS = set_Central_Wavelength(OS,lambda)
            
            if nargin < 2
                if isempty(F.lambda_array_) == 0
                    len = length(F.lambda_array_);
                    ind = ceil(len/2);
                    F.lambda0_ = F.lambda_array_(1,ind);
                else
                    error('Need to include lambda or set lambda_array');
                end
            else
                F.lambda0_ = lambda;
                if isempty(F.lambda_array_)
                    F.set_lambda_array(lambda);
                end
                
            end
        end % of set_Central_Wavelength
        
        function F = set_field(F,field)
            % elem = set_field(elem,field)
            
            F.field_ = field;
            F.set_datatype;
            
        end % of set_field
        
        function F = set_seed(F,seed)
            % elem = set_seed(elem,seed)
            
            F.seed = seed;
        end % of set_seed
        
        function F = set_gridsize(F,N,M)
            % elem = set_N(elem,N)
            if nargin < 3
                M = N;
            end
            
            F.gridsize_ = [N, M];
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
                        F.field_ = double(F.field_);
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
            
            if isreal(F.field_)
                for ii = 1:sz(3)
                    imagesc(F.field_(:,:,ii))
                    plotUtils(sprintf('Field for %s, wavelength %d',F.name,ii));
                    drawnow;
                end
            else
                F.ReIm2WF;
                for ii = 1:sz(3)
                    subplot(1,2,1);
                    imagesc(F.amp_(:,:,ii))
                    plotUtils(sprintf('Field Amplitude for %s, wavelength %d',F.name,ii));
                    
                    subplot(1,2,2);
                    imagesc(F.pha_(:,:,ii))
                    plotUtils(sprintf('Field Phase for %s, wavelength %d',F.name,ii));
                    drawnow;
                end
            end
            
        end % of show
        
        function [KX,KY,KR,kx] = Kcoords2D(F)
            % [KX,KY,KR,kx] = Kcoords(OS)
            
            N = F.gridsize_(1);
            dk = (2*pi) ./ (N*F.dx_);
            kx = ((1:N) - (N/2))*dk;
            [KX,KY] = meshgrid(kx);
            KR = sqrt(KX.^2 + KY.^2);
            
        end % of Kcoords2D
        
        function [X,Y,R,x] = Coords2D(F)
            % [X,Y,R,x] = Coords2D(F)
            N = F.gridsize_(1);
            dx = F.dx_;
            
            x = ((1:N) - (N/2))*dx;
            [X,Y] = meshgrid(x);
            R = sqrt(X.^2 + Y.^2);
            
        end % of Coords2D
        
        %% Propagation Utilities
        
        function Phi = phase(F,lambda,ind)
            % Phi = phase(elem,lambda,ind)
            
            if nargin < 2
                lambda = F.lambda0_;
            end
            if nargin < 3
                ind = 1;
            end
            
            k = (2*pi) / lambda(ind);
            
            
            Phi = k .* F.field_(:,:,ind);
        end % of phase
        
        function Psi = phasor(F,lambda,ind)
            % Psi = phasor(elem, lambda,ind)
            
            if nargin < 2
                lambda = F.lambda0_;
            end
            if nargin < 3
                ind = 1;
            end
            
            Psi = exp(-1i * F.phase(lambda,ind));
        end % of phasor
        
        function F = ReIm2WF(F)
            % F = ReIm2WF(F)
            % Converts F.field_ from real/imag to amp/phase, and stores that
            % in F.field_
            if isa(F.field_,'gpuArray')
                flag = true;
            else
                flag = false;
            end
            [F.amp_,F.pha_] = WFReIm2AmpPhase2(real(F.field_),imag(F.field_));
            if flag
                F.amp_ = gpuArray(F.amp_);
                F.pha_ = gpuArray(F.pha_);
            end
            F.AmpPhase2WF();
            
        end % of ReIm2WF
        
        function F = AmpPhase2WF(F,ind)
            % F = AmpPhase2WF(F)
            % Sets the field as the combination of the amplitude and phase
            % components
            if nargin < 2
                for ii = 1:size(F.amp_,3)
                    F.field_(:,:,ii) = F.amp_(:,:,ii) .* exp(1i * F.pha_(:,:,ii));
                end
            elseif nargin == 2
                F.field_(:,:,ind) = F.amp(:,:,ind) .* exp(1i * F.pha_(:,:,ind));
            end
            
            
        end % of AmpPhase2WF
        
        function F = FresnelProp( F, z, pscale, lambda )
            % F = FresnelProp(F,z,pscale,lambda)
            
            k = (2*pi)./lambda;
            [ny,nx,nlambda] = size(F.field_);
            
            Lx = pscale*nx;
            Ly = pscale*ny;
            
            dfx = 1./Lx;
            dfy = 1./Ly;
            
            u = ones(nx,1)*((1:nx)-nx/2)*dfx;
            v = ((1:ny)-ny/2)'*ones(1,ny)*dfy;
            
            sqrdist = u.^2 + v.^2;
            coeff = pi*z*lambda;
            
            G = fftshift(fft2(fftshift(F.field_)))* Lx * Ly;
            
            H = zeros(ny,nx,nlambda);
            for ii = 1:nlambda
                H(:,:,ii) = exp(1i.*k(ii)*z) .* exp((-1i .* coeff(ii)) .* sqrdist);
            end
            
            F.set_field(ifftshift(ifft2(ifftshift(G.*H))) / Lx / Ly);
            F.ReIm2WF;
            
        end % of FresnelProp
        
        
        function F = FraunhoferProp(F, propdist, pscale, lambda)
            % F = FraunhoferProp(F, propdist, pscale, lambda)
            % Monochromatic Calculation - doesn't not scale psf size with
            % wavelength
            
            sz = size(F.field_);
            if length(sz) == 2
                sz(3) = 1;
            end
            
            k = (2*pi)./lambda;
            
            if isa( F.field_,'gpuArray')
                flag = true;
                tmp = gather( F.field_);
                datatype = class(tmp);
                clear tmp;
            else
                datatype = class( F.field_);
                flag = false;
            end
            
            WFfocus = init_variable(sz(1),sz(2),sz(3),datatype,0);
            if flag
                WFfocus = gpuArray(WFfocus);
            end
            
            for ii = 1:length(lambda)
                if flag
                    coeff = gpuArray((exp(1i.*k(ii)*propdist))/(1i.*(lambda(ii))*propdist));
                else
                    coeff = ((exp(1i.*k(ii)*propdist))/(1i.*(lambda(ii))*propdist));
                end
                WFfocus(:,:,ii) = coeff * fftshift(fft2(fftshift(OptSys.halfpixelShift( F.field_(:,:,ii),1 ) ))) .* (sz(1).* sz(2) .* pscale .* pscale);
            end
            F.set_field(WFfocus);
            F.ReIm2WF;
            
        end % FraunhoferProp

        function [F,PSF] = computePSF_chirpzDFT(F,nld,D)
            % F = computePSF_chirpzDFT(F,WFin)
            
            sz = size(F.field_);
            if length(sz) == 2
                sz(3) = 1;
            end
            
            WFfocus = init_variable(sz(1),sz(2),sz(3),F.default_data_type,0);
            N = sz(1);
            M = N;
            
            lambda = F.lambda_array_;
            ld = F.lambda0_ / D;
            
            x_extent = N*F.dx_;
            
            theta = nld.*ld;
            xi = theta ./ lambda;
            df = N / (M*x_extent);
            q = (df*M) / 2 ./ xi;
            
            f0 = zeros(1,length(lambda));
            for ii = 1:length(f0)
                f0(ii) = (df*M) - ((df/q(ii)) * M / 2);
            end
            
            
            
            for ii = 1:sz(3)
                WFfocus(:,:,ii) = OptSys.chirp_zDFT2(F.field_(:,:,ii), 0, F.dx_,f0(ii),df/q(ii),M)*(N * N * F.dx_ * F.dx_);
            end
            
            F.set_field(WFfocus);
            F.ReIm2WF;
            if N ~= M
                F.planewave(0,sz(3));
            end
            F.AmpPhase2WF;
            PSF = abs(F.field_).^2;
            
            
        end % of computePSF_chirpzDFT
        
        
        function F = FP2PP_Fraunhofer(F, propdist, pscale, lambda)
            % F = FP2PP_Fraunhofer(F, propdist, pscale, lambda)
            
            WFin = F.field_;
            
            sz = size(WFin);
            if length(sz) == 2
                sz(3) = 1;
            end
            
            k = (2*pi)./lambda;
            
            WFpupil = init_variable(sz(1),sz(2),sz(3),F.default_data_type,0);
            
            for ii = 1:length(lambda)
                coeff = (((exp(1i.*k(ii)*propdist))/(1i.*(lambda(ii))*propdist)))^-1;
                WFpupil(:,:,ii) = coeff * ifftshift(ifft2(ifftshift(  WFin(:,:,ii) ))) .* (sz(1).* sz(2) .* (pscale) .* (pscale))^-1;
            end
            F.set_field(WFpupil);
            F.ReIm2WF;

            
        end % FP2PP_Fraunhofer
        
        %% Set a Field
        
        function F = planewave(F,amplitude,numLambdas,theta)
            % OS = planewave(OS,amplitude,numLambdas,theta)
            % Initializes a planewave as the stored wavefront
            %
            % TO DO: add off-axis planewave functionality
            
            % If no input, just use ones
            if nargin < 2
                amplitude = complex(1);
                if isempty(F.lambda_array_)
                    numLambdas = 1;
                else
                    numLambdas = size(F.lambda_array_,2);
                end
            end
            
            % If no number of lambdas, just use one
            if nargin < 3
                if isempty(F.lambda_array_)
                    numLambdas = 1;
                else
                    numLambdas = size(F.lambda_array_,2);
                end
            end
            
            % Theta is the off axis angle of the planewave (not supported
            % yet)
            if nargin < 4
                on_axis = true;
            else
                on_axis = false;
            end
            
            if on_axis
                amplitude = single(amplitude);
                WF = amplitude .* init_variable(F.gridsize_(1), F.gridsize_(2),numLambdas,F.default_data_type,1);
                if F.useGPU == 1
                    F.field_ = gpuArray(WF);
                else
                    F.field_ = WF;
                end
            else  % off-axis
                
                [X,Y] = F.Coords2D;
                k = (2*pi) ./ F.lambda_array_;
                thx = theta(1) / 206265;
                thy = theta(2) / 206265;
                kappax = k.*thx;
                kappay = k.*thy;
                
                alpha = exp(-1i.*(kappax.*X+kappay.*Y));
                WF = amplitude .* init_variable(F.gridsize_(1), F.gridsize_(2),numLambdas,F.default_data_type,1) .*alpha;
                if F.useGPU == 1
                    F.field_ = gpuArray(WF);
                else
                    F.field_ = WF;
                end
            end
        end % of planewave
        
        
        function F = ApplyPhaseScreen(F,PS,lambda_array)
            % F = ApplyPhaseScreen(F,PS,lambda_array)
            
            if isa(PS,'OptPhaseScreen')
                screen = PS.field_;
            else
                % Assume user input a matrix turbulence screen
                screen = PS;
            end
                
                
            WFin = F.field_;
            
            sz = size(WFin);
            if length(sz) == 2
                sz(3) = 1;
            end
            
            sz2 = size(screen);
            if length(sz2) == 2
                sz2(3) = 1;
            end
            
            if ~isequal(sz(3),sz2(3))
                error('Must have a phasescreen for each wavelength');
            end
                        
            if isa(WFin,'gpuArray')
                flag = true;
                tmp = gather(WFin);
                datatype = class(tmp);
                clear tmp;
            else
                datatype = class(WFin);
                flag = false;
            end
            
            Field = init_variable(sz(1),sz(2),sz(3),datatype,0);
            if flag
                Field = gpuArray(Field);
            end
            
            if isa(PS,'OptPhaseScreen')
                for ii = 1:sz(3)
                    if flag
                        Psi = gpuArray(PS.phasor(lambda_array,ii));
                    else
                        Psi = (PS.phasor(lambda_array,ii));
                    end
                    Field(:,:,ii) = (F.field_(:,:,ii) .* Psi);
                end
                F.set_field(Field);
                if F.verbose == 1
                    fprintf('Phase Screen %s Applied\n',F.name);
                end
            else
                for ii = 1:sz(3)
                    k = (2*pi) / lambda_array(ii);
                    if flag
                        Psi = gpuArray(exp(-1i * k * screen));
                    else
                        Psi = exp(-1i * k * screen);
                    end
                    Field(:,:,ii) = (F.field_(:,:,ii) .* Psi);
                end
                F.set_field(Field);
                if F.verbose == 1
                    fprintf('Phase Screen Applied\n');
                end
            end
        end % of ApplyPhaseScreen
        
        function F = ApplyElement(F,elem,n0,fnum)
            % F = ApplyElement(F,elem,n0,fnum)
            
            if nargin < 4
%                 fprintf('No F-number supplied\n');
                flag = 0;
            elseif nargin == 4
                flag = 1;
            end
            
            if isa(elem,'OptElement')
                if isa(elem,'OptFPM')
                    if flag
                        D = elem.diameter_;
                        z = fnum*D;
                    else
                        error('Need to Supply an F-number');
                    end
                    F.set_field(elem.ApplyElement(F.field_,F.lambda_array_,n0,F.dx_,z));
                    F.ReIm2WF;
                else
                    F.set_field(elem.ApplyElement(F.field_,F.lambda_array_,n0));
                    F.ReIm2WF;
                end
                
            else
                error('Matrix inputs unsupported. elem must a child of OptElement');
                
            end
            
            
        end % of ApplyElement
    end
    
end

