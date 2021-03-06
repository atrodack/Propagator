classdef OptSys < matlab.mixin.Copyable
    %OptSys Class for storing optical system parameters
    %   Detailed explanation goes here
    
    
    %% Properties
    properties(Constant=true, GetAccess='protected')
        
        
    end % of constant properties
    
    properties(GetAccess='public',SetAccess='public')
        name; % give a name to the structure/system
        savefile = 1;% save wavefronts as FITS
        verbose = 1; % print extra info
        
    end % of public properties
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        
        % System Parameters
        numElements_ = 0; % number of elements in system
        conjugations_ = zeros(1,1); % list of z locations
        f_number_; % F/# of system
        beam_radius_; % radius of beam in pixels
        pscale_; % size of pixels in meters
        gridsize_; %
        lambda0_; % central wavelength
        lambda_array_; % vector of lambdas to use
        lambda_scales_; % vector of scales for wavelengths
        
        % List of Elements
        ELEMENTS_; % cell array of element cards
        % Each element is of type OptElement
        % This includes the information:
        % 1) Element type
        % 2) Material (index of refraction)
        % 3) Element Position
        % 4) z-sag in meters
        % And possibly other pieces of information depending on the element
        % type
        
        % Wavefront
        WFamp;          % The amplitude of the complex wavefront
        WFphase;        % The phase of the complex wavefront
        WF_;            % The complex wavefront
        PSF_;           % Modulus Squared of complex Wavefront in Focal Plane
        
        % GPU properties
        useGPU = 0;         % use a GPU flag
        nGPUs;              % Number of GPUs MATLAB finds
        DEVICES;            % A list of those GPUs
        
        % Pool Characteristics (UNUSED)
        nWorkers; % number of workers for parpool
        
        
        % Misc
        interpolate_method = [];
        default_data_type = 'single';   % Default to type single for faster calculations
                                        % Change if higher precision is
                                        % required
        
        
    end % of protected properties
    
    %% Methods
    methods
        
        %% Constructor
        
        function OS = OptSys(sz)
            % OS = OptSys(sz)
            % Constructs OptSys object. Initializes WFCA to size sz
            
            if nargin < 1
                sz = [1024,1024];
            end
            
            OS.initOptSys(sz);
            
            
        end % of contructor
        
        
        %% Initializer for OptSys
        
        function OS = initOptSys(OS,sz)
            %OS = initOptSys(OS,sz)
            % Initialier for the object. Creates cell array in ELEMENTS_ to
            % store elements, and querys the computer if it has an NVIDIA
            % GPU. If it finds a GPU, it is initialized as well.
            
            OS.ELEMENTS_ = cell(OS.numElements_,1);
            OS.nGPUs = gpuDeviceCount;
            if OS.nGPUs > 0
                OS.initGPU;
            end
            OS.set_gridsize(sz);
            
        end % of initOptSys
        
        %% Utilities for adding and removing elements from the system
        
        function OS = addSequentialElement(OS,elem)
            % OS = addSequentialElement(OS,elem)
            % adds an Element object to the optical system located after
            % the previous stored element
            
            % Children of OptElement will return true to this
            if isa(elem,'OptElement')
                OS.numElements_ = OS.numElements_ + 1;
                OS.ELEMENTS_{OS.numElements_,1} = elem;
                %                 OS.set_conjugations; %maybe add in later
                
                if OS.verbose == 1
                    fprintf('Element %s Added to Optical System\n',elem.name);
                end
            else
                error('I do not understand anything but OptElement objects');
            end
        end % of addSequentialElement
        
        function OS = addSequentialElements(OS,elems)
            % OS = addSequentialElements(OS,elems)
            % adds list of elements "elems" to Optical system sequentially
            % elems is a cell array of OptElement objects
            
            numElems = length(elems);
            for ii = 1:numElems
                OS.addSequentialElement(elems{ii});
            end
        end % of addSequentialElements
            
        
        function OS = addElement(OS,elem,order_num)
            % OS = addElement(OS,elem,order_num)
            % Adds an Element to the optical system in the location given
            % in order_num
            
            if ~isa(elem,'OptElement')
                error('elem must be type OptElement');
            end
            
            element_list = OS.ELEMENTS_;
            tmpelement_list = cell(OS.numElements_ + 1,1);
            
            index = 1;
            for ii = 1:length(tmpelement_list)
                if isequal(ii,order_num) == 0
                    tmpelement_list{ii} = element_list{index};
                    index = index + 1;
                else
                    tmpelement_list{ii} = elem;
                    
                end
            end
            
            OS.numElements_ = OS.numElements_ + 1;
            OS.ELEMENTS_ = tmpelement_list;
            
            if OS.verbose == 1
                fprintf('Element %s Added to Optical System at position %d\n',elem.name,order_num);
            end;
            
        end % of addElement
        
        
        function OS = removeElement(OS,elem_num)
            % OS = removeElement(OS,elem_num)
            % removes and element from the Optical System
            
            element_list = OS.ELEMENTS_;
            tmpelement_list = cell(OS.numElements_ - 1,1);
            
            selected_elem = element_list{elem_num};
            prompt = sprintf('You have seleced element %s for removal\nIs this Correct?\nY/N [Y]: ',selected_elem.name);
            choice = input(prompt,'s');
            if isempty(choice)
                choice = 'Y';
            end
            
            if strcmpi(choice,'Y') % remove the element
                element_list{elem_num} = [];
                index = 1;
                
                for ii = 1:OS.numElements_
                    if isempty(element_list{ii}) == 0
                        tmpelement_list{index} = element_list{ii};
                        index = index + 1;
                    end
                end
                OS.numElements_ = OS.numElements_ - 1;
                OS.ELEMENTS_ = tmpelement_list;
                
                if OS.verbose == 1
                    fprintf('Element %s Removed from Optical System\n',selected_elem.name);
                end
                
            else
                fprintf('Element %s has not been removed\n',selected_elem.name);
            end
        end % of removeElement
        
        %% Property Setting Methods
        
        function OS = set_pscale(OS,val)
            % OS = set_pscale(OS,val)
            OS.pscale_ = val;
        end % of set_pscale
        
        function OS = set_beamrad(OS,val)
            % OS = set_beamrad(OS,val)
            
            OS.beam_radius_ = val;
            
        end % of set_beamrad
        
        function OS = set_lambda_array(OS,lambdalist)
            % OS = set_lambda_array(OS,lambdalist)
            
            if nargin < 2
                if isempty(OS.lambda0_) == 0
                    OS.lambda_array_ = OS.lambda0_;
                else
                    error('Need to include lambdalist or set central wavelength');
                end
            else
                OS.lambda_array_ = lambdalist;
                if isempty(OS.lambda0_)
                    OS.set_Central_Wavelength();
                end
            end
        end % of set_lambda_array
        
        
        function OS = set_Central_Wavelength(OS,lambda)
            % OS = set_Central_Wavelength(OS,lambda)
            
            if nargin < 2
                if isempty(OS.lambda_array_) == 0
                    len = length(OS.lambda_array_);
                    ind = ceil(len/2);
                    OS.lambda0_ = OS.lambda_array_(1,ind);
                else
                    error('Need to include lambda or set lambda_array');
                end
            else
                OS.lambda0_ = lambda;
                if isempty(OS.lambda_array_)
                    OS.set_lambda_array(lambda);
                end
                
            end
        end % of set_Central_Wavelength
        
%         
%         function OS = setWavelengthScales(OS)
%             % OS = setWavelengthScales(OS)
%             
%             tmp = OS.pscale_ / OS.lambda0_;
%             OS.lambda_scales_ = tmp * OS.lambda_array_;
%             
%         end % of setWavelengthScales

        
        function OS = set_field(OS,field)
            % OS = set_field(OS,field)
            if nargin < 2
                OS.planewave();
            elseif nargin == 2
                if(strcmp(OS.default_data_type, 'single'))
                    OS.WF_ = single(field);
                elseif(strcmp(OS.default_data_type, 'double'))
                    OS.WF_ = double(field);
                elseif(strcmp(OS.default_data_type, 'uint8'))
                    OS.WF_ = uint8(field);
                else
                    error('Data type not supported');
                end
            end
        end % of set_field
            
        function OS = setdatatype(OS,default_data_type)
            % elem = datatype(default_data_type)
            % Sets the datatype to use. Do not do anything if already of
            % the correct data type.
            % Currently supported:
            % single
            % double
            % uint8
            
            if nargin < 2
                default_data_type = OS.default_data_type;
            else
                OS.default_data_type = default_data_type;
            end
            
            switch default_data_type
                case 'single'
                    if ~isa(OS.WF_,'single')
                        OS.WF_ = single(OS.WF_);
                        if OS.verbose == 1
                            fprintf('Data Type set to single\n');
                        end
                    end
                    
                case 'double'
                    if ~isa(OS.WF_,'double')
                        OS.OS.WF_ = double(OS.WF_);
                        if OS.verbose == 1
                            fprintf('Data Type set to double\n');
                        end
                    end
                    
                case 'uint8'
                    if ~isa(OS.WF_,'uint8')
                        OS.WF_ = uint8(OS.WF_);
                        if OS.verbose == 1
                            fprintf('Data Type set to uint8\n');
                        end
                    end
                    
                otherwise
                    error('I do not understand that data type (yet)!');
            end
        end % of setdatatype
        
        function conjugations = set_conjugations(OS)
            % OS = set_conjugations(OS)
            % Adds the zpos_ of an added element to the conjugations list
            
            conjugations = zeros(OS.numElements_,1);
            for z = 1:OS.numElements_
                conjugations(z,1) = OS.ELEMENTS_{z}.z_position_;
            end
            OS.conjugations_ = conjugations;
            
        end % of set_conjugations
        
        function OS = set_Fnum(OS,val)
            % OS = set_Fnum(OS,val)
            
            OS.f_number_ = val;
        end % of set_Fnum
        
        function OS = get_element_Fnum(OS,element_num)
            % OS = get_element_Fnum(OS,element_num)
            % sets the f/# of the system to that of an element
            
            OS.f_number_ = OS.ELEMENTS_{element_num}.getFNumber;
        end % of get_element_Fnum
        
        
        function gridsize = set_gridsize(OS,sz)
            % gridsize = set_gridsize(OS,sz)
            
            if nargin < 2
                gridsize = size(OS.ELEMENTS_{1}.zsag_);
            else
                gridsize = sz;
            end
            
            OS.gridsize_ = gridsize;
            
        end % of set_gridsize
        
        function show(OS,fignum)
            % show(OS)
            % Plots matrix stored in zsag_
            numLambdas = size(OS.WF_,3);
            
            if isreal(OS.WF_) == 1
                if nargin < 2
                    figure;
                else
                    figure(fignum);
                end
                for ii = 1:numLambdas
                    imagesc(OS.WF_(:,:,ii))
                    plotUtils(sprintf('WF lambda %d',ii));
                    drawnow;
                end
            else
                if nargin < 2
                    figure;
                else
                    figure(fignum);
                end
                for ii = 1:numLambdas
                    subplot(1,2,1)
                    imagesc(OS.WFamp(:,:,ii));
                    plotUtils(sprintf('WF Amplitude, lambda %d',ii));
                    subplot(1,2,2)
                    imagesc(OS.WFphase(:,:,ii));
                    plotUtils(sprintf('WF Phase, lambda %d',ii));
                    drawnow;
                end
                
            end
            
        end % of show
        
        function show_PP(OS,fignum)
            % show(OS,fignum)
            numLambdas = size(OS.lambda_array_,2);
            
            [xx,yy] = OS.PPcoords;
            
            if isreal(OS.WF_) == 1
                if nargin < 2
                    figure;
                else
                    figure(fignum);
                end
                for ii = 1:numLambdas
                    imagesc(xx(:,:,ii),yy(:,:,ii),OS.WF_(:,:,ii))
%                     imagesc(OS.WF_(:,:,ii));
                    plotUtils(sprintf('WF\n lambda %g',OS.lambda_array_(ii)),'x','y');
                    drawnow;
                end
            else
                if nargin < 2
                    figure;
                else
                    figure(fignum);
                end
                for ii = 1:numLambdas
                    subplot(1,2,1)
                    imagesc(xx(:,:,ii),yy(:,:,ii),OS.WFamp(:,:,ii));
%                     imagesc(OS.WFamp(:,:,ii));
                    plotUtils(sprintf('WF Amplitude\n lambda %d',ii),'x','y');
                    subplot(1,2,2)
                    imagesc(xx(:,:,ii),yy(:,:,ii),OS.WFphase(:,:,ii));
%                     imagesc(OS.WFphase(:,:,ii));
                    plotUtils(sprintf('WF Phase\n lambda %d',ii),'x','y');
                    drawnow;
                end
                
            end
            
        end % of show_PP
        
%         function show_FP(OS,fignum)
%             % show(OS,fignum)
%             numLambdas = size(OS.lambda_array_,2);
%             
% %             [xx,yy,dTH] = OS.FPcoords;
%             
%             if isreal(OS.WF_) == 1
%                 if nargin < 2
%                     figure;
%                 else
%                     figure(fignum);
%                 end
%                 for ii = 1:numLambdas
%                     normalizer = max(max(OS.WF_));
%                     imagesc(xx(:,:,ii),yy(:,:,ii),OS.WF_(:,:,ii))
% %                     imagesc(OS.WF_(:,:,ii));
%                     plotUtils(sprintf('WF\n lambda %d',ii),'\lambda / D', '\lambda / D');
%                     drawnow;
%                 end
%             else
%                 if nargin < 2
%                     figure;
%                 else
%                     figure(fignum);
%                 end
%                 for ii = 1:numLambdas
%                     normalizer = max(max(OS.WFamp));
% %                     subplot(1,2,1)
%                     imagesc(xx(:,:,ii),yy(:,:,ii),log10(OS.WFamp(:,:,ii)/normalizer(ii)),[-6,0]);
% %                     imagesc(OS.WFamp(:,:,ii));
%                     plotUtils(sprintf('WF Amplitude\n lambda %g',OS.lambda_array_(ii)),'\lambda / D', '\lambda / D');
% %                     subplot(1,2,2)
% %                     imagesc(xx(:,:,ii),yy(:,:,ii),OS.WFphase(:,:,ii));
% % %                     imagesc(OS.WFphase(:,:,ii));
% %                     plotUtils(sprintf('WF Phase\n lambda %g',OS.lambda_array_(ii)),'\lambda / D', '\lambda / D');
%                     drawnow;
%                 end
%                 
%             end
%             
%         end % of show_FP
        
        
        %% Propagation Utility Methods
        
        function [numLambdas,propdists] = initPropagation(OS,starting_elem, ending_elem)
            % [numLambdas,propdist] = initPropagation(OS,starting_elem, ending_elem)
            % Runs checks on OS object to ensure properties needed to
            % propagate through system are set.
            %
            % Returns number of lambdas and list of propagation distances
            % for use in PropagateSystem method
            
            OK = false(4,1);
            % Check that given Elements exist
            if OS.numElements_ >= starting_elem
                if ending_elem <= OS.numElements_
                    OK(1,1) = true;
                else
                    fprintf('\n\n');
                    warning('Element %d does not exist',ending_elem);
                end
            else
                fprintf('\n\n');
                warning('Element %d does not exist',starting_elem);
            end
            
            % Check that pscale_ is set
            if isempty(OS.pscale_) == 0
                OK(2,1) = true;
            else
                fprintf('\n\n');
                warning('Pixel Scale is not set');
            end
            
            % Check Wavelength
            if isempty(OS.lambda_array_) == 0
                numLambdas = length(OS.lambda_array_);
                OK(3,1) = true;
            else
                if isempty(OS.lambda0_) == 0
                    numLambdas = 1;
                    OK(3,1) = true;
                else
                    fprintf('\n\n');
                    warning('Central Wavelength and Lambda_array are empty');
                end
            end
            
            % Check that lambda_scales is set
            if isempty(OS.lambda_scales_) == 1
                OS.setWavelengthScales;
            end
            
            % Check Element Gridsizes
            elements = OS.ELEMENTS_;
            numelems = ending_elem - starting_elem + 1;
            gridsizes = zeros(numelems,2);
            for ii = starting_elem:ending_elem
                gridsizes(ii,1) = elements{ii}.gridsize_(1,1);
                gridsizes(ii,2) = elements{ii}.gridsize_(1,2);
            end
            if prod(gridsizes(:,1)) == OS.gridsize_(1)^numelems
                if prod(gridsizes(:,2)) == OS.gridsize_(2)^numelems
                    OK(4,1) = true;
                else
                    fprintf('\n\n');
                    warning('2nd Dimensions are not equal');
                end
            else
                fprintf('\n\n');
                warning('1st Dimensions are not equal');
            end
                
                
            
            % Error out if a Check fails
            if prod(OK) == 0
                fprintf('\n\n');
                error('Check failed. See above warnings for help');
            end
            
            % Set the propagation distances
            propdists = OS.computePropdists();
            
            
            
        end % initPropagation
        
        
        function OS = ReIm2WF(OS)
            % OS = ReIm2WF(OS)
            % Converts OS.WF_ from real/imag to amp/phase, and stores that
            % in OS.WF_
            
            [OS.WFamp,OS.WFphase] = WFReIm2AmpPhase2(real(OS.WF_),imag(OS.WF_));
            OS.AmpPhase2WF();
            
        end % of ReIm2WF
        
        function OS = AmpPhase2WF(OS,ind)
            % OS = AmpPhase2WF(OS)
            % Sets the WF as the combination of the amplitude and phase
            % components
            if nargin < 2
                for ii = 1:size(OS.WFamp,3)
                    OS.WF_(:,:,ii) = OS.WFamp(:,:,ii) .* exp(1i * OS.WFphase(:,:,ii));
                end
            elseif nargin == 2
                OS.WF_(:,:,ind) = OS.WFamp(:,:,ind) .* exp(1i * OS.WFphase(:,:,ind));
            end
            
            
        end % of AmpPhase2WF
        
        
        function propdists = computePropdists(OS,z_list)
            % propdist = computePropdists(OS)
            
            if nargin < 2
                if length(OS.conjugations_) ~= OS.numElements_
                    z_list = OS.set_conjugations();
                else
                    z_list = OS.conjugations_;
                end
            end
            
            propdists = zeros(length(z_list),1);
            for ii = 2:length(z_list)
                propdists(ii) = z_list(ii) - z_list(ii-1);
            end
        
        end % of computePropdists
        
        
        
        function PSF = computePSF_FFT(OS,WFin)
            % OS = computePSF_FFT(OS,WFin)
            
            if nargin < 2
                WFin = OS.WF_;
            end
               
            sz = size(WFin);
            if length(sz) == 2
                sz(3) = 1;
            end

            WFfocus = fftshift(fft2(fftshift(OptSys.halfpixelShift(WFin,1)))) .* (sz(1).* sz(2) .* OS.pscale_ .* OS.pscale_);

            WFreal = real(WFfocus);
            WFimag = imag(WFfocus);
            [OS.WFamp,OS.WFphase] = WFReIm2AmpPhase2(WFreal,WFimag);
            OS.AmpPhase2WF();
            PSF = abs(OS.WF_).^2;
            OS.PSF_ = PSF;


        end % computePSF_FFT
        
        function output = initchirpzDFT(OS,ind)
            % output = initchirpzDFT(OS)
            
            D = 2*OS.beam_radius_ * OS.pscale_;
            lambda = OS.lambda_array_;
            ld = OS.lambda0_ / D;
            
            if isa(OS.ELEMENTS_{ind},'OptDetector')
                nld = OS.ELEMENTS_{ind}.FPregion_;
                M = ceil((2*nld) / OS.ELEMENTS_{ind}.pixel_size_);
            else
                nld = 20;
                M = size(OS.WF_,1);
            end
            
            N = size(OS.WF_,1);
            dx = OS.pscale_;
            x_extent = N*dx;
            
            theta = nld.*ld;
            xi = theta ./ lambda;
            df = N / (M*x_extent);
            q = (df*M) / 2 ./ xi;
            
            f0 = zeros(1,length(lambda));
            for ii = 1:length(f0)
                f0(ii) = (df*M) - ((df/q(ii)) * M / 2);
            end
            
            phaseshift = zeros(N,N,length(lambda));
            for ii = 1:length(lambda)
                shift = linspace(0,pi/(2*q(ii)),N+1);
                shift = shift(1:end-1);
                [SX,SY] = meshgrid(shift);
                Tip = exp(-1i.*SX);
                Tilt = exp(-1i.*SY);
                phaseshift(:,:,ii) = Tip.*Tilt;
            end
            
            output{1} = nld;
            output{2} = N;
            output{3} = dx;
            output{4} = f0;
            output{5} = df./q;
            output{6} = M;
            output{7} = phaseshift;
            
        end % initchirpzDFT
        
        function PSF = computePSF_chirpzDFT(OS,WFin,ind)
            % PSF = computePSF_chirpzDFT(OS,WFin)
            
            params = OS.initchirpzDFT(ind);
            
            sz = size(WFin);
            if length(sz) == 2
                sz(3) = 1;
            end
            
            % Too slow column by column on GPU -- bring this to the CPU
            if isa(WFin,'gpuArray')
                flag = false;
                WFin = gather(WFin);
                datatype = class(WFin);
            else
                datatype = class(WFin);
                flag = false;
            end
            
            if strcmp(datatype,'single')
                WFfocus = single(zeros(params{6},params{6},sz(3)));
                if flag
                    WFfocus = gpuArray(WFfocus);
                end
            elseif strcmp(datatype,'double')
                WFfocus = double(zeros(params{6},params{6},sz(3)));
            elseif strcmp(datatype,'uint8')
                WFfocus = uint8(zeros(params{6},params{6},sz(3)));
            else
                error('Unsupported data type');
            end
            
            for ii = 1:sz(3)
                WFfocus(:,:,ii) = OptSys.chirp_zDFT2(WFin(:,:,ii), 0, params{3},params{4}(ii),params{5}(ii),params{6})*(params{2} * params{2} * OS.pscale_ * OS.pscale_);
            end
            WFreal = real(WFfocus);
            WFimag = imag(WFfocus);
            [OS.WFamp,OS.WFphase] = WFReIm2AmpPhase2(WFreal,WFimag);
            if params{2} ~= params{6}
                OS.set_field(zeros(params{6},params{6},size(WFin,3)));
                OS.setdatatype(OS.default_data_type);
            end
            OS.AmpPhase2WF();
            PSF = abs(OS.WF_).^2;
            OS.PSF_ = PSF;
            
            
            
            
        end % of computePSF_chirpzDFT
        
        
        
        function [xx,yy] = PPcoords(OS)
            % [xx,yy] = PPcoords(OS)
            
            sz = size(OS.WF_);
            if length(sz) == 2
                sz(3) = 1;
            end
            
            lambda = OS.lambda_array_;
            
            xx = zeros(sz(2),1,sz(3));
            yy = xx;
            
            for jj = 1:length(lambda)
                xx(:,:,jj) = (-((sz(2)/2)):((sz(2))/2 -1))*(OS.pscale_);
                yy(:,:,jj) = (-((sz(1)/2)):((sz(1))/2 -1))*(OS.pscale_);
            end
        end % of PPcoords
        
        function [ldx,ldy,fx,fy] = FPcoords(OS,nld,N)
            % [xx,yy] = FPcoords(OS)
            % call after computing field in Focal Plane
            
            sz = size(OS.PSF_);
            if length(sz) == 2
                sz(3) = 1;
            end
            lambda = OS.lambda_array_;
            
            D = 2*OS.beam_radius_ * OS.pscale_;
            ld = lambda ./ D;
            theta = nld.*ld;
            xi = theta ./ lambda;
            x_extent = N*OS.pscale_;
            df = N / (sz(1)*x_extent);
            q = (df*sz(1)) / 2 ./ xi;
            
            fx = zeros(1,sz(1),sz(3));
            fy = fx;
            ldx = fx;
            ldy = fx;
           
            for ii = 1:length(lambda)
                fx(1,:,ii) = (-floor(sz(1)/2 +0.5):floor(sz(1)/2)-0.5) * (df/q(ii));
                fy(1,:,ii) = fx(1,:,ii);
                ldx(1,:,ii) = fx(1,:,ii) * D;
                ldy(1,:,ii) = fy(1,:,ii) * D;
            end
            
        end % of FPcoords
        
        function [KX,KY,KR,kx] = Kcoords2D(OS)
            % [KX,KY,KR,kx] = Kcoords(OS)
            
            N = OS.gridsize_(1);
            dk = (2*pi) ./ (N*OS.pscale_);
            kx = ((1:N) - (N/2))*dk;
            [KX,KY] = meshgrid(kx);
            KR = sqrt(KX.^2 + KY.^2);
            
        end % of Kcoords2D
        
        %% System Propagation Methods
        % See static methods for actual Fresnel propagation method
        
        function OS = propagate2Elem(OS, propdist, pscale, lambda)
            % OS = propagate2Elem(OS,propdist,pscale,lambda)
            
            % Old Propagation Method (pixel by pixel)
%             OS.set_field(OptSys.FresnelPropagateWF(OS.WF_,propdist,pscale,lambda));
%             OS.ReIm2WF;

            % Faster Propagation Method
            OS.set_field(OptSys.FresnelProp(OS.WF_,propdist,pscale,lambda));
            OS.ReIm2WF;
            
        end % of propagate2Elem
        
        function OS = ApplyElement(OS,elem_num,n0)
            % OS = ApplyElement(OS,elem_num,n0)
            
            if isa(OS.ELEMENTS_{elem_num},'OptFPM')
                % The FPM needs access to the pscale property and a
                % propdist
                D = 2*OS.beam_radius_ * OS.pscale_;
                z = OS.f_number_ * D;
                OS.set_field( OS.ELEMENTS_{elem_num}.ApplyElement(OS.WF_,OS.lambda_array_,n0,OS.pscale_,z));
                OS.ReIm2WF;
            else
                OS.set_field( OS.ELEMENTS_{elem_num}.ApplyElement(OS.WF_,OS.lambda_array_,n0));
                OS.ReIm2WF;
            end
        end % of ApplyElement
        
        function OS = ApplyPhaseScreen(OS,PS)
            % OS = ApplyPhaseScreen(OS,PS)
            
            
            sz = size(OS.WF_);
            if length(sz) == 2
                sz(3) = 1;
            end
            
            sz2 = size(PS.field_);
            if length(sz2) == 2
                sz2(3) = 1;
            end
            
            if ~isequal(sz(3),sz2(3))
                error('Must have a phasescreen for each wavelength');
            end
            
            WFin = OS.WF_;
            lambda_array = OS.lambda_array_;
            
            if isa(PS,'OptPhaseScreen')
                WFout = PS.ApplyPhaseScreen(WFin,lambda_array);
                OS.set_field(WFout);
                OS.ReIm2WF;
                if OS.verbose == 1
                    fprintf('Phase Screen %s Applied\n',PS.name);
                end
            else
                error('Phasescreen must be an OptPhaseScreen object');
            end
            
        end
        
        function OS = PropagateSystem1(OS,starting_elem, ending_elem,n0)
            % OS = PropagateSystem1(OS, starting_elem, ending_elem)
            %
            % Propagate from ELEMENTS_{starting_elem} to
            % ELEMENTS_{ending_elem} using static method below for Fresnel
            % Propagation
            
            fprintf('\n\n\n');
            [numLambdas, propdists] = OS.initPropagation(starting_elem, ending_elem);
            
            if numLambdas == 1
                lambda = OS.lambda0_;
            else
                lambda = OS.lambda_array_;
            end
            
            % Propagation Limit
            proplim = 1.0e-4;
            
            % Internal storage index of refraction
            n0_ = n0;
            
            % Initialze directory name to save files with current time
            dirname = datestr(now,'dd_mm_yy_HH:MM');
            
            % Loop over given elements
            for ii = starting_elem:ending_elem
                if ii == starting_elem
                    fprintf('******************************************\n*    Starting at Element %s    *\n******************************************\n',OS.ELEMENTS_{ii}.name);
                end
                
                % The last element is different (it will have no propdist)
                if ii ~= ending_elem
                    
                    % Propagate
                    if OS.ELEMENTS_{ii}.propagation_method_ == 0
                        if abs(propdists(ii)) > proplim
                            fprintf('\nPropagating %g meters to element %s\n\n',propdists(ii), OS.ELEMENTS_{ii}.name);
                            
                            
                            % Do the Fresnel Propagation
                            OS.propagate2Elem(propdists(ii),OS.pscale_,lambda);
                            
                        else
                            % If the distance isn't large enough,
                            % don't propagate: OS.WF_ is unchanged
                            fprintf('%g meters is less than proplim. Copying Field\n\n',propdists(ii));
                            OS.ReIm2WF;
                        end
                        
                    elseif OS.ELEMENTS_{ii}.propagation_method_ == 1
                        fprintf('Computing Focused Field using FFT of WF00%d\n',ii-1);
                        D = 2*OS.beam_radius_ * OS.pscale_;
                        z = OS.f_number_ * D;
                        OS.set_field( OptSys.FraunhoferPropWF(OS.WF_,z,OS.pscale_,OS.lambda_array_));
                        
                        
                    elseif OS.ELEMENTS_{ii}.propagation_method_ == 2
                        warning('Not Supporting Zoom-FFTs yet\n');
                        fprintf('Implementing Zoom-FFTs for FPM. This is done by Applying the FPM\n');
                        
                    elseif OS.ELEMENTS_{ii}.propagation_method_ == 3
                        warning('Not Supporting Convolution with FPM yet\n');
                        fprintf('Implementing Convolution with Babinets Principle for FPM. This is done by Applying the FPM\n');
                        
                    end
                    
                    % Save
                    if OS.savefile == 1
                        %                         [OS.WFamp,OS.WFphase] = WFReIm2AmpPhase2(real(OS.WF_),imag(OS.WF_));
                        OptSys.saveWFfits(dirname,OS.WFamp, OS.WFphase, ii-1,'postprop');
                    end
                    
                    descr = OS.ELEMENTS_{ii}.describe;
                    fprintf('\n====================================================\n');
                    fprintf('Applying %s',descr);
                    fprintf('\n====================================================\n\n');
                    
                    
                    %******************************************************
                    % I haven't checked this, probably needs some
                    % refinement based on what material the OptFPM actually
                    % is. For now I am going to comment it out until I have
                    % thought about it more.
                    %******************************************************
                    
                    if ii > 1
                        % Check if the field has just exited something with
                        % a different index than the original n0
                        if isa(OS.ELEMENTS_{ii-1},'OptLens')
                            n0 = OS.ELEMENTS_{ii-1}.material_.n_;
                            %                         elseif isa(OS.ELEMENTS_{ii-1}, 'OptFPM')
                            %                             n0 = OS.ELEMENTS_{ii-1}.material_.n_;
                        else
                            n0 = n0_;
                        end
                    end
                    
                    
                    % Apply the element
                    OS.ApplyElement(ii,n0);
                    if OS.verbose == true
                        OS.show_PP;
                    end
                    
                    % Save
                    if OS.savefile == 1
                        OptSys.saveWFfits(dirname,OS.WFamp, OS.WFphase, ii-1,'preprop');
                    end
                    
                else % last element
                    
                    descr = OS.ELEMENTS_{ii}.describe;
                    fprintf('\n====================================================\n');
                    fprintf('************* Reached the Last Element *************\n');
                    fprintf('%s',descr);
                    fprintf('\n====================================================\n\n');
                    
                    % Propagate
                    if OS.ELEMENTS_{ii}.propagation_method_ == 0
                        if abs(propdists(ii)) > proplim
                            % Do the Fresnel Propagation
                            OS.propagate2Elem(propdists(ii),OS.pscale_,lambda);
                            
                        else
                            % If the distance isn't large enough,
                            % don't propagate: OS.WF_ is unchanged
                            OS.ReIm2WF;
                        end
                        
                    elseif OS.ELEMENTS_{ii}.propagation_method_ == 1
%                         fprintf('Computing Focused Field using FFT of WF00%d\n',ii-1);
%                         PSFa0 = OS.computePSF_FFT(OS.WF_);
%                         [thx,thy] = OS.FPcoords;
                        fprintf('Computing Focused Field using Chirp-Z Transform of WF00%d\n',ii-1);
                        if isa(OS.ELEMENTS_{ii},'OptDetector')
                            nld = OS.ELEMENTS_{ii}.FPregion_;
                        else
                            nld = 20;
                        end
                        PSFa0 = computePSF_chirpzDFT(OS,OS.WF_,ii);
                        [thx,thy] = OS.FPcoords(nld,length(PSFa0));

                        
                        figure;
                        for jj = 1:length(OS.lambda_array_)
%                             imagesc(xx(:,:,jj),yy(:,:,jj),PSFa0(:,:,jj));
%                             imagesc(xx(:,:,jj),yy(:,:,jj),log10(PSFa0(:,:,jj) / max(max(PSFa0(:,:,jj)))),[-5,0]);
%                             imagesc(PSFa0(:,:,jj));
%                             imagesc(log10(PSFa0(:,:,jj) / max(max(PSFa0(:,:,floor(length(OS.lambda_array_)/2))))),[-6,0]);
%                             imagesc(log10(PSFa0(:,:,jj) / max(max(PSFa0(:,:,jj)))),[-6,0]);
                            imagesc(thx(:,:,jj),thy(:,:,jj),log10(PSFa0(:,:,jj) / max(max(PSFa0(:,:,jj)))),[-6,0]);


                            plotUtils(sprintf('PSFa0,\n lambda = %g',OS.lambda_array_(jj)),'\lambda / D ','\lambda / D ');
                            drawnow;
                        end
                        
                        if OS.savefile == 1
                            OptSys.savePSFfits(dirname,OS.WFamp, OS.WFphase, ii);
                        end
                        
                    elseif OS.ELEMENTS_{ii}.propagation_method_ == 2
                        fprintf('Not Supporting Zoom-FFTs yet\n');
                        
                        
                    elseif OS.ELEMENTS_{ii}.propagation_method_ == 3
                        fprintf('Not Supporting Convolution with FPM yet\n');
                        
                        
                    end
                    
                end % of check if last element
                
                
            end % of loop over elements
            
            if isa((OS.WF_),'gpuArray')
                OS.GPU2CPU;
%                 OS.resetGPU;
            end
            
            
        end % of PropagateSystem1
        
        
        
        %% Remnant of the old code. Here for debugging new code purposes
        function OS = PropagateSystem4(OS, WFin, starting_elem, ending_elem)
            %OS = PropagateSystem4(OS, WFin, starting_elem, ending_elem)
            
            fprintf('\n');
            [numLambdas, propdist] = OS.initPropagation(starting_elem, ending_elem);
            
            if numLambdas == 1
                lambda = OS.lambda0_;
            else
                lambda = OS.lambda_array_;
            end
            
            size1 = size(WFin,1);
            size2 = size(WFin,2);
            numpix = size1*size2;
            FPMnext = false;
            
            % Propagation Limit
            proplim = 1.0e-4;
            
            % Initialze directory name to save files with current time
            dirname = datestr(now,'dd_mm_yy_HH:MM');
            
            % Loop over given elements
            for ii = starting_elem:ending_elem
                if ii == starting_elem
                    fprintf('Starting at Element %s\n',OS.ELEMENTS_{ii}.name);
                end
                
                % The last element is different (it will have no propdist)
                if ii ~= ending_elem
                    
                    type = OS.ELEMENTS_{ii}.type_;
                    descr = OS.ELEMENTS_{ii}.describe;
                    fprintf('\n====================================================\n');
                    fprintf('%s',descr);
                    fprintf('====================================================\n\n');
                    
                    switch type % have this if needed
                        case 0 % System Pupil
                            % Apply the element to the Input Wavefront
                            WFtmp = OS.ApplyElement(ii,WFin,lambda);
                            % Store Real and Imaginary Parts of Wavefront before propagation
                            WFrealpre = real(WFtmp);
                            WFimagpre = imag(WFtmp);
                            
                            % Convert to Amplitude and Phase
                            [OS.WFamp,OS.WFphase] = WFReIm2AmpPhase2(WFrealpre,WFimagpre);
                            
                            % Plot Propagated Wavefront
                            if OS.verbose == 1
                                % for jj = 1:numLambdas
                                jj = floor(numLambdas / 2);
                                figure(9999+ii)
                                subplot(1,2,1);
                                imagesc(OS.WFamp(:,:,jj));
                                plotUtils(sprintf('WF00%damp^-,\n lambda = %0.3f',ii-1,lambda(jj)*1e6));
                                drawnow;
                                % end
                            end
                            
                            % Save the system pupil
                            if ii == starting_elem
                                ind = 0;
                                if OS.savefile == 1
                                    OptSys.saveWFfits(dirname,OS.WFamp, OS.WFphase, ind);
                                end
                            end
                            
                            if OS.ELEMENTS_{ii}.isFocal == 0
                                if abs(propdist(ii)) > proplim
                                    % Do the Fresnel Propagation
                                    WFout = OptSys.FresnelPropagateWF(WFtmp,propdist(ii),OS.pscale_,lambda);
                                else
                                    % If the distance isn't large enough,
                                    % copy the wavefront
                                    WFout = WFtmp;
                                end
                                
                            elseif OS.ELEMENTS_{ii}.isFocal == 1
                                fprintf('Computing Focused Field using FFT of WF00%d\n',ii-1);
                                OS.computePSF_normalize(WFout);
                                WFout = OS.WF;
                                
                                if OS.savefile == 1
                                    OptSys.saveWFfits(dirname,OS.WFamp, OS.WFphase, ii);
                                end
                                
                            elseif OS.ELEMENTS_{ii}.isFocal == 2
                                % Next Element is a FPM
                                if OS.ELEMENTS_{ii+1}.type_ == 6
                                    FPMnext = true;
                                    fprintf('Going to Focal Plane Mask\n');
                                else
                                    FPMnext = false;
                                end
                            end
                            
                            
                            % Store Real and Imaginary Parts of Wavefront
                            WFrealpost = real(WFout);
                            WFimagpost = imag(WFout);
                            
                            % Convert to Amplitude and Phase
                            [OS.WFamp,OS.WFphase] = WFReIm2AmpPhase2(WFrealpost,WFimagpost);
                            
                            % Plot Propagated Wavefront
                            if OS.verbose == 1
                                % for jj = 1:numLambdas
                                jj = floor(numLambdas / 2);
                                figure(9999+ii)
                                subplot(1,2,2);
                                imagesc(OS.WFamp(:,:,jj));
                                plotUtils(sprintf('WF00%damp^+,\n lambda = %0.3f',ii,lambda(jj)*1e6));
                                drawnow;
                                % end
                            end
                            
                            % Save
                            if OS.savefile == 1
                                OptSys.saveWFfits(dirname,OS.WFamp, OS.WFphase, ii);
                            end
                            
                            % Move Calculated Wavefront to WFin to keep propagating
                            WFin = WFout;
                            OS.WF = WFout;
                            
                        case 1 % Lens
                            
                        case 2 % Mirror
                            
                        case 3 % Aspheric Lens
                            
                        case 4 % Aspheric Mirror
                            % Apply the element to the Input Wavefront
                            WFtmp = OS.ApplyElement(ii,WFin,lambda);
                            % Store Real and Imaginary Parts of Wavefront before propagation
                            WFrealpre = real(WFtmp);
                            WFimagpre = imag(WFtmp);
                            
                            % Convert to Amplitude and Phase
                            [OS.WFamp,OS.WFphase] = WFReIm2AmpPhase2(WFrealpre,WFimagpre);
                            
                            % Plot Propagated Wavefront
                            if OS.verbose == 1
                                % for jj = 1:numLambdas
                                jj = floor(numLambdas / 2);
                                figure(9999+ii)
                                subplot(1,2,1);
                                imagesc(OS.WFamp(:,:,jj));
                                plotUtils(sprintf('WF00%damp^-,\n lambda = %0.3f',ii-1,lambda(jj)*1e6));
                                drawnow;
                                % end
                            end
                            
                            % Save the system pupil
                            if ii == starting_elem
                                ind = 0;
                                if OS.savefile == 1
                                    OptSys.saveWFfits(dirname,OS.WFamp, OS.WFphase, ind);
                                end
                            end
                            
                            if OS.ELEMENTS_{ii}.isFocal == 0
                                if abs(propdist(ii)) > proplim
                                    % Do the Fresnel Propagation
                                    WFout = OptSys.FresnelPropagateWF(WFtmp,propdist(ii),OS.pscale_,lambda);
                                else
                                    % If the distance isn't large enough,
                                    % copy the wavefront
                                    WFout = WFtmp;
                                end
                                
                            elseif OS.ELEMENTS_{ii}.isFocal == 1
                                fprintf('Computing Focused Field using FFT of WF00%d\n',ii-1);
                                OS.computePSF_normalize(WFout);
                                WFout = OS.WF;
                                
                                if OS.savefile == 1
                                    OptSys.saveWFfits(dirname,OS.WFamp, OS.WFphase, ii);
                                end
                                
                            elseif OS.ELEMENTS_{ii}.isFocal == 2
                                % Next Element is a FPM
                                if OS.ELEMENTS_{ii+1}.type_ == 6
                                    FPMnext = true;
                                    fprintf('Going to Focal Plane Mask\n');
                                else
                                    FPMnext = false;
                                end
                            end
                            
                            
                            % Store Real and Imaginary Parts of Wavefront
                            WFrealpost = real(WFout);
                            WFimagpost = imag(WFout);
                            
                            % Convert to Amplitude and Phase
                            [OS.WFamp,OS.WFphase] = WFReIm2AmpPhase2(WFrealpost,WFimagpost);
                            
                            % Plot Propagated Wavefront
                            if OS.verbose == 1
                                % for jj = 1:numLambdas
                                jj = floor(numLambdas / 2);
                                figure(9999+ii)
                                subplot(1,2,2);
                                imagesc(OS.WFamp(:,:,jj));
                                plotUtils(sprintf('WF00%damp^+,\n lambda = %0.3f',ii,lambda(jj)*1e6));
                                drawnow;
                                % end
                            end
                            
                            % Save
                            if OS.savefile == 1
                                OptSys.saveWFfits(dirname,OS.WFamp, OS.WFphase, ii);
                            end
                            
                            % Move Calculated Wavefront to WFin to keep propagating
                            WFin = WFout;
                            OS.WF = WFout;
                            
                            
                        case 5 % Pupil Mask
                            % Apply the element to the Input Wavefront
                            WFtmp = OS.ApplyElement(ii,WFin,lambda);
                            % Store Real and Imaginary Parts of Wavefront before propagation
                            WFrealpre = real(WFtmp);
                            WFimagpre = imag(WFtmp);
                            
                            % Convert to Amplitude and Phase
                            [OS.WFamp,OS.WFphase] = WFReIm2AmpPhase2(WFrealpre,WFimagpre);
                            
                            % Plot Propagated Wavefront
                            if OS.verbose == 1
                                % for jj = 1:numLambdas
                                jj = floor(numLambdas / 2);
                                figure(9999+ii)
                                subplot(1,2,1);
                                imagesc(OS.WFamp(:,:,jj));
                                plotUtils(sprintf('WF00%damp^-,\n lambda = %0.3f',ii-1,lambda(jj)*1e6));
                                drawnow;
                                % end
                            end
                            
                            % Save the system pupil
                            if ii == starting_elem
                                ind = 0;
                                if OS.savefile == 1
                                    OptSys.saveWFfits(dirname,OS.WFamp, OS.WFphase, ind);
                                end
                            end
                            
                            if OS.ELEMENTS_{ii}.isFocal == 0
                                if abs(propdist(ii)) > proplim
                                    % Do the Fresnel Propagation
                                    WFout = OptSys.FresnelPropagateWF(WFtmp,propdist(ii),OS.pscale_,lambda);
                                else
                                    % If the distance isn't large enough,
                                    % copy the wavefront
                                    WFout = WFtmp;
                                end
                                
                            elseif OS.ELEMENTS_{ii}.isFocal == 1
                                fprintf('Computing Focused Field using FFT of WF00%d\n',ii-1);
                                OS.computePSF_normalize(WFout);
                                WFout = OS.WF;
                                
                                if OS.savefile == 1
                                    OptSys.saveWFfits(dirname,OS.WFamp, OS.WFphase, ii);
                                end
                            end
                            
                            % Store Real and Imaginary Parts of Wavefront
                            WFrealpost = real(WFout);
                            WFimagpost = imag(WFout);
                            
                            % Convert to Amplitude and Phase
                            [OS.WFamp,OS.WFphase] = WFReIm2AmpPhase2(WFrealpost,WFimagpost);
                            
                            % Plot Propagated Wavefront
                            if OS.verbose == 1
                                % for jj = 1:numLambdas
                                jj = floor(numLambdas / 2);
                                figure(9999+ii)
                                subplot(1,2,2);
                                imagesc(OS.WFamp(:,:,jj));
                                plotUtils(sprintf('WF00%damp^+,\n lambda = %0.3f',ii,lambda(jj)*1e6));
                                drawnow;
                                % end
                            end
                            
                            % Save
                            if OS.savefile == 1
                                    OptSys.saveWFfits(dirname,OS.WFamp, OS.WFphase, ii);
                            end
                            
                            % Move Calculated Wavefront to WFin to keep propagating
                            WFin = WFout;
                            OS.WF = WFout;
                            
                        case 6 % Focal Plane Mask
                            DFTgridpad = 1;
                            % Store Real and Imaginary Parts of Wavefront before propagation
                            WFrealpre = real(WFin);
                            WFimagpre = imag(WFin);
                            gsize = 2*DFTgridpad+1;
                            offset = DFTgridpad;
                            
                            
                        case 7 % Detector
                            
                        otherwise
                            error('Unknown Element type');
                    end
                    
                    
                    
                    % Last element
                else
                    
                    
                end % of if not last element
                
                
                
                
            end % loop over elements
            
            
        end % of PropagateSystem4
        
        %% GPU Methods
        
        function OS = initGPU(OS)
            
            OS.DEVICES = cell(OS.nGPUs,1);
            
            % Initialize Detected GPUs
            if OS.nGPUs == 1
                OS.DEVICES{1} = gpuDevice(1);
                OS.useGPU = 1;  % initialize to use GPU if there is one
                
            elseif OS.nGPUs == 2
                OS.DEVICES{2} = gpuDevice(2);
                OS.DEVICES{1} = gpuDevice(1);
                OS.nWorkers = 2;
                OS.useGPU = 1;  % initialize to use GPU if there is one
                
            elseif OS.nGPUS > 2
                fprintf('Only supports up to 2 GPUs\n');
                n = input('Please Select a GPU to use:  ');
                m = input('Please Select a second GPU to use:  ');
                OS.DEVICES{2} = gpuDevice(m);
                OS.DEVICES{1} = gpuDevice(n);
                OS.nGPUS = 2;
                OS.nWorkers = 2;
                OS.useGPU = 1;  % initialize to use GPU if there is one
                
            else
                warning('GPU:noGPU','No GPU to use!');
                OS.useGPU = 0;
%                 c = parcluster('local'); % build the 'local' cluster object
%                 OS.nWorkers = c.NumWorkers; % get the number of CPU cores from it
            end
            
        end % of initGPU
        
        function OS = GPUify(OS)
            % OS = GPUify(OS)
            
            % Check if GPU use is allowed
            if(OS.useGPU ~=1)
                warning('GPU:CPUselected', 'GPU is not being used. One may not exist. If you know it does, run initGPU()');
                return;
            end
            % Get the number of Elements
            numElems = OS.numElements_;
            
            % Make sure the field is of type float
            OS.setdatatype('single');
            
            if OS.nGPUs	 == 1
            % Send field to gpu
            OS.set_field(gpuArray(OS.WF_));
            
            for ii = 1:numElems
                % Make sure the zsag_ is of type float
                OS.ELEMENTS_{ii}.setdatatype('single');
                % Send zsag_ to the active gpu
                OS.ELEMENTS_{ii}.set_zsag(gpuArray(OS.ELEMENTS_{ii}.zsag_));
            end
            
            fprintf('\n\n\n***************************************************\n');
            fprintf('*         Now Using GPU %s        *\n',OS.DEVICES{1}.Name);
            fprintf('***************************************************\n');
%             warning('GPU:PROPNS','Propagation is currently pixel by pixel, and not a matrix multiply. This will be incredibly slow on GPU. Consider using CPU');
            elseif OS.nGPUs == 2
                
                % Send field to gpu
                OS.set_field(gpuArray(OS.WF_));
                
                for ii = 1:numElems
                    % Make sure the zsag_ is of type float
                    OS.ELEMENTS_{ii}.setdatatype('single');
                    % Send zsag_ to the active gpu
                    OS.ELEMENTS_{ii}.set_zsag(gpuArray(OS.ELEMENTS_{ii}.zsag_));
                end
                
                fprintf('\n\n\n***************************************************\n');
                fprintf('*         Now Using GPU %s        *\n',OS.DEVICES{1}.Name);
                fprintf('***************************************************\n');
                
                
            end
                
        end % of GPUify
        
        function OS = GPU2CPU(OS)
            % OS = GPU2CPU(OS)
            % Gather any gpuArrays back to the CPU
            % Get the number of Elements
            numElems = OS.numElements_;
            
            % Send field to gpu
            OS.set_field(gather(OS.WF_));
            OS.PSF_ = gather(OS.PSF_);
            
            for ii = 1:numElems
                % Send zsag_ to the active gpu
                OS.ELEMENTS_{ii}.set_zsag(gather(OS.ELEMENTS_{ii}.zsag_));
            end
            
        end % of GPU2CPU
            
        function mem = GPUavailableMem(OS)
            % mem = GPUavailableMem(OS)
            % Returns the available memory on the active GPU
            
            mem = OS.DEVICES{1}.AvailableMemory;
            
            
        end % of GPUavailableMem
        
        function OS = useCPU(OS)
            % useCPU(OS)
            % Tell code to not use GPUs. Initialize CPU workers
            OS.useGPU = 0;
%             c = parcluster('local'); % build the 'local' cluster object
%             OS.nWorkers = c.NumWorkers; % get the number of CPU cores from it
            
            % Gather any gpuArrays back to the CPU
            % Get the number of Elements
            numElems = OS.numElements_;
            
            % Send field to gpu
            OS.set_field(gather(OS.WF_));
            OS.PSF_ = gather(OS.PSF_);
            
            for ii = 1:numElems
                % Send zsag_ to the active gpu
                OS.ELEMENTS_{ii}.set_zsag(gather(OS.ELEMENTS_{ii}.zsag_));
            end
        end % useCPU
        
        function OS = resetGPU(OS)
            % OS = resetGPU(OS)
            % Clears GPU memory
            
            gpuDevice(1);
            fprintf('GPU memory reset\n');
            
        end % of resetGPU
        
    end % of methods
    
    %% Static Methods
    methods(Static = true)
        
        function [Shifted_im] = pixelShift(input, sx, sy)
            % [Shifted_im] = pixelShift(input, sx, sy)
            % Static method for shifting matrices by subpixel distances
            % Only works for shifts of 0<=shift<=1
            % Damages the phase, so probably not the best to do for FPM
            %
            % input is 2D matrix
            % sx is pixel shift in x
            % sy is pixel shift in y
            
            sz = size(input);
            if length(sz) == 2
                sz(3) = 1;
            end
            Shifted_im = zeros(sz(1),sz(2),sz(3));
            
            
            if sx < 0
                input = circshift(input,[0,1]);
                sx = -sx;
            end
            if sy < 0
                input = circshift(input,[1,0]);
                sy = -sy;
            end
            for ii = 1:sz(3)
                tmp = conv2(input(:,:,ii),[sx, 1-sx],'same');
                Shifted_im(:,:,ii) = conv2(tmp,[sy; 1-sy],'same');
            end
        end % of pixelShift
        
        function [output] = halfpixelShift(input,direction)
            % [ouput] = halfpixelShift(input, direction)
            % Multiplies input by the appropriate phase to shift it by half
            % a pixel in x and y when using an FFT
            
            sz = size(input);
            if length(sz) == 2
                sz(3) = 1;
            end
            
            shift = linspace(0,pi,sz(2)+1);
            shift = shift(1:end-1);
            [SX,SY] = meshgrid(shift);
            
            if direction == 1
                Tip = exp(-1i.* SX);
                Tilt = exp(-1i .* SY);
            elseif direction == -1
                Tip = exp(1i .* SX);
                Tilt = exp(1i .* SY);
            else
                error('direction is 1 for forward fft, -1 for ifft');
            end
            
            if isa(input,'gpuArray')
                flag = true;
                tmp = gather(input);
                datatype = class(tmp);
                clear tmp;
            else
                datatype = class(input);
                flag = false;
            end
            
            if strcmp(datatype,'single')
                output = single(zeros(sz(1),sz(2),sz(3)));
                if flag
                    output = gpuArray(output);
                end
            elseif strcmp(datatype,'double')
                output = double(zeros(sz(1),sz(2),sz(3)));
            elseif strcmp(datatype,'uint8')
                output = uint8(zeros(sz(1),sz(2),sz(3)));
            else
                error('Unsupported data type');
            end
            
            for ii = 1:sz(3)
                output(:,:,ii) = input(:,:,ii) .* Tip .* Tilt;
            end
            
            
        end % of halfpixelShift
        
        function [ WFout ] = FresnelProp( WFin, z, pscale, lambda )
            
            k = (2*pi)./lambda;
            [ny,nx,nlambda] = size(WFin);
            
            Lx = pscale*nx;
            Ly = pscale*ny;
            
            dfx = 1./Lx;
            dfy = 1./Ly;
            
            u = ones(nx,1)*((1:nx)-nx/2)*dfx;
            v = ((1:ny)-ny/2)'*ones(1,ny)*dfy;
            
            sqrdist = u.^2 + v.^2;
            coeff = pi*z*lambda;
            
            G = fftshift(fft2(fftshift(WFin)))* Lx * Ly;
            
            H = zeros(ny,nx,nlambda);
            for ii = 1:nlambda
                H(:,:,ii) = exp(1i.*k(ii)*z) .* exp((-1i .* coeff(ii)) .* sqrdist);
            end
            
            WFout = ifftshift(ifft2(ifftshift(G.*H))) / Lx / Ly;
            
            
        end

        function output = chirp_zDFT2(input,x0,dx,f0,df,M)
            
            sz = size(input);
            if length(sz) == 2
                sz(3) = 1;
            end
            
            N = numel(input(1,:));
            P = 2^ceil(log2(M+N-1));
            
            a = exp(-1i*2*pi*((1:N)*dx*(f0-df)+dx*df*(1:N).^2/2));
            b = exp(-1i*2*pi*((x0-dx)*(f0+(0:M-1)*df)+dx*df*(1:M).^2/2));
            
            if isa(input,'gpuArray')
                flag = true;
                tmp = gather(input);
                datatype = class(tmp);
                clear tmp;
            else
                datatype = class(input);
                flag = false;
            end
            
            if strcmp(datatype,'single')
                X1 = single(zeros(N,M,sz(3)));
                output = single(zeros(M,M,sz(3)));
                c = single((exp(1i*2*pi*dx*df*(1-N:M-1).^2/2)));
                if flag
%                     X1 = gpuArray(X1);
%                     output = gpuArray(output);
%                     a = gpuArray(a);
%                     b = gpuArray(b);
%                     c = gpuArray(c);
                end
            elseif strcmp(datatype,'double')
                X1 = double(zeros(sz(1),sz(2),sz(3)));
                output = double(zeros(M,M,sz(3)));
                c = double((exp(1i*2*pi*dx*df*(1-N:M-1).^2/2)));
            elseif strcmp(datatype,'uint8')
                X1 = uint8(zeros(sz(1),sz(2),sz(3)));
                output = uint8(zeros(M,M,sz(3)));
                c = uint8((exp(1i*2*pi*dx*df*(1-N:M-1).^2/2)));
            else
                error('Unsupported data type');
            end
            
             phi = fft(c,P);   
            for ii = 1:sz(1)
                tmp = input(ii,:);
                tmpX = ifft( fft(tmp.*a,P) .* phi );
                tmpX = tmpX(N:M+N-1) .* b;
                X1(ii,:) = tmpX;
            end
            
            % Permute
            X1 = X1.';
            
            % Do second dimmension
            
            for ii = 1:size(X1,1)
                tmp = X1(ii,:);
                tmpX = ifft( fft(tmp.*a,P) .* phi );
                
                tmpX = tmpX(N:M+N-1) .* b;
                output(ii,:) = tmpX;
            end

        
        end % of chirp_zDFT2
        
        % Put here so it can be used without needing the class object
        function WFout = FresnelPropagateWF(WFin, propdist, pscale, lambda)
            % WFout = FresnelPropagateWF(WFin, propdist, pscale, lambda)
            % Fresnel Propagates a wavefront
            
            % Initializations
            naxes0 = size(WFin,1);
            naxes1 = size(WFin,2);
            n0h = naxes0 / 2;
            n1h = naxes1 / 2;
            
            numLambdas = size(WFin,3);
            if numLambdas > 1
                WFout = zeros(naxes0,naxes1,numLambdas);
            end
            
            for kk = 1:numLambdas
                WFintmp = WFin(:,:,kk);
                
                coeff = pi*propdist*lambda(kk) / (pscale * naxes0) / (pscale * naxes1);

                % Propagate
                tmp = fftshift(fft2(fftshift(WFintmp))) * pscale * pscale;
                tmpre = real(tmp);
                tmpim = imag(tmp);
                for jj = 0:naxes1-1
                    jj1 = naxes0*jj;
                    jj2 = (jj - n1h)^2;
                    for ii = 0:naxes0-1
                        ii1 = jj1+ii;
                        ii2 = ii-n0h;
                        sqdist = ii2^2 + jj2;
                        angle = -coeff*sqdist;
                        re = tmpre(ii1+1);
                        im = tmpim(ii1+1);
                        tmpre(ii1+1) = re*cos(angle) - im*sin(angle);
                        tmpim(ii1+1) = re*sin(angle) + im*cos(angle);
                    end
                end
                tmp = tmpre + 1i*tmpim;
                WFout(:,:,kk) = ifftshift(ifft2(ifftshift(tmp))) * (1/pscale) * (1/pscale);
            end
            
        end % of FresnelPropagateWF
        
        function WFout = FraunhoferPropWF(WFin, propdist, pscale, lambda)
            % OS = FraunhoferPropWF(OS,WFin, propdist, pscale, lambda)

            sz = size(WFin);
            if length(sz) == 2
                sz(3) = 1;
            end
            
            k = (2*pi)./lambda;
            
            if isa(WFin,'gpuArray')
                flag = true;
                tmp = gather(WFin);
                datatype = class(tmp);
                clear tmp;
            else
                datatype = class(WFin);
                flag = false;
            end
            
            if strcmp(datatype,'single')
                WFfocus = single(zeros(size(WFin,1),size(WFin,2),sz(3)));
                if flag
                    WFfocus = gpuArray(WFfocus);
                end
            elseif strcmp(datatype,'double')
                WFfocus = double(zeros(size(WFin,1),size(WFin,2),sz(3)));
            elseif strcmp(datatype,'uint8')
                WFfocus = uint8(zeros(size(WFin,1),size(WFin,2),sz(3)));
            else
                error('Unsupported data type');
            end
            
            for ii = 1:length(lambda)
                coeff = ((exp(1i.*k(ii)*propdist))/(1i.*(lambda(ii))*propdist));
                WFfocus(:,:,ii) = coeff * fftshift(fft2(fftshift(OptSys.halfpixelShift( WFin(:,:,ii),1 ) ))) .* (sz(1).* sz(2) .* pscale .* pscale);
            end
            WFreal = real(WFfocus);
            WFimag = imag(WFfocus);
            [WFamp,WFphase] = WFReIm2AmpPhase2(WFreal,WFimag);
            WFout = WFamp .* exp(1i .* WFphase);

        end % FraunhoferProp
        
        
        function WFout = move2PP(WFin, propdist, pscale, lambda)
            % WFout = static_move2PP(OS,WFin)

            sz = size(WFin);
            if length(sz) == 2
                sz(3) = 1;
            end
            
            k = (2*pi)./lambda;
            
            if isa(WFin,'gpuArray')
                flag = true;
                tmp = gather(WFin);
                datatype = class(tmp);
                clear tmp;
            else
                datatype = class(WFin);
                flag = false;
            end
            
            if strcmp(datatype,'single')
                WFpupil = single(zeros(size(WFin,1),size(WFin,2),sz(3)));
                if flag
                    WFpupil = gpuArray(WFpupil);
                end
            elseif strcmp(datatype,'double')
                WFpupil = double(zeros(size(WFin,1),size(WFin,2),sz(3)));
            elseif strcmp(datatype,'uint8')
                WFpupil = uint8(zeros(size(WFin,1),size(WFin,2),sz(3)));
            else
                error('Unsupported data type');
            end
            
            
            
            for ii = 1:length(lambda)
                coeff = (((exp(1i.*k(ii)*propdist))/(1i.*(lambda(ii))*propdist)))^-1;
                WFpupil(:,:,ii) = coeff * ifftshift(ifft2(ifftshift(  WFin(:,:,ii) ))) .* (sz(1).* sz(2) .* (pscale) .* (pscale))^-1;
            end
            WFreal = real(WFpupil);
            WFimag = imag(WFpupil);
            [WFamp,WFphase] = WFReIm2AmpPhase2(WFreal,WFimag);
            WFout = WFamp .* exp(1i .* WFphase);

            
        end % move2PP
        
        function saveWFfits(dirname,WFamp,WFphase,ind,location)
            % saveWFfits(dirname,WFamp,WFphase,ind,location)
            % saves WFamp and WFphase as fits files in Propfile directory.
            % If the directory doesn't exist, it is created.
            
            if nargin < 5
                location = '';
            end
            
            direc = sprintf('Propfiles/%s',dirname);
            val = exist(direc, 'dir');
            if val == 7
                fprintf('****************************************\nSaving files to %s\n****************************************\n\n',direc);
                filenameamp = sprintf('WFamp00%d%s.fits',ind,location);
                filenamephase = sprintf('WFphase00%d%s.fits',ind,location);
                current_dir = pwd;
                newdir = sprintf('%s/%s',current_dir,direc);
                cd(newdir);
                fitswrite(WFamp,filenameamp);
                fitswrite(WFphase,filenamephase);
                cd(current_dir);
            else
                mkdir(direc)
                val = exist(direc, 'dir');
                if val == 7
                    fprintf('****************************************\nSaving files to %s\n****************************************\n\n',direc);
                    filenameamp = sprintf('WFamp00%d%s.fits',ind,location);
                    filenamephase = sprintf('WFphase00%d%s.fits',ind,location);
                    current_dir = pwd;
                    newdir = sprintf('%s/%s',current_dir,direc);
                    cd(newdir);
                    fitswrite(WFamp,filenameamp);
                    fitswrite(WFphase,filenamephase);
                    cd(current_dir);
                else
                    error('Problem Creating %s Directory',direc);
                end
            end
        end % of saveWFfits
        
        function savePSFfits(dirname,PSFamp,PSFphase,ind)
            % saveWFfits(dirname,PSFamp,PSFphase,ind)
            % saves PSFamp and PSFphase as fits files in Propfile directory.
            % If the directory doesn't exist, it is created.

            direc = sprintf('Propfiles/%s',dirname);
            val = exist(direc, 'dir');
            if val == 7
                fprintf('****************************************\nSaving files to %s\n****************************************\n\n',direc);
                filenameamp = sprintf('PSFamp00%d.fits',ind);
                filenamephase = sprintf('PSFphase00%d.fits',ind);
                filenamepsf = sprintf('PSF00%d.fits',ind);
                current_dir = pwd;
                newdir = sprintf('%s/%s',current_dir,direc);
                cd(newdir);
                fitswrite(PSFamp,filenameamp);
                fitswrite(PSFphase,filenamephase);
                fitswrite(abs(PSFamp .* exp(1i*PSFphase)).^2,filenamepsf);
                cd(current_dir);
            else
                mkdir(direc)
                val = exist(direc, 'dir');
                if val == 7
                    fprintf('****************************************\nSaving files to %s\n****************************************\n\n',direc);
                    filenameamp = sprintf('PSFamp00%d.fits',ind);
                    filenamephase = sprintf('PSFphase00%d.fits',ind);
                    filenamepsf = sprintf('PSF00%d.fits',ind);
                    current_dir = pwd;
                    newdir = sprintf('%s/%s',current_dir,direc);
                    cd(newdir);
                    fitswrite(PSFamp,filenameamp);
                    fitswrite(PSFphase,filenamephase);
                    fitswrite(abs(PSFamp .* exp(-1i*PSFphase)).^2,filenamepsf);
                    cd(current_dir);
                else
                    error('Problem Creating %s Directory',direc);
                end
            end
        end % of savePSFfits
        
    end % of static methods
end
