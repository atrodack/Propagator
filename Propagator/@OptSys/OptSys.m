classdef OptSys < matlab.mixin.Copyable
    %OptSys Class for storing optical system parameters
    %   Detailed explanation goes here
    
    
    %% Properties
    properties(Constant=true, GetAccess='protected')
        
        
    end % of constant properties
    
    properties(GetAccess='public',SetAccess='public')
        name; % give a name to the structure/system
        savefile = 1;
        verbose = 1; % print extra info
        
    end % of public properties
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        
        % System Parameters
        numElements_ = 0; % number of elements in system
        f_number_; % F/# of system
        conjugations_ = zeros(1,1); % list of z locations
        pscale_; % size of pixels in meters
        gridsize_; %
        lambda0_; % central wavelength
        lambda_array_; % vector of lambdas to use
        apertureStop; % probably not needed
        
        % List of Elements
        ELEMENTS_; % cell array of element cards
        % Each element is of type OptElement
        % This includes the information:
        % 1) Element type
        % 2) Material (index of refraction)
        % 3) Element Position
        % 4) z-sag in meters
        % 5) Type of propagation to use (could change)
           
        % Wavefront
        WFamp;
        WFphase;
        WFreal;
        WFimag;
        WF;
        
        % GPU properties
        useGPU; % use a GPU flag
        nGPUs;
        DEVICES;
        
        % Pool Characteristics
        nWorkers; % number of workers for parpool
        
        
        % Misc
        interpolate_method = [];
        
        
    end % of protected properties
    
    %% Methods
    methods
        
        %% Constructor
        
        function OS = OptSys(sz)
            % OS = OptSys(sz)
            % Constructs OptSys object. Initializes WFCA to size sz
            
            if nargin < 1
                sz = 1024;
            end
            
            OS.initOptSys(sz);
            
            
        end % of contructor
        
        
        %% Set Properties
        
        function OS = initOptSys(OS,sz)
            %OS = initializeOptSys(OS,sz)
            
            OS.ELEMENTS_ = cell(OS.numElements_,1);
            OS.nGPUs = gpuDeviceCount;
            if OS.nGPUs > 0
                OS.initGPU;
            end
            OS.setGridsize(sz);
            
        end % of initOptSys
        
        
        function OS = addSequentialElement(OS,elem)
            % OS = addSequentialElement(OS,elem)
            % adds andElement object to the optical system located after
            % the previous stored element
            
            if isa(elem,'OptElement')
                OS.numElements_ = OS.numElements_ + 1;
                OS.ELEMENTS_{OS.numElements_,1} = elem;
%                 OS.setConjugations; %maybe add in later
                
                if OS.verbose == 1
                    fprintf('Element %s Added to Optical System\n',elem.name);
                end
                
            else
                error('I do not understand anything but OptElement objects');
            end
        end % of addSequentialElement
        
        
        function OS = addElement(OS,elem,order_num)
            % OS = addElement(OS,elem,order_num)
            % Adds an Element to the optical system in the location given
            % in order_num
            
            if isa(elem,'OptElement') == 0
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
                
        
        function OS = setPscale(OS,val)
            % OS = setPscale(OS,val)
            OS.pscale_ = val;
        end % of setPscale
        
        
        function OS = setLambdaarray(OS,lambdalist)
            % OS = setLambdaarray(OS,lambdalist)
            
            OS.lambda_array_ = lambdalist;
        end % of setLambdaarray
        
        
        function OS = setCentralWavelength(OS,lambda)
            % OS = setCentralWavelength(OS,lambda)
            
            OS.lambda0_ = lambda;
        end % of setCentralWavelength
        
        
        function conjugations = setConjugations(OS)
            % OS = setConjugations(OS)
            % Adds the zpos_ of an added element to the conjugations list
            
            conjugations = zeros(OS.numElements_,1);
            for z = 1:OS.numElements_
                conjugations(z,1) = OS.ELEMENTS_{z}.z_position_;
            end
            OS.conjugations_ = conjugations;
            
        end % of setConjugations
        
        
        function OS = getElementFNum(OS,element_num)
            % OS = getElementFNum(OS,element_num)
            % sets the f/# of the system to that of an element
            
            OS.f_number_ = OS.ELEMENTS_{element_num}.getFNumber;
        end % of getElementFNum
        
        
        function OS = setGridsize(OS,sz)
            % gridsize = setGridsize(OS,sz)
            
            if nargin < 2
                sz = size(OS.ELEMENTS_{1}.zsag_,1);
            end
            OS.gridsize_ = sz;
        end % of setGridsize
           
        
        %% Propagation Methods
        
        function [numLambdas,propdist] = initPropagation(OS,starting_elem, ending_elem)
            % [numLambdas,propdist] = initPropagation(OS,starting_elem, ending_elem)
            % Runs checks on OS object to ensure properties needed to
            % propagate through system are set.
            %
            % Returns number of lambdas and list of propagation distances
            % for use in PropagateSystem method
            
            OK = false(3,1);
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
            
            % Error out if a Check fails
            if prod(OK) == 0
                fprintf('\n\n');
                error('Check failed. See above warnings for help');
            end
            
            
            % Get the z_position list
            if isempty(OS.conjugations_) == 0
                if length(OS.conjugations_) == OS.numElements_
                    z_list = OS.conjugations_; % copy if already cached
                else
                    z_list = OS.setConjugations; % cache the values
                end
            end
            propdist = zeros(length(z_list)-1,1);
            for ii = 1:length(z_list)-1
                propdist(ii) = z_list(ii+1) - z_list(ii);
            end
            
        end % initPropagation
         
        
        function planewave = planewave(OS,amplitude,theta)
            % OS = planewave(OS,amplitude,theta)
            % Initializes a planewave as the stored wavefront
            %
            % TO DO: add off-axis planewave functionality
            
            % If no input, just use ones
            if nargin < 2
                amplitude = complex(1);
            end
            
            % Theta is the off axis angle of the planewave (not supported
            % yet)
            if nargin < 3
                on_axis = true;
            else
                on_axis = false;
            end
            
            if on_axis
                OS.WFamp = amplitude .* ones(OS.gridsize_);
                OS.WFphase = zeros(OS.gridsize_);
                OS.AmpPhase2WF;
                planewave = OS.WF;
            else  % off-axis
                error('Not supported yet'); 
            end
        end % of planewave
        
        
        function OS = AmpPhase2WF(OS,ind)
            % OS = AmpPhase2WF(OS)
            % Sets the WF as the combination of the amplitude and phase
            % components
            if nargin < 2
                ind = 1;
            end
            
            OS.WF(:,:,ind) = OS.WFamp(:,:,ind) .* exp(-1i * OS.WFphase(:,:,ind));
        end % of AmpPhase2WF
        
        
        function [WFreal,WFimag] = ReImWF(OS,WFin)
            % OS = ReImWF(OS,WFin)
            % Finds the real and imaginary components of a wavefront. If no
            % wavefront is given, acts on wavefront stored in OptSys
            % object, and stores results in object
            
            if nargin < 2
                if isempty(OS.WF) == 0
                    OS.WFreal = real(OS.WF);
                    OS.WFimag = imag(OS.WF);
                else
                    error('No wavefront stored in Object');
                end
            else
                WFreal = real(WFin);
                WFimag = imag(WFin);
            end
        end % of ReImWF
        
        
        function [WFreal,WFimag] = WFAmpPhase2ReIm(OS,WFamp,WFphase)
            % [WFreal,WFimag] = WFAmpPhase2ReIm(OS,WFamp,WFphase)
            % sets the real and imaginary parts of the wavefront from
            % amplitude and phase components. If no amplitude and phase
            % components are given, act on stored ones, and store resutls
            
            if nargin < 2
                sz = size(OS.WF);
                for ii = 1:sz(1)
                    for jj = 1:sz(2)
                        amp = OS.WFamp(ii,jj);
                        phase = OS.WFphase(ii,jj);
                        OS.WFreal(ii,jj) = amp*cos(phase);
                        OS.WFimag(ii,jj) = amp*sin(phase);
                    end
                end
            elseif nargin == 3
                sz = size(WFamp);
                WFreal = zeros(sz(1),sz(2));
                WFimag = zeros(sz(1),sz(2));
                for ii = 1:sz(1)
                    for jj = 1:sz(2)
                        amp = WFamp(ii,jj);
                        phase = WFphase(ii,jj);
                        WFreal(ii,jj) = amp*cos(phase);
                        WFimag(ii,jj) = amp*sin(phase);
                    end
                end
            end
        end % of WF2ReIm
        
        
        function [WFamp,WFphase] = WFReIm2AmpPhase(OS,WFreal,WFimag)
            % [WFamp,WFphase] = WFReIm2AmpPhase(OS,WFreal,WFimag)
            % Converts real and imaginary parts of WF to amplitude and
            % phase.
            
            if nargin < 2
                sz = size(OS.WF);
                for ii = 1:sz(1)
                    for jj = 1:sz(2)
                        re = OS.WFreal(ii,jj);
                        im = OS.WFimag(ii,jj);
                        amp = sqrt(re*re + im*im);
                        phase = atan2(im,re);
                        OS.WFamp(ii,jj) = amp;
                        OS.WFphase(ii,jj) = phase;
                    end
                end
            elseif nargin == 3
                sz = size(WFreal);
                WFamp = zeros(sz(1),sz(2));
                WFphase = zeros(sz(1),sz(2));
                for ii = 1:sz(1)
                    for jj = 1:sz(2)
                        re = WFreal(ii,jj);
                        im = WFimag(ii,jj);
                        amp = sqrt(re*re + im*im);
                        phase = atan2(im,re);
                        WFamp(ii,jj) = amp;
                        WFphase(ii,jj) = phase;
                    end
                end
            end
                
        end % of WFReIm2AmpPhase
        
        function WFout = ApplyElement(OS,elem_num,WFin,lambda)
            % OS = ApplyElement(OS,elem_num,WFin)
            % % Applys element to current wavefront
            
            if nargin < 3
                error('Must include the number of the element to apply (elem_num) and the complex wavefront to apply it to');
            elseif nargin < 4
                lambda = OS.lambda0_;
            end
            
            % Store element
            elem = OS.ELEMENTS_{elem_num};
            
            % Get the type of element
            type = elem.type_;
            
            switch type
                case 0 % System Pupil
                    if elem.amponly ~= 0 % amplitude mask
                        WFout = WFin .* elem.zsag_;
                    else % phase mask
                        if isempty(elem.phasefac_) == 1
                            elem.set_phasefactor(lambda);
                        end
                        WFout = WFin .* exp(-1i * elem.phasefac_ .* elem.zsag_);
                    end
                    
                case 1 % Lens
                    if isempty(elem.phasefac_) == 1
                        elem.set_phasefactor(lambda);
                    end
                    WFout = WFin .* exp(-1i * elem.phasefac_ .* elem.zsag_);
                    
                case 2 % Mirror
                    if isempty(elem.phasefac_) == 1
                        elem.set_phasefactor(lambda);
                    end
                    WFout = WFin .* exp(-1i * elem.phasefac_ .* elem.zsag_);
                    
                case 3 % Aspheric Lens
                    if isempty(elem.phasefac_) == 1
                        elem.set_phasefactor(lambda);
                    end
                    WFout = WFin .* exp(-1i * elem.phasefac_ .* elem.zsag_);
                    
                case 4 % Aspheric Mirror
                    if isempty(elem.phasefac_) == 1
                        elem.set_phasefactor(lambda);
                    end
                    WFout = WFin .* exp(-1i * elem.phasefac_ .* elem.zsag_);
                    
                case 5 % Pupil Mask
                    if elem.amponly ~= 0 % amplitude mask
                        WFout = WFin .* elem.zsag_;
                    else % phase mask
                        if isempty(elem.phasefac_) == 1
                            elem.set_phasefactor(lambda);
                        end
                        WFout = WFin .* exp(-1i * elem.phasefac_ .* elem.zsag_);
                    end
                    
                case 6 % Focal Plane Mask
                    if elem.amponly ~= 0 % amplitude mask
                        WFout = WFin .* elem.zsag_;
                    else % phase mask
                        if isempty(elem.phasefac_) == 1
                            elem.set_phasefactor(lambda);
                        end
                        WFout = WFin .* exp(-1i * elem.phasefac_ .* elem.zsag_);
                    end
                    
                case 7 % Detector
                    
                otherwise
                    error('Unknown Element type');
            end
            
            
        end % of ApplyElement
        
        
        function OS = computePSF(OS,WFin)
            % OS = computePSF(OS,WFin)
            
            sz = size(WFin);
            
            for jj = 1:sz(3)
                WFfocus(:,:,jj) = fftshift(fft2(fftshift(WFin(:,:,jj))));
                WFreal(:,:,jj) = real(WFfocus(:,:,jj));
                WFimag(:,:,jj) = imag(WFfocus(:,:,jj));
                [OS.WFamp(:,:,jj),OS.WFphase(:,:,jj)] = OS.WFReIm2AmpPhase(WFreal(:,:,jj),WFimag(:,:,jj));
                OS.AmpPhase2WF(jj);
                
                if OS.verbose == 1
                            figure(9999)
                            imagesc(OS.WFamp(:,:,jj));
                            axis square;
                            axis xy;
                            colorbar;
                            colormap(gray(256));
                            title(sprintf('WF00%damp, lambda = %g',jj,OS.lambda_array_(jj)));
                            drawnow;
                end
            end
            
        end % computePSF
        
        
        function OS = PropagateSystem(OS, WFin, starting_elem, ending_elem)
            % OS = PropagateSystem(OS, WFin, starting_elem, ending_elem)
            % Propagates WFin through the elements [starting_elem,
            % starting_elem + 1, starting_elem +2, ..., ending_elem].
            %
            % WFin should be a complex amplitude (amp .* exp(-i*phase))
 
            
            [numLambdas, propdist] = OS.initPropagation(starting_elem, ending_elem);
            
            if numLambdas == 1
                lambda = OS.lambda0_;
            else
                lambda = OS.lambda_array_;
            end
            
            WFout = zeros(size(WFin,1),size(WFin,2),numLambdas);
            WFtmp = WFout;
            
            for ii = starting_elem:ending_elem
                if ii ~= ending_elem
                    for jj = 1:numLambdas
                        WFtmp(:,:,jj) = OS.ApplyElement(ii,WFin(:,:,jj),lambda(jj));
                        WFout(:,:,jj) = OptSys.FresnelPropagateWF(WFtmp(:,:,jj),propdist(ii),OS.pscale_,lambda(jj));
                        
                        WFreal(:,:,jj) = real(WFout(:,:,jj));
                        WFimag(:,:,jj) = imag(WFout(:,:,jj));
                        [OS.WFamp(:,:,jj),OS.WFphase(:,:,jj)] = OS.WFReIm2AmpPhase(WFreal(:,:,jj),WFimag(:,:,jj));
                        
                        if OS.verbose == 1
                            figure(9999)
                            imagesc(OS.WFamp(:,:,jj));
                            axis square;
                            axis xy;
                            colorbar;
                            colormap(gray(256));
                            title(sprintf('WF00%damp, lambda = %g',ii,lambda(jj)));
                            drawnow;
                        end
                    end
                else % last element
                    if OS.ELEMENTS_{ii}.isFocal == true
                        fprintf('Computing PSF using FFT of WF0%d\n',ii);
                        OS.computePSF(WFout);
                    else
                        fprintf('At last element\n');
                    end
                end
                
                
                if OS.savefile == 1
                    OptSys.saveWFfits(OS.WFamp, OS.WFphase, ii);
                end
                
                % Move Calculated Wavefront to WFin to keep propagating
                WFin = WFout;
            end
            
            

        
                    
            
        end % of PropagateSystem
        
        
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
                warning('GPU:noGPU','No GPU to use!\n');
                OS.useGPU = 0;
                c = parcluster('local'); % build the 'local' cluster object
                OS.nWorkers = c.NumWorkers % get the number of CPU cores from it
            end
            
        end % of initGPU

        
        function mem = GPUavailableMem(OS)
            % mem = GPUavailableMem(OS)
            % Returns the available memory on the active GPU
            
            mem = OS.DEVICES{1}.AvailableMemory;
                
            
        end % of GPUavailableMem
        
        function OS = useCPU(OS)
            % useCPU(OS)
            % Tell code to not use GPUs. Initialize CPU workers
            
            OS.useGPU = 0;
            c = parcluster('local'); % build the 'local' cluster object
            OS.nWorkers = c.NumWorkers; % get the number of CPU cores from it
        end % useCPU    
        
    end % of methods
    
    %% Static Methods
    methods(Static = true)
        
        function WFout = FresnelPropagateWF(WFin, propdist, pscale, lambda)
            % WFout = FresnelPropagateWF(WFin, propdist, pscale, lambda)
            % Fresnel Propagates a wavefront
        
            % Initializations
            naxes0 = size(WFin,1);
            naxes1 = size(WFin,2);
            n0h = naxes0 / 2;
            n1h = naxes1 / 2;
            
            coeff = pi*propdist*lambda / (pscale * naxes0) / (pscale * naxes1);
            
            % Propagate
            tmp = fftshift(fft2(fftshift(WFin))) * pscale * pscale;
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
            WFout = ifftshift(ifft2(ifftshift(tmp))) * (1/pscale) * (1/pscale);
            
        end % of FresnelPropagateWF
        
        function saveWFfits(WFamp,WFphase,ind)
            % saveWFfits(WFamp,WFphase,ind)
            % saves WFamp and WFphase as fits files in Propfile directory.
            % If the directory doesn't exist, it is created.
            
            val = exist('Propfiles', 'dir');
            if val == 7
                filenameamp = sprintf('WF0%d_amp.fits',ind);
                filenamephase = sprintf('WF0%d_phase.fits',ind);
                current_dir = pwd;
                newdir = sprintf('%s/Propfiles',current_dir);
                cd(newdir);
                fitswrite(WFamp,filenameamp);
                fitswrite(WFphase,filenamephase);
                cd(current_dir);
            else
                mkdir Propfiles
                val = exist('Propfiles', 'dir');
                if val == 7
                    filenameamp = sprintf('WF0%d_amp.fits',ind);
                    filenamephase = sprintf('WF0%d_phase.fits',ind);
                    current_dir = pwd;
                    newdir = sprintf('%s/Propfiles',current_dir);
                    cd(newdir);
                    fitswrite(WFamp,filenameamp);
                    fitswrite(WFphase,filenamephase);
                    cd(current_dir);
                else
                    error('Problem Creating Propfiles Directory');
                end
            end
        end % of saveWFfits
            
        
    end % of static methods
end
