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
        f_number_; % F/# of system
        conjugations_ = zeros(1,1); % list of z locations
        pscale_; % size of pixels in meters
        gridsize_; %
        lambda0_; % central wavelength
        lambda_array_; % vector of lambdas to use
        
        % List of Elements
        ELEMENTS_; % cell array of element cards
        % Each element is of type OptElement
        % This includes the information:
        % 1) Element type
        % 2) Material (index of refraction)
        % 3) Element Position
        % 4) z-sag in meters
        
        % Wavefront
        WFamp;
        WFphase;
        WF_;
        
        % GPU properties
        useGPU; % use a GPU flag
        nGPUs;
        DEVICES;
        
        % Pool Characteristics
        nWorkers; % number of workers for parpool
        
        
        % Misc
        interpolate_method = [];
        default_data_type = 'single';
        
        
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
            %OS = initializeOptSys(OS,sz)
            
            OS.ELEMENTS_ = cell(OS.numElements_,1);
            OS.nGPUs = gpuDeviceCount;
            if OS.nGPUs > 0
                OS.initGPU;
            end
            OS.setGridsize(sz);
            
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
                %                 OS.setConjugations; %maybe add in later
                
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
        
        %% Property Setting Methods
        
        function OS = setPscale(OS,val)
            % OS = setPscale(OS,val)
            OS.pscale_ = val;
        end % of setPscale
        
        
        function OS = setLambdaarray(OS,lambdalist)
            % OS = setLambdaarray(OS,lambdalist)
            
            if nargin < 2
                if isempty(OS.lambda0_) == 0
                    OS.lambda_array_ = OS.lambda0_;
                else
                    error('Need to include lambdalist or set central wavelength');
                end
            else
                OS.lambda_array_ = lambdalist;
                if isempty(OS.lambda0_)
                    OS.setCentralWavelength();
                end
            end
        end % of setLambdaarray
        
        
        function OS = setCentralWavelength(OS,lambda)
            % OS = setCentralWavelength(OS,lambda)
            
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
                    OS.setLambdaarray(lambda);
                end
                
            end
        end % of setCentralWavelength
        
        function OS = setField(OS,field)
            % OS = setField(OS,field)
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
        end % of setField
            
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
        
        
        function gridsize = setGridsize(OS,sz)
            % gridsize = setGridsize(OS,sz)
            
            if nargin < 2
                gridsize = size(OS.ELEMENTS_{1}.zsag_);
            else
                gridsize = sz;
            end
            
            OS.gridsize_ = gridsize;
            
        end % of setGridsize
        
        function show(OS)
        % show(OS)
        % Plots matrix stored in zsag_
        numLambdas = size(OS.lambda_array_,2);
        
        if isreal(OS.WF_) == 1
            figure;
            for ii = 1:numLambdas
                imagesc(OS.WF_(:,:,ii))
                plotUtils(sprintf('WF lambda %d',ii));
                drawnow;
                pause(0.15);
            end
        else
            figure;
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
        
        
        function OS = planewave(OS,amplitude,numLambdas,theta)
            % OS = planewave(OS,amplitude,theta)
            % Initializes a planewave as the stored wavefront
            %
            % TO DO: add off-axis planewave functionality
            
            % If no input, just use ones
            if nargin < 2
                amplitude = complex(1);
                if isempty(OS.lambda_array_)
                    numLambdas = 1;
                else
                    numLambdas = size(OS.lambda_array_,2);
                end
            end
            
            % If no number of lambdas, just use one
            if nargin < 3
                if isempty(OS.lambda_array_)
                    numLambdas = 1;
                else
                    numLambdas = size(OS.lambda_array_,2);
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
                if strcmp( OS.default_data_type, 'single')
                    amplitude = single(amplitude);
                    WF = amplitude .* ones(OS.gridsize_(1),OS.gridsize_(2));
                    field = single(zeros(size(WF,1),size(WF,2),numLambdas));
                elseif strcmp( OS.default_data_type, 'double')
                    amplitude = double(amplitude);
                    WF = amplitude .* ones(OS.gridsize_(1),OS.gridsize_(2));
                    field = double(zeros(size(WF,1),size(WF,2),numLambdas));
                elseif strcmp( OS.default_data_type, 'uint8')
                    amplitude = uint8(amplitude);
                    WF = amplitude .* ones(OS.gridsize_(1),OS.gridsize_(2));
                    field = uint8(zeros(size(WF,1),size(WF,2),numLambdas));
                else
                    error('Data type not supported (yet)');
                end
                
                for ii = 1:numLambdas
                    field(:,:,ii) = WF;
                end;
                OS.WF_ = field;
                
            else  % off-axis
                error('Not supported yet');
            end
        end % of planewave
        
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
                    z_list = OS.setConjugations();
                else
                    z_list = OS.conjugations_;
                end
            end
            
            propdists = zeros(length(z_list),1);
            for ii = 2:length(z_list)
                propdists(ii) = z_list(ii) - z_list(ii-1);
            end
        
        end % of computePropdists
        
        
        % Switch these to using internal WF_, not input!
        function OS = computePSF(OS,WFin)
            % OS = computePSF(OS,WFin)
            
            sz = size(WFin);
            if length(sz) == 2
                sz(3) = 1;
            end;
            
            % Initialize
            WFfocus = zeros(sz(1),sz(2),sz(3));
            WFreal = WFfocus;
            WFimag = WFfocus;
            
            for jj = 1:sz(3)
                WFfocus(:,:,jj) = fftshift(fft2(fftshift(WFin(:,:,jj)))) .* (sz(1).* sz(2) .* OS.pscale_ .* OS.pscale_);
                WFreal(:,:,jj) = real(WFfocus(:,:,jj));
                WFimag(:,:,jj) = imag(WFfocus(:,:,jj));
                [OS.WFamp(:,:,jj),OS.WFphase(:,:,jj)] = WFReIm2AmpPhase(WFreal(:,:,jj),WFimag(:,:,jj));
                OS.AmpPhase2WF(jj);
                psfa0 = abs(OS.WF_).^2;
                
                % Plot
                if OS.verbose == 1
                    figure(99999)
                    imagesc(psfa0(:,:,jj));
                    plotUtils(sprintf('PSFa0,\n lambda = %g',OS.lambda_array_(jj)));
                    drawnow;
                end
            end
            
        end % computePSF
        
        function OS = computePSF_normalize(OS,WFin)
            % OS = computePSF(OS,WFin)
            
            if nargin < 2
                WFin = OS.WF_;
            end
               
            sz = size(WFin);
            if length(sz) == 2
                sz(3) = 1;
            end
            
            WFfocus = fftshift(fft2(fftshift(WFin))) .* (sz(1).* sz(2) .* OS.pscale_ .* OS.pscale_);
            WFreal = real(WFfocus);
            WFimag = imag(WFfocus);
            [OS.WFamp,OS.WFphase] = WFReIm2AmpPhase2(WFreal,WFimag);
            OS.AmpPhase2WF();
            psfa0 = abs(OS.WF_).^2;

            % Plot
            if OS.verbose == 1
                for jj = 1:sz(3)
                    normalizer = max(max(psfa0(:,:,jj)));
                    figure(99999)
                    imagesc(log10(psfa0(:,:,jj)/normalizer),[-5,0]);
                    plotUtils(sprintf('PSFa0,\n lambda = %g',OS.lambda_array_(jj)));
                    drawnow;
                end
            end
            
        end % computePSF_normalize
        
        
        
        %% System Propagation Methods
        % See static methods for actual Fresnel propagation method
        
        function OS = propagate2Elem(OS, propdist, pscale, lambda)
            % OS = propagate2Elem(OS,propdist,pscale,lambda)
            OS.setField(OptSys.FresnelPropagateWF(OS.WF_,propdist,pscale,lambda));
            OS.ReIm2WF;
            
        end % of propagate2Elem
        
        function OS = ApplyElement(OS,elem_num,n0)
            % OS = ApplyElement(OS,elem_num,n0)
            
            if isa(OS.ELEMENTS_{elem_num},'OptFPM')
                        % The FPM needs access to the pscale property
                        OS.setField( OS.ELEMENTS_{elem_num}.ApplyElement(OS.WF_,OS.lambda_array_,n0,OS.pscale_));
                        OS.ReIm2WF;
                    else
                        OS.setField( OS.ELEMENTS_{elem_num}.ApplyElement(OS.WF_,OS.lambda_array_,n0));
                        OS.ReIm2WF;
            end
        end % of ApplyElement
            
        
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
                    if OS.ELEMENTS_{ii}.isFocal_ == 0
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
                        
                    elseif OS.ELEMENTS_{ii}.isFocal_ == 1
                        fprintf('Computing Focused Field using FFT of WF00%d\n',ii-1);
                        OS.setField( OptSys.static_computePSF(OS.WF_,OS.pscale_));

                        
                    elseif OS.ELEMENTS_{ii}.isFocal_ == 2
                        warning('Not Supporting Zoom-FFTs yet\n');
                        fprintf('Implementing Zoom-FFTs for FPM. This is done by Applying the FPM\n');
                        
                    elseif OS.ELEMENTS_{ii}.isFocal_ == 3
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
                        OS.show;
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
                    if OS.ELEMENTS_{ii}.isFocal_ == 0
                        if abs(propdists(ii)) > proplim
                            % Do the Fresnel Propagation
                            OS.propagate2Elem(propdists(ii),OS.pscale_,lambda);

                        else
                            % If the distance isn't large enough,
                            % don't propagate: OS.WF_ is unchanged
                            OS.ReIm2WF;
                        end
                        
                    elseif OS.ELEMENTS_{ii}.isFocal_ == 1
                        fprintf('Computing Focused Field using FFT of WF00%d\n',ii-1);
                        OS.computePSF_normalize(OS.WF_);

                        
                        if OS.savefile == 1
                            OptSys.savePSFfits(dirname,OS.WFamp, OS.WFphase, ii);
                        end
                        
                    elseif OS.ELEMENTS_{ii}.isFocal_ == 2
                        fprintf('Not Supporting Zoom-FFTs yet\n');
                        
                        
                    elseif OS.ELEMENTS_{ii}.isFocal_ == 3
                        fprintf('Not Supporting Convolution with FPM yet\n');
                        
                        
                    end
                    
                end % of check if last element
                
                
            end % of loop over elements
            
            
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
        
        %% GPU Methods -- Nothing available yet...getting all functionality working before accelerating
        % Evaluations folder holds images of performance profiles. This can
        % be used to help. The likely candidates for acceleration are the
        % for loops handling multiple wavelengths, and 
        
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
                warning('GPU:CPUselected', 'GPU is not being used. Run initGPU()');
                return;
            end
            % Get the number of Elements
            numElems = OS.numElements_;
            
            % Make sure the field is of type float
            OS.setdatatype('single');
            % Send field to gpu
            OS.setField(gpuArray(OS.WF_));
            
            for ii = 1:numElems
                % Make sure the zsag_ is of type float
                OS.ELEMENTS_{ii}.setdatatype('single');
                % Send zsag_ to the active gpu
                OS.ELEMENTS_{ii}.set_zsag(gpuArray(OS.ELEMENTS_{ii}.zsag_));
            end
            
            warning('GPU:PROPNS','Propagation is currently pixel by pixel, and not a matrix multiply. This will be incredibly slow on GPU. Consider using CPU');
            
        end % of GPUify
        
        
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
            OS.setField(gather(OS.WF_));
            
            for ii = 1:numElems
                % Send zsag_ to the active gpu
                OS.ELEMENTS_{ii}.set_zsag(gather(OS.ELEMENTS_{ii}.zsag_));
            end
        end % useCPU
        
    end % of methods
    
    %% Static Methods
    methods(Static = true)
        
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
        
        
        function WFout = static_computePSF(WFin, pscale_)
            % WFout = static_computePSF(OS,WFin)
               
            sz = size(WFin);
            if length(sz) == 2
                sz(3) = 1;
            end
            
            WFout = WFin;
            
            WFfocus = fftshift(fft2(fftshift(WFin))) .* (sz(1).* sz(2) .* pscale_ .* pscale_);
            WFreal = real(WFfocus);
            WFimag = imag(WFfocus);
            [WFamp,WFphase] = WFReIm2AmpPhase2(WFreal,WFimag);
            for ii = 1:size(WFamp,3)
                    WFout(:,:,ii) = WFamp(:,:,ii) .* exp(-1i * WFphase(:,:,ii));
            end

            
        end % static_computePSF
        
        function WFout = static_move2PP(WFin, pscale_)
            % WFout = static_move2PP(OS,WFin)
               
            sz = size(WFin);
            if length(sz) == 2
                sz(3) = 1;
            end
            
            WFout = WFin;
            
            WFfocus = ifftshift(ifft2(ifftshift(WFin))) .* (sz(1).* sz(2) .* (pscale_) .* (pscale_))^-1;
            WFreal = real(WFfocus);
            WFimag = imag(WFfocus);
            [WFamp,WFphase] = WFReIm2AmpPhase2(WFreal,WFimag);
            for ii = 1:size(WFamp,3)
                    WFout(:,:,ii) = WFamp(:,:,ii) .* exp(-1i * WFphase(:,:,ii));
            end
            
            % Fix shift induced by MATLAB
            WFout = flipud(WFout);

            
        end % static_move2PP
        
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
