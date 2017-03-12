classdef OptElement < matlab.mixin.Copyable
    %Element Class for storing parameters of an optical element
    %   Detailed explanation goes here
    
    
    %% Properties
    properties(Constant=true, GetAccess='protected')
        
        
    end % of constant properties
    
    properties(GetAccess='public',SetAccess='public')
        name; % give a name to the structure/system
        verbose = 1; % print extra info
        interpolate_method = []; % selects a method of interpolation
    end % of public properties
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        
        % Data Type
        default_data_type = 'single'; % data_type of zsag_
        
        % Element Properties
        material_;
        z_position_; % position downstream from previous element (pupil is z = 0)
        diameter_; % diameter of element in meters
        isFocal_ = 0;
        gridsize_ = [1024, 1024]; % default gridsize
        
        % Map of Optic Surface
        zsag_; % (in meters)
        
        phasefac_; % factor used to convert zsag_ into phase
        
                
    end % of protected properties
    
    %% Methods
    methods

        
    %% Set Properties
     
    function elem = set_material(elem, code)
        % elem = set_material(elem,code)
        % This creates a material at a default wavelength set by the
        % OptMaterial class. 
        
        elem.material_ = OptMaterial(code);
    end
        
    function elem = set_phasefactor(elem,lambda,n0)
        % elem = set_phasefactor(elem)
        % computes and sets the phasefac_ property using the material
        elem.material_.setWavelength(lambda);
        elem.phasefac_ = elem.material_.ComputePhaseFactor(n0);
    end % of set_phasefactor
        
        
    
    
    function elem = set_z_position(elem,z)
        % elem = set_z_position(elem,z)
        % sets the z position of the element to z (in meters)
        
        elem.z_position_ = z;
        
    end % of set_z_position
    
    function elem = set_diameter(elem,D)
        % elem = set_diameter(elem,D)
        % sets the diameter of the element to D (in meters)
        
        elem.diameter_ = D;
    end % of set_diameter
    
    function elem = set_zsag(elem,A)
        % elem = set_zsag(elem,A)
        % sets the map of the optic surface
        % if A is a matrix, load it directly
        % if A is a file path to a fits file, load the fits file
        % sags should be in meters
        
        if ischar(A) == 1
            if elem.verbose == 1
                fprintf('Reading FITS file %s into zsag\n',A);
            end
            A = fitsread(A);
            elem.gridsize_ = size(A);
            
        elseif ismatrix(A)
            if numel(A) == 1
                A = ones(elem.gridsize_); % use size stored in gridsize_
            else
                elem.gridsize_ = size(A);
            end
            if elem.verbose == 1
                fprintf('Loading matrix into zsag\n');
            end
        end
        elem.zsag_ = A;
    end % of set_zsag
        
    function elem = set_name(elem,name)
        % elem = set_name(elem,name)
        % sets the name of the element
        
        elem.name = name;
    end % of set_name
    
    function elem = setdatatype(elem,default_data_type)
            % elem = datatype(default_data_type)
            % Sets the datatype to use. Do not do anything if already of
            % the correct data type.
            % Currently supported:
            % single
            % double
            % uint8
            
            if nargin < 2
                default_data_type = elem.default_data_type;
            else
                elem.default_data_type = default_data_type;
            end
            
            switch default_data_type
                case 'single'
                    if ~isa(elem.zsag_,'single')
                        elem.zsag_ = single(elem.zsag_);
                        if elem.verbose == 1
                            fprintf('Data Type set to single\n');
                        end
                    end
                    
                case 'double'
                    if ~isa(elem.zsag_,'double')
                        elem.zsag_ = double(elem.zsag_);
                        if elem.verbose == 1
                            fprintf('Data Type set to double\n');
                        end
                    end
                    
                case 'uint8'
                    if ~isa(elem.zsag_,'uint8')
                        elem.zsag_ = uint8(elem.zsag_);
                        if elem.verbose == 1
                            fprintf('Data Type set to uint8\n');
                        end
                    end
                    
                otherwise
                    error('I do not understand that data type (yet)!');
            end
        end % of setdatatype
    %% Read out properties
    
    function zpos = getZPosition(elem)
        %zpos = getZPosition(elem)
        % returns the value stored in z_position
        
        zpos = elem.z_position_;
    end % of getzposition
    
    function D = getDiameter(elem)
        % D = getDiameter(elem)
        % returns the value stored in diameter_
        
        D = elem.diameter_;
    end % of getDiameter
    
    function zsag = getZsag(elem)
        % zsag = getZsag(elem)
        % returns the array stored in zsag_
        
        zsag = elem.zsag_;
    end % of getZsag

    
    function show(elem)
        % show(elem)
        % Plots matrix stored in zsag_
        
        if isreal(elem.zsag_) == 1
            figure;
            imagesc(elem.zsag_)
            plotUtils(sprintf('Sag of Element %s',elem.name));
        else
            re = real(elem.zsag_);
            im = imag(elem.zsag_);
            [sagamp,sagphase] = WFReIm2AmpPhase(re,im);
            figure;
            subplot(1,2,1)
            imagesc(sagamp);
            plotUtils(sprintf('Sag Amplitude of Element %s',elem.name));
            subplot(1,2,2)
            imagesc(sagphase);
            plotUtils(sprintf('Sag Phase of Element %s',elem.name));
            
        end
        
    end % of show
    
    %% Utility Methods
    
    function elem = Cubify(elem,numLambdas)
        % elem = Cubify(elem,numLambdas)
        
        if strcmp( elem.default_data_type, 'single')
            tmp = single(zeros(elem.gridsize_(1),elem.gridsize_(2),numLambdas));
        elseif strcmp( elem.default_data_type, 'double')
            tmp = double(zeros(elem.gridsize_(1),elem.gridsize_(2),numLambdas));
        elseif strcmp( elem.default_data_type, 'uint8')
            tmp = uint8(zeros(elem.gridsize_(1),elem.gridsize_(2),numLambdas));
        end
        
        for ii = 1:numLambdas
            tmp(:,:,ii) = elem.zsag_;
        end
        elem.set_zsag(tmp);
        clear tmp;
    end % of Cubify
    
    
    
    end % of methods
    
    %% Static Methods
    methods(Static = true)
        
        function addnewline(numlines)
            % addnewline(numlines)
            for ii = 1:numlines
                fprintf('\n');
            end
        end % of addnewline
        
    end % of static methods
end
