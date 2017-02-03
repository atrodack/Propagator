classdef OptElement < matlab.mixin.Copyable
    %Element Class for storing parameters of an optical element
    %   Detailed explanation goes here
    
    
    %% Properties
    properties(Constant=true, GetAccess='protected')
        
        
    end % of constant properties
    
    properties(GetAccess='public',SetAccess='public')
        name; % give a name to the structure/system
        default_data_type = 'double'; % data_type of zsag_
        
        verbose = 1; % print extra info
        interpolate_method = []; % selects a method of interpolation
    end % of public properties
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        % Element Properties
        type_; % element type ---> see set_type method
        material_;
        focal_length_; % focal length in meters
        z_position_; % position along optical axis (pupil is z = 0)
        diameter_; % diameter of element in meters
        gridsize_ = [1024, 1024]; % default gridsize
        
        % Map of Optic Surface
        zsag_; % (in meters)
        amponly = 1; % flag for determining if Mask is phase or amplitude
                     % Defaults to amplitude, but is set by input
                     % focal_length_ for Mask types anyway. If not 0, then
                     % it is considered amplitude only (NOTE: This means
                     % mistakenly setting a focal length will make it
                     % amplitude mask. This follows the assumption that the
                     % user should take more care to ensure correctness if
                     % the mask is meant to be a phase mask.
        
        phasefac_; % factor used to convert zsag_ into phase
        
        % Propagation Type
        % 0 = Fresnel Propagation needed, 1 = Fourier Transform
        isFocal;
                
    end % of protected properties
    
    %% Methods
    methods
        %% Constructor
        function elem = OptElement(PROPERTIES)
            % elem = Element(A)
            % PROPERTIES is a 8x1 cell array containing element properties
            % If something is unknown or to be set later, leave that
            % position in A empty:
            %
%           PROPERTIES{1,1} = name
%           PROPERTIES{2,1} = type (0-7)
%           PROPERTIES{3,1} = material (0-10)
%           PROPERTIES{4,1} = focal length (m) *sets amponly for Mask types
%           PROPERTIES{5,1} = z_position (m)
%           PROPERTIES{6,1} = diameter (m)
%           PROPERTIES{7,1} = zsag (m) [file path or matrix]
%           PROPERTIES{8,1} = isFocal flag
            
            if size(PROPERTIES) == [8,1]
                if iscell(PROPERTIES) == 1
                    A = PROPERTIES;
                else
                    error('PROPERTIES must be 8x1 cell array!');
                end
            else
                error('PROPERTIES must be a 8x1 cell array!');
            end
            
            elem.set_name(A{1,1});
            elem.set_type(A{2,1});
            elem.set_material(A{3,1});
            elem.set_focal_length(A{4,1});
            elem.set_z_position(A{5,1});
            elem.set_diameter(A{6,1});
            elem.set_zsag(A{7,1});
            elem.set_isFocal(A{8,1});
            
            
            
            
            
            
        end % of contructor
        
    %% Set Properties
    
    function elem = set_type(elem,type)
        % elem = set_type(elem,type)
        % sets the type of optical element
        % Masks default to amplitude only.
        %
        % Types Supported
        % ===============
        % 0: System Pupil
        % 1: Lens
        % 2: Mirror
        % 3: Aspheric Lens
        % 4: Aspheric Mirror
        % 5: Pupil Mask
        % 6: Focal Plane Mask
        % 7: Detector
        
        elem.type_ = type;
        switch type % have this if needed
                case 0 % System Pupil
                    
                case 1 % Lens
                    
                case 2 % Mirror
                    
                case 3 % Aspheric Lens
                    
                case 4 % Aspheric Mirror
                    
                case 5 % Pupil Mask
                    
                case 6 % Focal Plane Mask
                    
                case 7 % Detector
                    
                otherwise
                    error('Unknown Element type');
        end
            
    end % of set_type
    
    function elem = set_material(elem,MATERIAL)
        % elem = set_material(elem,MATERIAL)
        % sets material to corresponding OptMaterial object
        % See OptMaterial.m for material more info
        
        if isa(MATERIAL,'OptMaterial')
            elem.material_ = MATERIAL;
            if elem.verbose == 1
                fprintf('Material %s loaded into Element %s\n',MATERIAL.material.name,elem.name);
            end
        elseif isempty(MATERIAL)
            if elem.verbose == 1
                fprintf('No Material for Element %s\n',elem.name);
            end
        else
            error('Material must be of class OptMaterial');
        end
        
    end % of set_material
        
    function elem = set_phasefactor(elem,lambda)
        % elem = set_phasefactor(elem)
        % computes and sets the phasefac_ property using the material
        
        elem.phasefac_ = elem.material_.ComputePhaseFactor(lambda);
    end % of set_phasefactor
        
        
    
    function elem = set_focal_length(elem,f)
        % elem = set_focal_length(elem,f)
        % sets the focal length of the element to f (in meters)
        % If the element is a mask, use this as a flag to set if the mask
        % is amplitude or phase mask (1 for amp, 0 for phase).
        
        switch elem.type_
                case 0 % System Pupil
                    elem.amponly = f;
                    
                case 1 % Lens
                    elem.focal_length_ = f;
                case 2 % Mirror
                    elem.focal_length_ = f;
                case 3 % Aspheric Lens
                    elem.focal_length_ = f;
                case 4 % Aspheric Mirror
                    elem.focal_length_ = f;
                case 5 % Pupil Mask
                    elem.amponly = f;
                case 6 % Focal Plane Mask
                    elem.amponly = f;
                case 7 % Detector
                    
                otherwise
                    error('Unknown Element type');
        end
        
    end % of set_focal_length
    
    function elem = set_z_position(elem,z)
        % elem = set_z_position(elem,z)
        % sets the z position of the element to z (in meters)
        if elem.type_ == 0
            elem.z_position_ = 0;
        else
            elem.z_position_ = z;
        end
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
                fprintf('Reading FITS file %s into zsag\n\n',A);
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
                fprintf('Loading matrix into zsag\n\n');
            end
        end
        elem.zsag_ = A;
    end % of set_zsag
    
    function elem = set_isFocal(elem,val)
        % elem = set_isFocal
        % sets the propagation type:
        % val = 0 --> Fresnel
        % val = 1 --> Fourier Transform
        elem.isFocal = val;
    end % of set_isFocal
    
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
    
    function type = getElementType(elem)
        % type = getElementType(elem)
        % returns the type of element 
        
        type = elem.type_;
        if type == 0
            descr = 'System Pupil';
        elseif type == 1
            descr = 'Lens';
        elseif type == 2
            descr = 'Mirror';
        elseif type == 3
            descr = 'Aspheric Lens';
        elseif type == 4
            descr = 'Aspheric Mirror';
        elseif type == 5
            descr = 'Pupil Plane Mask';
        elseif type == 6
            descr = 'Focal Plane Mask';
        elseif type == 7
            descr = 'Detector';
        end
        if elem.verbose == 1
            fprintf('Element is a %s \n',descr);
        end
    end % of getElementType
    
    
    function fl = getFocalLength(elem)
        % fl = getFocalLength(elem)
        % returns the value stored in the focal_length_
        
        fl = elem.focal_length_;
    end % of getFocalLength
    
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
    
    function descr = describe(elem)
        % descr = describe(elem)
        
        if elem.verbose == 1
            flag = true;
            elem.verbose = 0;
        end
        
        objtype = 'OptElement';
        
        type = elem.type_;
        if type == 0
            elemtype = 'System Pupil';
        elseif type == 1
            elemtype = 'Lens';
        elseif type == 2
            elemtype = 'Mirror';
        elseif type == 3
            elemtype = 'Aspheric Lens';
        elseif type == 4
            elemtype = 'Aspheric Mirror';
        elseif type == 5
            elemtype = 'Pupil Plane Mask';
        elseif type == 6
            elemtype = 'Focal Plane Mask';
        elseif type == 7
            elemtype = 'Detector';
        end
        
        
        zpos = abs(elem.z_position_);
        sz = size(elem.zsag_);
        
        descr = sprintf('%s:\nElement is a %s at %0.3f meters downstream\n',elem.name,elemtype,zpos);
        if flag == true
            elem.verbose = 1;
        end
        
    end % of describe
    
    function show(elem)
        % show(elem)
        % Plots matrix stored in zsag_
        
        if isreal(elem.zsag_) == 1
            figure()
            imagesc(elem.zsag_)
            plotUtils(sprintf('Sag of Element %s',elem.name));
        else
            re = real(elem.zsag_);
            im = imag(elem.zsag_);
            [sagamp,sagphase] = WFReIm2AmpPhase(re,im);
            figure();
            subplot(1,2,1)
            imagesc(sagamp);
            plotUtils(sprintf('Sag Amplitude of Element %s',elem.name));
            subplot(1,2,2)
            imagesc(sagphase);
            plotUtils(sprintf('Sag Phase of Element %s',elem.name));
            
        end
        
    end % of show
    
    
    end % of methods
    
    %% Static Methods
    methods(Static = true)
        
        
    end % of static methods
end
