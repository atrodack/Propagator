classdef OptElement < matlab.mixin.Copyable
    %Element Class for storing parameters of an optical element
    %   Detailed explanation goes here
    
    
    %% Properties
    properties(Constant=true, GetAccess='protected')
        
        
    end % of constant properties
    
    properties(GetAccess='public',SetAccess='public')
        name; % give a name to the structure/system
        
    end % of public properties
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        % Element Properties
        type_; % element type
        isMirror_; %flag for determining lens or mirror
        focal_length_; % focal length in meters
        z_position_; % position from previous element along optical axis
        diameter_; % diameter of element in meters
        
        % Map of Optic Surface
        zsag_; % (in meters)
        
        % Propagation Type
        % 0 = Fourier Optics, 1 = Fraunhofer, 2 = Fresnel
        propagation_type;
        
        
        % Misc
        verbose = 1; % print extra info
        interpolate_method;
        
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
            % PROPERTIES{1,1} = name
            % PROPERTIES{2,1} = type (0-8)
            % PROPERTIES{3,1} = isMirror (1 or 0)
            % PROPERTIES{4,1} = focal length (m)
            % PROPERTIES{5,1} = z_position (m)
            % PROPERTIES{6,1} = diameter (m)
            % PROPERTIES{7,1} = zsag (m) [file path or matrix]
            % PROPERTIES{8,1} = propagation code
            
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
            elem.set_isMirror(A{3,1});
            elem.set_focal_length(A{4,1});
            elem.set_z_position(A{5,1});
            elem.set_diameter(A{6,1});
            elem.set_zsag(A{7,1});
            elem.set_propagation_type(A{8,1});
            
            
            
            
            
            
        end % of contructor
        
    %% Set Properties
    
    function elem = set_type(elem,type)
        % elem = set_type(elem,type)
        % sets the type of optical element
        % 0: System Pupil
        % 1: Normally Behaving Lens/Mirror
        % 2: Detector
        % 3: PIAA
        % 4: FPM
        % 5: Lyot Stop
        % 6: Apdoizer
        % 7: DM
        % 8: WFS -- likely to break this into specific WFS
        
        elem.type_ = type;
    end % of set_type
    
    function elem = set_isMirror(elem,flag)
        % elem = set_isMirror(elem,flag)
        % sets the flag to determine if element is a mirror or lens
        
        elem.isMirror_ = flag;
    end % of set_isMirror
    
    function elem = set_focal_length(elem,f)
        % elem = set_focal_length(elem,f)
        % sets the focal length of the element to f (in meters)
        
        elem.focal_length_ = f;
    end % of set_focal_length
    
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
            
        elseif ismatrix(A)
            if elem.verbose == 1
                fprintf('Loading matrix into zsag\n');
            end
        end
        elem.zsag_ = A;
    end % of set_zsag
    
    function elem = set_propagation_type(elem,val)
        % elem = set_propagation_type
        % sets the propagation type:
        % val = 0 --> Fourier Optics
        % val = 1 --> Fraunhofer
        % val = 2 --> Fresnel
        
        elem.propagation_type = val;
    end % of set_propagation_type
    
    function elem = set_name(elem,name)
        % elem = set_name(elem,name)
        % sets the name of the element
        
        elem.name = name;
    end % of set_name
    
    
    %% Read out properties
    
    function type = getElementType(elem)
        % type = getElementType(elem)
        % returns the type of element 
        
        type = elem.type_;
        if type == 0
            descr = 'System Pupil';
        elseif type == 1
            descr = 'Normally Behaving Lens/Mirror';
        elseif type == 2
            descr = 'Detector';
        elseif type == 3
            descr = 'PIAA';
        elseif type == 4
            descr = 'FPM';
        elseif type == 5
            descr = 'Lyot Stop';
        elseif type == 6
            descr = 'Apdoizer';
        elseif type == 7
            descr = 'DM';
        elseif type == 8
            descr = 'WFS';
        end
        fprintf('Element is a %s \n',descr);
    end % of getElementType
    
    function morl = getisMirror(elem)
        %morl = getisMirror(elem)
        % returns if the element is a mirror or lens
        
        morl = elem.isMirror_;
        if morl == 0
            descr = 'Lens';
        elseif morl == 1
            descr = 'Mirror';
        end
       fprintf('Element is a %s\n',descr);
    end % of getisMirror(elem);
    
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
    
    function propagation_method = getPropagationMethod(elem)
        % propagation_method = getPropagationMethod(elem)
        % returns the method to be used for propagation
        
        val = elem.propagation_type;
        if val == 0
            propagation_method = 'Fourier Optics';
        elseif val == 1
            propagation_method = 'Fraunhofer';
        elseif val == 2
            propagation_method = 'Fresnel';
        end
    end % of getPropagationMethod
    
    
    %% Utilities
    
    function fnum = getFNumber(elem)
        % fnum = getFNumber(elem)
        % returns the f/# of the element
        
        fnum = elem.getFocalLength / elem.getDiameter;
    end % getFNumber
    
    
    end % of methods
    
    %% Static Methods
    methods(Static = true)
        
        
    end % of static methods
end
