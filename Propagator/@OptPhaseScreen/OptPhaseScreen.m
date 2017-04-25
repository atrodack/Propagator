classdef OptPhaseScreen < OptWF
    %OPTPHASESCREEN Class for Turbulent PhaseScreens
    %   Detailed explanation goes here
    %
    %
    % REFERENCES:
    % [1] Hardy, John W. Adaptive Optics for Astronomical Telescopes.
    % Chapter 3
    
    properties(GetAccess = 'public', SetAccess = 'public')

    end % of public properties
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        % 
        fixLF_;
        
        altitude_;          % Altitude of the phase screen [m]
        thickness_;         % Thickness of phase screen    [m]
        L0_ = 30;           % Outer scale                  [m]
        l0_ = 1e-3;         % Inner scale (1mm - Roddier)  [m]
        alpha_ = 11/3;      % PSD exponent
        
        Cn2_;               % Refractive index structure constant
        r0_;                % Fried length                 [m]
                
        Model_code_ = 3;    % Model type for PSD (see below)
        PSD_;               % Storage of the PSD
                
        KR_;
                
    end % of private properties
    
    methods
        %% Constructor
        function elem = OptPhaseScreen(PROPERTIES)
            % elem = OptPhaseScreen(PROPERTIES)
            % PROPERTIES is a 10x1 cell array containing element properties
            % If something is unknown or to be set later, leave that
            % position in A empty. Provide either Cn2 or r0 (4 or 5 should
            % be blank). The other will be computed
            %
%           PROPERTIES{1,1} = name                  [string]
%           PROPERTIES{2,1} = altitude              [float in m]
%           PROPERTIES{3,1} = thickness             [float in m]
%           PROPERTIES{4,1} = lambdaref             [float in m]
%           PROPERTIES{5,1} = Cn2                   [float in ]
%           PROPERTIES{6,1} = r0                    [float in m]
%           PROPERTIES{7,1} = Model Code            [int]
%           PROPERTIES{8,1} = KR                    [float in m]
%           PROPERTIES{9,1} = dx                    [float in m]
%           PROPERTIES{10,1} = N                    [int in pixels]

            elem = elem@OptWF(PROPERTIES);
            
            if size(PROPERTIES) == [10,1]
                if iscell(PROPERTIES) == 1
                    A = PROPERTIES;
                else
                    error('PROPERTIES must be 10x1 cell array!');
                end
            else
                error('PROPERTIES must be a 10x1 cell array!');
            end
            
            elem.set_name(A{1,1});
            elem.set_altitude(A{2,1});
            elem.set_thickness(A{3,1});
            elem.set_Central_Wavelength(A{4,1});
            if isempty(A{5,1})
                elem.set_r0(A{6,1});
            else
                elem.set_Cn2(A{5,1});
            end
            
            elem.set_Modelcode(A{7,1});
            elem.set_KR(A{8,1});
            elem.set_dx(A{9,1});
            elem.set_gridsize(A{10,1});
            elem.set_datatype;
            
        end % of constructor
        
        %% Utilities for setting Properties
        
        
        function elem = set_altitude(elem,h)
            % elem = set_altitude(elem,h)
            % Set the altitude of the screen
            % if h is a scalar, assume single layer screen
            % if h is a vector, assume Cn^2 varies through the thickness 
            %                   with samples at h(n)
            
            elem.altitude_ = h;
        end % of set_altitude
        
        function elem = set_thickness(elem,th)
            % elem = set_thickness(elem,th)
            elem.thickness_ = th;
        end % of set_thickness
        
        function elem = set_Cn2(elem, Cn2)
            % elem = set_Cn2(elem, Cn2)
            % Set Cn^2 profile:
            % if Cn2 is scalar, assume a constant Cn^2 for the screen
            % if Cn2 is a vector, Cn^2 values at altitudes in
            % elem.altitudes_
            
            n = length(elem.altitude_);
            if length(Cn2) ~= n
                error('Must provide a Cn^2 value for each provided altitude');
            else
                elem.Cn2_ = Cn2;
                elem.compute_r0(Cn2);
            end
            
        end % of set_Cn2
        
        function elem = set_r0(elem,r0)
            % elem = set_r0(elem,r0)
            
            elem.r0_ = r0;
            elem.compute_Cn2(r0);
            
        end % of set_r0
        
        function elem = set_Modelcode(elem,code)
            % elem = set_Modelcode(elem,code)
            % Tell the Screen what model for the PSD to use
            % See [Ref 1] for help
            % Valid Code Inputs:
            % 0: Kolmogorov - Standard 2/3 law
            % 1: Kolmogorov w/ Tatarski Adjustment for inner scale
            % 2: von Karman Spectrum - finite inner and outer scales
            
            if code >= 0
                if code < 3
                    elem.Model_code_ = code;
                else
                    error('Code must be 0,1, or 2');
                end
            else
                error('Code must be 0,1, or 2');
            end
        end % of set_Modelcode
        
        function elem = set_innerScale(elem,l0)
            % elem = set_innerscale(elem,l0)
            
            elem.l0_ = l0;
        end % of set_innerScale
        
        function elem = set_outerScale(elem,L0)
            % elem = set_outerScale(elem,L0)
            
            elem.L0_ = L0;
        end % of set_outerScale
        
        function elem = set_alpha(elem,alpha)
            % elem = set_alpha(elem,alpha)
            
            elem.alpha_ = alpha;
        end % of set_alpha
        
        function elem = set_PSD(elem,PSD)
            % elem = set_PSD(elem,PSD)
            
            elem.PSD_ = PSD;
        end % of set_PSD
        
%         function elem = set_screen(elem,screen)
%             % elem = set_screen(elem,screen)
%             
%             elem.field_ = screen;
%         end % of set_screen
        
        function elem = set_KR(elem,KR)
            % elem = set_KR(elem,KR)
            
            elem.KR_ = KR;
        end % of set_KR
        
        function elem = set_fixLF(elem,flag)
            % elem = set_fixLF(elem,flag)
            elem.fixLF_ = flag;
        end % of set_fixLF
        
        
        %% Compute Utilities
        
        function elem = compute_Cn2(elem,r0,lambda,thickness)
            % elem = compute_Cn2(elem,r0,lambda,thickness)
            %
            % Assumes zenith angle of 0 degrees, Cn^2 constant over
            % thickness
            %
            % r0 = Fried Length
            % lambda = wavelength
            % thickness = thickness of layer
            %            if a vector, this vector represents the thickness
            %            of each layer, and should be the difference
            %            between adjacent altitudes stored in the
            %            altitude_ property
            %
            % See Ref[1] Eq. (3.51)
            
            if nargin < 2
                r0 = elem.r0_;
            end
            if nargin < 3
                lambda = elem.lambda0_;
            end
            if nargin < 4
                thickness = elem.thickness_;
            end
                        
            k = (2*pi) ./ lambda;
            elem.Cn2_ = (r0.^(-5/3)) / 0.423 ./ (k.^2) ./ thickness;
        end % of compute_Cn2
        
        function elem = compute_r0(elem,Cn2,lambda,thickness)
            % elem = compute_r0(elem,Cn2,lambda,thickness)
            %
            % Assumes zenith angle of 0 degrees, Cn^2 constant over
            % thickness
            %
            % Cn2 = refractive index structure constant
            % lambda = wavelength
            % thicknes = thickness of layer
            %            if a vector, this vector represents the thickness
            %            of each layer, and should be the difference
            %            between adjacent altitudes stored in the
            %            altitude_ property
            %
            % See Ref[1] Eq. (3.51)
            
            if nargin < 2
                Cn2 = elem.Cn2_;
            end
            if nargin < 3
                lambda = elem.lambda0_;
            end
            if nargin < 4
                thickness = elem.thickness_;
            end
            
            k = (2*pi) ./ lambda;
            elem.r0_ = (0.423 .* (k.^2) .* Cn2 .* thickness).^(-3/5);
        end % of compute_r0
        
        function elem = compute_Cn2_zenith(elem,r0,lambda,thickness,zenith)
            % elem = compute_Cn2_zenith(elem,r0,lambda,thickness)
            %
            % Assumes zenith angle of 0 degrees, Cn^2 constant over
            % thickness
            %
            % r0 = Fried Length
            % lambda = wavelength
            % thickness = thickness of layer
            %            if a vector, this vector represents the thickness
            %            of each layer, and should be the difference
            %            between adjacent altitudes stored in the
            %            altitude_ property
            % zenith = angle from zenith in degrees
            %
            % See Ref[1] Eq. (3.51)
            
            if nargin < 2
                r0 = elem.r0_;
            end
            if nargin < 3
                lambda = elem.lambda0_;
            end
            if nargin < 4
                thickness = elem.thickness_;
            end
            if nargin < 5
                zenith = 0;
            end
            
            k = (2*pi) ./ lambda;
            elem.Cn2_ = (r0.^(-5/3)) / 0.423 ./ (k.^2) / secd(zenith) ./ thickness;
        end % of compute_Cn2_zenith
        
        function elem = compute_r0_zenith(elem,Cn2,lambda,thickness,zenith)
            % elem = compute_r0_zenith(elem,Cn2,lambda,thickness)
            %
            % Assumes zenith angle of 0 degrees, Cn^2 constant over
            % thickness
            %
            % Cn2 = refractive index structure constant
            % lambda = wavelength
            % thicknes = thickness of layer
            %            if a vector, this vector represents the thickness
            %            of each layer, and should be the difference
            %            between adjacent altitudes stored in the
            %            altitude_ property
            % zenith = angle from zenith in degrees
            %
            % See Ref[1] Eq. (3.51)
            
            if nargin < 2
                Cn2 = elem.Cn2_;
            end
            if nargin < 3
                lambda = elem.lambda0_;
            end
            if nargin < 4
                thickness = elem.thickness_;
            end
            if nargin < 5
                zenith = 0;
            end
            
            k = (2*pi) ./ lambda;
            elem.r0_ = (0.423 .* (k.^2) * secd(zenith) .* Cn2 .* thickness).^(-3/5);
        end % of compute_r0_zenith
        
        function r0 = compute_sensorLambdar0(elem, sensor_lambda, zenith)
            % r0 = compute_sensorLambdar0(elem,sensor_lambda,zenith)
            %
            % sensor_lambda must be in meters
            
            if nargin < 3
                zenith = 0;
            end
            
            % convert sensor_lambda to microns
            sensor_lambda = sensor_lambda * 1e6;
            ref_lambda = elem.lambda0_ * 1e6;
            r0 = (elem.r0_) .* (sensor_lambda / ref_lambda) .^(6/5) * (cos(zenith)).^(3/5);
            
            
        end % of compute_sensorLambdar0
        
        %% Utilities
        
        function Rf = Fresnel_Scale(elem,lambda,h)
            % Rf = Fresnel_Scale(lambda,h)
            
            if nargin < 2
                lambda = elem.lambda0_;
            end
            if nargin < 3
                h = elem.altitude_;
            end
            
            Rf = sqrt(lambda .* h);
        end % of Fresnel_Scale
        
        
        
        %% Overloaded Method
        function elem = set_datatype(elem,default_data_type)
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
                    if ~isa(elem.field_,'single')
                        elem.field_ = single(elem.field_);
                        elem.PSD_ = single(elem.PSD_);
                        elem.KR_ = single(elem.KR_);
                        if elem.verbose == 1
                            fprintf('Data Type set to single\n');
                        end
                    end
                    
                case 'double'
                    if ~isa(elem.field_,'double')
                        elem.OS.field_ = double(elem.field_);
                        elem.PSD_ = double(elem.PSD_);
                        elem.KR_ = double(elem.KR_);
                        if elem.verbose == 1
                            fprintf('Data Type set to double\n');
                        end
                    end
                    
                case 'uint8'
                    if ~isa(elem.field_,'uint8')
                        elem.field_ = uint8(elem.field_);
                        elem.PSD_ = uint8(elem.PSD_);
                        elem.KR_ = uint8(elem.KR_);
                        if elem.verbose == 1
                            fprintf('Data Type set to uint8\n');
                        end
                    end
                    
                otherwise
                    error('I do not understand that data type (yet)!');
            end
        end % of set_datatype
        
        
        %% Make the screen
        
        function elem = makeScreen(elem,N,KR,dx,fixLF,L0,l0,alpha,thickness,Cn2)
        % elem = makeScreen(elem,N,KR,dx,L0,l0,alpha,thickness,altitude,Cn2)
        
        if nargin < 2
            N = elem.gridsize_(1);
        end
        if nargin < 3
            KR = elem.KR_;
        end
        if nargin < 4
            dx = elem.dx_;
        end
        if nargin < 5
            fixLF = elem.fixLF_;
        end
        if nargin < 6
            L0 = elem.L0_;
        end
        if nargin < 7
            l0 = elem.l0_;
        end
        if nargin < 8
            alpha = elem.alpha_;
        end
        if nargin < 9
            thickness = elem.thickness_;
        end
        if nargin < 10
            Cn2 = elem.Cn2_;
        end
        
        K0 = (2*pi) / L0;
        Ki = 5.92 / l0;
        
        PSD_ = zeros(N,N,length(Cn2));
        screen_ = zeros(N,N,length(Cn2));
            
        for ii = 1:length(Cn2)
            switch elem.Model_code_
                case 0 % Kolmogorov
                    PSD = 0.033.*Cn2(ii) .* (KR.^2).^(-alpha/2);
                    
                case 1 % Kolmogorov + Tatarski
                    PSD = 0.033.*Cn2(ii) .* (KR.^2).^(-alpha/2);
                    PSD = PSD .* exp((-KR.^2) ./ (Ki^2));
                    
                    
                case 2 % von Karman
                    PSD = 0.033 .* Cn2(ii) .* (K0^2 + KR.^2).^(-alpha/2);
                    PSD = PSD .* exp((-KR.^2) ./ (Ki^2));
            end
            PSD(N/2,N/2) = 0;
            PSD = PSD * thickness * 1.25;
            
            rv = elem.complexNormalRV(1,N);
            grid = rv.* sqrt(PSD);
            grid_ft = fftshift(fft2(fftshift(grid)));
            tmp = grid_ft .* 24.6 / dx / N;
            
            
            if(fixLF)
                Nlf = 8;
                LFpart = imag(tmp(1:round(N/Nlf*1.1),1:round(N/Nlf*1.1)));
                LFpart = interp2(LFpart,3);
                LFpart = LFpart(1:N,1:N); % trim to the size of the screen.
                scaling = 0.8 * 8^(5/6);
                screen_(:,:,ii) = real(tmp) + LFpart*scaling;
                PSD_(:,:,ii) = PSD;
            else
                screen_(:,:,ii) = real(tmp);
                PSD_(:,:,ii) = PSD;
            end
        end
        
        elem.field_ = screen_;
        elem.PSD_ = PSD_;
        elem.set_datatype;
        
        end % of makeScreen
        
    end % of methods
    
end % of classdef

















