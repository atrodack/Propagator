classdef OptSHWFS < OptWFS
    %OPTSHWFS Class for modeling a Shack-Hartmann WFS
    %   Child of OptWFS
    %
    %   Adapted from JLCodona's AOShackHartmann class in AOSim2:
    %   https://github.com/jlcodona/AOSim2/blob/jlc1.4/AOSim2/%40AOShackHartmann/AOShackHartmann.m
    
    properties(GetAccess='public',SetAccess='public')
        
    end % of public properties
    
    properties(GetAccess='public',SetAccess='private')
        
        % Pupil Function
        A_;
        x_;
        y_;
        
        % Device Properties
        NsubAps_;
        NunMaksedsubAps_;
        subAp_spacing_;
        subApDiam_;
        Offset_ = [0,0];
        
        % Coordinates of subapertures
        XsubApCoords_ = [];
        YsubApCoords_ = [];
        
        MaskedAps_ = [];
        
        % Field Pieces
        subApFields_;
        
        % Measured Centroids
        Xslopes_ = [];
        Yslopes_ = [];
        
        % Initial/previous bias
        Xslopes0_ = [];
        Yslopes0_ = [];
        
        % Sensing Shape
        useCircMasks_ = 1;
        
    end % of protected properties
    
    methods
        %% Constructor
        function SHWFS = OptSHWFS(PROPERTIES)
            % SHWFS = OptSHWFS(PROPERTIES)
            % PROPERTIES is a 6x1 cell array containing element properties
            % If something is unknown or to be set later, leave that
            % position empty:
            %
%           PROPERTIES{1,1} = A                     [Pupil Function]
%           PROPERTIES{2,1} = subAp_spacing         [float in m]   
%           PROPERTIES{3,1} = Offset                [float in m]
%           PROPERTIES{4,1} = x                     [float in m]
%           PROPERTIES{5,1} = y                     [float in m]
%           PROPERTIES{6,1} = useCircMasks          [Boolean]
            
            if size(PROPERTIES) == [6,1]
                if iscell(PROPERTIES) == 1
                    A = PROPERTIES;
                else
                    error('PROPERTIES must be 6x1 cell array!');
                end
            else
                error('PROPERTIES must be a 6x1 cell array!');
            end
            
            N = size(PROPERTIES{1,1},1);
            SHWFS.setPupilFunction(PROPERTIES{1,1});
            SHWFS.setsubAp_spacing(PROPERTIES{2,1});
            SHWFS.setOffset(PROPERTIES{3,1});
            SHWFS.setCoords(PROPERTIES{4,1},PROPERTIES{5,1});
            SHWFS.setCircMaskFlag(PROPERTIES{6,1});
            SHWFS.setGridSize(N);
            SHWFS.interpolate_method = 'nearest';
            
        end % of constructor
        
        
        %% Property Setting Methods
        
        function SHWFS = setPupilFunction(SHWFS,A)
            SHWFS.A_ = A;
        end % of setPupilFunction
        
        function SHWFS = setCoords(SHWFS,x,y)
            if nargin < 2
                error('Need to include Coordinate Vector');
            elseif nargin < 3
                y = x;
            end
            SHWFS.x_ = x;
            SHWFS.y_ = y;
        end % of setCoords
        
        function SHWFS = setsubAp_spacing(SHWFS,spacing)
            % SHWFS = setsubAp_spacing(SHWFS,spacing)
            
            SHWFS.subAp_spacing_ = spacing;
        end % of setsubAp_spacing
        
        function SHWFS = setsubApDiam(SHWFS,diameter)
            % SHWFS = setsubApDiam(SHWFS,diameter)
            SHWFS.subApDiam_ = diameter;
        end % of setsubApDiam
        
        function SHWFS = setOffset(SHWFS,Coords)
            % SHWFS = setOffset(SHWFS, Coords)
            SHWFS.Offset_ = Coords;
        end % of setOffset
        
        function SHWFS = setCircMaskFlag(SHWFS,flag)
            % SHWFS = setFakeFlag(SHWFS,flag)
            SHWFS.useCircMasks_ = flag;
        end % of setCircMaskFlag
        
        %% SubAperture Methods
        
        function SHWFS = defineSubAps(SHWFS,D,thresh)
            % SHWFS = defineSUbAps(SHWFS,D,thresh)
            % Function for finding the coordinates of the centers of the
            % subapertures
            
            dx = SHWFS.x_(2) - SHWFS.x_(1);
            dy = SHWFS.y_(2) - SHWFS.y_(1);
            extent = [D/2,D/2];
            x = -extent(1):dx:extent(1);
            y = -extent(2):dy:extent(2);
            
            xmin = x(1)-dx/2;
            xmax = x(end)+dx/2;
            ymin = y(1)-dy/2;
            ymax = y(end)+dy/2;
            
            Nx = ceil(extent(1)/SHWFS.subAp_spacing_);
            Ny = ceil(extent(2)/SHWFS.subAp_spacing_);
            
            SHWFS.subApDiam_ = SHWFS.subAp_spacing_;
            
            xWFS = (-Nx:Nx)*SHWFS.subAp_spacing_ + SHWFS.Offset_(1);
            xWFS(xWFS<xmin | xWFS>xmax) = [];
            yWFS = (-Ny:Ny)*SHWFS.subAp_spacing_ + SHWFS.Offset_(2);
            yWFS(yWFS<ymin | yWFS>ymax) = [];
            
            [XWFS,YWFS] = meshgrid(xWFS,yWFS);
            SHWFS.XsubApCoords_ = XWFS(:);
            SHWFS.YsubApCoords_ = YWFS(:);
            
            SHWFS.NsubAps_ = length(SHWFS.XsubApCoords_);
            
            
            Acopy = SHWFS.A_;
            Acopy = padarray(Acopy,[ceil(SHWFS.subApDiam_/dx), ceil(SHWFS.subApDiam_/dy)]);
            Acopy = conv2(Acopy,fspecial('disk',SHWFS.subApDiam_/dx/2));
            Acopy_ext = [size(Acopy,1)*dx , size(Acopy,2)*dy];

            gx = linspace(-ceil(Acopy_ext(1)/2),ceil(Acopy_ext(1)/2),size(Acopy,1));
            gy = linspace(-ceil(Acopy_ext(2)/2),ceil(Acopy_ext(2)/2),size(Acopy,2));
            [GX,GY] = meshgrid(gx,gy);
            aper = interp2(GX,GY,Acopy,XWFS,YWFS,SHWFS.interpolate_method);
            
            SHWFS.MaskedAps_ = (aper < thresh);
            
        end % of defineSubAps
        
        function [XsubApCoords, YsubApCoords] = getMaskedsubApCoords(SHWFS)
            XsubApCoords = SHWFS.XsubApCoords_(~SHWFS.MaskedAps_);
            YsubApCoords = SHWFS.YsubApCoords_(~SHWFS.MaskedAps_);
            SHWFS.NunMaksedsubAps_ = length(XsubApCoords);
        end % of getMaskedsubApCoords
        
        %% Sensing Methods
        
        function SHWFS = setBias(SHWFS)
            
            SHWFS.Xslopes0_ = SHWFS.Xslopes_;
            SHWFS.Yslopes0_ = SHWFS.Yslopes_;
        end % of setBias
        
        
        function [subApFields,dx] = parseField(SHWFS,field)
            
            dx = SHWFS.x_(2) - SHWFS.x_(1);
            dy = SHWFS.y_(2) - SHWFS.y_(1);
            d = SHWFS.subApDiam_;
            Nf = ceil(d/dx/2)*2+1;
            x = linspace(-d/2,d/2,Nf);
            dx = x(2) - x(1);
            [Xf,Yf] = meshgrid(x);
            [XsubApCoords, YsubApCoords] = SHWFS.getMaskedsubApCoords;
            
            if ~SHWFS.useCircMasks_
                Rf = sqrt(Xf.^2 + Yf.^2);
                mask = (Rf<max(Xf(:)));
            else
                mask = (ones(size(Xf)));
            end
            
            subApFields = init_variable(size(Xf,1),size(Xf,2),SHWFS.NunMaksedsubAps_,SHWFS.default_data_type,0);
            for ii = 1:SHWFS.NunMaksedsubAps_
                subApFields(:,:,ii) = mask.*interp2(SHWFS.x_, SHWFS.y_, field,Xf+XsubApCoords(ii), Yf+YsubApCoords(ii),SHWFS.interpolate_method);
            end            
        end % of parseField
        
        
        function SHWFS = IdealSense(SHWFS,field,lambda)
            
            [subApFields,dx] = SHWFS.parseField(field);
            subApFields(isnan(subApFields)) = 0;
            
            dwf1 = angle(squeeze(mean(mean(subApFields(2:end,:,:).*conj(subApFields(1:end-1,:,:))))))*lambda/2/pi;
            dwf2 = angle(squeeze(mean(mean(subApFields(:,2:end,:).*conj(subApFields(:,1:end-1,:))))))*lambda/2/pi;
            
            SHWFS.Xslopes_ = 206265 * dwf2 / dx;
            SHWFS.Yslopes_ = 206265 * dwf1 / dx;
            
            if(length(SHWFS.Xslopes_) ~= length(SHWFS.Xslopes0_))
                SHWFS.Xslopes0_ = zeros(size(SHWFS.Xslopes_));
                SHWFS.Yslopes0_ = zeros(size(SHWFS.Yslopes_));
            end
        end % of IdealSense
            
        
        function [Sx, Sy] = getSlopes(SHWFS)
            Sx = SHWFS.Xslopes_ - SHWFS.Xslopes0_;
            Sy = SHWFS.Yslopes_ - SHWFS.Yslopes0_;
        end % of getSlopes
        
        
        %% Plotting Methods
        function plotsubApCenters(SHWFS,fignum)
            if nargin < 2
                fignum = 1;
            end
            figure(fignum)
            hold on
            plot(SHWFS.XsubApCoords_,SHWFS.YsubApCoords_,'g.');
            hold off
        end % of plotsubApCenters
        
        function plotSubAps(SHWFS,linespec)
            if nargin < 2
                linespec = 'b-';
            end
            N = numel(SHWFS.XsubApCoords_);
            for n = 1:N
                hold on;
                drawCircles(SHWFS.subApDiam_/2,[SHWFS.XsubApCoords_(n),SHWFS.YsubApCoords_(n)],1,linespec);
                if(SHWFS.MaskedAps_(n))
                    hold on;
                    plot(SHWFS.XsubApCoords_(n),SHWFS.YsubApCoords_(n),'rx');
                end
            end
            hold off;
            
        end % of plotSubAps
        
        function quiver(SHWFS,linespec)
            if nargin < 2
                linespec = 'g-';
            end
            
            [X,Y] = SHWFS.getMaskedsubApCoords;
            [Sx,Sy] = SHWFS.getSlopes;
            hold on
            quiver(X(:),Y(:),Sx(:),Sy(:),linespec);
        end % of quiver
        
    end % of methods
    
end

