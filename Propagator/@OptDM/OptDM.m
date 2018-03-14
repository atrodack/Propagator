classdef OptDM < OptMirror
    %OPTDM Class for modelling a Deformable Mirror
    %   
    
    properties(GetAccess='public',SetAccess='private')
        
        % Actuators
        actuators_; %[x,y,displacement,enabled]
        influencefxn_;
        ActMap_;
        nActuators_;
        
        % Coordinates
        X_;
        Y_;
        bconds_;
        
    end % protected properties
    
    methods
        %% Constructor
        function DM = OptDM(PROPERTIES,x,y)
            % DM = OptDM(PROPERTIES)
            % PROPERTIES is a 6x1 cell array containing element properties
            % If something is unknown or to be set later, leave that
            % position empty:
            %
%           PROPERTIES{1,1} = name                  [string]
%           PROPERTIES{2,1} = focalLength           [float in m]       
%           PROPERTIES{3,1} = isFocal_              [code {see below}]
%           PROPERTIES{4,1} = z_position            [float in m]
%           PROPERTIES{5,1} = diameter              [float in m]
%           PROPERTIES{6,1} = zsag                  [file path or matrix]

            DM = DM@OptMirror(PROPERTIES);
            DM.setCoords(x,y);
            DM.defineBC((PROPERTIES{5,1}/2)*1.05,4,'square');
        end % of contructor
        
        %% Init Methods
        
        function DM = setCoords(DM,x,y)
            if size(y) == size(x)
                sz = size(x);
                if(sz(1) == 1 || sz(2) == 1)
                    [X,Y] = meshgrid(x,y);
                elseif(sz(1) == sz(2))
                    X = x;
                    Y = y;
                end
            else
                error('Provide coordinate vector or meshgrid');
            end
            DM.X_ = X;
            DM.Y_ = Y;
        end % of setCoords
        
        function DM = defineBC(DM,radius,npoints,PATTERN,Offset)
            if nargin == 3
                PATTERN = 'square';
            end
            
            if nargin < 5
                Offset = [0 0];
            end
            
            if strcmp(PATTERN,'circle')
                theta = (1:npoints)'/npoints*2*pi;
                DM.bconds_ = [cos(theta) sin(theta)] * radius;
                
            elseif strcmp(PATTERN,'square')
                %change the meaning of npoints
                npoints_per_side = npoints;
                % Syntax requires this to be at least 2, which puts a BC on each corner
                if npoints_per_side == 1
                    npoints_per_side = 2;
                end
                
                points = linspace(-radius,radius,npoints_per_side);
                
                % Initialize vectors
                top = zeros(npoints_per_side,2);
                bottom = top;
                left = top;
                right = top;
                
                for ii = 1:npoints_per_side
                    top(ii,1) = points(ii);
                    top(ii,2) = radius;
                    bottom(ii,1) = points(ii);
                    bottom(ii,2) = -radius;
                    left(ii,1) = -radius;
                    left(ii,2) = points(ii);
                    right(ii,1) = radius;
                    right(ii,2) = points(ii);
                end
                
                % Remove duplicate points
                left = left(2:end-1,:);
                right = right(2:end-1,:);
                
                % Concatenate and place into object
                DM.bconds_ = vertcat(top,left,bottom,right);
            else
                error('PATTERN must be the string circle or square');

            end
			
			DM.bconds_(:,1) = DM.bconds_(:,1) + Offset(2);
			DM.bconds_(:,2) = DM.bconds_(:,2) + Offset(1);

			DM.bconds_(:,3) = 0;
        end % of defineBC
        
        %% Actuators
        
        function DM = addActuators(DM,COORDS,Offset)
        % DM = addActuators(DM,COORDS,Offset)
        %
        % COORDS is [x,y] coordinates of actuators
            if nargin < 3
                Offset = [0,0];
            end
            
            COORDS(:,1) = COORDS(:,1) + Offset(1);
            COORDS(:,2) = COORDS(:,2) + Offset(2);
            
            for ii = 1:size(COORDS,1)
                n = size(DM.actuators_,1)+1;
                DM.actuators_(n,1:2) = COORDS(ii,:);
                DM.actuators_(n,3) = 0;
                DM.actuators_(n,4) = true;
            end
            
            DM.countActuators;
        end % of addActuators
        
        
        function DM = countActuators(DM)
            DM.nActuators_ = size(DM.actuators_,1);
        end % of countActuators
        
        
        function DM = flatten(DM)
            DM.actuators_(:,3) = 0;
        end % of flatten
        
        
        function DM = poke(DM,actuator,displacement)
            DM.flatten;
            DM.actuators_(actuator,3) = displacement;
        end % of poke
        
        
        function DM = addPoke(DM,actuator,bump)
            DM.actuators_(actuator,3) = DM.actuators_(actuator,3) + bump;
        end % of addPoke
        
        
        function DM = updateDM(DM,updates)
            if length(updates)~=DM.nActuators_
                error('Need to provide update value for all actuators');
            end
            DM.actuators_(:,3) = DM.actuators_(:,3) + updates;
        end % of updateDM
        
        function DM = setDM(DM,displacements)
            if length(displacements) ~= DM.nActuators_
                error('Need to provide a displacement for all actuators');
            end
            DM.actuators_(:,3) = displacements;
        end % of setDM
        
        
      %% Influence Functions
      
      %MEMS
      function DM = MEMS_Ifxn(DM,V,d,r,r_a,r_m)
          eps0 = 8.85*10^-12;
          coeff = (eps0 * V^2) / (4*d^2);
          for ii = 1:length(r)
              if r(ii) < r_a
                  phi(ii) = coeff * (r_a^2 - r(ii)^2 + 2*(r_a^2) * log(r_a / r_m));
              elseif r(ii) < r_m
                  phi(ii) = coeff * (2*(r_a^2)*log(r_a / r_m));
              end
          end
          DM.influencefxn_ = phi;
      end % of MEMS_Ifxn
      
      % Gaussian
      function DM = Gaussian_Ifxn(DM,r,r_c,c_a)
          % DM = Gaussian_Ifxn(DM,r,r_c,c_a)
          % r = polar coordinate
          % r_c = interactuator spacing
          % c_a = actuator coupling - [0,1]
          
          coeff = log(c_a) / (r_c^2);
          DM.influencefxn_ = exp(coeff .* r.^2);
      end % of Gaussian_Ifxn
      
      %% Surface Model
      
      % Easy
      function DM = EZmodel(DM)
          g = griddata(double([DM.actuators_(:,1);DM.bconds_(:,1)]), double([DM.actuators_(:,2);DM.bconds_(:,2)]),double([DM.actuators_(:,3);DM.bconds_(:,3)]),DM.X_,DM.Y_,'natural');
          g(isnan(g)) = 0;
          
          if strcmpi(DM.default_data_type,'single')
              g = single(g);
          end
          DM.zsag_ = g;
      end % of EZmodel
      
      %% Plotting
      
      function plotActuators(DM,show_labels)
          if nargin < 2
              show_labels = false;
          end
          hold on
          OnActs = DM.actuators_(:,4)~=0;
          plot(DM.actuators_(OnActs,1), DM.actuators_(OnActs,2),'ko','MarkerSize',2);
          plot(DM.actuators_(~OnActs,1), DM.actuators_(~OnActs,2),'rx','MarkerSize',2);
          plot(DM.bconds_(:,1),DM.bconds_(:,2),'bs');
        if show_labels
            for ii = 1:DM.nActuators_
                text(DM.actuators_(ii,1), DM.actuators_(ii,2),sprintf('%d',ii),'FontSize',8);
            end
        end
        hold off;
      end % of plotActuators
        
    end % of methods
    
    methods(Static=true)
        
        function ActDisp = getDisplacements(DM,WF)
            ActDisp = interp2(DM.X_,DM.Y_,WF,DM.actuators_(:,1),DM.actuators_(:,2));
            
        end % of getDisplacements
        
    end % of static methods
end

