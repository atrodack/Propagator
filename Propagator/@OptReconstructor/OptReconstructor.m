classdef OptReconstructor < matlab.mixin.Copyable
    %OPTRECONSTRUCTOR Class for training and building AO reconstructors
    %   
    
    properties(GetAccess='public',SetAccess='public')
        verbose = 0;
    end % of public properties
    
    properties(GetAccess='public',SetAccess='private')
        A_;
        DM_;
        WFS_;
        lambda_;
        
        SLOPES_;
        Training_method_;
        
        RECONSTRUCTOR_;
        ACTS_;
        U_;
        s_;
        V_;
        
        Nmodes_;
        D_;
        
    end % of protected properties
    
    
    
    
    
    methods
        %% Constructor
        function recon = OptReconstructor(A,DM,WFS)
            % RECON = OptReconstructor(A,DM,WFS)
            %
            % A = OptPupil class object
            % DM = OptDM class object
            % WFS = OptWFS Child Class object
            % PS = empty OptPhaseScreen class object
            % lambda is wavelength
            
            recon.A_ = A;
            recon.DM_ = DM;
            recon.WFS_ = WFS;
            
        end % of constructor
        
        %% Build Reconstructor
        
        function recon = buildRecon(recon,cutoff)
            if nargin < 2
                cutoff = 1e-2;
            end
            
            if isempty(recon.s_)
                error('Train the Reconstructor before building!');
            end
            
            if cutoff < 1
                MODES = (recon.s_ / recon.s_(1) >= cutoff);
                recon.Nmodes_ = sum(MODES);
            else
                if(cutoff>length(recon.s_))
                    cutoff = length(recon.s_);
                end
                MODES = 1:round(cutoff);
                recon.Nmodes_ = MODES(end);
            end
            
            fprintf('Reconstructor Correcting %d modes. \n',recon.Nmodes_);
            
            MODES(recon.s_(MODES)<1e-16)=[];
            
            recon.RECONSTRUCTOR_ = ((recon.V_(:,MODES) * ...
                (diag(1./recon.s_(MODES)) * recon.U_(:,MODES)')) ...
                * recon.ACTS_')' ;
            
        end % of buildRecon
        
        %% Training Methods
        
        function recon = ZernikeTraining(recon,D,Zmax,field)
            % recon = ZernikeTraining(recon,D,Nmax,field,PS)
            %
            % D = diameter of pupil
            % Zmax is highest Zernike Noll Index
            % field is an OptWF class object
            
            % Initialize some variables
            recon.Training_method_ = 'zernike';
            recon.D_ = D;
            recon.lambda_ = field.lambda0_;
            k = (2*pi) / recon.lambda_;
            N = recon.WFS_.N_;
            x = recon.WFS_.x_;
            y = recon.WFS_.y_;
            dx = x(2) - x(1);
            dy = y(2) - y(1);
            flag = 1;
            
            % Prepare Zernike Basis
            fprintf('Setting up Zernike Basis set:\n');
            Noll_list = 2:Zmax;
            nZmodes = length(Noll_list);
            recon.Nmodes_ = nZmodes;
            recon.ACTS_ = zeros(recon.DM_.nActuators_,nZmodes);
            recon.SLOPES_ = zeros(2*recon.WFS_.nUnMaksedsubAps_,nZmodes);
            [Z_basis,~] = Zernike_Basis(Noll_list,ones(nZmodes,1),(D*1.125)/dx,N);
            Z = zeros(N,N,nZmodes);
            for ii = 1:nZmodes
                [nn,mm] = Noll(Noll_list(ii));
                weight = (nn^2 + mm^2) ^ (-5/6);
                Z(:,:,ii) = weight * reshape(Z_basis(ii,:),N,N) * recon.lambda_/4 / nn;
            end
            clear Z_basis
            
            % Set WFS bias
            field.planewave(1,1);
            field.ApplyElement(recon.A_,1); % Assuming in Air
            recon.WFS_.sense(field.field_,recon.lambda_).setBias;
            fprintf('Training with Zernike Polynomials: Nmodes = %d:\n',nZmodes);
            
            
            nn = 0;
            % Train
            for ii = 1:nZmodes
                nn_prev = nn;
                [nn,mm] = Noll(Noll_list(ii));
                if nn_prev ~= nn
                    flag = 1;
                else
                    flag = 0;
                end
                if flag
                    fprintf('\n');
                    fprintf('%d:',nn);
                end
                fprintf('\t %d',mm); 
                    
                ActDisp = recon.DM_.getDisplacements(recon.DM_,Z(:,:,ii));
                recon.DM_.setDM(ActDisp).EZmodel;
                
                field.planewave(1,1).ApplyElement(recon.A_,1).ApplyPhaseScreen(-recon.DM_.zsag_,recon.lambda_).ReIm2WF;
                recon.WFS_.sense(field.field_,recon.lambda_);
                
                recon.ACTS_(:,ii) = recon.DM_.actuators_(:,3);
                recon.SLOPES_(:,ii) = recon.WFS_.getSlopesArray;
                
                if(recon.verbose)
                    clf;
                    hold on;
                    subplot(1,2,1)
                    imagesc(x,x,field.pha_/k);
                    axis square;
                    axis xy;
                    title('SHWFS over Field Phase')
                    colormap(gray);
                    colorbar;
                    recon.WFS_.plotSubAps('b-');
                    recon.WFS_.quiver('g-');
                    
                    subplot(1,2,2)
                    imagesc(x,y,recon.A_.zsag_ .* recon.DM_.zsag_);
                    recon.DM_.plotActuators(false);
                    axis square;
                    axis xy;
                    colorbar;
                    title('DM Surface');
                    hold off;
                    drawnow;
%                     pause(0.1);
                end
            end % of Training for loop
            fprintf('\n');
            recon.SLOPES_(isnan(recon.SLOPES_)) = 0;
            recon.ACTS_(isnan(recon.ACTS_)) = 0;
            [UU,SS,VV] = svd(recon.SLOPES_','econ');
            recon.U_ = UU;
            recon.s_ = diag(SS);
            recon.V_ = VV;
            clear UU SS VV
            recon.buildRecon;
        end % of ZernikeTraining
        
        function recon = KLTraining(recon,D,Zmax,r0,field)
            % recon = KLTraining(recon,D,Nmax,field,PS)
            %
            % D = diameter of pupil
            % Zmax is highest Zernike Noll Index
            % field is an OptWF class object
            
            % Initialize some variables
            recon.Training_method_ = 'Karhunen–Loève';
            recon.D_ = D;
            recon.lambda_ = field.lambda0_;
            k = (2*pi) / recon.lambda_;
            N = recon.WFS_.N_;
            x = recon.WFS_.x_;
            y = recon.WFS_.y_;
            dx = x(2) - x(1);
            dy = y(2) - y(1);
            
            % Prepare KL Basis
            fprintf('Setting up Karhunen–Loève Basis Set\n');
            Noll_list = 2:Zmax;
            nKLmodes = length(Noll_list);
            recon.Nmodes_ = nKLmodes;
            recon.ACTS_ = zeros(recon.DM_.nActuators_,nKLmodes);
            recon.SLOPES_ = zeros(2*recon.WFS_.nUnMaksedsubAps_,nKLmodes);
            [KL,~,~] = KL_Basis(Noll_list, r0, recon.D_/2, dx, N);
            KL_Cube = init_variable(N,N,length(Noll_list),'single',0);
            for ii = 1:nKLmodes
                weight = 1;
                KL_Cube(:,:,ii) = weight * reshape(KL(ii,:),N,N) * recon.lambda_/4;
            end
            clear KL
            
            % Set WFS bias
            field.planewave(1,1);
            field.ApplyElement(recon.A_,1); % Assuming in Air
            recon.WFS_.sense(field.field_,recon.lambda_).setBias;
            fprintf('Training with Karhunen–Loève Polynomials: Nmodes = %d:\n',nKLmodes);
            % Train
            for ii = 1:nKLmodes
                fprintf('*');
                ActDisp = recon.DM_.getDisplacements(recon.DM_,KL_Cube(:,:,ii));
                recon.DM_.setDM(ActDisp).EZmodel;
                
                field.planewave(1,1).ApplyElement(recon.A_,1).ApplyPhaseScreen(-recon.DM_.zsag_,recon.lambda_).ReIm2WF;
                recon.WFS_.sense(field.field_,recon.lambda_);
                
                recon.ACTS_(:,ii) = recon.DM_.actuators_(:,3);
                recon.SLOPES_(:,ii) = recon.WFS_.getSlopesArray;
                
                if(recon.verbose)
                    clf;
                    hold on;
                    subplot(1,2,1)
                    imagesc(x,x,field.pha_/k);
                    axis square;
                    axis xy;
                    title('SHWFS over Field Phase')
                    colormap(gray);
                    colorbar;
                    recon.WFS_.plotSubAps('b-');
                    recon.WFS_.quiver('g-');
                    
                    subplot(1,2,2)
                    imagesc(x,y,recon.A_.zsag_ .* recon.DM_.zsag_);
                    recon.DM_.plotActuators(false);
                    axis square;
                    axis xy;
                    colorbar;
                    title('DM Surface');
                    hold off;
                    drawnow;
%                     pause(0.1);
                end
            end % of Training for loop
            fprintf('\n');
            recon.SLOPES_(isnan(recon.SLOPES_)) = 0;
            recon.ACTS_(isnan(recon.ACTS_)) = 0;
            [UU,SS,VV] = svd(recon.SLOPES_','econ');
            recon.U_ = UU;
            recon.s_ = diag(SS);
            recon.V_ = VV;
            clear UU SS VV
            recon.buildRecon;
        end % of KLTraining
        
        
        function recon = PokeTraining(recon,D,field)
            % recon = PokeTraining(recon,D,field,PS)
            %
            % D = diameter of pupil
            % Zmax is highest Zernike Noll Index
            % field is an OptWF class object
            
            % Initialize some variables
            recon.Training_method_ = 'zernike';
            recon.D_ = D;
            recon.lambda_ = field.lambda0_;
            k = 2*pi / recon.lambda_;
            x = recon.WFS_.x_;
            y = recon.WFS_.y_;

            
            % Prepare Poke Basis
            Poke_list = 1:recon.DM_.nActuators_;
            nmodes = length(Poke_list);
            recon.Nmodes_ = nmodes;
            recon.ACTS_ = zeros(recon.DM_.nActuators_,nmodes);
            recon.SLOPES_ = zeros(2*recon.WFS_.nUnMaksedsubAps_,nmodes);

            % Set WFS bias
            field.planewave(1,1);
            field.ApplyElement(recon.A_,1); % Assuming in Air
            recon.WFS_.sense(field.field_,recon.lambda_).setBias;
            fprintf('Training with DM pokes\n');
            
            % Train
            for ii = 1:nmodes
                fprintf('*');
                recon.DM_.poke(ii,1e-6);
                recon.DM_.EZmodel;
                
                field.planewave(1,1).ApplyPhaseScreen(-recon.DM_.zsag_,recon.lambda_).ReIm2WF;
                recon.WFS_.sense(field.field_,recon.lambda_);
                
                recon.ACTS_(:,ii) = recon.DM_.actuators_(:,3);
                recon.SLOPES_(:,ii) = recon.WFS_.getSlopesArray;
                
                if(recon.verbose)
                    clf;
                    hold on;
                    subplot(1,2,1)
                    imagesc(x,x,field.pha_/k);
                    axis square;
                    title('SHWFS over Field Phase')
                    colormap(gray);
                    colorbar;
                    recon.WFS_.plotSubAps('b-');
                    recon.WFS_.quiver('g-');
                    
                    subplot(1,2,2)
                    imagesc(x,y,recon.DM_.zsag_);
                    recon.DM_.plotActuators(false);
                    axis square;
                    colorbar;
                    title('DM Surface');
                    hold off;
                    drawnow;
%                     pause(0.1);
                end
            end % of Training for loop
            fprintf('\n');
            
            recon.SLOPES_(isnan(recon.SLOPES_)) = 0;
            recon.ACTS_(isnan(recon.ACTS_)) = 0;
            
            [UU,SS,VV] = svd(recon.SLOPES_','econ');
            recon.U_ = UU;
            recon.s_ = diag(SS);
            recon.V_ = VV;
            clear UU SS VV
            
            recon.buildRecon;
        end % of PokeTraining
        
    end
    
end

