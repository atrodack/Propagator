classdef OptMaterial < matlab.mixin.Copyable
    %OPTMATERIAL Class to Store Optical Materials
    %   Detailed explanation goes here
    
    
    %% Properties
    properties(Constant=true, GetAccess='protected')
        
        
    end % of constant properties
    
    properties(GetAccess = 'public', SetAccess = 'public')
        lambda; % must be in meters
        verbose = 0;
        
    end % of public properties
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        material = struct;
        material_code = 0;
        n_;
        k_;
        lambdaum_;

        
    end % of protected properties
    
    %% Methods
    methods
        %% Constructor
        
        function OM = OptMaterial(code,lambda)
            % OM = OptMaterial(code)
            % constructs the material structure for material given by code
            %
            %===========
            % Code List
            %===========
            % 0) Vacuum
            % 1) Air
            % 2) Mirror
            % 3) Si
            % 4) SiO2
            % 5) BK7
            % 6) SF10
            % 7) CaF2
            % 8) Aluminium
            % 9) Gold
            % 10) Silver
            
            
            if nargin < 2
                OM.lambda = 550e-9; % default to photopic (units of meters)
                OM.initMaterial(code);
                
            elseif nargin == 2
                OM.lambda = lambda;
                OM.initMaterial(code);
                
            else
                error('InputArgs:IncorrectNum','Incorrect Number of Input Arguments\nSee help for OptMaterial');
            end
            
            
            
        end % of constructor
        
        %% Material Properties
        
        function OM = initMaterial(OM,code)
            OM.material_code = code;
            switch code
                case 0
                    OM.material.name = 'Vacuum';
                    OM.n_ = ones(1,length(OM.lambda));
                    OM.lambdaum_ = OM.lambda * 1e6;
                    OM.material.band = [0, inf];
                    
                case 1 % Ciddor 1996: n 0.23-1.690 µm
                    OM.material.name = 'Air';
                    OM.material.Abbe = 89.30;
                    OM.material.band = [0.23, 1.69];
                    
                    OM.material.B1 = 0.05792105;
                    OM.material.B2 = 0.00167917;
                    OM.material.C1 = 238.0185;
                    OM.material.C2 = 57.362;
                    
                    OM.SellmeierDispersion_Air();
                    
                case 2
                    OM.material.name = 'Mirror';
                    OM.n_ = 2*ones(1,length(OM.lambda));
                    OM.lambdaum_ = OM.lambda * 1e6;
                    OM.material.band = [0, inf];
                    
                case 3 % Chandler-Horowitz and Amirtharaj 2005: n 2.5-22.222 µm
                    OM.material.name = 'Si';
                    OM.material.band = [2.5, 22.222];
                    
                    OM.material.A = 11.67316;
                    OM.material.B2 = 0.004482633;
                    OM.material.C2 = 1.108205^2;
                                        
                    OM.SellmeierDispersion_Si();
                    
                case 4 % Malitson 1965: n 0.21-3.71 µm
                    OM.material.name = 'SiO2';
                    OM.material.Abbe = 67.82;
                    OM.material.band = [0.21, 3.71];
                    
                    OM.material.A = 0;
                    OM.material.B1 = 0.6961663;
                    OM.material.B2 = 0.4079426;
                    OM.material.B3 = 0.8974794;
                    OM.material.C1 = 0.0684043^2;
                    OM.material.C2 = 0.1162414^2;
                    OM.material.C3 = 9.896161^2;
                    
                    OM.SellmeierDispersion();
                    
                case 5 % Schott Catalog: n 0.3-2.5 µm
                    OM.material.name = 'BK7';
                    OM.material.Abbe = 64.17;
                    OM.material.band = [0.3, 2.5];
                    
                    OM.material.A = 0;
                    OM.material.B1 = 1.03961212;
                    OM.material.B2 = 0.231792344;
                    OM.material.B3 = 1.01046945;
                    OM.material.C1 = 0.00600069867;
                    OM.material.C2 = 0.0200179144;
                    OM.material.C3 = 103.560653;
                    
                    OM.SellmeierDispersion();
                    
                case 6 % Schott Catalog: n 0.38-2.5 µm
                    OM.material.name = 'SF10';
                    OM.material.Abbe = 28.53;
                    OM.material.band = [0.38, 2.5];
                    
                    OM.material.A = 0;
                    OM.material.B1 = 1.62153902;
                    OM.material.B2 = 0.256287842;
                    OM.material.B3 = 1.64447552;
                    OM.material.C1 = 0.0122241457;
                    OM.material.C2 = 0.0595736775;
                    OM.material.C3 = 147.468793;
                    
                    OM.SellmeierDispersion();
                    
                case 7 % Li 1980: n 0.15-12.0 µm; 20 °C
                    OM.material.name = 'CaF2';
                    OM.material.Abbe = 95.31;
                    OM.material.band = [0.15, 12.0];
                    
                    OM.material.A = 0.33973;
                    OM.material.B1 = 0.69913;
                    OM.material.B2 = 0.11994;
                    OM.material.B3 = 4.35181;
                    OM.material.C1 = 0.09374^2;
                    OM.material.C2 = 21.18^2;
                    OM.material.C3 = 38.46^2;
                    
                    OM.SellmeierDispersion();
                    
                case 8 % Rakić 1998: n,k 0.2066-12.40 µm
                    OM.material.name = 'Aluminium';
                    OM.material.band = [0.2066, 12.40];
                    OM.getMetalData(8);
                    
                case 9 % Rakić 1998: n,k 0.2066-12.40 µm
                    OM.material.name = 'Gold';
                    OM.material.band = [0.2066, 12.40];
                    OM.getMetalData(9);
                    
                case 10 % Rakić 1998: n,k 0.2066-12.40 µm
                    OM.material.name = 'Silver';
                    OM.material.band = [0.2066, 12.40];
                    OM.getMetalData(10);
                    
                otherwise
                    error('MaterialError:MaterialnotSupported','That material is not (currently) supported');
                    
            end % of switch
            
        end % of initMaterial
        
        %% Utilities
        function OM = convert2Micron(OM)
            % OM = convert2Micron(OM)
            % Internal function for storing lambda property (meters) in
            % units of microns in order to use the Sellmeier dispersion
            % formula
            
            OM.lambdaum_ = OM.lambda * 1e6;
        end % of convert2Micron
        
        
        function PlotIndexvsWavelength(OM)
            % PlotIndexvsWavelength(OM)
            % Plots n vs lambda
            
            figure(8888);
            plot(OM.lambdaum_,OM.n_);
            axis square;
            xlabel('\lambda [\mum]');
            ylabel('n(\lambda)');
            title(sprintf('Refractive Index vs. Wavelength for %s',OM.material.name));
            
        end % of PlotIndexvsWavelength
        
        function PlotIndexPoint(OM,lambda,n)
            % PlotIndexPoint(OM,n)
            % Plots the found index on the curve of index vs. wavelength
            
            figure(8888)
            hold on
            plot(lambda,n,'kp');
            hold off
        end % of PlotIndexPoint
        
        function PlotIndexk(OM)
            % PlotIndexk(OM)
            % plots the imaginary term of the refractive index
            % if present in data
            
            figure(8888)
            hold on
            plot(OM.lambdaum_,OM.k_,'r');
            hold off
            legend('n','k','Location','Best');
        end % of PlotIndexk
        
        %% Index of Refraction Computations
        function n = SellmeierDispersion(OM,lambda_)
            % n = SellmeierDispersion(OM,)
            % Internal function to compute index for material at input
            % wavelength(s). Plots if verbose is on.
            
            if nargin < 2
                OM.convert2Micron;
                lambda_ = OM.lambdaum_;
            else
                lambda_ = lambda_ * 1e6;
            end
            
             band = [min(lambda_), max(lambda_)];
            if sum(isinf(band)) == 0
                if band(1) < OM.material.band(1) || band(2) > OM.material.band(2)
                    error('Asking for wavelengths not supported by implementation of Sellmeier Formula\nSupported band is [%g %g] µm\n',OM.material.band(1),OM.material.band(2));
                end
            else
                error('Why do you care what happens at infinite wavelength?');
            end
            
            m = OM.material.A + ((OM.material.B1 * lambda_.^2)./(lambda_.^2 - OM.material.C1)) + ((OM.material.B2 * lambda_.^2)./(lambda_.^2 - OM.material.C2)) + ((OM.material.B3 * lambda_.^2)./(lambda_.^2 - OM.material.C3));
            nn = m + 1;
            n = sqrt(nn);
            
            if nargin < 2
                OM.n_ = n;
                if OM.verbose == 1
                    OM.PlotIndexvsWavelength;
                end
            else
                if OM.verbose == 1
                    OM.PlotIndexvsWavelength;
                    OM.PlotIndexPoint(lambda_,n);
                end
            end

        end % of SellmeierDispersion
        
        
        function n = SellmeierDispersion_Si(OM,lambda_)
            % n = SellmeierDispersion_Si(OM)
            % Internal function to compute index for Silicon
            
            if nargin < 2
                OM.convert2Micron;
                lambda_ = OM.lambdaum_;
            else
                lambda_ = lambda_ * 1e6;
            end
            
            band = [min(lambda_), max(lambda_)];
            if sum(isinf(band)) == 0
                if band(1) < OM.material.band(1) || band(2) > OM.material.band(2)
                    error('Asking for wavelengths not supported by implementation of Sellmeier Formula\nSupported band is [%g %g] µm\n',OM.material.band(1),OM.material.band(2));
                end
            else
                error('Why do you care what happens at infinite wavelength?');
            end
            
            nn = OM.material.A + (1./lambda_.^2) + ((OM.material.B2) ./ (lambda_.^2 - OM.material.C2));
            n = sqrt(nn);
            
            if nargin < 2
                OM.n_ = n;
                if OM.verbose == 1
                    OM.PlotIndexvsWavelength;
                end
            else
                if OM.verbose == 1
                    OM.PlotIndexvsWavelength;
                    OM.PlotIndexPoint(lambda_,n);
                end
            end
        end % of SellmeierDispersion_Si
        
        function n = SellmeierDispersion_Air(OM,lambda_)
            % n = SellmeierDispersion_Si(OM)
            % Internal function to compute index for Air
            
            if nargin < 2
                OM.convert2Micron;
                lambda_ = OM.lambdaum_;
            else
                lambda_ = lambda_ * 1e6;
            end
            
             band = [min(lambda_), max(lambda_)];
            if sum(isinf(band)) == 0
                if band(1) < OM.material.band(1) || band(2) > OM.material.band(2)
                    error('Asking for wavelengths not supported by implementation of Sellmeier Formula\nSupported band is [%g %g] µm\n',OM.material.band(1),OM.material.band(2));
                end
            else
                error('Why do you care what happens at infinite wavelength?');
            end
            
            m = ((OM.material.B1) ./ (OM.material.C1 - lambda_.^-2)) + ((OM.material.B2) ./ (OM.material.C2 - lambda_.^-2));
            n = m + 1;
            
            if nargin < 2
                OM.n_ = n;
                if OM.verbose == 1
                    OM.PlotIndexvsWavelength;
                end
            else
                if OM.verbose == 1
                    OM.PlotIndexvsWavelength;
                    OM.PlotIndexPoint(lambda_,n);
                end
            end
            
        end % of SellmeierDispersion_Air
        
        function OM = getMetalData(OM,material_code)
            % OS = getMetalData(OS,material_code)
            % loads in data from Rakić 1998: n,k 0.2066-12.40 µm
            home_direc = pwd;

            switch material_code
                case 8 % Aluminium
                    filename = sprintf('%s/Propagator/Database/Aluminium.mat',home_direc);
                    load(filename);
                    lambdalist = Material_data(:,1);
                    n = Material_data(:,2);
                    k = Material_data(:,3);
                    
                case 9 % Gold
                    filename = sprintf('%s/Propagator/Database/Gold.mat',home_direc);
                    load(filename);
                    lambdalist = Material_data(:,1);
                    n = Material_data(:,2);
                    k = Material_data(:,3);

                case 10 % Silver
                    filename = sprintf('%s/Propagator/Database/Silver.mat',home_direc);
                    load(filename);
                    lambdalist = Material_data(:,1);
                    n = Material_data(:,2);
                    k = Material_data(:,3);

                otherwise
                    
            end
            OM.n_ = n;
            OM.k_ = k;
            OM.lambdaum_ = lambdalist;
            OM.lambda = lambdalist *1e-6;
            
            if OM.verbose == 1
                OM.PlotIndexvsWavelength;
                OM.PlotIndexk;
            end
            
        end % of getMetalData
        
        %% Phase
        function phasefac = ComputePhaseFactor(OM, lambda)
            % phasefac = ComputePhase(Om,lambda)
            % returns the phase factor (phase = k*OPL = (2pi/lambda)*n*z)
            
            material_code_ = OM.material_code;
            
            %             nair = 1; %Vacuum
            numLambdas = length(lambda);
            
            switch material_code_
                case 0
                    n = 1 * ones(1,numLambdas);
                case 1
                    n = OM.SellmeierDispersion_Air(lambda);
                case 2
                    n = 2 * ones(1,numLambdas);
                case 3
                    n = OM.SellmeierDispersion_Si(lambda);
                case 4
                    n = OM.SellmeierDispersion(lambda);
                case 5
                    n = OM.SellmeierDispersion(lambda);
                case 6
                    n = OM.SellmeierDispersion(lambda);
                case 7
                    n = OM.SellmeierDispersion(lambda);
                case 8
                    error('Not supported yet');
                case 9
                    error('Not supported yet');
                case 10
                    error('Not supported yet');
                otherwise
                    % Shouldn't be able to get this far, but just in case
                    error('MaterialError:MaterialnotSupported','That material is not (currently) supported');
            end
            
            % phasefac = (2*pi / lambda) * (n - nair);
            phasefac = (2*pi ./ lambda) .* n;
            
            
        end % of ComputePhaseFactor
        
        function descr = describe(OM)
            % descr = describe(OM)
            
            descr = sprintf('%s with Wavelength Band [%0.3f %0.3f]',OM.material.name, OM.material.band);
        end % of describe
        
    end % of methods
    
end

