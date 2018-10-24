classdef OptDetector < OptElement
    %OptDetector Class for Detector objects
    %   Placeholder for the end of a system no longer with chirp-z
    %   transform!
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        exposure_time_;
        frame_;
        centerFocal_;
        
        % Detector Pixel Specs
        pixel_Angsize_; % [lambda / D]
        lambdaRef_;
        flD_; % focal length times lambda/D
        FoVflD_; % Field of view of camera in lambda/D
        samplesPerflD_; % number of pixels in each flD
        N_;
        detectorLength_;
        
        
        % Noise Characteristics
        useNoise_;
        useReadNoise_ = 1;
        useShotNoise_ = 1;
        ReadNoise_ = 1;
        
        % Coordinates
        magnification_ = 1;
        dx_;
        dy_;
        x_;
        xlD_;
        y_;
        ylD_;
        X_;
        Y_;
        XlD_;
        YlD_;
        R_;
        RlD_;
        
    end
    
    methods
        %% Constructor
        function elem = OptDetector(PROPERTIES)
            % elem = OptDetector(PROPERTIES)
            % PROPERTIES is a 5x1 cell array containing element properties
            % If something is unknown or to be set later, leave that
            % position in A empty:
            %
            %           PROPERTIES{1,1} = name                  [string]
            %           PROPERTIES{2,1} = samplesPerflD_        [int]
            %           PROPERTIES{3,1} = lambdaRef             [float in m]
            %           PROPERTIES{4,1} = flD                   [float in m]
            %           PROPERTIES{5,1} = FoVflD                [float]
            %           PROPERTIES{6,1} = pixel_Angsize         [float in l/d]
            %           PROPERTIES{7,1} = exposure_time         [float in seconds]
            %           PROPERTIES{8,1} = useNoise              [Boolean]
            
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
            elem.set_samplesPerflD(A{2,1});
            elem.set_lambdaRef(A{3,1});
            elem.set_flD(A{4,1});
            elem.set_FoVflD(A{5,1});
            elem.set_pixel_Angsize(A{6,1});
            elem.set_noiseFlag(A{8,1});
            elem.set_exposureTime(A{7,1});
            elem.compute_pixel_Angsize;
            elem.compute_detector_px_spacing;
            elem.compute_detector_numPixels;
            elem.compute_detector_coords;
            elem.gridsize_ = size(elem.R_);
            elem.set_z_position(0); %isn't used, but needed to make code happy
            elem.setdatatype();
            elem.set_propagation_method(1);
            elem.addnewline(2);
            
        end % of contructor
        
        function elem = set_lambdaRef(elem,lambda)
            % elem = set_lambdaRef(elem,lambda)
            elem.lambdaRef_ = lambda;
        end % of set_lambdaRef
        
        function elem = set_flD(elem,flD)
            % elem = set_flD(elem,flD)
            elem.flD_ = flD;
        end % of set_flD
        
        function elem = set_FoVflD(elem, nld)
            % elem = set_FPregion(elem, nld)
            elem.FoVflD_ = nld;
            
        end % of set_FPregion
        
        function elem = set_pixel_Angsize(elem,M)
            % elem = set_pixel_size(elem,M)
            
            elem.pixel_Angsize_ = M;
        end % of set_pixel_size
        
        
        function elem = set_noiseFlag(elem,flag)
            % elem = set_noiseFlag(elem,flag)
            elem.useNoise_ = flag;
            
        end % of set_noiseFlag
        
        function elem = set_exposureTime(elem,dt)
            % elem = set_exposureTime(elem,dt)
            elem.exposure_time_ = dt;
            
        end % of set_exposureTime
        
        function elem = set_samplesPerflD(elem,samplesPerflD)
            % elem = set_samplesPerflD(elem,samplesPerflD)
            elem.samplesPerflD_ = samplesPerflD;
        end% of set_samplesPerflD
        
        function elem = set_magnification(elem, magnification)
            % elem = set_magnification(elem, magnification)
            elem.magnification_ = magnification;
        end % of set_magnification
        
        function elem = set_Noise_characteristics(elem, useShotNoise, useReadNoise, ReadNoise)
            % elem = set_Noise_characteristics(elem, useShotNoise, useReadNoise, ReadNoise)
            elem.useShotNoise_ = useShotNoise;
            elem.useReadNoise_ = useReadNoise;
            elem.ReadNoise_ = ReadNoise;
        end % of set_Noise_characteristics
        
        function elem = set_centerFocal(elem, centerFocal)
            % elem = set_centerFocal(elem, centerFocal)
            elem.centerFocal_ = centerFocal;
        end % of set_centerFocal
        
        function descr = describe(elem)
            % elem = describe(elem)
            
            elemtype = 'OptDetector';
            
            descr = sprintf('%s:\nElement is an %s',elem.name,elemtype);
            
        end
        
        %% Detector Specs
        
        function [elem,pxSz] = compute_pixel_Angsize(elem)
            % [elem, pxSz] = compute_pixel_Angsize(elem)
            
            pxSz = 1 / elem.samplesPerflD_;
            elem.set_pixel_Angsize(pxSz);
        end % of compute_pixel_Angsize
        
        function elem = compute_detector_px_spacing(elem)
            elem.dx_ = elem.flD_ / elem.samplesPerflD_;
            elem.dy_ = elem.flD_ / elem.samplesPerflD_;
        end % of compute_detector_px_spacing
        
        function elem = compute_detector_numPixels(elem)
            elem.N_ = ceil(elem.FoVflD_ * elem.samplesPerflD_);
            elem.detectorLength_ = elem.N_ * elem.dx_;
        end % of compute_detector_numPixels
        
        function elem = compute_detector_coords(elem)
            elem.x_ = -(elem.detectorLength_ - elem.dx_)/2 : elem.dx_ : (elem.detectorLength_ - elem.dx_)/2;
            elem.xlD_ = elem.x_ / elem.flD_;
            elem.y_ = -(elem.detectorLength_ - elem.dy_)/2 : elem.dy_ : (elem.detectorLength_ - elem.dy_)/2;
            elem.ylD_ = elem.y_ / elem.flD_;
            [elem.X_, elem.Y_] = meshgrid(elem.x_,elem.y_);
            elem.XlD_ = elem.X_ / elem.flD_;
            elem.YlD_ = elem.Y_ / elem.flD_;
            elem.R_ = sqrt(elem.X_.^2 + elem.Y_.^2);
            elem.RlD_ = elem.R_ / elem.flD_;
        end % of compute_detector_coords
        
        function elem = detector_add_noise(elem,ReadNoise, useShotNoise, useReadNoise)
            % elem = detector_add_noise(elem,ReadNoise, useShotNoise, useReadNoise)
            if elem.useNoise_
                frame = double(elem.frame_);
                if useShotNoise
                    frame = frame * (1e-12);
                    frame = (1e12) * imnoise(frame,'poisson');
                end
                if useReadNoise
                    frame = frame + ((randn(size(frame)).* (frame.^0.5)) * ReadNoise);
                end
                elem.frame_ = single(frame);
            end
        end % of detector_add_noise
        
        function [elem, field] = make_PSF(elem, WF, focal_length)
            dataType = 'double';
            
            info = struct;
            info.x1 = elem.x_(1);
            info.y1 = elem.y_(1);
            info.dx0 = WF.dx_;
            info.dy0 = WF.dx_;
            info.dx1 = elem.dx_;
            info.dy1 = elem.dy_;
            info.f = focal_length;
            info.lambdaList = WF.lambda_array_;
            
            [~,~,~,WFx] = WF.Coords2D;
            WFy = WFx;
            
            tmp = gpuArray(init_variable(size(WF.field_,1),size(WF.field_,2),size(WF.field_,3),dataType,0));
            t1 = czt_GPU_v2(WF.field_, elem.N_, info);
            
            % Apply the Multiplicative factor
            if length(info.lambdaList) > 1
                if info.lambdaList(2) ~= info.lambdaList(1) % if multiple wavelengths, use them all
                    for thisLambda = 1:length(lambda)
                        t2 = exp(-2*pi*1i * (WFx(1)*ones(elem.N_,1,'double')*elem.x_ + WFy(1)*elem.y_'*ones(1,elem.N_,'double'))./(info.lambdaList(thisLambda)*focal_length));
                        tmp(:,:,thisLambda) = t1.*t2 * (WF.dx_ * WF.dx_) / (1i*info.lambdaList(thisLambda)*focal_length);
                    end
                else % otherwise assume only 1 wavelength is needed and save time
                    t2 = exp(-2*pi*1i * (WFx(1)*ones(elem.N_,1,'double')*elem.x_ + WFy(1)*elem.y_'*ones(1,elem.N_,'double'))./(info.lambdaList*focal_length));
                    tmp = bsxfun(@times,t1, t2) * (WF.dx_ * WF.dx_) ./(1i*info.lambdaList(1)*focal_length);
                end
            else
                t2 = exp(-2*pi*1i * (WFx(1)*ones(elem.N_,1,'double')*elem.x_ + WFy(1)*elem.y_'*ones(1,elem.N_,'double'))./(info.lambdaList*focal_length));
%                 tmp = bsxfun(@times,t1, t2) * (WF.dx_ * WF.dx_) ./(1i*info.lambdaList(1)*focal_length);
                tmp = bsxfun(@times,t1, t2) * (WF.dx_ * WF.dx_) ./(1i*info.lambdaList(1)*focal_length);

            end
            tmp = gather(tmp);
%             WF.set_field(tmp);
            
%             elem.frame_ = abs(tmp).^2 * elem.dx_ * elem.dy_ * elem.exposure_time_;
%             elem.frame_ = abs(tmp).^2;
            pxsizeratio = elem.dx_ / WF.dx_;
            field = tmp .* pxsizeratio;
            elem.frame_ = abs(field).^2;
            elem.detector_add_noise(elem.ReadNoise_, elem.useShotNoise_, elem.useReadNoise_);
            
        end % of make_PSF
        
        end % of methods
        
end
    
