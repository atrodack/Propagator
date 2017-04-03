%% Start Clean
% Avoid problems with making source code changes to objects and then trying
% to test with an old object. Just destroy them!

clear all; clc; close all;

%% Initialize Optical System Parameters

% Define a central wavelength, in meters
lambda_0 = 565*1e-9;

% Define a bandwidth
bandwidth = 0.2;
lambda = linspace(lambda_0 - (lambda_0*(bandwidth/2)),lambda_0 + (lambda_0*(bandwidth/2)),10);

% Define index of refraction of medium Optical System is in
n0 = 1;

% Beam Diameter
D = 0.022 * 2;

% lambda / D
ld = lambda_0 / D;

%% Make the OptSys Object

% Construct
OS = OptSys();

% set the name
OS.name = 'Telescope Optical System';

% set the wavelength info
OS.setCentralWavelength(lambda_0);
OS.setLambdaarray(lambda);

% Set the beam radius (in pixels)
OS.setBeamrad(200);

% Set the pixel scale
OS.setPscale(0.00011);
OS.setWavelengthScales;

% Set the f-number
OS.setFnum(80);

% Set whether or not to save intermediate steps as FITS
OS.savefile = 0;

% Turn off verbose to suppress messages of adding elements in
OS.toggle_verbose('off');


%% Make an OptPupil Object
% Circular aperture
% 30% central obscuratoin
% 10% D spiders

% Make the properties cell array
props = cell(6,1);
props{1} = 'System Pupil';          % string of your choosing for name of object
props{2} = 0;                       % code for material
props{3} = true;                    % amplitude only flag
props{4} = 0;                       % z position
props{5} = D;                       % Physical Diameter [m]
props{6} = 'circ_pup1024.fits';     % mask [matrix or /path/filename.fits]

% Construct
PUPIL = OptPupil(props);
% PUPIL.describe
% PUPIL.show;

%% Make an OptDetector Object
props = cell(5,1);
props{1} = 'Detector';          % string of your choosing for name of object
props{2} = D;                   % diameter
props{3} = ones(1024);          % zsag
props{4} = 20;                  % nld
props{5} = 0.05;                % detector pixel size in fractions of lambda/D
DET = OptDetector(props);

%% Build the Optical System

% Make a list of the made elements
ELEMENTS = cell(2,1);
ELEMENTS{1} = PUPIL;
ELEMENTS{2} = DET;

% Add Elements to System
OS.addSequentialElements(ELEMENTS);

% Turn verbose back on to plot intermediate Propagator steps
OS.toggle_verbose('on');

%% Propagate through the Optical System

% Set the input field
OS.planewave(single(1),length(OS.lambda_array_));
% OS.show;

% Propagate
tic;
OS.PropagateSystem1(1,2,n0);
toc

% combined_exposure = 0;
% for ii = 1:length(lambda)
%     combined_exposure = combined_exposure + OS.PSF_(:,:,ii);
% end
% combined_exposure = combined_exposure / length(lambda);
% 
% [ldx,ldy] = OS.FPcoords(DET.FPregion_,1024);
% 
% 
% 
% figure(2);
% imagesc(ldx(:,:,length(lambda)/2),ldy(:,:,length(lambda)/2),log10(OS.PSF_(:,:,length(lambda)/2) / max(max(OS.PSF_(:,:,length(lambda)/2)))),[-4,0])
% axis xy; axis square;
% colorbar;
% colormap(gray(256));
% xlabel('\lambda / D');
% ylabel('\lambda / D');
% title('Central \lambda PSF');
% 
% figure(3);
% imagesc(ldx(:,:,length(lambda)/2),ldy(:,:,length(lambda)/2),log10(combined_exposure / max(max(combined_exposure))),[-4,0])
% axis xy; axis square;
% colorbar;
% colormap(gray(256));
% xlabel('\lambda / D');
% ylabel('\lambda / D');
% title('Full 20% bandwidth PSF');


%% Do it again on the GPU

% Set the input field
OS.planewave(single(1),length(OS.lambda_array_));

tic;
OS.GPUify;
OS.PropagateSystem1(1,2,n0);
toc
