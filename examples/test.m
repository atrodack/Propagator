
%% Start Clean
% Avoid problems with making source code changes to objects and then trying
% to test with an old object. Just destroy them!

clear all; clc; close all;

%% Initialize Optical System Parameters

% Define a central wavelength, in meters
lambda_0 = 565*1e-9;

% Define a bandwidth
bandwidth = 0.1;
lambda = linspace(lambda_0 - (lambda_0*(bandwidth/2)),lambda_0 + (lambda_0*(bandwidth/2)),10);

% Define index of refraction of medium Optical System is in
n0 = 1;

%% Make the OptSys Object

% Construct
OS = OptSys();

% set the name
OS.name = 'PIAACMC Optical System';

% set the wavelength info
OS.setCentralWavelength(lambda_0);
OS.setLambdaarray(lambda);

% Set the pixel scale
OS.setPscale(0.00011);

% Set whether or not to save intermediate steps as FITS
OS.savefile = 0;

% Turn off verbose to suppress messages of adding elements in
OS.toggle_verbose('off');


%% Make an OptPupil Object

% Make the properties cell array
props = cell(6,1);
props{1} = 'System Pupil';      % string of your choosing for name of object
props{2} = 0;                   % code for material
props{3} = true;                % amplitude only flag
props{4} = 0;                   % z position
props{5} = 1;                   % Physical Diameter [m]
props{6} = 'pup_1024.fits';     % mask [matrix or /path/filename.fits]

% Construct
PUPIL = OptPupil(props);
% PUPIL.describe
% PUPIL.show;


%% Make an OptMirror Object
props = cell(6,1);
props{1} = 'Mirror 1';          % string of your choosing for name of object
props{2} = 0;                   % focal length [m]
props{3} = 0;                   % isFocal code
props{4} = 1.19999997;          % z position
props{5} = 25.4e-3;             % Physical Diameter [m]
props{6} = 'piaam0z.fits';      % zsag [matrix or /path/filename.fits]

M1 = OptMirror(props);
% M1.describe
% M1.show;

%% Make an OptMirror Object
props = cell(6,1);
props{1} = 'Mirror 2';          % string of your choosing for name of object
props{2} = 0;                   % focal length [m]
props{3} = 0;                   % isFocal code
props{4} = -1.102609;           % z position
props{5} = 25.4e-3;             % Physical Diameter [m]
props{6} = 'piaam1z.fits';      % zsag [matrix or /path/filename.fits]

M2 = OptMirror(props);
% M2.describe
% M2.show;

%% Make an OptFPM Object

% Make a crude dark disk for mask
x = linspace(-511,512,1024)*OS.pscale_;
[X,Y] = meshgrid(x);
R = sqrt(X.^2 + Y.^2);
circ = single(R<=0.00035);

props = cell(7,1);
props{1} = 'FPM';               % string of your choosing for name of object
props{2} = 2;                   % code for material
props{3} = true;                % amplitude only flag
props{4} = 1;                   % isFocal_ [determines how the propagator applies the FPM]
props{5} = -1.102609;           % z position
props{6} = 1;                   % Physical Diameter [m]
props{7} = 1-circ;              % mask [matrix or /path/filename.fits]

FPM = OptFPM(props);
% FPM.describe
% FPM.show;


%% Make an OptLyot Object
props = cell(6,1);
props{1} = 'Lyot Stop';         % string of your choosing for name of object
props{2} = 0;                   % code for material
props{3} = true;                % amplitude only flag
props{4} = -0.185609;           % z position
props{5} = 25.4e-3;             % Physical Diameter [m]
props{6} = 'LyotStop0.fits';    % mask [matrix or /path/filename.fits]

LYOT = OptLyot(props);
% LYOT.describe
% LYOT.show;

%% Make an OptDetector Object
props = cell(3,1);
props{1} = 'Detector';          % string of your choosing for name of object
props{2} = 25.4e-3;             % diameter
props{3} = ones(1024);          % zsag
DET = OptDetector(props);


%% Build the Optical System

% Make a list of the made elements
ELEMENTS = cell(6,1);
ELEMENTS{1} = PUPIL;
ELEMENTS{2} = M1;
ELEMENTS{3} = M2;
ELEMENTS{4} = FPM;
ELEMENTS{5} = LYOT;
ELEMENTS{6} = DET;

% Add Elements to System
OS.addSequentialElements(ELEMENTS);

% Turn verbose back on to plot intermediate Propagator steps
OS.toggle_verbose('on');

%% Propagate through the Optical System

% Set the input field
OS.planewave(single(1),length(lambda));
% OS.show;

% Propagate
OS.PropagateSystem1(1,6,n0);
