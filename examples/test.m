
%% Start Clean
% Avoid problems with making source code changes to objects and then trying
% to test with an old object. Just destroy them!

clear all; clc; close all;

%% Make an OptPupil Object
props = cell(6,1);
props{1} = 'System Pupil';  % string of your choosing for name of object
props{2} = 0;               % code for material
props{3} = true;            % amplitude only flag
props{4} = 0;               % z position
props{5} = 1;               % Physical Diameter [m]
props{6} = 'pup_1024.fits'; % mask [matrix or /path/filename.fits]

PUPIL = OptPupil(props);
PUPIL.describe
PUPIL.show;


%% NOTE TO USER:
% I don't want the element to know the wavelength. It doesn't make sense 
% that the element know such a thing. This means the Optical system 
% (or wavefront if it is separated from OptSys) has to provide this 
% information. This makes sense, because a wavefront knowns what 
% wavelength it is.
%
% The implication of the system providing it is right now during testing,
% it needs to be supplied by hand to ensure things work right. That is done
% doing the following:

% Define a visible bandwidth, in meters
lambda = linspace(500,600,1000) * 1e-9;

% Inject it into the material stored in the object
% This automatically runs the micron conversion for future index
% calculations
PUPIL.material_.setWavelength(lambda);

% The only time the element needs to care about the wavelength is when the
% refractive index is used to find the "phase coefficient" that will turn
% the zsag map into an optical phase upon application of the element to the
% wavefront.....
% This first object is of class OptPupil, which as stated in its
% description, is intended to be a binary amplitude mask (although the
% hooks for making it more compilcated are present). 
%..... meaning this object, which follows the current intention of the
% code, and has vacuum material, has no need to care, and nothing
% interesting is going to happen.
%
% If one did want to make the element care, it would be done doing the 
% following:

PUPIL.material_.ComputePhaseFactor;

% Plot the computed refractive index vs. wavelength to see that it did
% something. As you will see, n = 1 for all wavelengths in Vacuum, as it
% should. I told you it wouldn't be interesting.
PUPIL.material_.PlotIndexvsWavelength;

% As stated at the beginning of this Note, all of these things will be
% handled internally by the Optical System when using the full code
% functionality. This toy example being used for debugging/teaching will
% show these steps in order to ensure their functionality, and provide
% users a peek under the hood of what is happening in a use of the code
% once it is completed.