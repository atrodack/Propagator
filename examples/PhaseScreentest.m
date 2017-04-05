%% Start Clean
% Avoid problems with making source code changes to objects and then trying
% to test with an old object. Just destroy them!

clear all; clc; close all;


%% Set Coordinate parameters

N = 1024;
dx = 0.00011;
dk = (2*pi) ./ (N*dx);
kx = ((1:N) - (N/2))*dk;
[KX,KY] = meshgrid(kx);
KR = sqrt(KX.^2 + KY.^2);

%% Set up Screen Properties

PROPERTIES = cell(10,1);
PROPERTIES{1,1} = 'PhaseScreen';
PROPERTIES{2,1} = 0;
PROPERTIES{3,1} = 100;
PROPERTIES{4,1} = 0.5e-6;
PROPERTIES{5,1} = [];
PROPERTIES{6,1} = 5e-2;
PROPERTIES{7,1} = 2;
PROPERTIES{8,1} = KR;
PROPERTIES{9,1} = dx;
PROPERTIES{10,1} = N;

PS = OptPhaseScreen(PROPERTIES);


%% Make the Screen

% PS.set_seed(0);
PS.set_fixLF(false);

for ii = 1:10
    PS.makeScreen();
    imagesc(PS.screen_);
    colorbar;
    axis square;
    drawnow;
    pause(0.15);
end