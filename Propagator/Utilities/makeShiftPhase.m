function [centerPupil, centerFocal] = makeShiftPhase(N,nLambdas,datatype)
% [centerPupil, centerFocal] = makeShiftPhase(N,nLambdas,datatype)

if nargin < 2
    nLambdas = 1;
    datatype = 'single';
elseif nargin < 3
    datatype = 'single';
end
M = N-1;

norm = 1/sqrt(2);
cnorm = complex(norm,norm);
arg = pi+(-2.0*pi*0.5*N/(N-3));


xcen = (0.5*(M));
ycen = (0.5*(M));
centerFocal = init_variable(N,N,nLambdas,datatype,0);
centerPupil = init_variable(N,N,nLambdas,datatype,0);
for ii = 1:N
    for jj = 1:N
        centerFocal(ii,jj) = cnorm*exp(complex(0.,arg*((ii-xcen)+(jj-ycen))));
        centerPupil(ii,jj) = cnorm*exp(complex(0., 0.5*pi - arg*((ii-xcen)+(jj-ycen))));
    end
end








end