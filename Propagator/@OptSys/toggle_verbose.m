function [ OS ] = toggle_verbose( OS, onoff )
%TOGGLE_VERBOSE 
% Toggles the verbose property of an OptSys and all its
% ELEMENTS

if ~isa(OS,'OptSys')
    error('INPUT:NOTOPTSYS','Input must be an OptSys object');
end

orig = zeros(1,OS.numElements_+1);
orig(1) = OS.verbose;
for ii = 1:length(orig)-1
    orig(ii+1)=OS.ELEMENTS_{ii}.verbose;
end



if nargin < 2
% Toggle the verbose flags to the other setting

    % OptSys verbose
    if orig(1) == 1
        % Turn off verbose
        OS.verbose = 0;
    else
        % Turn on verbose
        OS.verbose = 1;
    end
    
    %ELEMENT verbose
    for ii = 1:length(orig)-1
        if orig(ii+1) == 1
            % Turn off verbose
            OS.ELEMENTS_{ii}.verbose = 0;
        else
            % Turn on verbose
            OS.ELEMENTS_{ii}.verbose = 0;
        end
    end
elseif nargin == 2
    % Turn verbose flags to on/off
    if ~ischar(onoff)
        error('INPUT:notstring','onoff must be a string');
    end
    
    if strcmpi(onoff,'on')
        % OptSys verbose
        if orig(1) == 0
            % Turn on verbose
            OS.verbose = 1;
        end
        
        %ELEMENT verbose
        for ii = 1:length(orig)-1
            if orig(ii+1) == 0
                % Turn off verbose
                OS.ELEMENTS_{ii}.verbose = 1;
            end
        end
        
    elseif strcmpi(onoff,'off')
        % OptSys verbose
        if orig(1) == 1
            % Turn off verbose
            OS.verbose = 0;
        end
        
        %ELEMENT verbose
        for ii = 1:length(orig)-1
            if orig(ii+1) == 1
                % Turn off verbose
                OS.ELEMENTS_{ii}.verbose = 0;
            end
        end
    end
        



end

