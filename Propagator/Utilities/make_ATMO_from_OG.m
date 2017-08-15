function [Atmo_amp, Atmo_pha] = make_ATMO_from_OG(NB_TSPAN, WFTIME_STEP, TIME_SPAN,WFsize)


total_frames = (TIME_SPAN / WFTIME_STEP) * NB_TSPAN;
Atmo_amp = init_variable(WFsize,WFsize,total_frames,'single',0);
Atmo_pha = Atmo_amp;

filenames = dir;

counter_amp = 1;
counter_pha = 1;

for ii = 3:2+(2*NB_TSPAN)
    filename = filenames(ii).name;
    
    if mod(ii,2) == 0
        fprintf('Reading image %d\n',ii-3);
        tmp_amp = fitsread(filename);
        for jj = 1:size(tmp_amp,3)
            Atmo_amp(:,:,counter_amp) = tmp_amp(:,:,jj);
            counter_amp = counter_amp+1;
        end
    elseif mod(ii,2) == 1
        fprintf('Reading image %d\n',ii-3);
        tmp_pha = fitsread(filename);
        for jj = 1:size(tmp_pha,3)
            Atmo_pha(:,:,counter_pha) = tmp_pha(:,:,jj);
            counter_pha = counter_pha+1;
        end
    end
end









end