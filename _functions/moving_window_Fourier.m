% Given slopefield stacks sX and sY, computes moving-window wave dCP, the
% difference in phase speed between upwind and downwind-propagating waves
%
% Nathan Laxague 2020
%
function [out_t,out_k,out_f_down,out_f_up,out_dir_down,out_dir_up,out_omni_k,out_omni_Sk] = moving_window_Fourier(Sx,Sy,m_per_px,fps,window_length,window_res,num_substacks,ccw_rot_angle_degrees,run_name)

taper_width = 0.1;

w = circular_tukey(0*Sx+1,taper_width);
Sx = w.*Sx;
Sy = w.*Sy;

[s1,~,s3] = size(Sx);

frames_per_bit = fps*window_length;
frames_per_shift = floor(fps*window_res);
if frames_per_shift == 0
    frames_per_shift = 1;
end

subframes = s3 - frames_per_bit;

spect_struc = struct();

start_indices = 1:frames_per_shift:subframes;

spect_num = 1;
istart = start_indices(spect_num);
iend = start_indices(spect_num) + frames_per_bit - 1;

xbit = Sx(:,:,istart:iend);
ybit = Sy(:,:,istart:iend);

[dirspect,~] = compute_slope_spectrum(xbit,ybit,m_per_px,fps,s1);
Skw = dirspect.Skw;
Skw_rot = imrotate3(Skw,ccw_rot_angle_degrees,[0 0 1],'crop');

dopp_shift_struc = get_upwind_downwind_wave_frequencies(Skw_rot,m_per_px,fps);

nanvec = NaN*dopp_shift_struc.k_vec';

frames_per_stack = floor(s3/num_substacks);
start_indices_stack = 1:frames_per_shift:frames_per_stack-frames_per_bit;
num_spects_per_stack = length(start_indices_stack);

empty_checker = ones(s3,1);

for stack_num = 1:num_substacks
    
    stack_start = (stack_num-1)*frames_per_stack+1;
    stack_end = stack_num*frames_per_stack;
    
    Sx_stack = Sx(:,:,stack_start:stack_end);
    Sy_stack = Sy(:,:,stack_start:stack_end);
    
    stack_offset = stack_start - 1;
    
    parfor spect_num = 1:num_spects_per_stack
        
        e = [];
        
        try
            
            istart = start_indices(spect_num);
            iend = start_indices(spect_num) + frames_per_bit - 1;
            
            xbit = Sx_stack(:,:,istart:iend);
            ybit = Sy_stack(:,:,istart:iend);
            
            [dirspect,omnispect] = compute_slope_spectrum(xbit,ybit,m_per_px,fps,s1);
            Skw = dirspect.Skw;
            Skw_rot = imrotate3(Skw,ccw_rot_angle_degrees,[0 0 1],'crop');
            
            dopp_shift_struc = get_upwind_downwind_wave_frequencies(Skw_rot,m_per_px,fps);
            
            k_vec = dopp_shift_struc.k_vec';
            f_down = dopp_shift_struc.f_down';
            f_up = dopp_shift_struc.f_up';
            dir_down = mod(dopp_shift_struc.dir_down+ccw_rot_angle_degrees,360)';
            dir_up = mod(dopp_shift_struc.dir_up+ccw_rot_angle_degrees+180,360)';
            
            spect_struc(spect_num+stack_offset).t = (start_indices(spect_num)+stack_offset)/fps;
            spect_struc(spect_num+stack_offset).k_vec = k_vec;
            spect_struc(spect_num+stack_offset).f_down = f_down;
            spect_struc(spect_num+stack_offset).f_up = f_up;
            spect_struc(spect_num+stack_offset).dir_down = dir_down;
            spect_struc(spect_num+stack_offset).dir_up = dir_up;
            spect_struc(spect_num+stack_offset).omni_k = omnispect.k;
            spect_struc(spect_num+stack_offset).omni_Sk = omnispect.S;
            empty_checker(spect_num+stack_offset) = 0;
            
        catch e
            
            if ~isempty(e)
                
                disp(e)
                spect_struc(spect_num+stack_offset).t = NaN;
                spect_struc(spect_num+stack_offset).k_vec = nanvec;
                spect_struc(spect_num+stack_offset).f_down = nanvec;
                spect_struc(spect_num+stack_offset).f_up = nanvec;
                spect_struc(spect_num+stack_offset).dir_down = nanvec;
                spect_struc(spect_num+stack_offset).dir_up = nanvec;
                spect_struc(spect_num+stack_offset).omni_k = nanvec;
                spect_struc(spect_num+stack_offset).omni_Sk = nanvec;
                empty_checker(spect_num+stack_offset) = 0;
                
            end
            
        end
        
    end
    
    dstr = datestr(now,'mm/dd/yyyy HH:MM:SS');
    disp([dstr '... DONE WITH PART #' num2str(stack_num) '/' num2str(num_substacks) ', ' run_name])
    
end

empty_checker = logical(empty_checker);

inds_remaining = 1:s3;
inds_remaining = inds_remaining(empty_checker);
inds_remaining(inds_remaining>subframes) = [];

for i = 1:length(inds_remaining)
    
    e = [];
    
    try
        
        ind = inds_remaining(i);
        
        istart = ind;
        iend = ind + frames_per_bit - 1;
        
        xbit = Sx(:,:,istart:iend);
        ybit = Sy(:,:,istart:iend);
        
        [dirspect,omnispect] = compute_slope_spectrum(xbit,ybit,m_per_px,fps,s1);
        Skw = dirspect.Skw;
        Skw_rot = imrotate3(Skw,ccw_rot_angle_degrees,[0 0 1],'crop');
        
        dopp_shift_struc = get_upwind_downwind_wave_frequencies(Skw_rot,m_per_px,fps);
        
        k_vec = dopp_shift_struc.k_vec';
        f_down = dopp_shift_struc.f_down';
        f_up = dopp_shift_struc.f_up';
        dir_down = mod(dopp_shift_struc.dir_down+ccw_rot_angle_degrees,360)';
        dir_up = mod(dopp_shift_struc.dir_up+ccw_rot_angle_degrees+180,360)';
        
        spect_struc(ind).t = (ind)/fps;
        spect_struc(ind).k_vec = k_vec;
        spect_struc(ind).f_down = f_down;
        spect_struc(ind).f_up = f_up;
        spect_struc(ind).dir_down = dir_down;
        spect_struc(ind).dir_up = dir_up;
        spect_struc(ind).omni_k = omnispect.k;
        spect_struc(ind).omni_Sk = omnispect.S;
        
    catch e
        
        if ~isempty(e)
            
            disp(e)
            spect_struc(ind).t = NaN;
            spect_struc(ind).k_vec = nanvec;
            spect_struc(ind).f_down = nanvec;
            spect_struc(ind).f_up = nanvec;
            spect_struc(ind).dir_down = nanvec;
            spect_struc(ind).dir_up = nanvec;
            spect_struc(ind).omni_k = nanvec;
            spect_struc(ind).omni_Sk = nanvec;
            
        end
        
    end
    
end

dstr = datestr(now,'mm/dd/yyyy HH:MM:SS');
disp([dstr '... DONE WITH REMNANTS,  ' run_name])

out_t = [spect_struc.t];
out_k = [spect_struc.k_vec];
out_f_down = [spect_struc.f_down];
out_f_up = [spect_struc.f_up];
out_dir_down = [spect_struc.dir_down];
out_dir_up = [spect_struc.dir_up];
out_omni_k = [spect_struc.omni_k];
out_omni_Sk = [spect_struc.omni_Sk];