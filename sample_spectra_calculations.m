% Spectral Analysis applied to georectified sample wave image stack

addpath _data/
addpath _functions/

load _data/angle_stack_x.mat
load _data/angle_stack_y.mat
load _data/freeboard.mat

slopefield_stack_x = tand(double(ax_int)*1e-3);
slopefield_stack_y = tand(double(ay_int)*1e-3);

fps = 15;
m_per_px_mean = 2*0.01074; % found from the Motion-Correction-Georeferencing Library

alt_mean = mean(freeboard);

mss_x = var(slopefield_stack_x,[],'all','omitnan');
mss_y = var(slopefield_stack_y,[],'all','omitnan');
mss = mss_x + mss_y;

slopefield_stack_x(isnan(slopefield_stack_x)) = 0;
slopefield_stack_y(isnan(slopefield_stack_y)) = 0;

Sx = double(slopefield_stack_x(:,:,1));
Sy = double(slopefield_stack_y(:,:,1));

% Select a region
[s1,s2] = size(Sx);
N_spect = 2^floor(log(min([s1 s2]))/log(2));
N_spect_2 = N_spect/2;


center_rows = floor(s1/2)-N_spect_2:floor(s1/2)+N_spect_2-1;
center_cols = floor(s2/2)-N_spect_2:floor(s2/2)+N_spect_2-1;

sx = slopefield_stack_x(center_rows,center_cols,:);
sy = slopefield_stack_y(center_rows,center_cols,:);

% COMPUTE SPECTRA FROM 3D STACK

w = circular_tukey(1+0*sx,0.2);
[dirspect,~] = compute_slope_spectrum(w.*sx,w.*sy,m_per_px_mean,fps,N_spect);
Skf = dirspect.Skf;
mean_heading =  0; % absolute heading is already considered in the rectification process, so the input is zero.
nframes = size(slopefield_stack_y,3);
df = fps/nframes; %fs/Nt
dk = 2*pi/(m_per_px_mean*length(center_rows)); % 2*pi/L_max
spectra_struc = compute_sub_spectra(Skf,dk,df,mean_heading,6);


% Figures
plot_sub_spectra(spectra_struc.f,spectra_struc.k,spectra_struc.nu,spectra_struc.theta,spectra_struc,1)

% Compare F_f and F_f obtained by F_k 
k = spectra_struc.k;
f = spectra_struc.f;
F_k = spectra_struc.F_k;
F_f = spectra_struc.F_f;
f_disp = sqrt(9.81*k)/(2*pi);
c = 2*pi*f_disp./k;
cg = c/2;
F_f_disp = F_k.*k./cg/(2*pi);
figure(101);loglog(f,F_f,'-',f_disp,F_f_disp,'-','linewidth',2)


% Plot k-f Spectrum both for Cross-Slope and Along-Slope
Skf_center = N_spect_2;

kvec = 0:dk:dk*2*Skf_center;
kvec = kvec - dk*Skf_center;
figure(998);imagesc(kvec,f,log10(squeeze(Skf(:,Skf_center,:)))');colorbar;view([0 -90]);xlim([-1 1]*100);hold on;plot(k,sqrt(9.81*k+0.072/1000*k.^3)/(2*pi),'-r',-k,sqrt(9.81*k+0.072/1000*k.^3)/(2*pi),'-r','linewidth',2);hold off
figure(999);imagesc(kvec,f,log10(squeeze(Skf(Skf_center,:,:)))');colorbar;view([0 -90]);xlim([-1 1]*100);hold on;plot(k,sqrt(9.81*k+0.072/1000*k.^3)/(2*pi),'-r',-k,sqrt(9.81*k+0.072/1000*k.^3)/(2*pi),'-r','linewidth',2);hold off

