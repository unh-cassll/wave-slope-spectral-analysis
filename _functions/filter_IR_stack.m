% Lowpass and ocean wave dispersive bandstop filter for IR data
% Nathan J. M. Laxague, 3/2018
% Originally based on an unnamed script created by Arete Associates
function filt_stack = filter_IR_stack(raw_stack,dx,dt,param1,param2)

g = 9.81; %m/s^2

% Measure input image stack dimensions
[Ny,Nx,Nt] = size(raw_stack);

% Produce frequency and wavenumber vectors
fmax = 1/dt;
f = (-1/2:1/Nt:1/2)*fmax;
df = median(diff(f));
f(Nt/2+1) = [];
kmin_x = 2*pi/(Nx*dx);
kmin_y = 2*pi/(Ny*dx);
kmax = pi/dx;
kx = -kmax:kmin_x:kmax;
ky = -kmax:kmin_y:kmax;
kx(1) = [];
ky(1) = [];

% Produce frequency and wavenumber 3D arrays
[Kx,Ky,F] = meshgrid(kx,ky,f);
K = sqrt(Kx.^2+Ky.^2);
maskblock = 0*raw_stack + 1;
%Fa = abs(F);

% Compute 3D FFT of input image stack
Mf = fftshift(fftn(raw_stack,[Ny Nx Nt]));

% Define wave dispersive bandstop mask
f_low = (1/(2*pi))*sqrt(g*K)-(param1*Kx+param2*Ky) - 0.1;
f_high = (1/(2*pi))*sqrt(g*K)-(param1*Kx+param2*Ky) + 0.1;
maskblock(F>f_low-df/2 & F<f_high+df/2) = 0;

% Create the lowpass filter
%Kmax = max(K(:,1,1));
%Kcut = Kmax/10;
%maskblock(K>Kcut) = 0;

% Ensure that energy remains the same post-filtration
coeff_scale = squeeze(Nx*Ny*Nt/nansum(nansum(nansum(maskblock))));

% Export the filtered image stack
filt_stack = real(ifftn(fftshift(coeff_scale*maskblock.*Mf)));

% Old way
%
% % Compute 3D FFT of input image stack
% Mf = fftn(raw_stack,[Ny Nx Nt]);
% 
% % Define wave dispersive bandstop mask
% f_low = (1/(2*pi))*sqrt(g*K)-(param1*Kx+param2*Ky) - 1.0;
% f_high = (1/(2*pi))*sqrt(g*K)-(param1*Kx+param2*Ky) + 1.0;
% bandstop_mask = ~((Fa>f_low-df/2) & (Fa<f_high+df/2));
% 
% % Combine that mask with the positive frequency condition
% Tmask = ~bandstop_mask & Fa>0;
% 
% % Create the lowpass filter
% Kmax = max(K(:,1,1));
% Kcut = Kmax/100;
% Kpass = ones(size(K));
% Kfilt = 1 - ((K-Kcut)/(Kmax-Kcut)).^4;
% Kfilt(Kfilt<0) = 0;
% Kpass(K>Kcut) = Kfilt(K>Kcut);
% coeff_pass = squeeze(Nx*Ny*Nt/nansum(nansum(nansum(Kpass))));
% 
% % Export the filtered image stack
% filt_stack = real(ifftn(coeff_pass*Mf.*(Tmask).*Kpass,[Ny Nx Nt]));

