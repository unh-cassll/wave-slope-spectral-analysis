% Given 3D spectrum with dimensions kx,ky,omega, plots wavenumber-frequency
% slice in direction 'angle' (clockwise from y-axis)
function plot_wavenumber_frequency_spectrum(S,dx,dt,rot_angle)

% Check input array size
[s1,~,s3] = size(S);
mid_ind = s3/2+1;

f = (-1/2:1/s3:1/2)'/dt;
f(1:mid_ind) = [];

% Compute wavenumber limits
kmin = 2*pi/(s1*dx);
kmax = pi/dx;

% Compute wavenumber arrays
[kx,ky] = meshgrid(-kmax:kmin:kmax,-kmax:kmin:kmax);
kx(:,1) = [];
ky(:,1) = [];
kx(end,:) = [];
ky(end,:) = [];
ky = flipud(ky);

fmat = repmat(f',s1,1);
kmat = repmat(ky(:,1),1,s3/2);

S_rot = imrotate3(S,-rot_angle,[0 0 1]);
spect_slice = squeeze(S_rot(:,mid_ind,mid_ind:end));

figure(1)
pcolor(kmat,fmat,log10(abs(kmat.^2).*spect_slice))
shading('flat')
colorbar
xlabel('k [rad m^{-1}]')
ylabel('f [Hz]')