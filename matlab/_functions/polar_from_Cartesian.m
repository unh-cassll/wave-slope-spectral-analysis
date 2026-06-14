% Conversion of directional wavenumber slope spectrum from Cartesian to
% polar coordinates
%
% Original function from Juhnze (Jerry) Liu [2023]
% Annotation and variance conservation by Nathan Laxague [2023]
%

function [Sktheta_polar,k_vec,theta_vec] = polar_from_Cartesian(Skxy,kx,ky,dtheta,camera_heading_deg)

% Suppress warning involved in interpolation over NaNs
warning('off','MATLAB:interp1:NaNstrip')

% Preallocate wavenumber array
[Nx, Ny] = size(Skxy);
num_bins = floor(min([Nx Ny])/2);
kmin = mean(abs(diff(kx(1,:))));
kmax = num_bins*kmin;
kmat_xy = sqrt(kx.^2 + ky.^2);
k_vec = linspace(kmin,kmax,num_bins)';

% Preallocate direction array
theta_mat_xy = atan2(kx,ky);
theta_mat_xy = (theta_mat_xy - pi)*-1;
theta_mat_xy = flipud(theta_mat_xy);
theta_mat_xy = pi/180*mod(theta_mat_xy*180/pi+camera_heading_deg,360);
theta_vec = (0:dtheta*180/pi:359)*pi/180;
theta_mat = repmat(theta_vec, [num_bins, 1]);

% Loop over all wavenumber bins
Sktheta_polar = NaN(size(theta_mat));
for n = 1:num_bins

    % Grab indices within a particular wavenumber range
    k_val = n*kmin;
    k_low = k_val - kmin*0.5;
    k_high = k_val + kmin*0.5;
    inds = find( kmat_xy > k_low & kmat_xy <= k_high);
    inds(isinf(Skxy(inds))) = [];

    % Grab angle and spectral energy density values at those indices
    thetavec = theta_mat_xy(inds)';
    Svec = Skxy(inds)';
    [thetavec, order] = sort(thetavec);
    Svec = Svec(order);
    Svar = sum(Svec)*kmin*kmin;
    [thetavec, inds_uni] = unique(thetavec);
    Svec = Svec(inds_uni);

    % Interpolate to full theta vector
    [thetainterp_uni, inds_uni] = unique([thetavec thetavec+2*pi thetavec+4*pi]);
    Svec_uni0 = [Svec Svec Svec];
    Svec_uni = Svec_uni0(inds_uni);
    Svec_interp = interp1( thetainterp_uni, Svec_uni, theta_vec+2*pi, 'pchip');

    % Conserve variance in each wavenumber bin
    C = Svar/(k_val*sum(Svec_interp)*kmin*dtheta);
    Sktheta_polar(n,:) = C*Svec_interp;

end

% Enforce variance conservation
Svar_Cartesian = sum(Skxy,'all','omitnan')*kmin*kmin;
Svar_polar = sum(k_vec.*Sktheta_polar,'all','omitnan')*kmin*dtheta;
Sktheta_polar = Sktheta_polar*Svar_Cartesian/Svar_polar;
