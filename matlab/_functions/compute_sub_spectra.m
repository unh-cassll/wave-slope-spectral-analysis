%
% Given wavenumber-frequency directional slope spectrum Skf, computes
%
% * wavenumber directional spectrum
% ** slope S(k,theta) & elevation F(k,theta)
% * frequency directional spectrum
% ** slope F(f,theta) & elevation F(f,theta)
% * inverse phase speed directional spectrum
% ** slope Q_s(nu,theta) & elevation Q_eta(nu,theta)
% * phase speed difference about +/-180 deg.
%
% Function by Nathan Laxague (2023)
%
% Polar coordinate transformation from Jerry Liu
% Inverse phase speed spectrum concept from Jan-Victor Bjorkqvist
% Phase speed difference concept from Plant & Wright [1980]
%
function spectra_struc = compute_sub_spectra(Skf,dk,df,camera_heading_deg,parallel_size)

% Produce frequency and wavenumber arrays
[Nx,Ny,Nt] = size(Skf);
m = mean([Nx Ny]);
Kx = (dk) * (-ceil((Nx-1)/2): floor((Nx-1)/2));
kx = repmat(Kx, [Ny,1]);
Ky = (dk) * (ceil((Ny-1)/2):-1: -floor((Ny-1)/2));
ky = repmat(Ky', [1, Nx]);
f = df:df:Nt*df;

% Convert from Cartesian to polar coordinates
k = (dk:dk:m/2*dk)';
dtheta = 5*pi/180;
theta = (0:dtheta*180/pi:359)*pi/180;
Skftheta = NaN*ones(m/2,length(theta),length(f));

if parallel_size == 0
    
    for fnum = 1:length(f)
        
        Skxy = Skf(:,:,fnum);
        
        [Sktheta_polar,~,~] = polar_from_Cartesian(Skxy,kx,ky,dtheta,camera_heading_deg);
        
        Skftheta(:,:,fnum) = Sktheta_polar;
        
    end
    
else
    
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
        parpool(parallel_size)
%     else
%         disp('USING EXISTING PARALLEL POOL...')
    end
    
    parfor fnum = 1:length(f)
        
        Skxy = Skf(:,:,fnum);
        
        [Sktheta_polar,~,~] = polar_from_Cartesian(Skxy,kx,ky,dtheta,camera_heading_deg);
        
        Skftheta(:,:,fnum) = Sktheta_polar;
        
    end
    
end

% Enforce variance conservation after conversion to polar coordinates
Skftheta_norm = Skftheta*sum(Skf,'all','omitnan')*dk^2*df/(sum(k.*Skftheta,'all','omitnan')*dk*dtheta*df);
Fkftheta_norm = k.^-2.*Skftheta_norm;
Bkftheta_norm = k.^2.*Skftheta_norm;

clear Skf

% Trim off portion of spectra for which frequency spectrum tends towards
% noise-afflicted range (or does not resolve variation)
f_low = 0; %0.5;
f_high = 15; % 10;
inds_f = find(f>f_low,1,'first'):find(f<f_high,1,'last');

% Compute location of dispersion shell in (f,theta) space
kp = NaN*ones(size(Skftheta_norm,2),size(Skftheta_norm,3));
Sp = kp;
Skftheta_filt = medfilt3(Skftheta_norm,[3 1 1]);
for i = 1:length(f)
    for j = 1:length(theta)
    Skbit = Skftheta_filt(:,j,i);
    Smax = max(Skbit,[],'all','omitnan');
    ind_p = find(Skbit==Smax,1,'first');
    kp(j,i) = k(ind_p);
    Sp(j,i) = Smax*k(ind_p)*dk;       % slope^2/df/dtheta
    end
end

% Convert dispersion shell data to (k,theta) space
f_k_theta = NaN*ones(length(k),length(theta));
Smax_k_theta = f_k_theta;
for i = 2:length(k)
    for j = 1:length(theta)
        kp_slice = kp(j,:);
        Sp_slice = Sp(j,:);
        fbit = f;
        fbit(abs(kp_slice-k(i))>=dk|f>f_high) = NaN;
        f_k_theta(i,j) = mean(fbit,'all','omitnan');
        mask = fbit;
        mask(~isnan(mask)) = 1;
        Smax_k_theta(i,j) = sum(mask.*Sp_slice,'all','omitnan')*df/(k(i)*dk);    % slope^2/(k*dk*dtheta)
    end
end

% Compute the difference in phase speed between a wave and its 180
% degree-flipped pair
f_k_theta_opposed = [f_k_theta(:,length(theta)/2+1:end) f_k_theta(:,1:length(theta)/2)];
delta_c = 2*pi*(f_k_theta-f_k_theta_opposed)./k;

% Use the fraction of good phase speed difference to find the cutoff k
frac_good = sum(~isnan(delta_c),2)/length(theta);
inds_k = 1:find(frac_good>0.1,1,'last')+1;

% Integrate with respect to frequency
S_k_theta = squeeze(sum(Skftheta_norm,3,'omitnan'))*df;
F_k_theta = squeeze(sum(Fkftheta_norm,3,'omitnan'))*df;
B_k_theta = squeeze(sum(Bkftheta_norm,3,'omitnan'))*df;

% Integrate with respect to wavenumber
S_f_theta = squeeze(sum(k.*Skftheta_norm,1,'omitnan'))*dk;
F_f_theta = squeeze(sum(k.*Fkftheta_norm,1,'omitnan'))*dk;
B_f_theta = squeeze(sum(k.*Bkftheta_norm,1,'omitnan'))*dk;

% Integrate with respect to direction
S_k = sum(k.*S_k_theta,2)*dtheta;
F_k = sum(k.*F_k_theta,2)*dtheta;
B_k = sum(k.*B_k_theta,2)*dtheta;
S_f = sum(S_f_theta)*dtheta;
F_f = sum(F_f_theta)*dtheta;
B_f = sum(B_f_theta)*dtheta;

% Compute slope and elevation variances
slope_var = sum(k.*S_k_theta,'all','omitnan')*dk*dtheta;
elev_var = sum(k.*F_k_theta,'all','omitnan')*dk*dtheta;

% Make cuts
Skf_theta_norm_cut = Skftheta_norm(inds_k,:,inds_f);
Fkf_theta_norm_cut = Fkftheta_norm(inds_k,:,inds_f);
Bkf_theta_norm_cut = Bkftheta_norm(inds_k,:,inds_f);

% Define inverse phase speed vector and 2D array
dnu_base = 0.01;
nu = [0.05:dnu_base:0.5 0.5+dnu_base*2:dnu_base*2:1 1+dnu_base*4:dnu_base*4:2 2+dnu_base*8:dnu_base*8:4];
dnu_vec = [dnu_base diff(nu)];
k_mat = repmat(k(inds_k),[1 length(f(inds_f))]);
f_mat = repmat(f(inds_f),[length(k(inds_k)) 1]);
nu_mat = k_mat./(2*pi*f_mat);

% Calculate inverse phase speed slope and elevation spectra
Qs_nu_theta = NaN*ones(length(theta),length(nu));
Qeta_nu_theta = Qs_nu_theta;
Qb_nu_theta = Qs_nu_theta;
% f_s_weighted = Qs_nu_theta;
% k_s_weighted = Qs_nu_theta;
% f_eta_weighted = Qs_nu_theta;
% k_eta_weighted = Qs_nu_theta;
for theta_num = 1:length(theta)

    % Compute spectral variance density per direction
    Skwtheta_slice = k_mat.*squeeze(Skf_theta_norm_cut(:,theta_num,:))*dk*df; % slope^2/rad
    Fkwtheta_slice = k_mat.*squeeze(Fkf_theta_norm_cut(:,theta_num,:))*dk*df; % m^2/rad
    Bkwtheta_slice = k_mat.*squeeze(Bkf_theta_norm_cut(:,theta_num,:))*dk*df; % m^-2/rad

    for nu_num = 1:length(nu)

        mask_nu = NaN*Skwtheta_slice;
        mask_nu(nu_mat>=nu(nu_num)-dnu_vec(nu_num)/2 & nu_mat<nu(nu_num)+dnu_vec(nu_num)/2) = 1;

        Qs_nu_theta(theta_num,nu_num) = sum(mask_nu.*Skwtheta_slice,'all','omitnan')/(nu(nu_num)*dnu_vec(nu_num));   % slope^2/rad/(s m^-1)^2   = m^2/s^2/rad
        Qeta_nu_theta(theta_num,nu_num) = sum(mask_nu.*Fkwtheta_slice,'all','omitnan')/(nu(nu_num)*dnu_vec(nu_num)); % m^2/rad/(s m^-1)^2       = m^4/s^2/rad                                
        Qb_nu_theta(theta_num,nu_num) = sum(mask_nu.*Bkwtheta_slice,'all','omitnan')/(nu(nu_num)*dnu_vec(nu_num));   % m^-2/rad/(s m^-1)^2      = s^-2/rad                                
        
        % f_s_weighted(theta_num,nu_num) = sum(mask_nu.*f_mat.*Skwtheta_slice,'all','omitnan')/sum(mask_nu.*Skwtheta_slice,'all','omitnan');
        % k_s_weighted(theta_num,nu_num) = sum(mask_nu.*k_mat.*Skwtheta_slice,'all','omitnan')/sum(mask_nu.*Skwtheta_slice,'all','omitnan');
        % f_eta_weighted(theta_num,nu_num) = sum(mask_nu.*f_mat.*Fkwtheta_slice,'all','omitnan')/sum(mask_nu.*Skwtheta_slice,'all','omitnan');
        % k_eta_weighted(theta_num,nu_num) = sum(mask_nu.*k_mat.*Fkwtheta_slice,'all','omitnan')/sum(mask_nu.*Skwtheta_slice,'all','omitnan');

    end

end

% Integrate with respect to direction
Qs = sum(nu.*Qs_nu_theta)*dtheta;       % s/m
Qeta = sum(nu.*Qeta_nu_theta)*dtheta;   % m^3/s
Qb = sum(nu.*Qb_nu_theta)*dtheta;       % 1/(s*m)

% Output
spectra_struc.k_cutoff = max(k(inds_k));
spectra_struc.total_slope_variance = slope_var;
spectra_struc.total_elevation_variance = elev_var;
% independent variables
spectra_struc.k = k(:);
spectra_struc.f = f(:);
spectra_struc.nu = nu(:);
spectra_struc.theta = theta(:)';
% wavenumber spectra
spectra_struc.S_k = S_k;
spectra_struc.S_k_theta = S_k_theta;
spectra_struc.F_k = F_k;
spectra_struc.F_k_theta = F_k_theta;
spectra_struc.B_k = B_k;
spectra_struc.B_k_theta = B_k_theta;
% frequency spectra
spectra_struc.S_f = S_f';
spectra_struc.S_f_theta = S_f_theta';
spectra_struc.F_f = F_f';
spectra_struc.F_f_theta = F_f_theta';
spectra_struc.B_f = B_f';
spectra_struc.B_f_theta = B_f_theta';
% inverse phase speed spectra
spectra_struc.Qs = Qs';
spectra_struc.Qs_nu_theta = Qs_nu_theta';
spectra_struc.Qeta = Qeta';
spectra_struc.Qeta_nu_theta = Qeta_nu_theta';
spectra_struc.Qb = Qb';
spectra_struc.Qb_nu_theta = Qb_nu_theta';
% phase speed difference
spectra_struc.f_down_k_theta = f_k_theta;
spectra_struc.f_up_k_theta = f_k_theta_opposed;
spectra_struc.delta_c_k_theta = delta_c;
spectra_struc.Sp_k_theta = Smax_k_theta;
spectra_struc.Sp_f_theta = Sp';
spectra_struc.kp_f_theta = kp';
