% Retrieves upwind and downwind wave frequency/wavenumber pairs
% from 3D wave spectrum
%
% Nathan Laxague 2020
% Technique borrowed from Plant & Wright (1980), JGR: Oceans
%
function out_struc = get_upwind_downwind_wave_frequencies(kw_spect,scale,fps)

% Set median filter size and wavenumber cutoffs
k_filtnum = 3;
f_filtnum = 1;
high_wavenumber_cutoff = 50;
low_wavenumber_cutoff = 0;

% Retrieve input array sizes
[s1,~,s3] = size(kw_spect);

% Compute frequency vector and create repeated array block
f_obs = (-s3:1:s3)'*fps/2/s3;
f_obs(1:s3+1) = [];
fmat = repmat(f_obs,[1 s1 s1]);
fmat = permute(fmat,[2 3 1]);

% Create wavenumber arrays
kmin = 2*pi/(s1*scale);
kmax = pi/scale;
[kx,ky] = meshgrid(-kmax:kmin:kmax,-kmax:kmin:kmax);
kx(:,1) = [];
ky(:,1) = [];
kx(end,:) = [];
ky(end,:) = [];
ky = flipud(ky);

% Trim wavenumber and frequency arrays to cutoff value
kslice = kx(1,:);
inds_trim = abs(kslice) < high_wavenumber_cutoff;
kx = kx(inds_trim,inds_trim);
ky = ky(inds_trim,inds_trim);
fmat = fmat(inds_trim,inds_trim,:);
k = sqrt(kx.^2+ky.^2);
kmat = repmat(k,[1 1 s3]);
Cp_mat = 2*pi*fmat./kmat;

% Find indices of separation between downwind/upwind halves
[s1,~,s3] = size(kmat);
dummy_vec = 1:s1;
inds_downwind = dummy_vec(ky(:,1) > 0);
inds_upwind = dummy_vec(ky(:,1) < 0);

% Create mask block to exclude wavenumbers outside of cutoffs
mask_inds = k < high_wavenumber_cutoff;
mask = ones(s1,s1)*NaN;
mask(mask_inds) = 1;
mask(k<low_wavenumber_cutoff) = NaN;
mask_block = repmat(mask,[1 1 s3]);
mask_block(Cp_mat>5 & fmat>1) = NaN;
mask_block(fmat>10) = NaN;

% Compute angle arrays from Cartesian wavenumber arrays
angle_mat = 180/pi*atan2(kx,ky);
angle_mat_down = repmat(angle_mat,[1 1 s3]);
angle_mat = flipud(angle_mat);
angle_mat_up = repmat(angle_mat,[1 1 s3]);

mask_block(abs(angle_mat_down)>45 & abs(angle_mat_up)>45) = NaN;

% Pre-allocate vectors
k_vec = kx(1,:);
k_vec = k_vec(k_vec>0);
kL = length(k_vec);
f_down_vec = NaN*k_vec;
f_up_vec = f_down_vec;
theta_down_vec = f_down_vec;
theta_up_vec = f_down_vec;
S_down_vec = f_down_vec;
S_up_vec = f_down_vec;

% Prepare spectrum for analysis
%S_filt_masked = mask_block.*medfilt3(log10(kmat.^2.*kw_spect(inds_trim,inds_trim,:)),[k_filtnum k_filtnum f_filtnum]);
%S_filt_masked = mask_block.*log10(kmat.^0.*kw_spect(inds_trim,inds_trim,:));
S_filt_masked = mask_block.*kmat.^0.*kw_spect(inds_trim,inds_trim,:);
S_downwind = S_filt_masked;
S_upwind = S_filt_masked;
S_downwind(inds_upwind,:,:) = NaN;
S_upwind(inds_downwind,:,:) = NaN;

r = k/nanmax(k_vec);
r = floor(r*kL);
rblock = repmat(r,[1 1 s3]);

for i = 1:kL
    
    try
    
    % Separate spectrum into upwind/downwind halves
    S_k_downwind = S_downwind;
    S_k_upwind = S_upwind;
    
    % Keep only values corresponding to this iteration's wavenumber
    %kdiff = abs(kmat-k_vec(i));
    %inds_trim = kdiff > kmin/100;
    inds_trim = rblock ~= i;
    S_k_downwind(inds_trim) = NaN;
    S_k_upwind(inds_trim) = NaN;
    
    % Compute spectral energy cutoff
    
    % Percentage of spectral maximum
%     S_k_downwind_max = nanmax(nanmax(nanmax(S_k_downwind)));
%     S_k_upwind_max = nanmax(nanmax(nanmax(S_k_upwind)));
%     S_k_downwind_cutoff = S_k_downwind_max + S_k_downwind_max/100;
%     S_k_upwind_cutoff = S_k_upwind_max + S_k_upwind_max/100;
    
    % Nth percentile
    percentile_val = 99.9;
    S_k_downwind_vec = sort(reshape(S_k_downwind,[],1));
    S_k_upwind_vec = sort(reshape(S_k_upwind,[],1));
    S_k_downwind_vec(isnan(S_k_downwind_vec)) = [];
    S_k_upwind_vec(isnan(S_k_upwind_vec)) = [];
    S_k_downwind_cutoff = S_k_downwind_vec(floor(percentile_val/100*length(S_k_downwind_vec)));
    S_k_upwind_cutoff = S_k_upwind_vec(floor(percentile_val/100*length(S_k_upwind_vec)));
    
    % Find indices above cutoff
    inds_downwind = S_k_downwind > S_k_downwind_cutoff;
    inds_upwind = S_k_upwind > S_k_upwind_cutoff;
    
    % Energy
    S_d = S_k_downwind(inds_downwind);
    S_u = S_k_upwind(inds_upwind);
    S_down_vec(i) = nanmean(S_d);
    S_up_vec(i) = nanmean(S_u);
    
    % Frequency
    f_d = fmat(inds_downwind);
    f_u = fmat(inds_upwind);
    f_down_vec(i) = nanmean(f_d.*S_d)/nanmean(S_d);
    f_up_vec(i) = nanmean(f_u.*S_u)/nanmean(S_u);
    
    % Direction
    t_d = angle_mat_down(inds_downwind);
    t_u = angle_mat_up(inds_upwind);
    theta_down_vec(i) = nanmean(t_d.*S_d)/nanmean(S_d);
    theta_up_vec(i) = nanmean(t_u.*S_u)/nanmean(S_u);
    
    end
    
end

% inds_clean = ~isnan(f_down_vec) & ~isnan(f_up_vec);
% 
% f_down_vec = interp1(k_vec(inds_clean),f_down_vec(inds_clean),k_vec,'linear');
% f_up_vec = interp1(k_vec(inds_clean),f_up_vec(inds_clean),k_vec,'linear');
% S_down_vec = interp1(k_vec(inds_clean),S_down_vec(inds_clean),k_vec,'linear');
% S_up_vec = interp1(k_vec(inds_clean),S_up_vec(inds_clean),k_vec,'linear');

% Output
s.k_vec = k_vec;
s.f_down = f_down_vec;
s.f_up = f_up_vec;
s.dir_down = theta_down_vec;
s.dir_up = theta_up_vec;
s.S_down = S_down_vec;
s.S_up = S_up_vec;

out_struc = s;