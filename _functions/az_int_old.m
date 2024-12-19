% Azimuthally integrate slope spectrum to retrieve omnidirectional spectrum
% Additionally produces angle of maximum degree of saturation for each K
% UPDATE 9/2017: produces directional spreading function
function [out_waven,out_spect,out_theta,out_delta] = az_int_old(in_spect,m_per_px,wind_dir)

in_spect = imrotate(in_spect,wind_dir,'crop');

[m,~] = size(in_spect);

kmax = pi/m_per_px;
kmin = 2*pi/(m*m_per_px);

[kx,ky] = meshgrid(-kmax:kmin:kmax,-kmax:kmin:kmax);
kx(:,1) = [];
ky(:,1) = [];
kx(end,:) = [];
ky(end,:) = [];
ky = flipud(ky);
[x,y] = meshgrid(-m/2:m/2,-m/2:m/2);
x(:,1) = [];
y(:,1) = [];
x(end,:) = [];
y(end,:) = [];
y = flipud(y);

k = sqrt(kx.^2+ky.^2);
wavenvec = k(m/2:m,m/2);

r = sqrt(x.*x+y.*y);

theta = atan2(x,y);
theta = (theta - pi)*-1;
theta = flipud(theta);

mask = ones(m,m)*NaN;
mask(r<m/2) = 1;

r = r.*mask;
r_int = floor(r);

integrated_spect = NaN*ones(m/2,1);

r_vec = reshape(r_int,m*m,1);
t_vec = reshape(theta,m*m,1);
p_vec = reshape(in_spect,m*m,1);
k_vec = reshape(k,m*m,1);

[r_vec,order] = sort(r_vec);
t_vec = t_vec(order);
p_vec = p_vec(order);
k_vec = k_vec(order);
bad_inds = isnan(r_vec);
r_vec(bad_inds) = [];
t_vec(bad_inds) = [];
p_vec(bad_inds) = [];
k_vec(bad_inds) = [];

angle_max = NaN*integrated_spect;
delta = angle_max;

r_half = r_vec;
r_half(t_vec>pi) = [];
t_half = t_vec;
t_half(t_vec>pi) = [];
k_half = k_vec;
k_half(t_vec>pi) = [];
p_half = p_vec;
p_half(t_vec>pi) = [];

for j = 1:m/2+1
    
    e = [];
    
    try
        
        inds = j == r_vec;
        
        holder_theta = t_vec(inds);
        holder_theta_d = holder_theta*180/pi;
        holder_spect = p_vec(inds);
        
        inds_nan = isnan(holder_spect);
        holder_theta(inds_nan) = [];
        holder_theta_d(inds_nan) = [];
        holder_spect(inds_nan) = [];
        
        inds_half = j == r_half;
        theta_half = t_half(inds_half);
        B_half = k_half(inds_half).*k_half(inds_half).*p_half(inds_half);
        
        [holder_theta,order] = sort(holder_theta);
        holder_spect = holder_spect(order);
        
        %     start_ang = pi/2;
        
        %     inds = holder_theta > start_ang*pi/180 & holder_theta < (start_ang+180)*pi/180;
        %     holder_theta = holder_theta(inds);
        %     holder_spect = holder_spect(inds);
        
        %     ind_half = floor(length(holder_spect)/2)+1;
        %     holder_theta = holder_theta(ind_half:end);
        %     holder_spect = holder_spect(ind_half:end);
        
        integrated_spect(j) = trapz(holder_theta,holder_spect);
        
        ind_max = find(B_half == nanmax(B_half),1,'first');
        angle_max(j) = theta_half(ind_max);
        
        zero_inds = abs(holder_theta_d - 0) < 5;
        ninety_inds = abs(holder_theta_d - 90) < 5;
        oneeighty_inds = abs(holder_theta_d - 180) < 5;
        
        D = pi/2*(nanmedian(holder_spect(zero_inds))+nanmedian(holder_spect(oneeighty_inds))-2*nanmedian(holder_spect(ninety_inds)))./integrated_spect(j);
        
        %D = (nanmedian(holder_spect(zero_inds)) - nanmedian(holder_spect(ninety_inds)))/(nanmedian(holder_spect(zero_inds)) + nanmedian(holder_spect(ninety_inds)));
        
        delta(j) = D;
        
    catch e
        
    end
    
    if ~isempty(e)
        
        %disp(e)
        
    end
    
end

out_waven = wavenvec(2:end);
out_spect = wavenvec(2:end).*integrated_spect(2:end);
out_theta = angle_max(2:end);
out_delta = delta(2:end);
