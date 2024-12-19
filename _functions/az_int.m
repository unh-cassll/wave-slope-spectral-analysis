% Azimuthally integrate slope spectrum to retrieve omnidirectional spectrum
%
% Nathan Laxague, 2014-2020
%
function [out_k,out_S] = az_int(in_S,m_per_px)

[m,~] = size(in_S);

kmax = pi/m_per_px;
kmin = 2*pi/(m*m_per_px);

[kx,ky] = meshgrid(-kmax:kmin:kmax,-kmax:kmin:kmax);
kx(:,end) = [];
ky(:,end) = [];
kx(1,:) = [];
ky(1,:) = [];
ky = flipud(ky);
[x,y] = meshgrid(-m/2:m/2,-m/2:m/2);
x(:,end) = [];
y(:,end) = [];
x(1,:) = [];
y(1,:) = [];
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
p_vec = reshape(in_S,m*m,1);

[r_vec,order] = sort(r_vec);
t_vec = t_vec(order);
p_vec = p_vec(order);
bad_inds = isnan(r_vec);
r_vec(bad_inds) = [];
t_vec(bad_inds) = [];
p_vec(bad_inds) = [];

for j = 1:m/2+1
        
        inds = j == r_vec;
        
        holder_theta = t_vec(inds);
        holder_spect = p_vec(inds);
        
        inds_nan = isnan(holder_spect);
        holder_theta(inds_nan) = [];
        holder_spect(inds_nan) = [];
        
        [holder_theta,order] = sort(holder_theta);
        holder_spect = holder_spect(order);
                
        integrated_spect(j) = 2*pi*trapz(holder_theta,holder_spect);
            
end

out_k = wavenvec(2:end);
out_S = wavenvec(2:end).*integrated_spect(2:end);

