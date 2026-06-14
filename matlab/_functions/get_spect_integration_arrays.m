%
function out_struc = get_spect_integration_arrays(n_k,m_per_px)

% Compute wavenumber limits
kmin = 2*pi/(n_k*m_per_px);
kmax = pi/m_per_px;

[kx,ky] = meshgrid(-kmax:kmin:kmax,-kmax:kmin:kmax);
kx(:,1) = [];
ky(:,1) = [];
kx(end,:) = [];
ky(end,:) = [];
ky = flipud(ky);
[x,y] = meshgrid(-n_k/2:n_k/2,-n_k/2:n_k/2);
x(:,1) = [];
y(:,1) = [];
x(end,:) = [];
y(end,:) = [];
y = flipud(y);

k = sqrt(kx.^2+ky.^2);
wavenvec = k(n_k/2+1:n_k,n_k/2+1);

r = sqrt(x.*x+y.*y);

theta = atan2(x,y);
theta = (theta - pi)*-1;

mask = ones(n_k,n_k)*NaN;
mask(r<n_k/2) = 1;

r = r.*mask;
r_int = floor(r);

r_vec = reshape(r_int,n_k*n_k,1);
t_vec = reshape(theta,n_k*n_k,1);
p_vec = reshape(in_spect,n_k*n_k,1);
k_vec = reshape(k,n_k*n_k,1);

[r_vec,order] = sort(r_vec);
t_vec = t_vec(order);
p_vec = p_vec(order);
k_vec = k_vec(order);
bad_inds = isnan(r_vec);
r_vec(bad_inds) = [];
t_vec(bad_inds) = [];
p_vec(bad_inds) = [];
k_vec(bad_inds) = [];