% Projects Cartesian spectrum S(kx,ky) onto S(k,theta) grid

function [out_k,out_d,out_S] = spect_cart_to_polar(in_spect,m_per_px,dnum)

[m,~] = size(in_spect);

kmax = pi/m_per_px;
kmin = 2*pi/(m*m_per_px);

[kx,ky] = meshgrid(-kmax:kmin:kmax,-kmax:kmin:kmax);
[x,y] = meshgrid(-m/2:m/2,-m/2:m/2);

x(:,end) = [];
x(end,:) = [];
y(:,end) = [];
y(end,:) = [];

k = sqrt(kx.^2+ky.^2);
wavenvec = k(m/2+2:m,m/2+1);

r = sqrt(x.*x+y.*y);

theta = atan2(x,y);
theta = (theta - pi)*-1;

mask = ones(m,m)*NaN;
mask(r<m/2) = 1;

r = r.*mask;
r_int = floor(r);

r_vec = reshape(r_int,m*m,1);
t_vec = reshape(theta,m*m,1);
p_vec = reshape(in_spect,m*m,1);

[r_vec,order] = sort(r_vec);
t_vec = t_vec(order);
p_vec = p_vec(order);
bad_inds = isnan(r_vec);
r_vec(bad_inds) = [];
t_vec(bad_inds) = [];
p_vec(bad_inds) = [];

knum = (m/2-1);

kmat = NaN*ones(knum,dnum);
dmat = kmat;
Smat = kmat;

dvec = linspace(0,360,dnum);
onevec = ones(1,dnum);

for j = 1:knum
    
    e = [];
    
    try
        
        inds = j == r_vec;
        
        holder_theta = t_vec(inds);
        holder_theta_d = holder_theta*180/pi;
        holder_spect = p_vec(inds);
        
        inds_nan = isnan(holder_spect);
        holder_theta_d(inds_nan) = [];
        holder_spect(inds_nan) = [];
        
        dmat(j,:) = dvec;
        kmat(j,:) = wavenvec(j)*onevec;
        Smat(j,:) = interp1(holder_theta_d,holder_spect,dvec,'pchip');
        
    catch e
        
    end
    
    if ~isempty(e)
        
        disp(e)
        
    end
    
end

out_k = kmat;
out_d = dmat;
out_S = abs(Smat);
