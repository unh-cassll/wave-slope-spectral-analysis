
function [spect_out,k_out,theta_out] = get_polar_spectrum(spect_cart,kx,ky)

[s1,s2] = size(spect_cart);
num_bins = floor(min([s1 s2])/2);

k_pol = sqrt(kx.^2+ky.^2);
theta_pol = atan2(kx,ky);
theta_pol = (theta_pol - pi)*-1;
theta_pol = flipud(theta_pol);
theta_out = linspace(0,2*pi,s2);
theta_vec_triple = [theta_out theta_out+2*pi theta_out+4*pi];

kmin = median(diff(kx(1,:)));
kmax = num_bins*kmin;
k_out = linspace(kmin,kmax,num_bins)';

spect_out = NaN*zeros(num_bins,s2*3);

for n = 1:num_bins
    
    k_low = (n-1)*kmin;
    k_high = n*kmin + 0.1;
    
    inds = k_pol > k_low & k_pol <= k_high;
    
    thetavec = theta_pol(inds)';
    spectvec = spect_cart(inds)';
    
    [thetavec,unique_inds] = unique(thetavec);
    spectvec = spectvec(unique_inds);
    
    [thetavec,order] = sort(thetavec);
    spectvec = spectvec(order);
    
    thetavec_triple = [thetavec thetavec+2*pi thetavec+4*pi];
    spectvec_triple = [spectvec spectvec spectvec];
    
    spect_out(n,:) = interp1(thetavec_triple,spectvec_triple,theta_vec_triple,'linear');
    
end

spect_out(isnan(spect_out)) = 0;

spect_out = spect_out(:,s2+1:2*s2);

%k_out = repmat(k_out,1,s2);
%theta_out = repmat(theta_vec,num_bins,1);
