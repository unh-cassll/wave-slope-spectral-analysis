% Given water surface slope joint probability density function, compute
% least-squares Gram-Charlier fit and skewness/kurtosis coefficients
%
% Procedure following Cox & Munk [1954]
% 
% N. Laxague 2023
%
function out_struc = fit_GramCharlier_slope_PDF(slope_centers,P_slope_c_u,mss_u,mss_c)

[xi,zeta] = meshgrid(slope_centers/sqrt(mss_c),slope_centers/sqrt(mss_u));

% Function to fit
fit = @(b,x,y) (2*pi*sqrt(mss_u)*sqrt(mss_c))^-1*exp(-(x.^2+y.^2)/2).*...
      (1 + ...
        - 1/2*b(1)*(x.^2-1).*y...
        - 1/6*b(2)*(y.^3-3*y)...
       + 1/24*b(3)*(x.^4-6*x.^2+3)...
       + 1/24*b(4)*(y.^4-6*y.^2+3)...
        + 1/4*b(5)*(x.^2-1).*(y.^2-1));

% Least-Squares cost function
fcn = @(b) sum((fit(b,xi,zeta) - P_slope_c_u).^2,'all');

% Minimize Least-Squares
s = fminsearch(fcn,[0 0 0 0 0]);

% Send G-C fit and skewness/kurtosis coefficients to output structure
out_struc.P_fit = fit(s,xi,zeta);   % Gram-Charlier expansion fit
out_struc.c21 = s(1);               % skewness, upwind
out_struc.c03 = s(2);               % skewness, upwind
out_struc.c40 = s(3);               % kurtosis
out_struc.c04 = s(4);               % kurtosis
out_struc.c22 = s(5);               % kurtosis
