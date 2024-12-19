% Given slopefield stacks Sx and Sy, scale in meters/pixel, and fps,
% provides directional wavenumber-frequency spectrum Skw, directional
% wavenumber spectrum Skk, omnidirectional spectrum omnispect, and the
% wavenumber and frequency arrays kx, ky, and f
%
% Nathan Laxague, built over the years 2014 - 2023
%
function [dirspect,omnispect] = compute_slope_spectrum(Sx,Sy,scale,fps,framesize)

% Check input array size
[s1,s2,s3] = size(Sx);

% Compute wavenumber limits
kmin = 2*pi/(framesize*scale);
kmax = pi/scale;

% Compute wavenumber arrays
m = mean([s1 s2]);

Kx = (kmin) * (-ceil((m-1)/2): floor((m-1)/2));
kx = repmat(Kx, [m,1]);
Ky = (kmin) * (ceil((m-1)/2):-1: -floor((m-1)/2));
ky = repmat(Ky', [1, m]);

% Create mask to exclude wavenumbers beyond Nyquist
k = sqrt(kx.^2+ky.^2);
mask = 0*k + 1;
mask(k>kmax) = NaN;

% Compute frequency vector
f = (-fps/2:fps/s3:fps/2);
f = f(s3/2+2:end);

% Compute discrete wavenumber/frequency elements
dk = kmin;
df = fps/s3;

% Compute FFT
Ax = fftshift(fftn(Sx,[framesize framesize s3]));
Ay = fftshift(fftn(Sy,[framesize framesize s3]));
clear Sx Sy

if s3 >1

    % Grab positive frequencies from spectrum
    Ax = flip(Ax(:,:,1:s3/2),3);
    Ay = flip(Ay(:,:,1:s3/2),3);
    
    % Compute scaling coefficient
    N = framesize*framesize*s3;
    C = 2*(N^2*dk*dk*df)^-1;
    
    % Compute spectrum
    A = (Ax.*conj(Ax)+Ay.*conj(Ay));
    clear Ax Ay
    maskblock = repmat(mask,[1 1 s3/2]);
    Skf = maskblock.*A*C;
    
    % Integrate directional spectrum with respect to f
    Skxky = sum(Skf,3,'omitnan')*df;
    
else
    
    N = s1*s2;
    C = (N^2*dk*dk)^-1;
    Skxky = mask.*(Ax.*conj(Ax)+Ay.*conj(Ay))*C;
    
end

% Integrate directional wavenumber spectrum with respect to k
[out_k,out_Sk] = az_int(Skxky,scale);

% Output
omnispect.k = out_k;
omnispect.S = out_Sk;

dirspect.kx = kx;
dirspect.ky = ky;

if s3 > 1
    
    dirspect.f = f;
    dirspect.Skf = Skf;
    
end

dirspect.Skxky = Skxky;
