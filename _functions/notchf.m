function [sireal,siimag] = notchf( si, dt, lowp , highp )
%        [sireal,siimag] = notch( si, dt, lowp , highp )
%
% This function fourier transforms a signal   and removes
%  energy between periods 'lowp' and 'highp'. Mean and trend are preserved.
%
%   si      = input signal;
%   dt      = sampling rate [s];
%   lowp    = low period cutoff [s];
%   highp    = high period cutoff [s];

% Frequency range estimation.
%
si = si(:) ; lsi = length(si) ; 
if rem(lsi,2)==1
  si=[si;si(lsi)];
end
corrsi=si-detrend(si);
clsi = ceil( lsi/2 ) ;
f = ( -clsi : 1 : clsi-1 )' / lsi / dt ;
%i = sqrt( -1 ) ;

% Signal transformation to the frequency domain.
%
SI = fft( detrend( si ) );
SIshift = fftshift( SI ) ;

% Double integration.
%
SI(1:clsi,1)    = SIshift(1:clsi) ;
SI(clsi+1,1)    = 0 ;
SI(clsi+2:lsi,1) = SIshift(clsi+2:lsi);

% Energy elimination at frequencies smaller than 1/cut-off.
%
ind =  (f >= -1/lowp & f <= -1/highp ) | (f <= 1/lowp & f >= 1/highp ); 
%if length(ind) == 1, SI(ind) = 0; else, SI(ind) = zeros( ind ) ; end
SI(ind) = 0;

% Signal transformation back to the temporal domain.
%
si = ifft( fftshift( SI ) ) ;
sireal = real( si ) + corrsi ;
siimag = imag( si ) ;
if rem(lsi,2)==1
   sireal=sireal(1:lsi);
   siimag=siimag(1:lsi);
end