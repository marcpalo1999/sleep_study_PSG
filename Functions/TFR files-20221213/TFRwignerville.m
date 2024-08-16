function TFR=TFRwignerville(x0,fs,fignum,analytical,freqband)

% Function that computes the Wigner-Ville distribution:
%
% TFR=TFRwignerville(x0,fs,fignum,analytical,freqband)
%
% Inputs:
%
% x0: input signal
% fs: sampling frequency
% fignum: figure where we want to fignum the result (0, does not draw, <0, creates new figure)
% analytical: indicate if it is desired to calculate WV from the analytical form of 
% the signal (to avoid aliasing and reduce interferences) or calculate it directly from the original signal
% (0: original signal, 1  (default): analytical form)
% frecband: range of frecbands to be plotted 
%
% output: TFR structure with the following fields:
%
% 
% type: 'WV'.
% NFFT: Number of points used to compute the fft
% signal: Structure containing 2 fields:
%       -x (signal amplitude vector Nx1)
%       -fs (sample frequency).
% TFR: Time-frequency representation matrix (NFFTxN)
% f: frequency vector (1xNFFT), from 0 to fs
% t: Time vector (1xN), from 0 to (N-1)/fs
% freqband: lowest and highest frequencies  (1x2 vector) to be plotted with the plotRTF function 
%
% Example:
%
% fs=1;
% t=0:1/fs:400;                 
% y=chirp(t,0.1,t(end),0.4);   
% TFR=wignerville(y,fs,1,1,[0 fs/2])
%
%
% Abel Torres (abel.torres@upc.edu), IBEC-ESAII-UPC
%


if nargin<1, % Example signal
       t=0:400;                  
       x0=chirp(t,0.1,t(end),0.4); 
end
if nargin<2,fs=1;end % Default sample frequency
if nargin<3,fignum=-1;end % Default: plot in a new figure 
if nargin<4, analytical=1;end % Default analytical
if nargin<5,freqband=[0 fs/2];end % Default whole frequency range

x0=x0(:); 
x0=detrend(x0); 
L=length(x0);

if analytical>0,
    % The analytical form of the signal is used in order to avoid aliasing
    % and reduce cross-terms
    x0=real(x0);
    x0 = hilbert(x0);
end

% The number of display frequency points is set to the next power of 2 greater than L
NFFT=2^nextpow2(L);

x=zeros(NFFT,1);
x(1:L,1)=x0;
% Zero padding
x(NFFT,1)=0;
cx=conj(x);
wx=zeros(NFFT,L);

% compute the product  cx(t-tau/2)*x(t+tau/2)
for i=1:L,
    taumaximo=min(i-1,L-i);
    tau2=-taumaximo:taumaximo;
    indices=1+tau2+NFFT*(tau2<0);
    wx(indices,i)=cx(i-tau2).*x(i+tau2);
end

% DFT of cx(t-tau/2)*x(t+tau/2) regarding delay tau (rows):
wx=fft(wx);

% The coefficients of the WVD have only a real part.
wx=real(wx);

% Time and frequency vectors
t=(0:L-1)/fs;
f=(0:NFFT-1)*(fs/2)/(NFFT);

% Output struct TFR:
TFR.type='WV';
TFR.NFFT=NFFT;
TFR.signal.x=real(x0);
TFR.signal.fs=fs;
TFR.TFR=wx;
TFR.f=f;
TFR.t=t;
TFR.freqband=freqband;    


if fignum,
  imagecontour=1;
  putcolorbar=0;
  if fignum>0,
      plotTFR(TFR,fignum,imagecontour,putcolorbar)
  else,
      plotTFR(TFR,0,imagecontour,putcolorbar)
  end
end
