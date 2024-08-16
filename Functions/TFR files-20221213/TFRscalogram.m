
function TFR=TFRscalogram(x,fs,fl,fh,N,k,envelope,freqband,fignum)

% TFR=TFRscalogram(x,fs,fl,fh,N,k,envelope,freqband,fignum)
%
% OUTPUT:
%
% TFR: Structure with these fields:
% - .type: character string indicating the type of TFR (SC or SC2 for
%    Morlet or Gaussian Derivative Scalogram)
% - .TFR: Matrix of the time-frequency representation
% - .f: Frequency axis
% - .t: Time axis
% - .fl: Minimum frequency computed
% - .fh: Maximum frequency computed
% - .freqband: 2-component vector indicating  minimum and maximum frequencies
%    to be plotted with the plotTFR function
% - .N: Number of frequency points used
% - .k: Number of oscillations of the wavelet inside the
%    envelope (in case of Morlet type wavelet) or order of the derivative
%    (in case of Gaussian derivative wavelet type).
% - .envelope: Type of envelope used
% - .signal: Structure with 2 components:
%   - .x: input time signal
%   - .fs: sampling frequency
%
% If envelope=-1 it generates the scalogram with a wavelet
% of the type k-th derivative (1<=k<=6) of Gaussian (exp(-t^2/2).
%
% If  envelope>0 (1,2,4 or 4) it generates the scalogram 
% with a wavelet of the morlet type with k periods in the time support
% of the envelope used (1: exponential envelope, 2: hanning envelope, 
%  3: blackmanning envelope, 4: blackmanning envelope),
%
% For the Morlet wavelet k=20/pi and envelope=1 (exponential).
% 
% INPUTS:
%
% x: input time signal
% fs: sampling frequency
% fl: Minimum frequency that should be calculated
% fh: Maximum frequency that should be calculated
% N: Number of frequency points
% k: Number of wavelet oscillations within the envelope 
% envelope: Number of wavelet oscillations within the envelope
% envelope: Envelope used (-1, 1,2,3 or 4)
% freqband: Range of drawn frequencies (vector 1x2) 
% fignum: Indicates if the figure where the result is to be plotted 
%
%   Example:
%
%       fm=1;
%       t=0:1/fm:400;                 
%       y=chirp(t,0.1,t(end),0.4);        
%       TFR=escalograma(y,fm,0.01,50,1024,20/pi,1,1[0 fm/2])
%
% Biomedical Signal Analysis Toolbox
% Abel Torres (abel.torres@upc.edu), IBEC-ESAII-UPC

if nargin<1, 
       t=0:400;                  
       x=chirp(t,0.1,t(end),0.4); 
end
x=detrend(x);
x=x(:)';
L=length(x);

if nargin<2,fs=1;end
if nargin<3,fl=0;end
if nargin<4,fh=fs/2;end
if nargin<5,N=1024;end
if nargin<6,k=20/pi;end
if nargin<7,envelope=1;end
if nargin<8 | length(freqband)<2,freqband=[fl fh];end
if nargin<9,fignum=-1;end

signal.x=x;
signal.fs=fs;
TFR.type='SC';
TFR.N=N;
TFR.k=k;
TFR.envelope=envelope;
TFR.signal=signal;
TFR.freqband=freqband;



Tmax=L/fs;
if envelope==-1,
      
else,
  fl=max(fl,k/Tmax);
end
fh=min(fs/2,fh);
f=linspace(fl,fh,N);

if nargin<2,fs=1;end


for i=1:N,
    if envelope==-1, % Gaussian derivative Wavelet 
      [h,T]=dergauss(fs,f(i),k,0);        
    else, % Morlet type Wavelet
      [h,T]=morletwavelet(fs,f(i),envelope,k,0);
    end
    C(i,:)=wkeep(conv(x,conj(h)),L);
end

t=(1:L)/fs;
TFR.TFR=abs(C).^2;
TFR.f=f;
TFR.t=t;
TFR.freqband=freqband;
TFR.fl=fl;
TFR.fh=fh;

if fignum<0,
    plotTFR(TFR)
elseif fignum>0,
    plotTFR(TFR,fignum)
end
