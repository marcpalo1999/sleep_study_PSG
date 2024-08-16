function TFR=TFRspectrogram(x,fs,NFFT,Nwindow,Noverlap,window,freqband,fignum)

% OUTPUTS:
%
% TFR: Structure with 8 fields:
% - .type: string of characters indicating the type of TFR (SP in this case)
% - .TFR: Matrix of the time-frequency representation
% - .f: Frequency axis
% - .t: Time axis
% - .freqband: 2-component vector indicating minimum and maximum frequencies to be plotted with the plotTFR function
% - .NFFT: Number of points used in the FFT calculation 
% - .Nwindow: Number of samples of the time window
% - .Noverlap: Number of overlapping samples between windows
% - .window: Type of window used (1: Hamming, 2: Hanning, 3: rectwin)
% - .signal: 2-component structure:
%    - .x: Input signal
%    - .fs: Sampling frequency
%
% INPUTS:
%
% x: Input signal
% fs: sampling frequency
% NFFT: Number of frequency points
% Nwindow: Number of window samples
% Noverlap: Number of overlapping samples between windows (<Nvent)
% (if we put Noverlap=-1: Noverlap=Nwindow-1)
% window: Type of window (1.Hamming, 2.Hanning,3. Rectangular)
% fignum: Indicates the figure number where the TFR is plotted
% freqband: Range of plotted frequencies (vector 1x2) 
%
% Example:
%
% fs=1;
% t=0:1/fs:400;                 
% y=chirp(t,0.1,t(end),0.4);        
% Nwindow=round(L/20);
% Noverlap=-1;
% TFR=TFRspectrogram(y,fs,Nwindow,Noverlap,1,1,[0 fs/2]);
%
% Biomedical Signal Analysis Toolbox
% Abel Torres (abel.torres@upc.edu), IBEC-ESAII-UPC

if nargin<1, 
       t=0:400;                  
       x=chirp(t,0.1,t(end),0.4); 
end

L=length(x);
if nargin<2,fs=1;end
if nargin<4,Nwindow=round(L/20);end
if nargin<5,Noverlap=Nwindow-1;end
if nargin<6,window=1;end
if nargin<7,freqband=[0 fs/2];end
if nargin<8,fignum=-1;end

if Nwindow>=L,Nwindow=L;end 
if nargin<3,
NFFT=max(2^nextpow2(Nwindow),256);
end

signal.x=x;
signal.fs=fs;
TFR.type='SP';
TFR.NFFT=NFFT;
TFR.Nwindow=Nwindow;
TFR.Noverlap=Noverlap;
TFR.window=window;
TFR.signal=signal;


   
if Noverlap>=Nwindow | Noverlap<0,Noverlap=Nwindow-1;end
if window==1,
  window=hamming(Nwindow);
elseif window==2,
  window=hanning(Nwindow);
elseif window==3,
  window=rectwin(Nwindow);  
end

x=detrend(x);

[STFT,f,t,TFR.TFR]=spectrogram(x,window,Noverlap,NFFT,fs);
TFR.f=f(:)';
TFR.t=t(:)';
TFR.freqband=freqband;


if fignum<0,
    plotTFR(TFR)
elseif fignum>0,
    plotTFR(TFR,fignum)
end