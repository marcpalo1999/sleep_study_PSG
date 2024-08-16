function A=TFRambiguity(x0,fs,fignum,zoomtaunu)

% Function that computes the ambiguity function of the signal x0:
%
% A=ambiguity(x0,fs,fignum)
%
% Inputs:
%
% x0: signal
% fs: sampling frequency
% fignum: figure where you want to draw the result (0, no drawing, <0, create new figure)
% zoomtaunu: Maximum displayed value of the time and frequency delays (tau and nu, respectively)
%
% Output: A structure with the following fields:
%
%
% A = 
% 
% A: Ambiguity function in a format suitable for drawing with imagesc
% tau: Time delay vector (1xNFFT)
% nu: Frequency delay vector (1xNFFT)
%
% Example:
%
%       t=0:400;                  
%       x0=chirp(t,0.1,t(end),0.2); 
%       A=ambiguity(x0,1,-1)
%
% Abel Torres (abel.torres@upc.edu), IBEC-ESAII-UPC

if nargin<1, % Example signal
       t=0:400;                  
       x0=chirp(t,0.1,t(end),0.2); 
end       
if nargin<2,fs=1;end % Default sample frequency
if nargin<3,fignum=-1;end % Default figure (new)
L=length(x0);
if nargin<4,zoomtaunu=[L/fs 0.5*fs];end % Default range of time and frequency delays displayed

%x0=resample(x0,2,1);
%L=length(x0);
%fs=2;
% The number of display frequency points is set to the next power of 2 greater than L
NFFT=2^(nextpow2(L));


x0=x0(:); 
x0=detrend(x0); 

% The analytical form of the signal is used in order to avoid aliasing
% and reduce cross-terms
x0=real(x0);
x0 = hilbert(x0);


x=zeros(NFFT,1);
x(1:L,1)=x0;
% Zero padding
x(NFFT,1)=0;
cx=conj(x);
wx=zeros(NFFT,NFFT);

% compute the product  cx(t-tau/2)*x(t+tau/2)
for i=1:L,
    taumaximo=min(i-1,L-i);
    tau2=-taumaximo:taumaximo;
    indices=1+tau2+NFFT*(tau2<0);
    wx(indices,i)=cx(i-tau2).*x(i+tau2);
end


% Inverse DFT of cx(t-tau/2)*x(t+tau/2) regarding time  (columns)
wx=ifft(wx.').';

tau=(-NFFT/2:NFFT/2-1)/(fs/2);
nu=(-NFFT/2:NFFT/2-1)/(NFFT/2)*(fs/2);


A.A=fftshift(wx);
A.tau=tau;
A.nu=nu;
A.signal.x=real(x0);
A.signal.fs=fs;
A.NFFT=NFFT;
A.zoomtaunu=zoomtaunu;

if fignum,

if fignum>0,figure(fignum),else,figure,end
clf
SC=get(0,'ScreenSize');
set(gcf,'Position',[SC(3)*0.05 SC(4)*0.05 SC(3)*0.9 SC(4)*0.8]);
subplot(2,2,2)
imagesc(A.tau,A.nu,abs(A.A'))
set(gca,'XLim',A.tau([1 end]))
set(gca,'YLim',A.nu([1 end]))
set(gca,'XLim',zoomtaunu(1)*[-1 1])
set(gca,'YLim',zoomtaunu(2)*[-1 1])
XL=get(gca,'XLim');
YL=get(gca,'YLim');
title('Ambiguity function A_x_x(tau,nu)')	
xlabel('Tau (s)');
ylabel('Nu (Hz)');


subplot(2,2,3)
t=(0:length(x0)-1)/fs;
plot(t,real(x0))
axis tight
title('Time signal x(t)')	
xlabel('t(s)');
ylabel('Amplitud (V)');
set(gca,'XLim',t([1 end]))



subplot(2,2,1)
plot(abs(A.A(:,NFFT/2+1)),A.nu)
set(gca,'XDir','reverse')
set(gca,'YLim',YL);
title('A(tau=0,nu) (frecuency autocorrelation)')	
ylabel('Nu (Hz)');


subplot(2,2,4)
plot(A.tau,real(A.A(:,NFFT/2+1)))
set(gca,'XLim',XL);
title('A(tau,nu=0) (time autocorrelation)')	
xlabel('Tau (s)');



end
