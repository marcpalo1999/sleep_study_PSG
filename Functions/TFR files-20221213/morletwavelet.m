function [psi,t] = morletwavelet(fs,f0,envelopetype,k,fignumber)

% [psi,t] = morletwavelet(fs,f0,envelopetype,k,fignumber)
%
% Generates a Morlet type wavelet:
%
% OUTPUTS:
%
% psi: complex wavelet
% t: Time vector %
%
% INPUTS:
%
% fs: sampling frequency
% f0: wavelet center frequency
% envelopetype: defines the type of envelope used
%   =1: damped exponential envelope: exp((-t.^2)/2) with time support
%        between -4 and 4
%   =2: Hanning window type envelope: (1+cos(pi/4*t))/2, with time support
%       between -4 and 4
%   =3: Blackman window with time support between -4 and 4
%   =4: Rectangular window with time support between -4 and 4
% k: Defines the number of periods of the complex oscillating signal (btween -4 and 4)
%       that is multiplied by the envelope (in the Morlet wavelet, k=20/pi, approximately 6.3662)
% fignumber: draws the wavelets generated with the 4 types of envelope
%
%  Biomedical Signal Analysis Toolbox
%  Abel Torres (abel.torres@upc.edu)

if nargin<1,fs=1;end
if nargin<2,f0=0.1;end
if nargin<3,envelopetype=1;end
if nargin<4,k=20/pi;end
if nargin<5,fignumber=1;end

N=round(k*fs/f0)+1; % Number of samples of the wavelet
t=linspace(-4,4,N); % Time support of the wavelet
%=================================================================
% Generation of the envelope function:
% type 1: Exponential (Morlet):
envelopes(1,:)=exp(-(t.^2)/2);
% type 2: Hanning window:
envelopes(2,:)=hanning(N)'; % =(1+cos(pi/4*t))/2;
% type 3: Blackman window:
envelopes(3,:)=blackman(N)';
% type 4: Rectangular window:
envelopes(4,:)=rectwin(N)';
%=================================================================
% Generation of the complex oscillatory function with k oscillation in the
% time support (in the Mortlet wavelet k=20/pi)
complexoscillation=exp(i*2*pi*(k/8)*t);

%=================================================================
% Generation of the complex wavelet
psi=envelopes(envelopetype,:).*complexoscillation;

% Energy normalization:
T=N/fs/2; % time support
c=sqrt(N/(2*T*sum(abs(psi).^2)));
psi=c*psi;

% Time vector:
t=linspace(-T,T,N);

if fignumber,
    
    figure(fignumber)
    SC=get(0,'ScreenSize');
    set(gcf,'Position',[SC(3)*0.05 SC(4)*0.05 SC(3)*0.9 SC(4)*0.8]);
    
    
    NFFT=2^15;NFFT2=NFFT/2;
    f=(0:NFFT2)/(NFFT2)*(fs/2);
    texto={'Exponential','Hanning','Blackman','Rectangular'};
    for n=1:4,
        psi_n=envelopes(n,:).*complexoscillation;
        c=sqrt(N/(2*T*sum(abs(psi_n).^2)));
        psi_n=c*psi_n;
        PSI_N=fft(psi_n,NFFT);
        PSI_N=PSI_N(1:NFFT2+1);
        aPSI_N2=abs(PSI_N).^2;
        [~,posmax]=max(aPSI_N2);
        fmax=f(posmax);
        
        subplot(4,3,1+(n-1)*3)
        plot(t,c*envelopes(n,:),'r:',t,real(c*complexoscillation),'b:',t,real(psi_n))
        title([texto{n} ' window (real part)'])
        if n==4,xlabel('t(s)'),end
        subplot(4,3,2+(n-1)*3)
        plot(t,c*envelopes(n,:),'r:',t,imag(c*complexoscillation),'b:',t,imag(psi_n))
        title([texto{n} ' window (imaginary part)'])
        if n==4,xlabel('t(s)'),end
        subplot(4,3,3+(n-1)*3)
        plot(f,aPSI_N2),set(gca,'XTick',[0 fmax min([fs/2 2*f0])],'XLim',[0 min([fs/2 4*f0])])
        title('PSD (V^2/Hz)')
        if n==4,xlabel('f(Hz)'),end
        
    end
    
end