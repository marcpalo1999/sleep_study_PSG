function [s,t,fs]=gensig(parameters,fignum)


% Function to generate the Time-Frequency Representations test signal
% composed of a convex chirp, several sinusoidal, two attenuated 
% oscillatory components and background random noise.
% 
% The parameters input should be a struct with the following fields:
%
% parameters.fsin: Vector 1x4 with the frequency of the stationary sinusoidals
% parameters.tsin: Vector 1x2 with the times of change of the sinusoidals
% parameters.fchirp: Vector 1x2 with the initial and final frequencies of the chirp
% parameters.tchirp: Vector 1x2 with the initial and final time of the chirp
% parameters.fwave: Vector 1x2 with the wavelets frequencies
% parameters.twave: Vector 1x2 with the times of the wavelets
%
% The outputs are the signal samples (s), the time vector (t) and the 
% sampling frequency (fs).
% 
% If we do not indicate input or we make parameters='default', some default
% values are loaded.
%
% The second input (fignum) indicates the figure number where the resulting 
% signal is to be plotted (if it is equal to zero, the signal is not plotted).
% 
% Biomedical Signal Analysis 
% Barcelona East School of Engineering
% abel.torres@upc.edu
% 


if nargin<1,
    parameters='default';
end
if strcmp(parameters,'default'), % Default parameters
    
    % COMPONENT1
    tsin=[3 6];
    fsin= [10 7 5 15];
    % COMPONENT2
    fchirp=[3 5];
    tchirp=[3 7];
    % COMPONENT4
    twave=[5 7.5];
    fwave=[10 15];
    
else,
    tsin=parameters.tsin;
    fsin=parameters.fsin;
    tchirp=parameters.tchirp;
    fchirp=parameters.fchirp;
    twave=parameters.twave;
    fwave=parameters.fwave;

end
if nargin<2,fignum=1;end


% Signal duration (in s)
TF=10;
% Sampling frequency (Hz)
fs = 100;
% Time vector
t = (0:1/fs:TF);
% Signal length
N=length(t);



% COMPONENT 1: Stationary sinusoidals
ts1=[0 tsin(1);tsin(1) TF;0 tsin(2);tsin(2) TF];
frequency=fsin;
s1 = zeros(1,N);
amplitude=0.3;
for j=1:4
    s1(t>=ts1(j,1) & t<=ts1(j,2)) = s1(t>=ts1(j,1) & t<=ts1(j,2)) + amplitude*sin(2*pi*frequency(j)*t(t>=ts1(j,1) & t<=ts1(j,2)));
end


% COMPONENT 2: Chirp betweem times t0 and t1, and frequencies between f0 and f1
t0=tchirp(1);t1=tchirp(2);
f0=fchirp(1);f1=fchirp(2);
s2=zeros(1,N);
for i=1:length(f0),
    t2=t(t0(i)*fs:t1(i)*fs);
    y2=chirp(t2-t2(1),f0(i),t1(i)-t0(i),f1(i),'q',[],'convex');
    s2(t>=t2(1) & t<=t2(end))=s2(t>=t2(1) & t<=t2(end))+y2;
end
amplitude=1;
s2=amplitude*s2;

% COMPONENT 3: Random noise
randn('state',1);
ruido = randn(1,N);
amplitude=0.1;
s3=ruido*amplitude;


%  COMPONENT 4:  Wavelets at times t1 and t2 and frequencies f1 and f2
t1=twave(1);f1=fwave(1);
t2=twave(2);f2=fwave(2);
s4=zeros(1,N);
psi=morletwavelet(fs,f1,1,30/pi,0);
N=round(t1*fs-length(psi)/2);
s4(N+1:N+length(psi))=psi;
psi=morletwavelet(fs,f2,1,30/pi,0);
N=round(t2*fs-length(psi)/2);
s4(N+1:N+length(psi))=psi;
amplitude=1.6;
s4=amplitude*real(s4);
s=s1+s2+s3+s4;




if fignum,
    figure(fignum),clf(fignum,'reset')
    SC=get(0,'ScreenSize');
    set(gcf,'Position',[SC(3)*0.05 SC(4)*0.05 SC(3)*0.9 SC(4)*0.8]);
    subplot(5,1,1)
    plot(t,s1)
    title(['Stationary sinusoidal signals (frequencies: ' num2str(fsin) ' Hz, times: ' num2str(tsin) ' s'])
    subplot(5,1,2)
    plot(t,s2)
    title(['Chirp signal (frequencies: ' num2str(fchirp) ' Hz, times: ' num2str(tchirp) ' s'])
    subplot(5,1,3)
    plot(t,s3)
    title('Random noise')
    subplot(5,1,4)
    plot(t,s4)
    title(['Wavelets (frequencies: ' num2str(fwave) ' Hz, times: ' num2str(twave) ' s'])
    subplot(5,1,5)
    plot(t,s)
    title('xsint')
    xlabel('t(s)')
end





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