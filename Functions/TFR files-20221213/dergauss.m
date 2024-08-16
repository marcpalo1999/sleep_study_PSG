function [psi,T,N]=dergauss(fm,f_equiv,Nderiv,dibujar)

% First we generate the wavelet in the time support of the exponential envelope (from -5 s to 5s).
fs=100;
t=-5:1/fs:5;
N=length(t);
T=range(t);
x0=exp(-1/2*t.^2);
if Nderiv==1,
    psi=(-t).*x0;
elseif Nderiv==2,
    psi=(-1+t.^2).*x0;
elseif Nderiv==3,
    psi=(3*t-t.^3).*x0;
elseif Nderiv==4,
    psi=(3-6*t.^2+t.^4).*x0;
elseif Nderiv==5,
    psi=(-15*t+10*t.^3-t.^5).*x0;
elseif Nderiv==6,
    psi=(-15+45*t.^2-15*t.^4+t.^6).*x0;
end
% Normalisation of wavelets to unit energy [E=(T/N)*sum(x.^2)]:
psi=psi/sqrt(N/(T*sum(psi.^2)));

% Calculation of the peak frequency of the wavelet
NFFT=1024;
f=(0:NFFT/2)/(NFFT/2)*(fs/2);
P=fft(psi,NFFT);
P=P(1:NFFT/2+1);
P=abs(P);
[a,b]=max(P);
fmax=f(b);



%Calculamos ya la wavelet del tipo que nos piden y con la fmax=f_equiv:
% Calculamos soporte temporal correspondiente a f_equiv [Ti=T*(fmaxi/f_equiv)]:
%if Nderiv<0,Nderiv=abs(Nderiv);signo=-1;else,signo=1;end

T=T*(fmax/f_equiv);
N=floor(T*fm);

t=linspace(-5,5,N);
x0=exp(-1/2*t.^2);
if Nderiv==1,
    psi=(-t).*x0;
elseif Nderiv==2,
    psi=(-1+t.^2).*x0;
elseif Nderiv==3,
    psi=(3*t-t.^3).*x0;
elseif Nderiv==4,
    psi=(3-6*t.^2+t.^4).*x0;
elseif Nderiv==5,
    psi=(-15*t+10*t.^3-t.^5).*x0;
elseif Nderiv==6,
    psi=(-15+45*t.^2-15*t.^4+t.^6).*x0;
end
%psi=psi*signo;
% Normalizamos la wavelet:
c=sqrt(N/(T*sum(psi.^2)));
%psi=c*psi;


P=fft(psi,NFFT);
f=(0:NFFT/2)/(NFFT/2)*(fm/2);
t=(1:N)/fm;
P=P(1:NFFT/2+1);P=abs(P);
[a,b]=max(P);fmax=f(b);

if dibujar,
    figure(dibujar+1)
    set(gcf,'Name','Wavelet utilizada')
    subplot(1,2,1),plot(t,psi)
    title(['T=' num2str(T) ', N=' num2str(N)])
    subplot(1,2,2),plot(f,abs(P)),
    set(gca,'XLim',[0 2.5*fmax]),
    set(gca,'XTick',[0 fmax floor(2.5*fmax)])
end