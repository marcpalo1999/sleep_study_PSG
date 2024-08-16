function TFR=TFRCohen(A,kernel,fignum,freqband)

% Función que computa una TFR de la clase Cohen generalizada a partir de la
% función de ambigüedad A y del kernel (window) utilizado:
%
%   TFR=TFRCohen(Ashift,window,NFFT,L,fs,type)
%
% Entradas:
%   Ashift: Función de ambigüedad de la señal (ver función ambiguedad.m)
%   ventanta: kernel utilizado (ver función genvent.m)
%   NFFT: dimensión de A
%   L: Número de muestras de la señal original
%   fs: frecuencia de muestro
%   type: string donde indicamos el type de distribución generado
%   fignum: figura donde se desea fignum el resultado (0, no dibuja, <0, crea figura nueva)
%   freqband: rango de freqband que se fignuman
%   zoomtauni: retardos temporal (tau) y frecuencial (ni) máximos de la
%   función de ambigüedad que se representarán
%
% Salida: Estructura TFR con los siguientes campos:
%
% TFR = 
% 
%            type: string donde se indica el type de distribución generado
%            NFFT: Número de puntos de la fft
%          senyal: Estructura que contiene 2 campos: x (señal) y t (vector
%          tiempo)
%             TFR: Representación tiempo frecuencia (NFFT*N)
%               f: Vector de freqband (1*NFFT), de 0 a fs/2
%               t: Vector de tiempo (1*N), de 0 a (N-1)/fs
%     freqband: Rango de freqband a fignum por plotTFR (1*2) 
%
%   Ejemplo:
%
%       fs=100;
%       t=0:1/fs:4;                   % +/-2 secs @ 1kHz sample rate
%       y=chirp(t,5,4,10);       % Start @ 100Hz, cross 200Hz at t=1sec 
%       [As,A]=ambiguedad(y,fs,1);
%       NFFT=length(A.tau);
%       L=length(y);
%       par.type=1;
%       par.sigma=0.1;
%       par.NFFT=NFFT;
%       window=genvent(1,par);
%       TFR=RID2(As,window,NFFT,L,fs,'Dist. Choi-Williams, sigma=0.1',1,[0 15])
%
% Abel Torres (abel.torres@upc.edu), IBEC-ESAII-UPC


if nargin<1, % Si no ponemos entrada se ejecuta un ejemplo
       t=0:400;                  
       x0=chirp(t,0.1,t(end),0.4); 
       A=TFRambiguity(x0,1,0);
end
if nargin<2,
    par.type=3;
    par.NFFT=A.NFFT;
    kernel=TFRkernel(par,0);
end
if nargin<3,
    fignum=-1;
end
if nargin<4,freqband=[0 A.signal.fs/2];end % Default whole frequency range


% Amplicamos kernel a la función de ambigüedad:
Ak=fftshift(A.A).*fftshift(kernel.kernel);

% Computamos la FFT en tiempo y frecuencia (TFR)
wx=fft2(Ak);
wx=(real(wx.')');
wx=wx(:,1:length(A.signal.x));
TFR.type='Cohen';
TFR.TFR=wx;
TFR.f=(1:A.NFFT)/A.NFFT*(A.signal.fs/2);
TFR.t=(0:length(A.signal.x)-1)/A.signal.fs;
TFR.A=A;
TFR.kernel=kernel;
TFR.signal=A.signal;
TFR.NFFT=A.NFFT;
TFR.freqband=freqband;



if fignum<0,
    plotTFR(TFR)
elseif fignum>0,
    plotTFR(TFR,fignum)
end


  

