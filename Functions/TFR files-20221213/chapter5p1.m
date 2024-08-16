% 1. Energy and power: time domain
APPLICATION='Chapter5.1: Introduction (1/3)';
[existFlag,figNumber1]=figflag(APPLICATION);
if ~existFlag,  
     figNumber1=figure('Name',APPLICATION,'NumberTitle','off');
else,
    figure(figNumber1),clf
end


w1=2*pi*0.05;
w2=2*pi*0.1;
Fs=1;
t=0:1999;
x_A=[sin(w1*t) cos(w2*t)];
x_B=[cos(w2*t) sin(w1*t) ];
L=length(x_A);
t=0:length(x_A)-1;
x_C=0.5*sin(w1*t)+0.5*cos(w2*t);

NFFT=2^(nextpow2(L)+5);
X_A=fft(x_A,NFFT);
X_B=fft(x_B,NFFT);
X_C=fft(x_C,NFFT);
aX2_A=fftshift(abs(X_A).^2);
aX2_B=fftshift(abs(X_B).^2);
aX2_C=fftshift(abs(X_C).^2);
f=(0:NFFT-1)/NFFT*Fs-0.5;

letrasI='ABC';
letrasJ1={'x','aX2'};
letrasJ2='tf';
titulos={'x_A(t)=[sin(w1*t)  cos(w2*t)]','x_B(t)=[cos(w2*t) sin(w1*t)]','x_C(t)=[0.5*sin(w1*t)+0.5*cos(w2*t)]'};

subplot 211
plot(t,x_A)
xlim(2000+160*[-1 1])
title(titulos{1}([1 4:end]))
xlabel('t(s)')
ylabel('x (V)')

subplot 212
plot(f,aX2_A)
title(['|X_' letrasI(1) '(f)|^2'])
xlabel('f(Hz)')
ylabel(['|X_' letrasI(1) '|^2 (V^2)'])


APPLICATION='Chapter5.1: Introduction (2/3)';
[existFlag,figNumber2]=figflag(APPLICATION);
if ~existFlag,  
     figNumber2=figure('Name',APPLICATION,'NumberTitle','off');
else,
    figure(figNumber2),clf
end


for i=1:3
    for j=1:2,
        eval(['ejesX=' letrasJ2(j) ';'])
        eval(['ejesY=' letrasJ1{j} '_' letrasI(i) ';'])
        subplot(2,3,i+3*(j-1))
        plot(ejesX,ejesY)
        if j==1,xlim(2000+80*[-1 1]),title(titulos{i}),xlabel('t(s)'),ylabel(['x_' letrasI(i) '(V)'])
        else,title(['|X_' letrasI(i) '(f)|^2']),xlabel('f(Hz)'),ylabel(['|X_' letrasI(i) '|^2 (V^2)'])
        end
    end
end
        

APPLICATION='Chapter5.1: Introduction (3/3)';
[existFlag,figNumber3]=figflag(APPLICATION);
if ~existFlag,  
     figNumber3=figure('Name',APPLICATION,'NumberTitle','off');
else,
    figure(figNumber3),clf
end
for i=1:6,
    subplots(i)=subplot(2,3,i);
end

Twindow=100;
Nwindow=Twindow*Fs;
Noverlap=Nwindow-1;
NFFT=2^max([nextpow2(Nwindow) 12]);
[STFT_A,f2,t2] = spectrogram(x_A,hamming(Nwindow),Noverlap,NFFT,Fs);
[STFT_B,f2,t2] = spectrogram(x_B,hamming(Nwindow),Noverlap,NFFT,Fs);
[STFT_C,f2,t2] = spectrogram(x_C,hamming(Nwindow),Noverlap,NFFT,Fs);


for i=1:3
    
        eval(['TFR.signal.x=x_' letrasI(i) ';'])
        TFR.signal.fs=Fs;
        eval(['TFR.TFR=abs(STFT_' letrasI(i) ').^2;'])
        TFR.t=t2(:)';
        TFR.f=f2(:)';
        TFR.freqband=[0 Fs/2];
        TFR.type='SP';
        TFR.NFFT=NFFT;
        TFR.Nwindow=Nwindow;
        TFR.Noverlap=Noverlap;
        TFR.window=1;

        plotTFR(TFR,subplots([i+3 i]))

end
set(subplots,'XLim',2000+200*[-1 1])
for i=1:3,
    set(gcf,'CurrentAxes',subplots(i)),
    title(['x_' letrasI(i) ' (t)'])
    set(gcf,'CurrentAxes',subplots(i+3))
    title(['Spectrogram of x_' letrasI(i)])
end

