%--------------------------------------------------
% GENERACIÓN SEÑAL Y DENSIDAD ESPECTRAL DE POTENCIA

dibujar=0;
[s,t,fs]=gensig('default',dibujar);

dibujar=100;
Nplots=4;
figure(dibujar),clf
subplot(Nplots,1,1)
plot(t,s)
xlabel('t[s]'),ylabel('[V]'),title('xsint(t)')

%--------------------------------------------------------------
% ESPECTROGRAMA:
factor=[0.25 .5 2];
for i=1:3,

    subplot(Nplots,1,i+1)
    D=round(fs*factor(i));S=D-1;
    NFFT=2^max([10 nextpow2(D)]);
    window=hamming(D);
    [STFT,f,tSTFT,SPEC]=spectrogram (s,window,S,NFFT,fs);
    imagesc(tSTFT,f,SPEC)
    xlabel('t[s]'),ylabel('f[Hz]')
    title(['Spectrogram of xsint(t) with a ' num2str(D/fs) ' s Hamming window [V^2/Hz]'])
    set(gca,'YDir','normal')
    set(gca,'XLim',[t(1) t(end)])
    YL=[0 20];
    set(gca,'YLim',YL)

end


