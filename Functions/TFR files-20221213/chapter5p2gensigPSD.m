%--------------------------------------------------
% SIGNAL GENERATION AND POWER SPECTRAL DENSITY

numfig=1;
[s,t,fs]=gensig('default',numfig);
t=(0:length(s)-1)/fs;
% Welch's periodogram with higher frequency resolution
[pxx1,f] = pwelch(s,hanning(4*fs),2*fs,2048,fs);
% Welch's periodogram with lower frequency resolution
[pxx2,f] = pwelch(s,hanning(fs/2),fs/4,2048,fs);

numfig=numfig+1;
figure(numfig),clf
SC=get(0,'ScreenSize');
set(gcf,'Position',[SC(3)*0.05 SC(4)*0.05 SC(3)*0.9 SC(4)*0.8]);
plot(f,pxx1,f,pxx2)
set(gca, 'Xlim',[0 20])
xlabel('Hz')
ylabel('V^2/Hz')
title('Welch periodogram with Hamming window and 50% of overlap ')
legend('Window size: 4 s' , 'Window size: 0.5 s')

