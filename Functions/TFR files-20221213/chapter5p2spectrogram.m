[x,tx,fs]=gensig;
Nw=round(0.2*fs)*5;
window=hamming(Nw);
Noverlap=fix(Nw/2);
Noverlap=Nw-1;
NFFT=pow2(nextpow2(Nw));
fs=100;

[s,f,t,p]=spectrogram (detrend(x),window,Noverlap,NFFT,fs);
%U=window(:)'*window(:)/Nw;
%p=(abs(s).^2)/fs/Nw/U;
%p(2:end-1,:)=2*p2(2:end-1,:);
figure
subplot 211
plot(tx,x)
xlabel('Time (s)')
Ex=sum(x.^2)/fs;
title(['Test signal, Energy: ' num2str(Ex) ' V^2*s'])
ylabel('Amplitude (V)')

subplot 212
imagesc(t,f,p)
set(gca,'YDir','normal')
set(gca,'XLim',[tx(1) tx(end)])
xlabel('Time (s)');
ylabel('Frequency (Hz)');
Ep=sum(sum(p))*(fs/NFFT)*((Nw-Noverlap)/fs);
title(['Spectrogram, Hamming window (' num2str(Nw/fs) ' s) with ' num2str(Noverlap/fs) ' s of overlap, Energy: ' num2str(Ep) 'V^2*s'])
