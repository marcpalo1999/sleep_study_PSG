%--------------------------------------------------
% SIGNAL GENERATION 
numfig=0;
[s,t,fs]=gensig('default',numfig);
s=detrend(s);
%--------------------------------------------------------------
% SPECTROGRAM:
numfig=numfig+1;
NFFT=1024;
[TFRsp1]=TFRspectrogram(s,fs,NFFT,fs/2,-1,1,[0 20],numfig);
set(gca, 'Xlim',[0 10])

numfig=numfig+1;
[TFRsp2]=TFRspectrogram(s,fs,NFFT,fs*2,-1,1,[0 20],numfig);
set(gca, 'Xlim',[0 10])

numfig=numfig+1;
[TFRsp3]=TFRspectrogram(s,fs,NFFT,fs/4,-1,1,[0 20],numfig);
set(gca, 'Xlim',[0 10])

