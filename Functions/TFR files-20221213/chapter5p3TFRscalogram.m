%--------------------------------------------------
% SIGNAL GENERATION 
numfig=0;
[s,t,fs]=gensig('default',numfig);
s=detrend(s);
%--------------------------------------------------------------
% SCALOGRAM:

numfig=6;
TFR=TFRscalogram(detrend(s),fs,.05,20,2048,20/pi,1, [0 20],numfig);
axis([0 10 0 20])

numfig=numfig+1;
TFR=TFRscalogram(detrend(s),fs,.05,20,2048,10/pi,1, [0 20],numfig);
axis([0 10 0 20])

numfig=numfig+1;
TFR=TFRscalogram(detrend(s),fs,.05,20,2048,40/pi,1,[0 20],numfig);
axis([0 10 0 20])


