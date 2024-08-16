function [Pmband,Pmband2Pm,Pm]=powerband(Pxx,f,fband,showfigure)

%  function [Pmband,Pmband2Pm,Pm]=powerband(Pxx,f,fband)
%
% Function that given the PSD of a signal characterized by its vectors
% Pxx and f, calculates 'Pmband', the mean power in the frequency
% band (Pmband) delimited by the frequencies entered in the two-element
% vector 'fband'.

% It also calculates the total power (Pm) and the relative power 
% in the band (Pmband2Pm)
%

if nargin<1,
    NFFT=pow2(12);
    f=(0:NFFT)/NFFT*0.5;
    Pxx=f.*exp(-20*f);
    fband=f(round(length(f)*[0.1 0.25]));
elseif nargin<2,
    NFFT=length(Pxx);
    f=(0:NFF^T)/NFFT*0.5;
    fband=f(round(length(f)*[0.1 0.25]));
elseif nargin<3,
    fband=f(round(length(f)*[0.1 0.25]));
end
if nargin<4,showfigure=1;end

deltaf=mean(diff(f));
Pm=deltaf*sum(Pxx);
Pmband=deltaf*sum(Pxx(f>=fband(1) & f<=fband(2)));
Pmband2Pm=Pmband/Pm;

if showfigure>0,
    APPLICATION='Chapter 4.3: PSD power in a band';
    [existFlag,figNumber]=figflag(APPLICATION);
    if ~existFlag,
        figNumber=figure('Name',APPLICATION,'NumberTitle','off');
    else,
        figure(figNumber),clf
    end
    
    %plot(f,Pxx,'b')
    ffill=[f(1) f f(end)];
    Pxxfill=[0 Pxx 0];
    h=fill(ffill,Pxxfill,[1 0 0],'FaceAlpha',.5,'LineWidth',5,'EdgeColor',[1 0 0]);
    set(gca,'XLim',f([1 end]));
    title(['PSD: Power in the band ' num2str(fband(1)) '-' num2str(fband(2)) ' Hz']);
    xlabel('f(Hz)');
    ylabel('Pxx(f) [V^2s]');
    
    
    
    hold on
    ffill=f(f>=fband(1) & f<=fband(2));
    ffill=[ffill(1) ffill ffill(end)];
    Pxxfill=Pxx(f>=fband(1) & f<=fband(2));
    Pxxfill=[0 Pxxfill 0];
    fill(ffill,Pxxfill,[0 1 0.5],'FaceAlpha',1,'LineWidth',0.5,'EdgeColor',[0 1 0]);
    hold off
    
    legendtext=['Pm: ' num2str(Pm) ' V^2'];
    addtext=['Pm in band: ' num2str(Pmband) ' V^2 (' num2str(round(Pmband2Pm*10000)/100) ' %)'];
    legendtext=str2mat(legendtext,addtext);
    legend(legendtext,'FontSize',12);
end

