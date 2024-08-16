function plotTFR(TFR,fignum,imagecontour,putcolorbar)

% plotTFR(TFR,fignum,imagecontour,putcolorbar)
%
% - TFR: Structure with time-frequency representation (TFR) information
% generated with one of the TFR toolbox functions
% - fignum: figure number (<=0 for newplot) or axis labels where the TFR is
% to be displayed (default: newplot)
% - imagecontour: Scalar indicating wheter to display TFR with the MATLAB
% functions imagesc (imagecontour=1), contour (imagecontour=2) or meshc
% (imagecontour=3) (default: imagecontour=1)
% - putcolorbar: adds (putcolorbar=1) or does not add (putcolorbar=0)
% a colour bar to the TFR axes (default: putcolorbar=1);
%
% This structure must have the following fields:
%
%   .TFR: Time-frequency representation matrix
%   .f: Frequency vector
%   .t: Time vector
%
% and also could have the following fields to improve the representation:
%
%   .freqband: vector 1x2 with the range of frequencies to be represented
%   .type: string with the type of TFR (SP, AR, SC,SC2, WV, Cohen, CHW, GED, BJD or CKD)
%
%
% Time-Frequency Representation Toolbox
% abel.torres@upc.edu
%
%


if nargin<1,
    t=0:400;
    x0=chirp(t,0.1,t(end),0.4);
%     A=TFRambiguity(x0,1,0);
%     par.type=3;
%     par.NFFT=A.NFFT;
%     kernel=TFRkernel(par,0);
%     TFR=TFRCohen(A,kernel,0);
TFR=TFRspectrogram;
end
if nargin<2,fignum=-1;end
if nargin<3, imagecontour=1; end
if nargin<4, putcolorbar=0;end

if ~isstruct(TFR),
    TFR2=TFR;
    clear TFR
    TFR.TFR=TFR2;
    clear TFR2
end


if ~isfield(TFR,'type'),
    TFR.type='-';
end
if ~isfield(TFR,'t'),
    [nf,nt]=size(TFR.TFR);
    TFR.t=(0:nt-1);
end
if ~isfield(TFR,'f'),
    [nf,nt]=size(TFR.TFR);
    TFR.f=(0:nf-1)/nf*0.5;
end
if ~isfield(TFR,'freqband'),
    freqband=TFR.f([1 end]);
else,
    freqband=TFR.freqband;
end


if length(fignum)==1,
    if ishandle(fignum),
        if strcmp(get(fignum,'Type'),'figure'),
            set(0,'CurrentFigure',fignum),clf
            nejes=0;
        elseif  strcmp(get(fignum,'Type'),'axes'),
            set(gcf,get(fignum,'Parent'))
            nejes=1;
        else,
            fignum=figure;
            nejes=0;
        end
    else,
        if fignum>0,fignum=figure(ceil(fignum));clf
        else, fignum=figure;
        end
        nejes=0;
    end
else,
    if all(ishandle(fignum)) 
        if all(strcmp(get(fignum,'Type'),'axes'))
            ejes=fignum;
            nejes=length(fignum);
        else,
            fignum=figure;
            nejes=0;
        end
    else,
        fignum=figure;
        nejes=0;
    end
end
SC=get(0,'ScreenSize');
set(gcf,'Position',[SC(3)*0.05 SC(4)*0.05 SC(3)*0.9 SC(4)*0.8]);

if nejes==0,
    if isfield(TFR,'signal')
        SignalAxes=axes;
    end
TFRAxes=axes;
else,
    if nejes>=1,
        TFRAxes=fignum(1);
    end
    if nejes>=2,
        SignalAxes=fignum(2);
    end
end


if strcmp('SC2',TFR.type) & (nejes==0 | nejes>=3),
    if nejes==0,PsiAxes=axes;else,PsiAxes=fignum(3);end
    xlong=0.775*length(TFR.psi)/length(TFR.signal.x);
    set(PsiAxes,'Position',[.13 .85 xlong .1])
    plot((1:length(TFR.psi))/TFR.signal.fs,TFR.psi)
    title('Mother Wavalet (psi)')
    set(gca,'XLim',[1 length(TFR.psi)]/TFR.signal.fs)
end
if  ~strcmp('Cohen',TFR.type),
    if nejes==0 & isfield(TFR,'signal'),
        set(SignalAxes,'Position',[.13 .71 0.775 .2])
        set(TFRAxes,'Position',[.13 .06 0.775 .55])
    end
else,
    if nejes==0,
        kernelAxes=axes;
        ambAxes=axes;
        ambfAxes=axes;
        set(kernelAxes,'Position',[.675 .38 0.3 .23])
        set(ambAxes,'Position',[.675 .7 0.3 .23])
        set(ambfAxes,'Position',[.675 .06 0.3 .23])
        set(SignalAxes,'Position',[.13 .71 0.5 .2])
        set(TFRAxes,'Position',[.13 .06 0.5 .55])
    else,
        if nejes>=3,ambfAxes=fignum(3);end
        if nejes>=4,ambAxes=fignum(4);end
        if nejes>=5,kernelAxes=fignum(5);end
    end
end
if (nejes==0 | nejes>=2) & isfield(TFR,'signal'),
% Signal axes
set(gcf,'CurrentAxes',SignalAxes)
t=(1:length(TFR.signal.x))/TFR.signal.fs;
plot(t,TFR.signal.x)
set(gca,'XLim',[t(1) t(end)])

titulo=title(['Analysed Signal; fs=' num2str(TFR.signal.fs) ' Hz']);
ejex=xlabel('Time (s)');
ejey=ylabel('Amplitude (V)');
ticksX=get(gca,'XTick');
end

% TFR axes
set(gcf,'CurrentAxes',TFRAxes)

% We plot only the frequencies indicated in freqband

selectf=(TFR.f>=freqband(1) & TFR.f<=freqband(2));
TFR.f=TFR.f(selectf);
TFR.TFR=TFR.TFR(selectf,:);

if imagecontour==1,
    imagesc(TFR.t,TFR.f,TFR.TFR)
elseif imagecontour==2,
    contour(TFR.t,TFR.f,TFR.TFR)
elseif imagecontour==3,
    meshc(TFR.t,TFR.f,TFR.TFR)
end
if putcolorbar==1,
    colorbar
end
set(gca,'YDir','normal')
if isfield(TFR,'signal')
set(gca,'XLim',[t(1) t(end)])
end

if strcmp('-',TFR.type),
    title('Time-Frequency representation (TFR)')
elseif strcmp('SP',TFR.type),
    if ~isfield(TFR,'window'),TFR.window=0;end
    if TFR.window==1,window='Hanning';
    elseif TFR.window==2,window='rectangular';
    elseif TFR.window==3,window='Hanning';
    elseif TFR.window==0,window='Unknown';
    end
    title(['Spectrogram (V^2/Hz), NFFT=' num2str(TFR.NFFT) ', Window ' window ', Nwindow=' num2str(TFR.Nwindow) ', Noverlap=' num2str(TFR.Noverlap)])
elseif strcmp('AR',TFR.type),
    title(['AR spectrogram (V^2/Hz), Order=' num2str(TFR.Order)  ', Nwindow=' num2str(TFR.Nwindow) ])
elseif strcmp('SC',TFR.type),
    if TFR.envelope==1,window='Exponential';
    elseif TFR.envelope==2,window='Hanning';
    elseif TFR.envelope==3,window='Blackman';
    elseif TFR.envelope==4,window='Rectangular';
    end
    if TFR.envelope==-1,
        title(['Scalogram (V^2/Hz), N=' num2str(TFR.N) ', Gaussian derivative wavelet of order ' num2str(TFR.k) ])
    else,
        title(['Scalogram (V^2/Hz), N=' num2str(TFR.N) ', ' window ' Envelope, k=' num2str(TFR.k)])
    end
elseif strcmp('SC2',TFR.type),
    title(['Scalogram (V^2/Hz), N=' num2str(TFR.N) ', wavelet psi'])
elseif strcmp('WV',TFR.type),
    title(['Wigner-Ville Dist., NFFT=' num2str(TFR.NFFT) ])
elseif strcmp('Cohen',TFR.type),
    texto=[strtrim(TFR.kernel.parameters.name) ' distribution, NFFT=' num2str(TFR.NFFT)];
    if isfield(TFR.kernel.parameters,'sigma'),
        texto=[texto ', sigma=' num2str(TFR.kernel.parameters.sigma)];
    end
    if isfield(TFR.kernel.parameters,'N'),
        texto=[texto ', N=' num2str(TFR.kernel.parameters.N)];
    end
    if isfield(TFR.kernel.parameters,'M'),
        texto=[texto ', M=' num2str(TFR.kernel.parameters.M)];
    end
    if isfield(TFR.kernel.parameters,'beta'),
        texto=[texto ', beta=' num2str(TFR.kernel.parameters.beta)];
    end
    if isfield(TFR.kernel.parameters,'factor_nu'),
        texto=[texto ', factor_n_u=' num2str(TFR.kernel.parameters.factor_nu)];
    end
    title(texto)
elseif strcmp('CHW',TFR.type),
    title(['Choi-Williams Dist, NFFT=' num2str(TFR.NFFT) ', Sigma=' num2str(TFR.sigma) ])
elseif strcmp('GED',TFR.type),
    title(['Generalized Exponential Dist., NFFT=' num2str(TFR.NFFT) ', Sigma=' num2str(TFR.sigma) ', N=' num2str(TFR.N) ', M=' num2str(TFR.M) ])
elseif strcmp('BJD',TFR.type),
    title([TFR.type ', NFFT=' num2str(TFR.NFFT)])
elseif strcmp('CKD',TFR.type),
    title([TFR.type ', NFFT=' num2str(TFR.NFFT) ', Beta=' num2str(TFR.beta) ])
end
ejeTFRx=xlabel('Time (s)');
ejeTFRy=ylabel('Frequency (Hz)');

if  strcmp('Cohen',TFR.type),

    if nejes==0 | nejes>=5,
    % Kernel axes
    set(gcf,'CurrentAxes',kernelAxes)
    if imagecontour==1,
        imagesc(TFR.A.tau,TFR.A.nu,TFR.kernel.kernel')
    elseif imagecontour==2,
        contour(TFR.A.tau,TFR.A.nu,TFR.kernel.kernel')
    elseif imagecontour==3,
        meshc(TFR.A.tau,TFR.A.nu,TFR.kernel.kernel')
    end
    if putcolorbar==1,
        colorbar
    end
    title('Filtering kernel Phi(tau,nu)')
    xlabel('Tau (s)');
    ylabel('Nu (Hz)');
    %set(gca,'XLim',[-t(end) t(end)])
    %set(gca,'YLim',[-TFR.freqband(2) TFR.freqband(2)])
    set(gca,'XLim',TFR.A.zoomtaunu(1)*[-1 1])
    set(gca,'YLim',TFR.A.zoomtaunu(2)*[-1 1])
    end

    if nejes==0 | nejes>=4,
    % Ambiguity axes
    set(gcf,'CurrentAxes',ambAxes)
    if imagecontour==1,
        imagesc(TFR.A.tau,TFR.A.nu,abs(TFR.A.A'))
    elseif imagecontour==2,
        contour(TFR.A.tau,TFR.A.nu,abs(TFR.A.A'))
    elseif imagecontour==3,
        meshc(TFR.A.tau,TFR.A.nu,abs(TFR.A.A'))
    end
    if putcolorbar==1,
        colorbar
    end
    title('Ambiguity function A_x_x(tau,nu)')
    xlabel('Tau (s)');
    ylabel('Nu (Hz)');
    %set(gca,'XLim',[-t(end) t(end)])
    %set(gca,'YLim',[-TFR.freqband(2) TFR.freqband(2)])
    set(gca,'XLim',TFR.A.zoomtaunu(1)*[-1 1])
    set(gca,'YLim',TFR.A.zoomtaunu(2)*[-1 1])
    end

    if nejes==0 | nejes>=3,
    % Filtered ambiguity axes
    set(gcf,'CurrentAxes',ambfAxes)
    if imagecontour==1,
        imagesc(TFR.A.tau,TFR.A.nu,abs(TFR.kernel.kernel'.*TFR.A.A'))
    elseif imagecontour==2,
        contour(TFR.A.tau,TFR.A.nu,abs(TFR.kernel.kernel'.*TFR.A.A'))
    elseif imagecontour==3,
        meshc(TFR.A.tau,TFR.A.nu,abs(TFR.kernel.kernel'.*TFR.A.A'))
    end
    if putcolorbar==1,
        colorbar
    end
    title('Filtered ambiguity function Phi(tau,nu)*A_x_x(tau,nu)')
    xlabel('Tau (s)');
    ylabel('Nu (Hz)');
    %set(gca,'XLim',[-t(end) t(end)])
    %set(gca,'YLim',[-TFR.freqband(2) TFR.freqband(2)])
    set(gca,'XLim',TFR.A.zoomtaunu(1)*[-1 1])
    set(gca,'YLim',TFR.A.zoomtaunu(2)*[-1 1])
    end
end