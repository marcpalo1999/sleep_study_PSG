function kernel=TFRkernel(parameters,fignum)

% Function that computes kernels suitable for multiplying the 
% ambiguity function and compute a Cohen class TFR.
%
% kernel=TFRkernel(parameters,fignum)
%
% Inputs:
%
% parameters: Structure with the following fields:
% - type: Type of kernel (0: WV, 1: Choi-Williams, 2: Generalized exponential
% generalized, 3: Born-Jordan distribution, 4: Cone Kernel Distribution,
%  5: exponential)
% - NFFT: Dimension of the kernel (NFFT*NFFT)
% - sigma: Used in kernels type 1, 2 or 5
% - N: Used in kernel type 2
% - M: Used in kernel type 2
% - beta: Used in kernel type 4
% - factor_nu: Used in kernel of type 5
%
% fignum: figure where is displayed the result (0, do not display, <0, create new figure)
%
% Output: 
% kernel: Structure with the following fields:
% - kernel: generated kernel (NFFT*NFFT)
% - type: Type of kernel (0: WV, 1: Choi-Williams, 2: Generalized exponential
% generalized, 3: Born-Jordan distribution, 4: Cone Kernel Distribution,
%  5: exponential)
% - NFFT: Dimension of the kernel (NFFT*NFFT)
% - sigma: Used in kernels type 1, 2 or 5
% - N: Used in kernel type 2
% - M: Used in kernel type 2
% - beta: Used in kernel type 4
% - factor_nu: Used in kernel of type 5
%
% Example:
%
% parameters.type=1;
% parameters.sigma=0.1;
% parameters.NFFT=1024;
% kernel=TFRkernel(-1,par);
%
%
% Abel Torres (abel.torres@upc.edu), IBEC-ESAII-UPC

if nargin<1, % By default a Choi-Williams kernel with NFFT=1024 and sigma=0.1 is generated and displayed
    parameters.type=1;
    parameters.sigma=0.1;
    parameters.NFFT=1024;
    kernel=TFRkernel(parameters,-1);
end
if nargin<2,
    fignum=-1;
end
name=str2mat('Wigner-Ville',...
    'Choi-Williams',...
    'Generalized exponential',...
    'Born-Jordan',...
    'Cone-Kernel', ...
    'Exponential');
parameters.name=name(parameters.type+1,:);

kernel.kernel=-parameters.NFFT/2:(parameters.NFFT/2-1);

if parameters.type==0, % Dist WV
  kernel.kernel=ones(parameters.NFFT,parameters.NFFT);
elseif parameters.type==1,
  % Choi Williams 
  kernel.kernel=[kernel.kernel'*kernel.kernel]/parameters.NFFT;
  kernel.kernel=exp(-(kernel.kernel.^2)/(parameters.sigma));
elseif parameters.type==2,
  % Generalyzed exponencial
  tau=(kernel.kernel/parameters.NFFT).^(2*parameters.M);
  nu=(kernel.kernel).^(2*parameters.N);
  kernel.kernel=exp(-tau'*nu/parameters.sigma);
elseif parameters.type==3 | parameters.type==4,
  if parameters.type==4,
    kernel_CKD=abs(kernel.kernel).*exp(-parameters.beta*(kernel.kernel/(parameters.NFFT/2)).^2);
  end 
  % Born-Jordan:
  kernel.kernel=sinc([kernel.kernel'*kernel.kernel]/parameters.NFFT);
  if parameters.type==4,
    % Cone Kernel
    kernel_CKD=kernel_CKD'*ones(size(kernel_CKD));
    kernel.kernel=kernel_CKD.*kernel.kernel;
  end
elseif type==5,
  % Exponential
  kernel.kernel=kernel.kernel/(parameters.NFFT/2);
  kernel=exp(-kernel.^2/parameters.sigma)'*exp(-kernel.^2/parameters.nu_factor/parameters.sigma);
end

kernel.parameters=parameters;

if fignum,
  if fignum>0,figure(fignum),else,figure,end
  
    fs=1;
  
    tau=(-parameters.NFFT/2:parameters.NFFT/2-1)/(fs/2);
    nu=(-parameters.NFFT/2:parameters.NFFT/2-1)/(parameters.NFFT/2)*(fs/2);
  
    imagesc(tau,nu,kernel.kernel')

    parametersnames=fieldnames(parameters);
    parametersnames=parametersnames(~(strcmpi('type',parametersnames) | strcmpi('name',parametersnames)));
    texto=[strtrim(parameters.name) ' distribution' ];
    for i=1:length(parametersnames)
        texto=[texto ', ' char(parametersnames(i)) '=' num2str(getfield(parameters,char(parametersnames(i))))];
    end
    texto=[texto ' (normalized frequency:edi fs=1 Hz)'];
    title(texto)	
    xlabel('Tau (s)');
    ylabel('Nu (Hz)');
  
end


