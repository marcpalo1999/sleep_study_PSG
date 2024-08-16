function [Pxx_rel,Pxx_M]=relative_power_epoch(eeg,t_eeg,fwave, fs, epoch_time)
    %Property: Marc Palomer
    N = length(eeg)-1;
    Pxx_M=zeros(1025,1);
    Pxx_rel=zeros();
    
    
    
    for n_epoch = 1:N/fs/epoch_time %360 epochs
        eeg_epoch = eeg(t_eeg>=((n_epoch-1)*epoch_time+t_eeg(1)) & t_eeg<=(n_epoch*epoch_time)+t_eeg(1));
    
        % Calculate PSD parameters
        L = length(eeg_epoch);
        NFFT=2^max([10 nextpow2(L)]);
        D=1*fs;
        S=0.5*D;
    
        [Pxx,f]=pwelch(detrend(eeg_epoch),D,S,NFFT,fs,'onesided'); %Power per frequency
        Pxx_M(:,n_epoch)=Pxx;
    
        [~,Pmband2Pm_delta,~]=powerband(Pxx',f',fwave,0);
    
    
        Pxx_rel(:,n_epoch)=Pmband2Pm_delta;
    end

end