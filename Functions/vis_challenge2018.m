function vis_challenge2018(data)

n_samples=data.n_samples;
stages=data.sleepstages;
stages_names=data.sleepstages_names;
subject=data.subject;
fs=data.fs;

APPLICATION=['Sleep Stages and Arousals (subject ' subject ')'];
[existFlag,figNumber]=figflag(APPLICATION);
if ~existFlag,  
     figNumber=figure('Name',APPLICATION,'NumberTitle','off');
else,
    figure(figNumber),clf
end
screensize = get(groot, 'ScreenSize');
sizeX=fix(.95*screensize(3));
sizeY=fix(.45*screensize(4));
position = [ceil((screensize(3)-sizeX)/2) ceil((screensize(4)-sizeY)/2), sizeX, sizeY];
set(gcf,'Position',position)


subplot 211
imagesc((0:n_samples-1)/fs/3600,1:length(stages_names),stages)
set(gca,'YTick',1:length(stages_names),'YTickLabel',stages_names,'FontSize',14)
xlabel('Time (hours)')
title(['Sleep stages of subject ' subject])

% pause
% 
% length(stages_names)

mapcolor=zeros(256,3);
colortrue=[.35 0.85 0.25];
colorfalse=[1 0.93 0.73];
mapcolor(1:128,:)=ones(128,1)*colorfalse;
mapcolor(129:256,:)=ones(128,1)*colortrue;
colormap(mapcolor)
handle_colorbar=colorbar;
set(handle_colorbar,'Ticks',[0.25 0.75])
set(handle_colorbar,'Limits',[0 1])
set(handle_colorbar,'TickLabels',{'False','True'})

stages=data.arousals;
stages_names=data.arousals_names;
for i=1:length(stages_names),
stages_names{i}(stages_names{i}=='_')=' ';
end



subplot 212
imagesc((0:n_samples-1)/fs/3600,1:length(stages_names),stages)
set(gca,'YTick',1:length(stages_names),'YTickLabel',stages_names,'FontSize',14)
xlabel('Time (hours)')
title(['Arousals of subject ' subject])
mapcolor=zeros(256,3);
colortrue=[.35 0.85 0.25];
colorfalse=[1 0.93 0.73];
mapcolor(1:128,:)=ones(128,1)*colorfalse;
mapcolor(129:256,:)=ones(128,1)*colortrue;
colormap(mapcolor)
handle_colorbar=colorbar;
set(handle_colorbar,'Ticks',[0.25 0.75])
set(handle_colorbar,'Limits',[0 1])
set(handle_colorbar,'TickLabels',{'False','True'})
end
