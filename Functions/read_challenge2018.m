function data=read_challenge2018(subject)

if nargin<1,subject='tr03-0005';end

% Extract the info from the header file
infosubject=headerinfo([subject '.hea']);

% load all the *.mat data associated with this subject
signals=load(infosubject.signal_location);
signals=signals.val;
arousal=load(infosubject.arousal_location); 
%   undefined: [5147000×1 logical]
%       nonrem3: [5147000×1 logical]
%       nonrem2: [5147000×1 logical]
%       nonrem1: [5147000×1 logical]
%           rem: [5147000×1 logical]
%          wake: [5147000×1 logical]
sleep_stages=arousal.data.sleep_stages;

target_arousals = arousal.data.arousals;


stages_names={'wake','nonrem1','nonrem2','nonrem3','rem','undefined'};
for i=1:length(stages_names),
    eval(['stages(i,:)=sleep_stages.' stages_names{i} ''';'])
end

fs           = str2num(infosubject.fs);
n_samples    = str2num(infosubject.n_samples);
sid          = infosubject.subject_id;
signal_names = infosubject.signal_names;

data.header_text=read_header([subject '.hea']);
data.fs=fs;
data.n_samples=n_samples;
data.subject=sid;
data.signals=signals;
data.signal_names=signal_names;
data.sleepstages=stages;
data.sleepstages_names=stages_names;
data.arousals=[];
data.arousals_names=[];

if exist('rdann')==2, % the rdann funtion belongs to the WFDB Toolbox 
    % Read annotation file:
    [sleep_stages,arousals]=readann_challenge2018(subject);
    data.sleepstages_annotations=sleep_stages;
    data.arousals_annotations=arousals;
    [data.arousals,data.arousals_names]=readarousals(arousals,n_samples);
else,
    disp('Install the WFDB Toolbox for MATLAB to read the *.arousal annotation file')
    disp('(https://archive.physionet.org/physiotools/matlab/wfdb-app-matlab/)')
    disp('Original sleep_stages and arousals annotation  info from this file will not be generated')
    disp('The information in the target_arousals,sleepstages and sleepstages_names fields comes from the *-arousal.mat file generated for the Challenge18')
end

data.target_challenge18_arousals=target_arousals;


eval(['save ' subject 'data data'])

vis_challenge2018(data)


end

function [arousals,arousals_names]=readarousals(arousalsann,n_samples);
    [labels,ia,ic]=unique(arousalsann.label_id);
    n_arousals=length(labels)/2;
    arousals_names=arousalsann.label(ia(1:2:end));
    arousals=logical(zeros(n_arousals,n_samples));
    for i=1:n_arousals,
        arousals_names{i}=arousals_names{i}(2:end);
        start=arousalsann.sample(arousalsann.label_id==labels(2*i-1));
        last=arousalsann.sample(arousalsann.label_id==labels(2*i));
        if start(1)>last(1),start=[1 start];end
        if last(end)<start(end),last=[last n_samples];end
        for k=1:length(start),
            arousals(i,start(k):last(k))=ones(size(start(k):last(k)));
        end
    end
    
end

function infosubject=headerinfo(header_file_name)

% Import information of header file.
% The output of headerinfo is a struct with the following fields:
%     - subject_id: subject identification name
%     - fs: signal sampling rate
%     - n_samples: number of samples
%     - signal_names: cell with the names of the signals
%     - signal_location: location subdirectory of the signals file
%     - arousal_location: location of the arousal annotations (for training set)

% Import the header file
fid = fopen(header_file_name,'rt');
raw_header = textscan(fid,'%s','Delimiter','\n');
raw_header = raw_header{1};
fclose(fid);

% Process the first row of the header file
header_first_row = strsplit(raw_header{1}, ' ');
infosubject.subject_id = header_first_row{1};
infosubjects.n_channels= header_first_row{2};
infosubject.fs = header_first_row{3};
infosubject.n_samples = header_first_row{4};

% Extract the signal names from the remainder of the file.
for j = 2:length(raw_header)
    header_row = strsplit(raw_header{j}, ' ');
    signal_names{j-1} = header_row{end};
end
infosubject.signal_names = signal_names;

% Extract the signal location
[rec_dir, rec_name, ~] = fileparts(header_file_name);
infosubject.signal_location = [rec_name '.mat'];
infosubject.arousal_location = [rec_name '-arousal.mat'];

end




function [sleep_stages,arousals]=readann_challenge2018(subject)

persistent javaWfdbExec config
if(isempty(javaWfdbExec))
    javaWfdbExec=getWfdbClass('rdann');
    [~,config]=wfdbloadlib;
end

if nargin<1,subject='tr03-0005';end
%----------------------------------------------------------
% Reading annotations files (*.arousal) with the wfdb toolbox:
annotation='arousal';
wfdb_argument={'-r',subject,'-a',annotation};
dataJava=javaWfdbExec.execToStringList(wfdb_argument);
data=dataJava.toArray(); % Java array with Nann annotations info

Nann=length(data); % Number of annotations
time_ann0='00:00:00.000';
for i=1:Nann,
    textannotation=data(i,:); %'    6:00.000    72000     "    5    0    0→W'
    aux=textscan(textannotation,'%s %d %s %d %d %d %s',1);
    aux1=char(aux{1});
    timetext(i,1:12)=[time_ann0(1:12-length(aux1)) aux1];
    [~,~,~,HH,MM,SS]=datevec(timetext(i,:),'HH:MM:SS.FFF');
    time(i,1:3)=[HH MM SS]; % Annotation time (hour, minutes, seconds.miliseconds)
    sample(i)=aux{2}+1;
    label_id(i)=aux{4};
    arousal(i)=aux{5};
    label{i}=char(aux{7});
end

% Sleep stages annotations:
index=arousal==0;
sleep_stages.Nann=sum(index);
sleep_stages.time=time(index,1:3);
sleep_stages.timetext=timetext(index,1:12);
sleep_stages.sample=sample(index);
sleep_stages.label_id=label_id(index);
sleep_stages.label=label(index);
sleep_stages.summary=['N.',char(9),'time',char(9),'sample',char(9),'label_id',char(9),'label'];
for i=1:sleep_stages.Nann,
    addtext=[num2str(i),'.',char(9),sleep_stages.timetext(i,:),char(9),num2str(sleep_stages.sample(i)),char(9),num2str(sleep_stages.label_id(i)),char(9),char(sleep_stages.label(i))];
    sleep_stages.summary=str2mat(sleep_stages.summary,addtext);
end


% Arousals annotations
index=arousal==1;
arousals.Nann=sum(index);
arousals.time=time(index,1:3);
arousals.timetext=timetext(index,1:12);
arousals.sample=sample(index);
arousals.label_id=label_id(index);
arousals.label=label(index);
arousals.summary=['N.',char(9),'time',char(9),'sample',char(9),'label_id',char(9),'label'];
for i=1:arousals.Nann,
    addtext=[num2str(i),'.',char(9),arousals.timetext(i,:),char(9),num2str(arousals.sample(i)),char(9),num2str(arousals.label_id(i)),char(9),char(arousals.label(i))];
    arousals.summary=str2mat(arousals.summary,addtext);
end



end


function header_text = read_header(header_file_name)

% This function reads the header file information file of the
% data for the 2018 PhysioNet Challenge. 

if nargin>0,
    fid = fopen(header_file_name,'rt');
else
    fid=-1;
end
if fid>=3,
    header_text = textscan(fid,'%s','Delimiter','\n');
    header_text = char(header_text{1});
    fclose(fid);
    
else,
    disp('Header file not found')
    header_text='';
end

end


