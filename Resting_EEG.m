%% Settings
I_GenSettings;
Instructions;

%Trials
NrOfTrials = 7;   % How many Cycles to run (8 if  you want to run 6 cycles)
eyeO = 3:60:303; % Audio cues
eyeC = 23:60:323;

% Setting the Audiofiles
wavfilename_probe1 = '/home/stimuluspc/Neurometric/resting/open_eyes.wav'; %Open Eyes
wavfilename_probe2 = '/home/stimuluspc/Neurometric/resting/close_eyes.wav';%Close Eyes

% Setting the Trigger codes
par.CD_START  = 90;
par.CD_eyeO = 20;
par.CD_eyeC = 30;
par.CD_END  = 80;

%% Setting Par
if testmode == 1 % Testmode, modify to your liking
    par.runID = 1;
    par.recordEEG = 0;
    par.useEL = 0;
    par.useEL_Calib = 0;
else
    I_SetPar; % Setting the par file according to call
end;

par.runID= app.run{1};
par.ExaminationDate=app.run{2};
savePath = app.savePath;

%% Screen Calculations
[scresw, scresh]=Screen('WindowSize',whichScreen);  % Get screen resolution
center = [scresw scresh]/2;     % useful to have the pixel coordinates of the very center of the screen (usually where you have someone fixate)
fixRect = [center-2 center+2];  % fixation dot
hz=Screen('FrameRate', whichScreen,1);
cm2px = scresw/monitorwidth_cm;     % multiplication factor to convert cm to pixels
deg2px = dist_cm*cm2px*pi/180;      % multiplication factor to convert degrees to pixels (uses aproximation tanT ~= T).
load gammafnCRT;   % load the gamma function parameters for this monitor - or some other CRT and hope they're similar! (none of our questions rely on precise quantification of physical contrast)
maxLum = GrayLevel2Lum(255,Cg,gam,b0);
par.BGcolor=Lum2GrayLevel(maxLum/2,Cg,gam,b0);

%% Connect Eyetracker & Calibrate
if par.useEL
    window=Screen('OpenWindow', whichScreen, par.BGcolor);
    EL_Connect; %Connect the Eytracker, it needs a window
    try % open file to record data to
        edfFile=[num2str(par.runID),'_RES.edf'];
        Eyelink('Openfile', edfFile);
    catch
        disp('Error creating the file on Tracker');
    end;
    
    if par.useEL_Calib
        EL_Calibrate
    end; %If needed, run Calibration
    Eyelink('command', 'record_status_message "RESTING EEG"');
else
    window=Screen('OpenWindow', whichScreen, par.BGcolor);
    disp('No Eyetracker');
end;

%% Initiate NetStation Connection, Synchronization, and Recording
if par.recordEEG
    %[status,info] = NetStation('Connect','100.1.1.3',55513);
    NetStation('Synchronize');
    WaitSecs(1);
%     if status ~= 0
%         error(info);
%     end
else
    disp('No EEG');
end

%% Sound Stuff
%dev = PsychPortAudio('GetDevices')
%count = PsychPortAudio('GetOpenDeviceCount')
try
    PsychPortAudio('Close');
catch
end
try
    [y_probe1, freq1] = audioread(wavfilename_probe1);
    [y_probe2, freq2] = audioread(wavfilename_probe2);
    wavedata_probe1 = y_probe1';
    wavedata_probe2 = y_probe2';
    nrchannels = size(wavedata_probe1,1); % Number of rows == number of channels.
    % Add 15 msecs latency on Windows, to protect against shoddy drivers:
    sugLat = [];
    if IsWin
        sugLat = 0.015;
    end
    try
        InitializePsychSound;
        pahandle = PsychPortAudio('Open', [], [], 0, freq1, nrchannels, [], sugLat);
        duration_probe1 = size(wavedata_probe1,2)/freq1;
        duration_probe2 = size(wavedata_probe2,2)/freq1;
    catch
        error('Sound Initialisation Error');
    end
catch
    error('Sound Error');
end

i = 1;
t = 1;
tt = 1;

%% Experiment Window
clc;
if ~Screen(window, 'WindowKind') %Check if Window is open
    window=Screen('OpenWindow', whichScreen, par.BGcolor);
end

Screen('TextSize', window, tSize2);
DrawFormattedText(window, ins.resting.inst, 0.08*scresw, 0.25*scresh, colorText);
DrawFormattedText(window, ins.misc.mouse,'center', 0.9*scresh, colorText);
Screen('Flip', window);
HideCursor(whichScreen);

clc;
disp('THE SUBJECT IS READING THE INSTRUCTIONS');
if par.useEL 
    Eyelink('StartRecording');
    Eyelink('command', 'record_status_message "READING INSTRUCTIONS"');
end;
if par.recordEEG, NetStation('StartRecording'); end;
disp('DID THE EEG START RECORDING DATA? IF NOT, PRESS THE RECORD BUTTON!');
[clicks,x,y,whichButton] = GetClicks(whichScreen,0); % Waits for the user to press a button before starting

%% Experiment Block
time = GetSecs;
if par.useEL
    Eyelink('Message',['TR',num2str(par.CD_START)]);
    Eyelink('command', 'record_status_message "EXPERIMENT STARTED"');
end; %Sending Starttrigger
if par.recordEEG,  NetStation('Event', num2str(par.CD_START)); end

%while ~KbCheck
if testmode == 1
    trials = 4; % Test Mode
else
    trials = NrOfTrials;
end
fprintf('Running Trials\n');
while t < trials
    Screen('DrawLine', window,[0 0 0],center(1)-7,center(2), center(1)+7,center(2));
    Screen('DrawLine', window,[0 0 0],center(1),center(2)-7, center(1),center(2)+7);
    vbl = Screen('Flip',window); % clc
    if vbl >=time+eyeO(t) %Tests if a second has passed
        
        if par.recordEEG,  NetStation('Event', num2str(par.CD_eyeO(1))); end
        if par.useEL
            Eyelink('Message',['TR',num2str(par.CD_eyeO)]);
            Eyelink('command', 'record_status_message "EYES OPEN"');
        end; %Send Eytracker Trigger
        disp('Eyes Open');
        
        PsychPortAudio('FillBuffer', pahandle,wavedata_probe1);
        PsychPortAudio('Start', pahandle, 1, 0, 1);
        t = t+1;
    end
    
    if vbl >=time+eyeC(tt) %Tests if a second has passed
        if par.recordEEG,  NetStation('Event', num2str(par.CD_eyeC(1))); end
        if par.useEL
            Eyelink('Message',['TR',num2str(par.CD_eyeC)]);
            Eyelink('command', 'record_status_message "EYES CLOSED"');
        end; %Send Exetracker Trigger
        disp('Eyes Closed');
        PsychPortAudio('FillBuffer', pahandle,wavedata_probe2);
        PsychPortAudio('Start', pahandle, 1, 0, 1);
        tt = tt+1;
    end
end

if par.useEL, Eyelink('Message',['TR',num2str(par.CD_END)]); end; %Send Exetracker END Trigger
if par.recordEEG,  NetStation('Event', num2str(par.CD_END)); end

disp('Done');
Screen('TextSize', window, tSize3);
DrawFormattedText(window, ins.misc.finished,'center', 0.4*scresh, colorText);
Screen('TextSize', window, tSize2);
DrawFormattedText(window, ins.resting.end,'center', 0.5*scresh, colorText);
Screen('Flip', window);
ShowCursor(whichScreen);
WaitSecs(5);

if par.recordEEG
    fprintf('Stop Recording EEG\n');
    NetStation('StopRecording'); %Stop Recording
    %NetStation('Disconnect');
end

if par.useEL
    fprintf('Stop Recording Track\n');
    Eyelink('StopRecording'); %Stop Recording
    Eyelink('CloseFile');
    fprintf('Downloading File\n');
    EL_DownloadDataFile % Downloading the file
    EL_Cleanup %Shutdown Eyetracker and close all Screens
end
save([savePath,num2str(par.runID),'_RES'],'par','eyeO','eyeC');
clearvars -except select subj_ID metafile subj_Name app event
sca; %If Eyetracker wasn't used, close the Screens now
try
    PsychPortAudio('Close');
catch
end

