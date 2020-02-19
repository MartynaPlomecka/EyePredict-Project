%-------------------------------------------------------------------------
% Script to run multiple experimental scripts in a row
% Programmed by Kirsten Adam, June 2014
%-------------------------------------------------------------------------
%clear all;  % clear everything out!
%close all;  % close existing figures

tic;

% ONLY FOR DEBUGGING PURPOSES - DELETE OR COMMENT ME:
testmode = 0;

if testmode
    app.run{1}='MT2';
    app.run{2}='today';
    app.savePath=[pwd filesep];
    addpath(genpath(pwd));
end
% END DEBUGGING STUFF

I_GenSettings;

testmode=0;
Instructions;
prefs = getPreferences();  % function that grabs all of our preferences (at the bottom of this script)



if testmode
    warning('You are running the task in testmode!');
end

if testmode
    par.block=1;
else
    %par.block = app.run{6};
    par.block=1;
end

if testmode == 1 % Testmode, modify to your liking
    par.runID = 1;
    par.recordEEG = 0;
    par.useEL = 1;
    par.useEL_Calib = 1;
else
    I_SetPar; % Setting the par file according to call
end;

%--------------------------------------------------------------------------
%Trigger definition
%--------------------------------------------------------------------------

switch par.block
    case 1
        par.CD_START1  = 951;
        par.CD_START2  = 952;
        par.CD_START3  = 953;
        par.CD_START4  = 954;
end

par.CD_TRIAL  = 200+[1:prefs.numTrials]; %geht das`100er triggers auch als taskStart verwendet?
%par.CD_MEMSET_START=60;
%CD_MEMSET_START is dynamically!
% can have values 110,111,120,121,610,611,620,621
% first digit encodes setseize (1/6)
% second digit encodes cue position(1=left, 2 =right)
% third digit encodes change (0=nochange, 1=change);
par.CD_MEMSET_END=61;
par.CD_PROBE=80;
par.CD_REINST_START=70;
par.CD_REINST_END=71;
par.CD_BLOCK_START=20;
par.CD_BLOCK_END=21;

par.CD_RESP  = 30;
par.CD_END  = 50;

% -------------------------------------------------------------------------



par.runID= app.run{1};
par.ExaminationDate=app.run{2};
if testmode
    par.savePath=[app.savePath filesep 'testResultsVisualChange' filesep];
else
    par.savePath=[app.savePath];
end

warning('off','MATLAB:dispatcher:InexactMatch');  % turn off the case mismatch warning (it's annoying)
%dbstop if error  % tell us what the error is if there is one
AssertOpenGL;    % make sure openGL rendering is working (aka psychtoolbox is on the path)
%-------------------------------------------------------------------------
% Build a GUI to get subject number
%-------------------------------------------------------------------------
% prompt = {'Subject Number'};            % what information do we want from the subject?
% defAns = {''};                                           % fill in some stock answers - here the fields are left blank
% box = inputdlg(prompt,'Enter Subject Info');       % build the GUI

p.clockOutput = clock; % record time and date!!
p.rndSeed = round(sum(100*p.clockOutput));

p.subNum = app.run{1};
rand('state',p.rndSeed);

%-------------------------------------------------------------------------
% Important options
%-------------------------------------------------------------------------
p.is_PC = ispc; % detects whether this is a PC or windows machine.
p.windowed = 0; % 1 = smaller window for easy debugging!
%-------------------------------------------------------------------------
% Build an output directory & check to make sure it doesn't already exist
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
% Build psychtoolbox window & hide the task bar
%-------------------------------------------------------------------------



%% Connect and Calibrate Eyetracker
load gammafnCRT;   % load the gamma function parameters for this monitor - or some other CRT and hope they're similar! (none of our questions rely on precise quantification of physical contrast)
maxLum = GrayLevel2Lum(255,Cg,gam,b0);
par.BGcolor=Lum2GrayLevel(maxLum/2,Cg,gam,b0);

%delete me later:
%edfFileCell{1}=[num2str(par.runID),'_VM',num2str((par.block-1)*3+1),'.edf'];
disp('connecting ET for the first time');
if par.useEL
    window = Screen('OpenWindow', whichScreen, par.BGcolor);
    EL_Connect; %Connect the Eytracker, it needs a window
    
    try % open file to record data to
        % change here, how many "subplocks per block?? here it
        % works if its 2 or more experimental blocks containng 2
        % subblocks each:
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        edfFileCell{1}=[num2str(par.runID),'_VM',num2str((par.block-1)*2+1),'.edf'];
        Eyelink('Openfile', edfFileCell{1});
    catch
        fprintf('Error creating the file on Tracker\n');
        EL_Cleanup;
    end;
    if par.useEL_Calib
        EL_Calibrate
    end; %If needed, run Calibration
    Eyelink('command', 'record_status_message "VisualChangeDetection"');
else
    edfFileCell=[];
end
win= openWindow(p,whichScreen);
par.whichScreen=whichScreen;

%% Initiate NetStation Connection, Synchronization, and Recording
% if par.recordEEG
%     %try and set up connection to eeg
%     try
%         i = NetStation('Synchronize');
%         if i == 0
%             disp('already connected');
%         else
%             disp('need to connect');
%            [status,info] = NetStation('Connect','100.1.1.3',55513);
%         end
%     catch
%     end
%    
%     WaitSecs(1);
%     if status ~= 0
%         error(info);
%     end
%     NetStation('Synchronize');
% else
%     disp('No EEG');
% end

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

%-------------------------------------------------------------------------
% Run Experiment 1
%-------------------------------------------------------------------------
edfFileCell=ChangeDetection_Color_Function_9colors(p,win,prefs,par,ins,edfFileCell,par.savePath);

%-------------------------------------------------------------------------
%Stop recording and download ET Data
%-------------------------------------------------------------------------
if par.recordEEG,  NetStation('Event', num2str(par.CD_END)); end
if par.useEL, Eyelink('Message',['TR',num2str(par.CD_END)]); end;

if par.recordEEG
    fprintf('Stop Recording EEG\n');
    NetStation('StopRecording'); %Stop Recording
end


savePath=par.savePath;
if par.useEL
    fprintf('Stop Recording\n');
    Eyelink('StopRecording'); %Stop Recording
    Eyelink('CloseFile');
    fprintf('Downloading File\n');
    edfFile=edfFileCell{end};
    EL_DownloadDataFile % Downloading the file
    EL_Cleanup %Shutdown Eyetracker and close all Screens
end



% addpath(genpath('/home/stimuluspc/Tools/Tooboxes/eeglab'))
% pathEdf2Asc = '/home/stimuluspc/Tools/Tools/edf2asc';
% if par.useEL
%     for i=1:length(edfFileCell)
%         system([pathEdf2Asc ' "' [par.savePath, par.runID '_VM' num2str(((par.block-1)*2)+i) ] '" -y'])
%         parseeyelink(strrep([par.savePath, par.runID '_VM' num2str(((par.block-1)*2)+i) '.edf' ],'.edf','.asc'),[strrep([par.savePath, par.runID '_VM' num2str(((par.block-1)*2)+i)],'.edf','') '_ET.mat'],'TR');
%     end
% end
toc;
%-------------------------------------------------------------------------
% Close psychtoolbox window and clear it all out!
%-------------------------------------------------------------------------
sca;
ListenChar(0);
if p.is_PC
    ShowHideWinTaskbarMex(1);
end
close all;
clearvars -except select subj_ID metafile subj_Name app event








%-------------------------------------------------------------------------
%  CHANGE PREFERENCES!
%-------------------------------------------------------------------------
function prefs = getPreferences
%%%% Design conditions
prefs.numBlocks = 4; % was 6 / 9 before!
prefs.nTrialsPerCondition = 12;%12;
prefs.setSizes =[2,6]; %[1,6,8];
prefs.change = [0,1]; % 0 = no change, 1 = change!

%%%%% timing
prefs.retentionInterval =  [1.000]; % win.refRate;% 1 sec  (or, if we don't do this we can jitter .... )
prefs.stimulusDuration = [.500]; %win.refRate/2;% 500 ms
prefs.ITI = 900.000;  %prefs.retentionInterval;
prefs.breakLength = .1*5; % number of minutes for block
prefs.cueDuration=0.200;
prefs.fixDuration=0.400;
prefs.realITI=0.400;
prefs.firstFixationDuration=1;


%%%%% stimulus size & positions
prefs.stimSize = 51/1.57;
prefs.minDist = prefs.stimSize*1.5;
prefs.fixationSize = 6/1.57;

%%%%% randomize trial order of full factorial design order
prefs.fullFactorialDesign = fullfact([length(prefs.setSizes), ...
    length(prefs.change), ...
    length(prefs.retentionInterval), ...
    length(prefs.stimulusDuration), ...
    prefs.nTrialsPerCondition]);  %add prefs.numBlocks? No, because we are using fully counterbalanced blocks.

%%%%% total number of trials in each fully-crossed block.
prefs.numTrials = size(prefs.fullFactorialDesign,1);
end


