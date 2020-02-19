%-------------------------------------------------------------------------
% Script to run multiple experimental scripts in a row
% Programmed by Kirsten Adam, June 2014
%-------------------------------------------------------------------------
%clear all;  % clear everything out!
%close all;  % close existing figures

tic;

% ONLY FOR DEBUGGING PURPOSES - DELETE OR COMMENT ME:
testmode=0;
if testmode
app.run{1}='M01';
app.run{2}='today';
app.savePath=[pwd filesep];
addpath(genpath(pwd));
end
% END DEBUGGING STUFF

I_GenSettings;
Instructions;
prefs = getPreferences();  % function that grabs all of our preferences (at the bottom of this script)

par.runID= app.ID{1};
par.ExaminationDate =date;
par.block=prefs.numBlocks;

% -------------------------------------------------------------------------
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

p.subNum = par.runID;
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
win= openWindow(p,whichScreen);
par.whichScreen=whichScreen;

%% Initiate NetStation Connection, Synchronization, and Recording
%-------------------------------------------------------------------------
% Run Experiment 1
%-------------------------------------------------------------------------
ChangeDetection_Color_Function_9colors_Example(p,win,prefs,par,ins,par.savePath);

%-------------------------------------------------------------------------
%Stop recording and download ET Data
%-------------------------------------------------------------------------

savePath=par.savePath;

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
prefs.numBlocks = 2; % was 6 / 9 before!
prefs.nTrialsPerCondition = 4 %10;
prefs.setSizes =[2,6]; %[1,6,8]; 
prefs.change = [0,1]; % 0 = no change, 1 = change!

%%%%% timing
prefs.retentionInterval =  [1.000]; % win.refRate;% 1 sec  (or, if we don't do this we can jitter .... )
prefs.stimulusDuration = [.250]; %win.refRate/2;% 500 ms
prefs.ITI = 900.000;  %prefs.retentionInterval;
prefs.breakLength = .1*5; % number of minutes for block
prefs.cueDuration=0.200;
prefs.fixDuration=0.400;
prefs.realITI=1.00;
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


