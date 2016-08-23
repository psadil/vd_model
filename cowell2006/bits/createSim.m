function [] = createSim(firstRat, lastRat)
% listLength_createSim -- first function to
% NOTE: each iteration of the simulation is refered to as a 'rat.' One rat
% is run for every session, and is run as though it were a lesioned rat
% (which is to say it doesn't have a PRC grid) and a control rat (the rat
% has BOTH a PRC and caudal grid)

% overall, simulation replicates a Kohonen-style, self-organizing map (SOM)


% called by: YOU
% calls: runSim

% input:
%    firstRat: first rat to run (usually 1)
%    lastRat: last rat to run. defines how many rats are being run. Each
%        rat will run for 8 sessions, 4 lesion and 4 control

% output:
%   NA -- however, commented line at bottom (call to 'plotRecognition.m')
%      will generate summary graphs of the simulation


%% define output folder for session info

% request for experiment name
expt = input('\nEnter experiment name: \n1 - delay, \n2 - listLength, \n3 - tUnique \n');


% for cleanliness, expts get stored in this folder
outerDir = strcat(pwd, '/graphsAndSession');
if ~exist(outerDir, 'dir'),
    mkdir(outerDir);
end

if expt == 1
    exptFolder = 'delay';
elseif expt == 2
    exptFolder = 'listLength';
elseif expt == 3
    exptFolder = 'tUnique';
else
    error('not valid experiment');
end
fprintf('\nBeginning experiment, %s.\n', exptFolder);


%%
fprintf('\nSetting new seed by clock.\n');
rng('shuffle');


%% begign running rats
parfor rat = firstRat:lastRat
    
    % initialzie main p, to be used by all expts
    p = init_exptParms(expt);
    p.ratNum = rat;
    
    
    
    %% data directory information
    
    % tweak as necessary to reflect peculiarities in any given simulation
    p.nameOfFolder = ['eta', num2str(p.etaExp), '_g', num2str(p.G_exp), ...
        '_K', num2str(p.k_expt), '_A', num2str(p.A) ,'_B', num2str(p.B),...
        '_',num2str(p.nTrainCycles),'nTrnCyc_', num2str(p.nEncodingCycles),'enc_',...
        num2str(p.nSess),'sess_',num2str(p.nGrids_Caudal),'cGrds',...
        num2str(p.components),'nCmpts_', ...
        num2str(p.nRows),'rws_', 'wCor'];
    
    p.dataDir = strcat(outerDir, '\', exptFolder, '\' ,p.nameOfFolder);
    if ~exist(p.dataDir, 'dir'),
        mkdir(p.dataDir);
    end
    
       
    
    fprintf('\n Creating a new experiment...\n %s \n', p.nameOfFolder);
    
    
    %% enter into runSim, which controls flow of sessions
    
    [~] = runSim(p);
    
end

% call following to make graphs
% plotFamilDiffs(1, lastRat, p.nameOfFolder);

end