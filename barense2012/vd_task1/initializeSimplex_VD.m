% initialize fminsearchbnd (simplex with bounded parameter estimation)
% search of parameters for VD model
startTime = GetSecs;

global consts;

% set seed based on internal clock time
consts.seed = fix(1e6*sum(clock));

% need to make sure am estimating enough rats on each iteration of the
% model to acocunt for noise of output.
consts.nRats = 50;
consts.nIterations = 1000;

% name of folder containing stimuli to use
consts.exptName = '2jul2016';
% consts.nameOfFolder = 'fminsearchbnd_tempSessionInfo';

% data is in DIFFERENCE of dPrime, from first to second half of trials.
% order of [caudal_LA, caudal_HA, PRC_LA, PRC_HA]
% want data: caudal LA => no diff
%            caudal HA => dPrime at 2, then 0 (diff of 2?)
%            PRC LA => no diff
%            PRC HA => no diff     
data = [0, -2, 0, 0];

% initialize starting value of parameters
% etaExp
% k 
%%% numTrainCycles
%%% numEncodingCycles
startParms = [.51, .25];
consts.minParms = [.01, .01,];
consts.maxParms = [1, 2];
% NOTE: to deal with the need for integer parameters, 'nums' are rounded
% inside create_sim

% call wrapper function of fminsearchbnd
[finalParms, fVal] = wrapper4fminbnd(startParms, data);

totalTime = GetSecs - startTime;
fprintf('\n\n %d. \r', totalTime)

%print final parms
save finalParms.mat finalParms
fprintf('\neta: %f, k: %d\n',finalParms(1), finalParms(2))