function [ x,fVal ] = wrapper4fminbnd( pArray, data, nOfFolder )
%wrapper4min out wrapper for call to simplex algorithm on VD model

% pArray = parameter array, given by initializeSimplex_VD
% data will be the dPrime differenes between first and second half of
%   trials, as given by Barnese et al 2012.

% structure consts initialized mainly inside initializeSimplex_VD
global consts;

% rand('state', consts.seed)

defOpts = optimset('fminsearch');
options = optimset(defOpts, 'Display', 'iter', 'MaxFunEvals', consts.nIterations);

% call fminsearchbnd with min and max allowed parameters
[x, fVal] = fminsearchbnd(@bof,pArray,consts.minParms, consts.maxParms, options);

    function rmsd=bof(parms)
        
        % predictions will be related to dPrime, data will be desired
        % dPrime diffs to model.
        predictions = getdPrimePred(parms, nOfFolder)
        sqdist = (predictions - data).^2;
        rmsd=sqrt(mean(sqdist));
        fprintf('\n\n rmsd:%d., eta:%d, \r', rmsd, parms)
    end

end

