function [weights, initial_selec, p] = present_stimulus(stimuli, weights, p, trial)
% present_stimulus -- present entire network (rat) with a single stimulus,
% and update weights accordingly

% called by: model
% calls: findWinningNode, calc_selectivity, calc_neigh

% input
%   stim: stimulus to be trained upon
%   weights: weight matrix for network
%   p: experimental structure
%   trial: current trial (either experimental or test)

% output:
%   weights: updated weights
%   initial_selec: selectivity of network (divided into each layer and each
%      grid) of model for stim. To be given to calc_recognition. Calculated
%      with logisitic activation
%   p: output because it contains, now, p.winning.


%% prelim

% flag to indicate whether weights will be updated with experimental eta or
% with delay cycle eta
delay = 0;

% stim is just a vector of each feature. inp_mat tiles that stim into the
% same size as the overall set of weights
grandStim = repmat(reshape(stimuli,[1 1 length(stimuli)]), [p.nRows p.nRows 1]);

% storage for selectivity. To be fed into calc_recognition
initial_selec = zeros(p.nLayers,max(p.nGrids));

acts = zeros(p.nLayers,p.nRows,p.nRows,max(p.nGrids));

%% Expose network to stimuli and update weights

for layer=1:p.layer
    
    % defines which feature are attended to by each grid in a layer
    firstFeatureToCheck=(1:p.nInputDims(layer):p.nGrids(layer)*p.nInputDims(layer));
    lastFeatureToCheck = firstFeatureToCheck+p.nInputDims(layer)-1;
    
    % very probable that I don't need to be calculating each grid
    % individually like this. Could probably update them all at once...
    for grid = 1:p.nGrids(layer),
        
        % broken down by features attended to, this is the actual aspect of
        % the stim that is looked at by a given grid. For PRC layer, this
        % will be entire stimulus.
        inp = grandStim(:,:,(firstFeatureToCheck(grid):lastFeatureToCheck(grid)));
        
        % pull out just the weights of 1 single grid.
        w = squeeze(weights(layer,:,:,1:p.nInputDims(layer),grid));
        
        %% update weights on new stimuli for initial calc of selectivity
        
        
        for cycle=1:p.nEncodingCycles
            
            
            
            % for this stim, find the given winning node, and calculate
            % every node's distance away from the input in mse
            [win_row, win_col, dist_mse] = findWinningNode(w, inp);
            
            
            % for subsequent analysis, making sure that different nodes are
            % being declared as winners.
            p.winning(layer,grid,trial,:) = [win_row, win_col];
            
            %--------------------------------------------------------------
            % grab selectivity of this grid to stim, and calc the amount
            % that needs to be updated.
            %--------------------------------------------------------------
            if cycle == p.nEncodingCycles
                [selectivity, ~, acts_grid] = ...
                    calc_selectivity(win_row, win_col, dist_mse, p);
                acts(layer,:,:,grid) = acts_grid;
            else
                [selectivity, ~, ~] = ...
                    calc_selectivity(win_row, win_col, dist_mse, p);
            end
            
            % because we don't update weights during the choice phase, grab
            % the selectivity only with the fresh weights.
            if cycle==1,
                initial_selec(layer,grid) = selectivity;
            end
            
            %--------------------------------------------------------------
            % Calculate each unit's distance from winner and resultant activation
            % f is gaussian neighborhood.
            %--------------------------------------------------------------
            f = calc_neigh(win_row, win_col, layer, p, delay);
            
            % Update Weights
            w = w + f.*(inp-w);
            
        end
        
        % load updated weights back into full weight matrix.
        weights(layer,:,:,1:p.nInputDims(layer),grid) = w;
        
        
    end
end

end