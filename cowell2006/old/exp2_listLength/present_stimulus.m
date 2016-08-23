function [W, initial_selec, p] = present_stimulus(stim, W, p, trial)
% present_stimulus -- present entire network (rat) with a single stimulus,
% and update weights accordingly

% called by: model
% calls: findWinningNode, calc_selectivity, calc_neigh

% input
%   stim: stimulus to be trained upon
%   W: weight matrix for network
%   p: experimental structure
%   trial: current trial (either experimental or test)

% output:
%   W: updated weights
%   initial_selec: selectivity of network (divided into each layer and each
%      grid) of model for stim. To be given to calc_recognition. Calculated
%      with logisitic activation
%   p: output because it contains, now, p.winning.


%% prelim

% stim is just a vector of each feature. inp_mat tiles that stim into the
% same size as the overall set of weights
inp_mat = repmat(reshape(stim,[1 1 length(stim)]), [p.nRows p.nRows 1]);

% storage for selectivity. To be fed into calc_recognition
initial_selec = zeros(p.numLayers,max(p.numGrids));
initial_selec_gauss = zeros(size(initial_selec));


%% Expose network to stimuli and update weights

for layer=1:p.layer
    
    % defines which feature are attended to by each grid in a layer
    firstFeatureToCheck=(1:p.numInputDims(layer):p.nGrids(layer)*p.numInputDims(layer));
    lastFeatureToCheck = firstFeatureToCheck+p.numInputDims(layer)-1;
    
    % very probable that I don't need to be calculating each grid
    % individually like this. Could probably update them all at once...
    for grid = 1:p.nGrids(layer),
        
        % broken down by features attended to, this is the actual aspect of
        % the stim that is looked at by a given grid. For PRC layer, this
        % will be entire stimulus.
        input_mat=inp_mat(:,:,(firstFeatureToCheck(grid):lastFeatureToCheck(grid)));
        
        % pull out just the weights of 1 single grid. 
        weights = squeeze(W(layer,:,:,1:p.numInputDims(layer),grid));
        
        %% update weights on new stimuli for initial calc of selectivity
        
        
        for cycle=1:p.numEncodingCycles
            
            % for this stim, find the given winning node, and calculate
            % every node's distance away from the input in mse
            [win_row, win_col, dist_mat] = findWinningNode(weights, input_mat);
            
            % for subsequent analysis, making sure that different nodes are
            % being declared as winners.
            p.winning(layer,grid,trial,1:2) = [win_row, win_col];
            
            %--------------------------------------------------------------
            % grab selectivity of this grid to stim, and calc the amount
            % that needs to be updated.
            %--------------------------------------------------------------
            [selectivity, ~, ~] = ...
                calc_selectivity(win_row, win_col, dist_mat, p);
            
            % because we don't update weights during the choice phase, grab
            % the selectivity only with the fresh weights.
            if cycle==1,  
                initial_selec(layer,grid) = selectivity;
            end
            
            
            %--------------------------------------------------------------
            % Calculate each unit's distance from winner and resultant activation
            % f is gaussian neighborhood.
            %--------------------------------------------------------------
            f = calc_neigh(win_row, win_col,layer,p);
            
            % Update Weights
            weights = weights + f.*(input_mat-weights); 
            
        end 
        
        % load updated weights back into full weight matrix.
        W(layer,:,:,1:p.numInputDims(layer),grid) = weights;
        
        
    end
end

end