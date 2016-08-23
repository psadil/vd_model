function [ p, weights ] = interfere( p, weights )
%delay_interfere simulates intereference during delay cycles
%   Detailed explanation goes here

% expt 1 is delay expt
if p.expt == 1
    fprintf('\n%d interference cycles being executed...', p.delayCycles(p.stimCond));
end

delay = 1;

%%
for layer = 1:max(p.nLayers)
    for grid = 1:p.nGrids(layer),
        
        % put the variable 'weights' into the format previously accepted by the model.
        w = squeeze(weights(layer,:,:,1:p.nInputDims(layer),grid));
        
        for cycle=1:p.delayCycles(p.stimCond),
            
            
            inp = gen_randInp(p, p.nInputDims(layer)); %generate an input vector
            
            %--------------------------------------------------------------
            % Find winning node
            %--------------------------------------------------------------
            [win_row, win_col, ~] = findWinningNode(w, inp);
            
            
            %--------------------------------------------------------------
            % Calculate each unit's distance from winner, for use in
            % updating
            f = calc_neigh(win_row, win_col, layer, p, delay);
            
            %%% Update Weights
            w = w + f.*(inp-w);
            
            
        end
        weights(layer,:,:,1:p.nInputDims(layer),grid) = w;
    end
    
    
end

end