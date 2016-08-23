function [p,weights] = pretrain(p)

fprintf('\ndelay_pretrain being executed...');
interfere = 0;

%%%%%%%%%%%%%%%%%%%%% Perform pre-training of network %%%%%%%%%%%%%%%%
weights = zeros(max(p.numLayers),p.numRows,p.numRows,p.numInputDims(p.numLayers),p.nGrids(1));
for layer = 1:max(p.numLayers)
    
    nInpDims = p.numInputDims(layer);
    
    for grid = 1:p.nGrids(layer),
        fprintf('\n\nLayer number %d, Grid no. %d...\n(1==Caudal, 2==PRC)', layer, grid);
        
        
        % initialize random weights for grids in this layer
        w = init_weights(p, layer);
        
        %------------------------------------------------------------------
        % begin training cycle of newly generated grid
        %------------------------------------------------------------------
        
        for cycle=1:p.numTrainCycles(layer),
            
            p.eta = cycle^(-p.A);		% Learning rate (how quickly weights are adapted)
            p.G = 0.5 + 10*cycle^(-p.B);		% Gaussian width parameter
            %%% Note: G<0.5 is boring because the Gaussian only covers one node
            
            if cycle == 1 || cycle == p.numTrainCycles(layer),
                fprintf('\nWithin pretrain, Cycle %d, G = %f, ETA = %f', cycle, p.G, p.eta);
            end
            
            [inp_mat] = gen_limited_input(nInpDims/p.nDimReps,p); %generate an input vector
            
            %--------------------------------------------------------------
            % Find winning node
            %--------------------------------------------------------------
            [win_row, win_col, dist_mat] = findWinningNode(w, inp_mat, nInpDims);
            
            %------------------------------------------------------------------
            % Calculate each unit's distance from winner and resultant activation
            [f, ~] = calc_act_fast(win_row, win_col, dist_mat,layer,p,interfere);
            
            
            %%% Update Weights
            w = w + f.*(inp_mat-w);
            
            
        end  % end of training cycles loop
        
        %% Add random noise to the weights

        % add uniform noise [-1,1], then add 1, and divide it all by 3
        % ending weights are again distributed on [0,1]
%         w = ((w + (1 - 2*(rand(p.numRows,p.numRows,nInpDims)))) + 1)./3 ;
        
        % Is there a problem with this, in that only some weight get the maximum possible
        % increment or decrement, and only a subset of these were high or low to start off with, therefore very few end up anywhere
        % near the extremes of the distribution (0 or 1) and thus stimuli that have high or low
        % values tend to have lower matches. Perhaps not as bad as flat_noise_decay version of squidging, since at
        % least the noise increment/decrements will be normally distributed.
        
        
        weights(layer,:,:,1:p.numInputDims(layer),grid) = w;
    end % end of grid loop
    
    
end % end of layer loop

end