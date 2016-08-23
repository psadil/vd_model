function [p,weights] = pretrain(p)

fprintf('\ndelay_pretrain being executed...');
interefere = 0;

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
        
        for cycle=1:p.numTrainCycles,
            
            p.eta = cycle^(-p.A);		% Learning rate 
            p.G = 0.5 + 10*cycle^(-p.B);		% Gaussian width parameter
            
            if cycle == 1 || cycle == p.numTrainCycles,
                fprintf('\nWithin pretrain, Cycle %d, G = %f, ETA = %f', cycle, p.G, p.eta);
            end
            
            inp_mat = gen_limited_input(nInpDims,p); %generate an input vector
            
            %--------------------------------------------------------------
            % Find winning node
            %--------------------------------------------------------------
            [win_row, win_col, dist_mat] = findWinningNode(w, inp_mat, nInpDims);
            
            %------------------------------------------------------------------
            % Calculate each unit's distance from winner and resultant activation
            [f, ~] = calc_act_fast(win_row, win_col, dist_mat,layer,p,interefere);
            
            
            %%% Update Weights
            w = w + f.*(inp_mat-w);
            
            
        end  % end of training cycles loop
        
        
        weights(layer,:,:,1:p.numInputDims(layer),grid) = w;
    end % end of grid loop
    
    
end % end of layer loop

end