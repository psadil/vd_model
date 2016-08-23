function [p,weights] = pretrain(p)
% pretrain -- initializes the weights of a rat and trains it to stimuli
% randomly generated from all possible stimuli

% called by: runSim
% calls: init_weights, gen_limited_input, findWinningNode, calc_act_fast

% input
%    p: experimental structure

% output
%   p: experimental structure, now with p.eta and p.G set (NOTE: these
%      should match p.etaExp and p.G_exp
%   weights: weights associated with each grid. Every node in a grid will
%      have a vector of associated weights. Vector length is determined by
%      how many dimensions that grid is tuned to. For example, nodes in the
%      caudal grid currently attend to 2 dims, nodes in PRC layer attend to
%      8 dims.


%%
fprintf('\ndelay_pretrain being executed...');

% dims of weights: layers x rows x rows x inputDims x grid
% NOTE: this scheme means that much of this array will remain at 0. That
%    is, since the PRC layer has only 1 grid, only its first grid dim will
%    be used.
weights = zeros(max(p.numLayers),p.numRows,p.numRows,max(p.numInputDims),max(p.nGrids));

%% Perform pre-training of network

for layer = 1:max(p.numLayers)
    
    
    for grid = 1:p.nGrids(layer),
        
        fprintf('\n\nLayer number %d, Grid no. %d...\n(1==Caudal, 2==PRC)', layer, grid);
        
        
        % initialize random weights for grids in this layer
        % w indicates temporary weight matrix, to be loaded into weights at
        % end of training
        w = init_weights(p, layer);
        
        %------------------------------------------------------------------
        % begin training cycle of newly generated grid
        %------------------------------------------------------------------
        
        for cycle=1:p.numTrainCycles,
            
            % Learning rate
            p.eta = cycle^(-p.A);		 
            % Gaussian width
            p.G = 0.5 + 10*cycle^(-p.B);		
            
            if cycle == 1 || cycle == p.numTrainCycles,
                fprintf('\nWithin pretrain, Cycle %d, G = %f, ETA = %f', cycle, p.G, p.eta);
            end
            
            %generate a sample input to train weights on.
            inp_mat = gen_limited_input(p.numInputDims(layer),p); 
            
            %--------------------------------------------------------------
            % Find winning node
            %--------------------------------------------------------------
            [win_row, win_col, ~] = findWinningNode(w, inp_mat);
            
            %--------------------------------------------------------------
            % Calculate each unit's distance from winner and resultant activation
            % f is gaussian neighborhood.
            %--------------------------------------------------------------
            f = calc_neigh(win_row, win_col,layer,p);
            
            
            % Update Weights
            w = w + f.*(inp_mat-w); 
        end 
        
        % load pretrained weights into grand weight matrix
        weights(layer,:,:,1:p.numInputDims(layer),grid) = w;
    end
    
    
end

end