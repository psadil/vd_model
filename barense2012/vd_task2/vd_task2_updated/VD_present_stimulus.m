function [W, selec, initial_selec, p, pktot, usePRC, act_out, initial_acts] = VD_present_stimulus(stim, W, p, features_sampled, trial,pktot)

% Function called by VD_present_stimulus.m. Presents the network with the two stimuli on this trial.
% Choose whichever layer demonstrates higher selectivity

% present chosen features of stimulus to all layers. pick the layer with
% the highest selectivity, and output that selectivity

% 5/18/15 ps
% commented out all assignments of activations (in an effort to speed up
% processing), since not looking at them directly right now.

inp_mat = repmat(reshape(stim,[1 1 length(stim)]), [p.nRows p.nRows 1]);

act_out = zeros(p.layer,p.nRows,p.nRows,p.maxNumGrids); % 5 grids (only 1 stimulus)
initial_acts = zeros(p.layer,p.nRows,p.nRows,p.maxNumGrids);
selec = zeros(p.numLayers,max(p.numGrids));
initial_selec = zeros(p.numLayers,max(p.numGrids));
% pktot.fin_act_peak = zeros(p.numLayers,max(p.numGrids));
% pktot.fin_act_total = zeros(p.numLayers,max(p.numGrids));
% pktot.init_act_peak = zeros(p.numLayers,max(p.numGrids));
% pktot.init_act_total = zeros(p.numLayers,max(p.numGrids));
% winact = zeros(p.numLayers,p.maxNumGrids,2);  % activation in each grid in x,y coordinates?

% only use PRC layer if available and enough features were sampled
if (p.layer == 2) && (length(features_sampled) == p.numGrids_Caudal)
    2;
    usePRC = 1;
else
    usePRC = 0;
end



%% Expose network to stimuli and update weights
% switchRatio = 0;

% fix = 0;
% while switchRatio < p.fixn_ratio_lowHigh(p.stimCond)
%
%     % if fixation wasn't outside of stimulus, incriment fixations
%     % for each encoding cycle
%     p.fixations(trial) = p.fixations(trial) + 1; %% total fixations across both stimuli
%
for layer=1:p.layer
    
    firstFeatureToCheck=(1:p.numInputDims(layer):p.nGrids(layer)*p.numInputDims(layer));
    lastFeatureToCheck = firstFeatureToCheck+p.numInputDims(layer)-1;
    
    
    for grid = 1:p.nGrids(layer),
        
        input_mat=zeros(p.nRows,p.nRows,p.numInputDims(layer));
        
        % should end with a (:,:,15) input_mat for PRC layer that checks features 1:15, and a (:,:,3) input_mat for caudal, where each 4th-dim checks three features
        input_mat(:,:,:)=inp_mat(:,:,(firstFeatureToCheck(grid):lastFeatureToCheck(grid)));
        
        % put the variable 'weights' into the format previously accepted by the model.
        weights = W(layer,:,:,1:p.numInputDims(layer),grid);
        weights = squeeze(weights);
        
        % if no features that this grid pays attention to were sampled,
        % or if dealing with PRC but 5 features weren't sampled, skip
        % presenting
        %         if (layer==1 && (~any(features_sampled==grid) || usePRC)) || (layer==2 && ~usePRC)
        if ~p.fives
            if (layer==1 && (~any(features_sampled==grid))) || (layer==2 && ~usePRC)
%                 if layer ==2
%                     2;
%                 end
                
                continue
            end
        end
        
        nEncodCycles = p.nEncodCycles;
        if p.variableEncode
            stick_switch = rand;
            while stick_switch < p.fixn_ratio_lowHigh(p.stimCond)
                nEncodCycles = nEncodCycles+p.nEncodCycles;
                p.fixations(trial) = p.fixations(trial) + 1; %% total fixations across both stimuli
                stick_switch = rand;
            end
        end
        
        
        %% inside the encoding cycle?
        %------------------------------------------------------------------
        % Find winning node
        %--------------------------------------------------------------
%         [win_row, win_col, dist_mat] = findWinningNode(weights, input_mat, p.numInputDims(layer));
        
%         if p.layer==2
%             if layer ==2
%                 2;
%             end
%         end
        
%         [~, acts, selectivity, p, act_peak, act_total] = VD_calc_selectivity_fast(win_row, win_col, dist_mat, p, p.numInputDims(layer));
%         
%         initial_selec(layer,grid) = selectivity;
%         pktot.init_act_peak(layer,grid) = act_peak;
%         pktot.init_act_total(layer,grid) = act_total;
%         
%         % initial weights of new grid.
%         initial_acts(layer,:,:,grid) = acts;
        
        
        for cycle=1:nEncodCycles
            
            [win_row, win_col, dist_mat] = findWinningNode(weights, input_mat, p.numInputDims(layer));
            
            
            p.winning(layer,grid,trial,1:2) = [win_row, win_col];
            
            %--------------------------------------------------------------
            %Calculate each unit's distance from winner and activation
            %--------------------------------------------------------------
            [f, acts, selectivity, p, act_peak, act_total] = VD_calc_selectivity_fast(win_row, win_col, dist_mat, p, p.numInputDims(layer));
            
            if cycle==1,  %need to compare last set of weights of old grid to
                initial_selec(layer,grid) = selectivity;
                pktot.init_act_peak(layer,grid) = act_peak;
                pktot.init_act_total(layer,grid) = act_total;
                
                % initial weights of new grid.
                initial_acts(layer,:,:,grid) = acts;
                
                if layer == 2
                   2; 
                end
                
            end
            
            
            %%% Update Weights
            weights = weights + f.*(input_mat-weights);  % update based on spire around winning node
            
            
        end %%% Go to next cycle (if switchRatio is low enough)
        [win_row, win_col, dist_mat] = findWinningNode(weights, input_mat, p.numInputDims(layer));
        
        [~, acts, selectivity, p, act_peak, act_total] = VD_calc_selectivity_fast(win_row, win_col, dist_mat, p, p.numInputDims(layer));
        
        
        W(layer,:,:,1:p.numInputDims(layer),grid) = weights; %% Overwrite weight values in the orignial W matrix, for this grid.
        act_out(layer,:,:,grid) = acts; % use activations as calculated AFTER final weight update
        selec(layer,grid) = selectivity;
        pktot.fin_act_peak(layer,grid) = act_peak;
        pktot.fin_act_total(layer,grid) = act_total;
        
    end  % end of loop for each grid
    
end % end of layer loop

%     nothingRatio = rand;
%     if nothingRatio > p.outsideRatio(p.stimCond)
%         p.fixations(trial) = p.fixations(trial) + 1; %% total fixations across both stimuli
%     end

%     fix = fix+1;
%     % determine whether to stay or switch to other stimulus
%     switchRatio = rand;
% % end


