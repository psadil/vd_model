function[p] = save_peakTotActs(pktot, p, trial, whichCaudal)

%% save up 

if p.layer == 2
    2;
end

if (p.layer == 2) && (all(p.usePRC(:,trial)))
        % note: caudal is layer 1 (first index), PRC layer 2 (first index)
    p.prevStimFin_act_peak(2,trial) = pktot.prevStimFin_act_peak(2,1);
    p.prevStimFin_act_total(2,trial) = pktot.prevStimFin_act_total(2,1);
    
    % don't know how to care about previous stim's initial acts yet
    p.prevStimInit_act_peak(2,trial) = pktot.prevStimInit_act_peak(2,1);
    p.prevStimInit_act_total(2,trial) = pktot.prevStimInit_act_total(2,1);
    
    % look at new stimuli (which has final different from initial only
    % in trials that were determined to be matching
    p.newStimInit_act_peak(2,trial) = pktot.init_act_peak(2,1);
    p.newStimInit_act_total(2,trial) = pktot.init_act_total(2,1);
    
    % look at new stimuli (which has final different from initial only
    % in trials that were determined to be matching
    p.newStimFin_act_peak(2,trial) = pktot.fin_act_peak(2,1);
    p.newStimFin_act_total(2,trial) = pktot.fin_act_total(2,1);

else
    % note: caudal is layer 1 (first index), PRC layer 2 (first index)
    p.prevStimFin_act_peak(1,trial) = pktot.prevStimFin_act_peak(1,whichCaudal);
    p.prevStimFin_act_total(1,trial) = pktot.prevStimFin_act_total(1,whichCaudal);
    
    % don't know how to care about previous stim's initial acts yet
    p.prevStimInit_act_peak(1,trial) = pktot.prevStimInit_act_peak(1,whichCaudal);
    p.prevStimInit_act_total(1,trial) = pktot.prevStimInit_act_total(1,whichCaudal);
    
    % look at new stimuli (which has final different from initial only
    % in trials that were determined to be matching
    p.newStimInit_act_peak(1,trial) = pktot.init_act_peak(1,whichCaudal);
    p.newStimInit_act_total(1,trial) = pktot.init_act_total(1,whichCaudal);
    
    % look at new stimuli (which has final different from initial only
    % in trials that were determined to be matching
    p.newStimFin_act_peak(1,trial) = pktot.fin_act_peak(1,whichCaudal);
    p.newStimFin_act_total(1,trial) = pktot.fin_act_total(1,whichCaudal);

end