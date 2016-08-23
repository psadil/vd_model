function [ p ] = calc_recognition( p, selec_forComp, trial)
%% calc_recognition: calcs recognition socres


%%

if p.layer == 2
    
    tmp_prc = squeeze(selec_forComp(2,1,:));
    tmp_caudal = squeeze(mean(selec_forComp(1,:,:),2));
    p.recogByLayer(trial,2) = (tmp_prc(1) - tmp_prc(2)) / (tmp_prc(1) + tmp_prc(2));
    selec = mean([tmp_prc,tmp_caudal],2);
else
    selec = squeeze(mean(selec_forComp(1,:,:),2));
end


tmp_caudal = squeeze(mean(selec_forComp(1,:,:),2));
p.recogByLayer(trial,1) = (tmp_caudal(1) - tmp_caudal(2)) / (tmp_caudal(1) + tmp_caudal(2));
p.recognition(trial) = (selec(1) - selec(2)) / (selec(1) + selec(2));

end