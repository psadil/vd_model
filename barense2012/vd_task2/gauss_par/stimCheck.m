
tmp = 0;
for checker = 1:p.nMismatch
    checking = fid.stimuli1(checker,:,1);
    for row = checker+1:p.nMismatch
        stim1 = sum(checking==fid.stimuli1(row,:,1));
        stim2 = sum(checking==fid.stimuli2(row,:,1));
        
        tmp = max(tmp,stim1);
        tmp = max(tmp,stim2);
    end

end
