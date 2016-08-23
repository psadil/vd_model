function m_out = VD_init_weights(p, layer)

% avail_values = [.05 .333 .667 .95];

% (column X row X inputDim)

m_out = rand(p.numRows,p.numRows,p.numInputDims(layer));


end


