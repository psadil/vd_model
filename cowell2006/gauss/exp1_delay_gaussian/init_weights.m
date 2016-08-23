function m_out = init_weights(p, layer)


% (column X row X inputDim)

m_out = rand(p.numRows,p.numRows,p.numInputDims(layer));


end


