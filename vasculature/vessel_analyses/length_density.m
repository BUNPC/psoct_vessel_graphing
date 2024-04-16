function ld = length_density(data, t_mask)
%%
    vox = data.Graph.vox;
    vox_vol = vox(1) .* vox(2) .* vox(3);
    t_vol = sum(t_mask(:)) .* vox_vol;
    ld = sum(data.Graph.segInfo.segLen_um) ./ t_vol;
end

    