function branchden = branch_density(data, t_mask)
% volume is t_vol

    vox = data.Graph.vox;
    vox_vol = vox(1) .* vox(2) .* vox(3);   % Creat function to calculate tissue volume
    volume = sum(t_mask(:)) .* vox_vol;
    nb = data.Graph.nB;
    branchden = sum(nb) / volume;
end