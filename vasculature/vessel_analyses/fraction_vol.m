function fv = fraction_vol(data, t_mask)
    seg = data.angio;
    t_vol = sum(t_mask(:));
    fv = sum(seg(:)) ./ t_vol;
end
