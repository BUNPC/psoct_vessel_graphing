function fv = fraction_vol(data, t_mask)
%FRACTION_VOL calculate the volume fraction of the vasculature.
% Numerator = segmentation. Denominator = tissue mask
% The quotient is a unitless metric.
%
% INPUTS:
%   data (struct): the data structure including the segmentation
%           (data.angio) and the graph
%   t_mask (logical): the tissue mask for the region
% OUTPUTS:
%   fv (double): volume fraction (unitless)

%% Extract values
% Extract the segmentation from the data struct
seg = data.angio;
% Verify that segmentation is a logical
if ~isa(seg,'logical')
    seg = logical(seg);
    sprintf('\nConverted segmentation to logical\n')
end
% Verify that tissue mask is a logical
if ~isa(t_mask,'logical')
    t_mask = logical(t_mask);
    sprintf('\nConverted mask to logical\n')
end

%% Calculations
% Calculate the sum of the tissue mask (including segmentation)
t_vol = sum(t_mask(:));
% Calculate the ratio of vessel:tissue ratio
fv = sum(seg(:)) ./ t_vol;

end
