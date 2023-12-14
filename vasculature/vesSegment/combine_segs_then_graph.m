function [segfinal, graph_final] = combine_segs_then_graph(ov)
% Combine segmentations for the same volume from multiple seg params.
%   The Frangi vessel enhancement algorithm uses a Gaussian filter to
%   enhance vessels of a specific diameter. It is necessary to run the
%   Frangi with multiple Gaussian kernels to segment a range of vessel
%   diameters. This function will combine the segmentations from each of
%   these segmentations, to generate a segmentation that includes small to
%   large vessels.

%%% INPUTS:
%       ov (struct): struct containing filepath
%           ov(#).subid: subject ID
%           ov(#).basepath: base filepath to segmentation files
%           ov(#).fpath: cell array of strings of full file path to the
%                       .TIFs for each segmentation
%%% OUTPUTS:
%       segfinal (matrix): combined segmentation
%       graphfinal (struct): graph from combined segmentation

%% Overlay all segmentations
for ii = 1:length()

%% Initialize Graph, remove loops, and save output

%% Generate quality assurance overlays

end