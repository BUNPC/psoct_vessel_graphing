function [seg, skel] = overlay_seg_skel_graph(sz, graph, res, ds_flag, fig_title, seg)
% overlay the skeleton and segmentation
%   This code generates several overlaid 3D figures:
%       - segmentation + skeleton
%       - segmentation + graph
%       - skeleton + graph
% INPUTS:
%   sz (array): size of output volume. The order is [y x z]
%   graph (struct): struct containing graph nodes and edges. The nodes'
%                   coordinates have the order [x y z]
%   res (): resolution of the nifti image. Since we are not saving a nifti,
%           this is set at [0.01, 0.01, 0.01], but it's unused.
%            i.e. [0.01, 0.01, 0.01] 0.01cm (10um) isotropics
%   ds_flag (): 1 = the graph was downsampled
%               0 = the graph was not downsampled
%   fig_title (string): name of figure
%   skel (matrix): skeleton of graph (binary)
%   seg (matrix): segmentation (binary)
% OUTPUTS:
%   

%% Graph to skeleton logical
% Convert a graph from nodes + edges to a logical matrix
skel = sk3D(sz, graph, 'foo', res, ds_flag, 'false');

%% Verify dimensions of skeleton and segmentation matrices are equivalent
% Skeleton dimensions
skel_sz = size(skel);
% Segmentation dimensions
seg_sz = size(seg);
% Assert the two dimensions match, otherwise raise an error
assert(all(skel_sz == seg_sz), ['The dimensions of the skeleton and ' ...
    'the segmentation do not match.']);

%% If skeleton matrix dim exceeds 2048 -> downsample skeleton & segmentation
maxdim = 2048;
tmp = sz > maxdim;
if any(tmp)
    %%% Calculate scaling factor
    % Find largest dimension exceeding 2048
    m = max(sz(tmp));
    % Calculate smallest scaling factor to reach 2048
    scalar = ceil(m ./ maxdim);
    
    %%% Perform scaling according to scaling factor.
    % Maintain z-axis dimension, so long as it doesn't exceed 2048
    if sz(3) < 2048
        %%% Segmentation reshape and squeeze (fresh lemonade)
        seg_re = reshape(seg, scalar, sz(1)/scalar, scalar,...
                            sz(2)/scalar, 1, sz(3));
        seg = squeeze(sum(seg_re, [1,3,5])) / scalar.^2;
        %%% Skeleton reshape and squeeze (fresh lemonade)
        skel_re = reshape(skel, scalar, sz(1)/scalar, scalar,...
                            sz(2)/scalar, 1, sz(3));
        skel = squeeze(sum(skel_re, [1,3,5])) / scalar.^2;
    % Downsample z-axis dimension, if its => 2048
    else
        %%% Skeleton reshape and squeeze (fresh lemonade)
        seg_re = reshape(seg, scalar, sz(1)/scalar, scalar,...
                            sz(2)/scalar, scalar, sz(3)/scalar);        
        seg = squeeze(sum(seg_re, [1,3,5])) / scalar.^3;
        %%% Skeleton reshape and squeeze (fresh lemonade)
        skel_re = reshape(skel, scalar, sz(1)/scalar, scalar,...
                            sz(2)/scalar, scalar, sz(3)/scalar);        
        skel = squeeze(sum(skel_re, [1,3,5])) / scalar.^3;
    end
end

%% Crop to get a large vessel for poster
% Create figure for a poster the SPIE Photonics West Bios conference (2024)
%{
seg_crop = seg(498:525, 450:650, 100:end-75);
skel_crop = skel(498:525, 450:650, 100:end-75);
volshow(seg_crop);
seg = seg_crop;
skel = skel_crop;
%}

%% Initialize the 3D figure properties
view_panel = uifigure(figure,'Name',fig_title); close;
v = viewer3d(view_panel);
v.BackgroundColor = 'w';
v.BackgroundGradient = 'off';

%%% Display volume of skeleton (from graph)
h = volshow(skel,'Parent',v);
h.Parent.BackgroundColor = 'w';
% Make skeleton red
skelcolor = repmat([1, 0, 0], [256,1]);
h.Colormap = skelcolor;

%%% Overlay the segmentation
h.OverlayData = seg;
h.OverlayAlphamap = 0.1;
% Make overlay gray
gray = repmat([0.5, 0.5, 0.5], [256,1]);
h.OverlayColormap = gray;

pause(0.1)

end