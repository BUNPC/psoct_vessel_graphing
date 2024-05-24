function volshow_overlay(inner, outer, figtitle)
%volshow_overlay overlay two segmentations or skeletons
%   Initialize the 3D viewer, show the inner volume or skeleton, then
%   overlay with the outer volume or skeleton.
%
%   INPUTS:
%       inner (logical matrix): the segmentation that will be displayed as
%           the inner region of the overlay
%       inner (logical matrix): the segmentation that will be displayed as
%           the outer region of the overlay
%       figtitle (string): title of figure


%% Verify the dimensions are correct
if any(size(inner)>2000) || any(size(outer)>2000)
    error('One dimension of the inner exceeds the matlab limit of 2000.')
elseif size(outer) ~= size(inner)
    error('The dimensions of inner and outer do not match.')
end

%% Make 3D figure
%%% Initialize the figure for displaying 3D volume
view_panel = uifigure('Name',figtitle); close;
v = viewer3d(view_panel);
v.BackgroundColor = 'w';
v.BackgroundGradient = 'off';

%%% Display volume of skeleton (from graph)
h = volshow(inner,'Parent',v);
h.Parent.BackgroundColor = 'w';
% Make the inner red
skelcolor = repmat([1, 0, 0], [256,1]);
h.Colormap = skelcolor;

%%% Overlay the segmentation
h.OverlayData = outer;
h.OverlayAlphamap = 0.4;
% Make outer overlay gray
gray = repmat([0.5, 0.5, 0.5], [256,1]);
h.OverlayColormap = gray;


end