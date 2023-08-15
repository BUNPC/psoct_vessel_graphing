function SaveData( obj, eventdata, fpath, fname )
global im


try
    handles = guidata(im.gcf);
catch
    im.gcf = gcf;
    handles = guidata(im.gcf);
end


if ~exist('fpath')
    fpath = '';
    fname = 'autosave';
    hwait = waitbar(0,'Autosaving...');
else
    hwait = waitbar(0,'Saving...');
end
if length(fname)<=5
    fname = sprintf('%s.seed',fname);
elseif ~strcmp(fname(end-4:end),'.seed')
    fname = sprintf('%s.seed',fname);
end

im2 = rmfield(im,'III');
if isfield(im2,'vImarisApplication');
    im2 = rmfield(im2,'vImarisApplication');
end

%im2.gui. = get(handles.,'enable');
im2.gui.pushbuttonRunSeeds = get(handles.pushbuttonRunSeeds,'enable');
im2.gui.pushbuttonProcSeeds = get(handles.pushbuttonProcSeeds,'enable');
im2.gui.pushbuttonSeedConnectivity = get(handles.pushbuttonSeedConnectivity,'enable');
im2.gui.pushbuttonGraph = get(handles.pushbuttonGraph,'enable');
%im2.gui.pushbuttonImaris = get(handles.pushbuttonImaris,'enable');
im2.gui.pushbuttonRegraphNodes = get(handles.pushbuttonRegraphNodes,'enable');
im2.gui.pushbuttonImageGraph = get(handles.pushbuttonImageGraph,'enable');
im2.gui.pushbuttonUpdateVesselOverlay = get(handles.pushbuttonUpdateVesselOverlay,'enable');
im2.gui.pushbuttonUpdateVesselMask = get(handles.pushbuttonUpdateVesselMask,'enable');


im2.gui.uipanelRunAll = get(handles.uipanelRunAll,'visible');

%im2.gui. = get(handles.,'string');
im2.gui.editAutoSeedHxy = get(handles.editAutoSeedHxy,'string');
im2.gui.editAutoSeedHz = get(handles.editAutoSeedHz,'string');
im2.gui.editAutoSeedThresh = get(handles.editAutoSeedThresh,'string');
im2.gui.editAutoSeedNtrks = get(handles.editAutoSeedNtrks,'string');
im2.gui.editGraphHxy = get(handles.editGraphHxy,'string');
im2.gui.editGraphHz = get(handles.editGraphHz,'string');
im2.gui.editSeedConnHxy = get(handles.editSeedConnHxy,'string');
im2.gui.editSeedConnHz = get(handles.editSeedConnHz,'string');
im2.gui.editSeedConnThresh = get(handles.editSeedConnThresh,'string');
im2.gui.editCapillaryMaxDiam = get(handles.editCapillaryMaxDiam,'string');
im2.gui.editImageThresh = get(handles.editImageThresh,'string');

im2.gui.checkboxCenterStep = get(handles.checkboxCenterStep,'value');

im2.gui.menu_Flow_useSegments = get(handles.menu_Flow_useSegments,'checked');
im2.gui.menu_ViewSeedPanels = get(handles.menu_ViewSeedPanels,'checked');
im2.gui.menu_ViewSeeds = get(handles.menu_ViewSeeds,'checked');
im2.gui.menu_VisualizeCentering = get(handles.menu_VisualizeCentering,'checked');

save([fpath fname],'im2');

close(hwait);
