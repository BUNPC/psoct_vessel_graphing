%filenm_seed = 'HandlesTotalStack_466x466x351_Nov_7_2008';
filenm_seed = 'HandlesTotalStack_466x466x411_Nov_17_2008';

%%
% Load xml file
tg = SuperellipseXMLReader(sprintf('%s.xml',filenm_seed));

%%
% create vertices and edges and diam
v = reshape([tg.position],[3 length(tg)])';
e = [];
diam = zeros(length(tg),1);
for iN = 1:length(tg)
    lst = tg(iN).nbr+1;
    if tg(iN).ID~=iN-1
        disp( sprintf('ID ~= In-1 for iN=%d',iN) )
        break
    end
    lst(find(lst<iN)) = [];
    if ~isempty(lst)
        e(end+[1:length(lst)],:) = [iN*ones(length(lst),1) lst'];
    end
    
    diam(iN) = max(tg(iN).scale(1:2));
end

%%
figure(10);
trisurf([e e(:,1)],v(:,1),v(:,2),v(:,3))
axis image

%%
% create im2 structure
info = imfinfo(filenm_seed,'tiff');

im2.curSeed = 0;

im2.Hvox = [1 1 1];
im2.maskMaxRad = 20;
im2.maskSphereSpacing = 0.7000;
im2.maskCircleFlag = 1;

im2.filenm = filenm_seed;
im2.filetype = 'tiff';
im2.nZ = length(info);
im2.nodePos = v;
im2.nodeEdges = e;
im2.nodeDiam = diam;

im2.gui.menu_Flow_useSegments = 'on';
im2.gui.menu_ViewSeedPanels = 'off';
im2.gui.menu_ViewSeeds = 'off';
im2.gui.menu_viewPO2pts = 'off';
im2.gui.editImageThresh = '0.1';
im2.gui.uipanelRunAll = 'off';

save(sprintf('%s.seed',filenm_seed),'im2');

