function varargout = imView3d(varargin)
% IMVIEW3D M-file for imView3d.fig
%      IMVIEW3D, by itself, creates a new IMVIEW3D or raises the existing
%      singleton*.
%
%      H = IMVIEW3D returns the handle to a new IMVIEW3D or the handle to
%      the existing singleton*.
%
%      IMVIEW3D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMVIEW3D.M with the given input arguments.
%
%      IMVIEW3D('Property','Value',...) creates a new IMVIEW3D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imView3d_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to imView3d_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help imView3d

% Last Modified by GUIDE v2.5 03-Jun-2009 15:00:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @imView3d_OpeningFcn, ...
                   'gui_OutputFcn',  @imView3d_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before imView3d is made visible.
function imView3d_OpeningFcn(hObject, eventdata, handles, varargin)
global im

im = [];

% Choose default command line output for imView3d
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

if length(varargin)==1
    im.filenm = varargin{1};
    error('Does this work?')
else
%    [fname fpath fidx] = uigetfile( {'*.img','Image files';'*.tif;*.tiff','Tiff files';'*.seed','Seed Files'},'Load file');
    [fname fpath fidx] = uigetfile( '*.img;*.tif;*.tiff;*.mat;*.seed','Load file');
    if fname==0
        return
    end
    cd(fpath)
    if strcmpi(fname(end-2:end),'img')
        im.filenm = fname(1:end-4);
        im.filetype = 'img';
    elseif strcmpi(fname(end-2:end),'tif')
        im.filenm = fname(1:end-4);
        im.filetype = 'tif';
    elseif strcmpi(fname(end-3:end),'tiff')
        im.filenm = fname(1:end-5);
        im.filetype = 'tiff';
    elseif strcmpi(fname(end-2:end),'mat')
        im.filenm = fname(1:end-4);
        im.filetype = 'mat';        
    else
        im.filenm = fname(1:end-5);
        im.filetype = '';
    end
end

if strcmp(im.filetype,'img')
    filenm = sprintf( '%s.img',im.filenm );
    fid = fopen( filenm, 'rb' );
    n = fread( fid, 3, 'int');
    nx = n(1); ny = n(2); nz = n(3);
    I = fread( fid, nx*ny*nz, 'float' );
    I = reshape(I, [ny nx nz] );
    fclose(fid);

    im.III = uint8(zeros(size(I)));
    im.III = uint8(32*(I-0.0)/(1-0.0));
    [im.nY,im.nX,im.nZ] = size(im.III);

    im.nSeed = 0;
    im.seedDirectionFlag = 0;
    im.curSeed = 0;
elseif strcmp(im.filetype,'tif') | strcmp(im.filetype,'tiff')
    finfo = imfinfo(sprintf('%s.%s',im.filenm,im.filetype));
    ny = finfo(1).Height;
    nx = finfo(1).Width;
    nz = length(finfo);
    I = zeros(ny,nx,nz);
    for ii=1:nz
        I(:,:,ii) = imread(sprintf('%s.%s',im.filenm,im.filetype),ii);
    end
    I = I / max(I(:));
        
    im.III = uint8(zeros(size(I)));
    im.III = uint8(32*(I-0.0)/(1-0.0));
    [im.nY,im.nX,im.nZ] = size(im.III);


    im.nSeed = 0;
    im.seedDirectionFlag = 0;
    im.curSeed = 0;
elseif strcmp(im.filetype,'mat')
    load(im.filenm,'I');
    I = I / max(I(:));
    
    im.III = uint8(zeros(size(I)));
    im.III = uint8(32*(I-0.0)/(1-0.0));
    [im.nY,im.nX,im.nZ] = size(im.III);


    im.nSeed = 0;
    im.seedDirectionFlag = 0;
    im.curSeed = 0;
else
    load([fpath fname],'-mat','im2');
    im = im2;

    updateImageOverlay(0, handles);
    UpdateSeeds( handles )

    if isfield(im,'gui')
        guis = fieldnames(im.gui);
        for ii=1:length(guis)
            foos = guis{ii};
            try
                if strcmpi(foos(1:4),'push')
                    eval( sprintf('set(handles.%s,''enable'',''%s'');',foos,getfield(im.gui,foos)) );
                elseif strcmpi(foos(1:2),'ui')
                    set(handles.uipanelRunAll,'visible',im.gui.uipanelRunAll);
                elseif strcmpi(foos(1:4),'menu')
                    eval( sprintf('set(handles.%s,''checked'',''%s'');',foos,getfield(im.gui,foos)) );
                elseif strcmpi(foos(1:4),'chec')
                    eval( sprintf('set(handles.%s,''value'',%d);',foos,getfield(im.gui,foos)) );
                else
                    eval( sprintf('set(handles.%s,''string'',''%s'');',foos,getfield(im.gui,foos)) );
                end
            catch
            end
        end
    end

    set(handles.pushbuttonImageGraph,'enable','on')
    if isfield(im,'nodeDiam')
        set(handles.radiobuttonSelSegment,'value',1)
    end
    if isfield(im,'RunPath')
        set(handles.pushbuttonRunPath,'string',sprintf('Path: %s',im.RunPath));
    end
    
    pushbuttonImageGraph_Callback([], [], handles);

end

% set flagUseSegments
if strcmpi(get(handles.menu_Flow_useSegments,'checked'),'on')
    im.flagUseSegments = 1;
else
    im.flagUseSegments = 0;
end


if ~isfield(im,'Hvox')
    ch = menu('Please enter voxel dimensions Hx, Hy, Hz in the command window','okay');
    im.Hvox = str2num(input('Enter Hx Hy Hz?','s'));
    
    while length(im.Hvox)~=3
        im.Hvox = str2num(input('Must be 3 numbers... Enter Hx Hy Hz?','s'));
    end
end

% this is temporary and can be removed soon DAB 3/26/08
if isfield(im,'nodeFlag')
    im.edgeFlag = im.nodeFlag;
    im =rmfield(im,'nodeFlag');
end
if ~isfield(im,'nodeType') & isfield(im,'nodePos')
    im.nodeType = zeros(size(im.nodePos,1),1);
end
if ~isfield(im,'nodeBCType') & isfield(im,'nodePos')
    im.nodeBCType = zeros(size(im.nodePos,1),1);
end


set(handles.vesselGraph,'Name',sprintf('vesselGraph %s',im.filenm))

set(handles.slider1,'max',im.nZ)
set(handles.slider1,'value',1)
set(handles.slider1,'min',1)
set(handles.slider1,'sliderstep',[1/im.nZ .1])

axes( handles.axes1 );
xlim([1 size(im.III,2)]);
ylim([1 size(im.III,1)]);

updateAxes( handles );

set( handles.axes1,'ButtonDownFcn','imView3d(''axes1ButtonDown_Callback'',gcbo,[],guidata(gcbo))' );
set(get(handles.axes1,'children'), ...
    'ButtonDownFcn', 'imView3d(''axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))' );

if strcmpi(get(handles.menu_ViewSeedPanels,'checked'),'off')
    set(handles.uipanel_AutoSeed,'visible','off')
    set(handles.uipanelSeed,'visible','off')
    set(handles.uipanel_SeedConnectivity,'visible','off')
    set(handles.radiobuttonSeedPoint,'visible','off')
    set(handles.uipanelRunAll,'visible','off')
else
    set(handles.uipanel_AutoSeed,'visible','on')
    set(handles.uipanelSeed,'visible','on')
    set(handles.uipanel_SeedConnectivity,'visible','on')
    set(handles.radiobuttonSeedPoint,'visible','on')
    set(handles.uipanelRunAll,'visible','on')
end
%set(handles.menu_ViewSeedPanels,'checked','off')
%set(handles.menu_VisualizeCentering,'checked','off')

im.timer = timer('TimerFcn',@imView3d_SaveData, 'Period', 3600, 'StartDelay', 3600, 'ExecutionMode','fixedSpacing','TasksToExecute',10000);
im.gcf = gcf;
start(im.timer);


% UIWAIT makes imView3d wait for user response (see UIRESUME)
% uiwait(handles.vesselGraph);


% --- Outputs from this function are returned to the command line.
function varargout = imView3d_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
global im

z = round(get(handles.slider1,'value'));
nS = str2num(get(handles.editNslices,'string'));
if im.ImageView == 1
    if (z+nS-1)>size(im.III,3)
        z = size(im.III,3)-nS+1;
        set(handles.slider1,'value',z)
    end
elseif im.ImageView == 2
    if (z+nS-1)>size(im.III,1)
        z = size(im.III,1)-nS+1;
        set(handles.slider1,'value',z)
    end
else
    if (z+nS-1)>size(im.III,2)
        z = size(im.III,2)-nS+1;
        set(handles.slider1,'value',z)
    end
end
set(handles.textSlices,'string',sprintf('Slice %d to %d',z,z+nS-1))

updateAxes( handles );

if get(handles.radiobuttonZoom,'value')
    pan off
    zoom on
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateAxes( handles )
global im

axes( handles.axes1 );
xrng = xlim();
yrng = ylim();


zzMin = round(get(handles.slider1,'value'));
nS = str2num(get(handles.editNslices,'string'));
zzMax = zzMin+nS-1;

cm = gray(32);
cm(33:64,:) = hsv2rgb([ones(32,1)*.1 ones(32,1) [1:32]'/32]);
cm(65:255,:) = 0;
cm(65:96,:) = hsv2rgb([ones(32,1)*0 ones(32,1) [1:32]'/32]);
cm(97:128,:) = hsv2rgb([ones(32,1)*0.35 ones(32,1) [1:32]'/32]);
cm(129:160,:) = hsv2rgb([ones(32,1)*0.68 ones(32,1) [1:32]'/32]);
cm(161:192,:) = jet(32);
cm(255,:) = [1 0 0];
cm(251,:) = [1 1 0];
cm(250,:) = [1 0 0];
if get(handles.checkboxDisplayGrps,'value') && isfield(im,'nodeGrp')
    if ~get(handles.checkboxHighlightGrp,'value')
        nGrps = length(unique(im.nodeGrp));
        cm(250+[1:nGrps],:) = jet(nGrps);
    else
        nGrps = 2;
        cm(251,:) = [1 1 0];
        cm(252,:) = [0 1 1];
    end
else
    nGrps = 1;
end
thresh = ceil(str2num(get(handles.editImageThresh,'string'))*32);
cm(1:thresh,2:3) = 0;
colormap(cm)

if ~isfield(im,'ImageView');
    im.ImageView = 1;
end
if im.ImageView == 1
    imagesc( max(im.III(:,:,zzMin:zzMax),[],3), [1 255])%65+nGrps] )
elseif im.ImageView ==2
    imagesc( squeeze(max(im.III(zzMin:zzMax,:,:),[],1))', [1 255])%65+nGrps] )
else
    imagesc( squeeze(max(im.III(:,zzMin:zzMax,:),[],2))', [1 255])%65+nGrps] )
end
xlim(xrng)
ylim(yrng)

set( handles.axes1,'ButtonDownFcn','imView3d(''axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))' );
set(get(handles.axes1,'children'), ...
    'ButtonDownFcn', 'imView3d(''axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))' );

if strcmpi(get(handles.menu_ViewSeeds,'checked'),'on')
    for ii=1:im.nSeed
        if im.seed(ii).pos(3)>=zzMin & im.seed(ii).pos(3)<=zzMax
            ht=text(im.seed(ii).pos(1),im.seed(ii).pos(2),num2str(ii));
            set(ht,'color',[1 0 0])%[0.1667 1 1])
            set(ht,'clipping','on')
            set(ht,'hittest','off')
            if ~isempty(im.seed(ii).pos2)
                hl = line([im.seed(ii).pos(1) im.seed(ii).pos2(1)],[im.seed(ii).pos(2) im.seed(ii).pos2(2)]);
                set(hl,'color',[0.1667 1 1])
            end
        end
    end
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function editNslices_Callback(hObject, eventdata, handles)
global im

if im.ImageView==1
    nZ = im.nZ;
elseif im.ImageView==2
    nZ = im.nY;
else
    nZ = im.nX;
end

z = get(handles.slider1,'value');
nS = str2num(get(handles.editNslices,'string'));

if (z+nS) > nZ
    nS = floor(nZ - z - 2);
    set(handles.editNslices,'string',num2str(nS))
end

set(handles.slider1,'value',min(z,nZ-nS+1));
set(handles.slider1,'max',nZ-nS+1);

z=round(z);
set(handles.textSlices,'string',sprintf('Slice %d to %d',z,z+nS-1))

updateAxes( handles );





% --- Executes on button press in radiobuttonZoom.
function radiobuttonZoomPan_Callback(hObject, eventdata, handles)
set(handles.pushbuttonNodeUndo,'visible','off');
switch eventdata
    case 1
        pan off
        zoom on
%        set(handles.radiobuttonZoom,'value',1);
%        set(handles.radiobuttonPan,'value',0);
%        set(handles.radiobuttonSeedPoint,'value',0);
    case 2
        zoom off
        pan on
%        set(handles.radiobuttonZoom,'value',0);
%        set(handles.radiobuttonPan,'value',1);
%        set(handles.radiobuttonSeedPoint,'value',0);
    case {3,4,5}
        zoom off
        pan off
%        set(handles.radiobuttonZoom,'value',0);
%        set(handles.radiobuttonPan,'value',0);
%        set(handles.radiobuttonSeedPoint,'value',1);
        set( handles.axes1,'ButtonDownFcn','imView3d(''axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))' );
        set(get(handles.axes1,'children'), ...
            'ButtonDownFcn', 'imView3d(''axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))' );
end




% --- Executes on key press .
function axes1_KeyPressFcn(hObject, eventdata, handles )
global im

ch=get(handles.vesselGraph,'currentcharacter');
if ch=='a'
    dz = 1;
elseif ch=='z'
    dz = -1;
else
    return;
end

if im.ImageView==1
    nZ = im.nZ;
elseif im.ImageView==2
    nZ = im.nY;
else
    nZ = im.nX;
end

z = round(get(handles.slider1,'value'))+dz;
nS = str2num(get(handles.editNslices,'string'));

if (z+nS) > nZ 
    z = z -1;
elseif z<1
    z = 1;
end
set(handles.slider1,'value',round(z));
set(handles.textSlices,'string',sprintf('Slice %d to %d',z,z+nS-1))

updateAxes( handles );

if get(handles.radiobuttonZoom,'value')
    pan off
    zoom on
end


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
global im

pos = get(gca,'CurrentPoint');

nSeed = im.nSeed;

xx = round(pos(1,1));
yy = round(pos(1,2));

z = round(get(handles.slider1,'value'));
nS = str2num(get(handles.editNslices,'string'));
zzMin = z;
zzMax = z+nS-1;
zz = (zzMin + zzMax)/2;

if im.ImageView==1
    pIdx = [1 2 3];
    pIdx2 = [1 2 3];
elseif im.ImageView==2
    pIdx = [3 1 2];
    pIdx2 = [2 3 1];
    foo = xx;
    xx = yy;
    yy = foo;
else
    pIdx = [3 2 1];
    pIdx2 = [3 2 1];
    foo = xx;
    xx = yy;
    yy = foo;
end

if get(handles.radiobuttonSelSegment,'value')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Select a segment
    %
    % find graph point closest to clicked point
    h = 20;
    if im.ImageView==1
        [Ig,idxZ]=max(im.III(:,:,zzMin:zzMax),[],3);
    elseif im.ImageView==2
        [Ig,idxZ]=max(im.III(zzMin:zzMax,:,:),[],1);
    else
        [Ig,idxZ]=max(im.III(:,zzMin:zzMax,:),[],2);
    end
    [igr,igc] = find(squeeze(Ig)>=250);
    lst = find(xx>igc-h & xx<igc+h & ...
               yy>igr-h & yy<igr+h );
    for ii=1:length(lst)
        rsep(ii) = norm([igc(lst(ii))-xx  igr(lst(ii))-yy]);
    end    
    [foo,nIdx] = min(rsep);
    nIdx = lst(nIdx);
    yy = igr(nIdx);
    xx = igc(nIdx);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find edge coinciding with choosen graph point
    nodeEdges = im.nodeEdges;
    nodePos = im.nodePos;

    h=20;
    %%% modified by FangQ 2008/07/08  %%%
%     lst=[];
%     for i=1:length(nodeEdges)
%         if(i==99) 
%             disp('here');
%         end
%         bbxmin=min([nodePos(nodeEdges(i,1),:);nodePos(nodeEdges(i,2),:)]);
%         bbxmax=max([nodePos(nodeEdges(i,1),:);nodePos(nodeEdges(i,2),:)]);
%         if(xx>=bbxmin(1) & xx<=bbxmax(1) & yy>=bbxmin(2) & yy<=bbxmax(2) & zz>=zzMin & zz<=zzMax )
%            lst=[lst i];
%         end
%     end

%     for i=1:length(nodeEdges)
%         if(sum((nodePos(nodeEdges(i,1),:)-nodePos(nodeEdges(i,2),:)).^2)>40000)
%            fprintf(1,'bad segment: %d',i);
%         end
%     end
    % David's original code
    lst = find(xx>im.nodePos(:,pIdx(1))-h & xx<im.nodePos(:,pIdx(1))+h & ...
        yy>im.nodePos(:,pIdx(2))-h & yy<im.nodePos(:,pIdx(2))+h & ...
        im.nodePos(:,pIdx(3))>=zzMin & im.nodePos(:,pIdx(3))<=zzMax);

   %%% end of FangQ's modification %%%

    if isempty(lst)
        return;
    end
    clear rsep
    for ii=1:length(lst)
        rsep(ii) = norm([im.nodePos(lst(ii),pIdx(1))-xx  im.nodePos(lst(ii),pIdx(2))-yy]);
    end
% first method
%    [foo,nIdx] = min(rsep);
% second method
%    [foo,sidx] = sort(rsep,2,'ascend');
%    nIdx1b = lst(sidx(1));
%    nIdx2b = lst(sidx(2));
% third method
    minDiff = 999;
    for n1 = 1:length(lst)-1
        for n2 = n1+1:length(lst)
            if ~isempty(find( (nodeEdges(:,1)==lst(n1) & nodeEdges(:,2)==lst(n2)) | (nodeEdges(:,1)==lst(n2) & nodeEdges(:,2)==lst(n1)) ))
                foo = (round(nodePos(lst(n1),pIdx(1:2))) -[xx yy]) ./ (round(nodePos(lst(n1),pIdx(1:2)))-round(nodePos(lst(n2),pIdx(1:2)))+eps);
%                foo = foo / abs(foo(1));
                if abs(foo(1)-foo(2)) < minDiff & foo(1)>=0 & foo(2)>=0 & (foo(1)-eps)<=1 & (foo(2)-eps)<=1
                    minDiff = abs(foo(1)-foo(2));
                    nIdx1 = lst(n1);
                    nIdx2 = lst(n2);
                end            
            end
        end
    end
    
    if minDiff == 999
        error( 'This should not happen!' )
    end

    flag = 0;
    if length(find(nodeEdges(:,1)==nIdx1 | nodeEdges(:,2)==nIdx1))<=2
        nIdx = nIdx1;
        flag = 1;
    elseif length(find(nodeEdges(:,1)==nIdx2 | nodeEdges(:,2)==nIdx2))<=2
        nIdx = nIdx2;
        flag = 1;
    else
        nIdx = nIdx1;
        lst = find( (nodeEdges(:,1)==nIdx1 & nodeEdges(:,2)==nIdx2) | (nodeEdges(:,1)==nIdx2 & nodeEdges(:,2)==nIdx1) );
        if length(lst)==1;
            lstEdgesSel = lst;
        else
            flag = 1;
        end
    end
    
    
    im.nodeSelected = nIdx;

    % find all edges in the segment
    if flag
        lstNodes = [];
        nLst = 0;
        lstEdgesSel = [];

        lstEdges = find(nodeEdges(:,1)==nIdx | nodeEdges(:,2)==nIdx);
        if ~isempty(find(nodeEdges(:,1)==nIdx2 & nodeEdges(:,2)==nIdx1))
            lst1 = find(nodeEdges(:,1)==nIdx2 & nodeEdges(:,2)==nIdx1);
        else
            lst1 = find(nodeEdges(:,1)==nIdx1 & nodeEdges(:,2)==nIdx2);
        end
        lst2 = find(lstEdges~=lst1);
        lst3 = find(lstEdges==lst1);
        lstEdges = lstEdges([lst2 lst3]);
        %    if length(lstEdges)<=2
        while ~isempty(lstEdges)
            if ~isempty(nIdx)
                nLst = nLst + 1;
                lstNodes(nLst) = nIdx;
            end
            eIdx = lstEdges(end);
            lstEdgesSel(end+1) = eIdx;
            lstEdges = lstEdges(1:end-1);
            nIdx = setdiff( nodeEdges(eIdx,:), lstNodes );
            if ~isempty(nIdx)
                lstEtmp = find(nodeEdges(:,1)==nIdx | nodeEdges(:,2)==nIdx);
                if length(lstEtmp)<=2
                    lstEdges = [lstEdges; setdiff(lstEtmp,eIdx)];
                end
            end
        end
        %    end
    end
    

    % select the whole segment if 'extend' or just the edge
    if strcmp(get(gcf,'Selectiontype'),'normal') | strcmp(get(gcf,'Selectiontype'),'open')
        lstEdgesSel = lstEdgesSel(1);
    end
    im.edgeFlag(lstEdgesSel) = xor(im.edgeFlag(lstEdgesSel),1);
    im.segmentSelectedEdges = lstEdgesSel;
    im.segmentSelectedNodes = unique(nodeEdges(lstEdgesSel,:));

    for ii=lstEdgesSel
        pos0 = nodePos(nodeEdges(ii,1),:);
        pos1 = nodePos(nodeEdges(ii,2),:);
        rsep = norm(pos1-pos0);
        if rsep>0
            cxyz = (pos1-pos0) / rsep;
            rstep = 0;
            pos = pos0;
            while rstep<rsep
                im.III(max(round(pos(2)),1),max(round(pos(1)),1),max(round(pos(3)),1)) = uint8(251 - im.edgeFlag(ii));
                pos = pos + cxyz*0.5;
                rstep = rstep + 0.5;
            end
        end
        im.III(round(pos0(2)),round(pos0(1)),round(pos0(3))) = 255;
        im.III(round(pos1(2)),round(pos1(1)),round(pos1(3))) = 255;
    end

    updateAxes( handles );
    
    set(handles.uipanelSegments,'visible','on');
    
    vSel = unique(nodeEdges(lstEdgesSel,:));
    if ~isempty(vSel)
        foos = '';
        if isfield(im,'nodeDiam')
            foos = sprintf( 'Diam (min avg max), nE\n%.1f  %.1f  %.1f, %d', ...
                min(im.nodeDiam(vSel)), ...
                mean(im.nodeDiam(vSel)), ...
                max(im.nodeDiam(vSel)), length(lstEdgesSel) );
        end
        if isfield( im,'nodePressure' )
            foos = sprintf('%s\nVel %.1f, Pres %.1f', foos,...
                mean(im.nodeVel(vSel))/1e3, mean(im.nodePressure(vSel)) );
        end
        if isfield( im, 'segDiam' )
            if length(vSel)>2
                iSeg = median(im.nodeSegN(vSel));
            else
                iSeg = im.nodeSegN(vSel(1));
            end
            foos = sprintf( '%s\nSeg %d: d%.1f', foos, iSeg, im.segDiam(iSeg) );
        end
        set(handles.textSegmentInfo,'string',foos);

        if length(unique(im.nodeType(im.segmentSelectedNodes)))==1
            if unique(im.nodeType(im.segmentSelectedNodes))==1
                set(handles.radiobuttonSegmentCapillary,'value',0);
                set(handles.radiobuttonSegmentVein,'value',0);
                set(handles.radiobuttonSegmentArterial,'value',1);
            elseif unique(im.nodeType(im.segmentSelectedNodes))==2
                set(handles.radiobuttonSegmentArterial,'value',0);
                set(handles.radiobuttonSegmentVein,'value',0);
                set(handles.radiobuttonSegmentCapillary,'value',1);
            elseif unique(im.nodeType(im.segmentSelectedNodes))==3
                set(handles.radiobuttonSegmentArterial,'value',0);
                set(handles.radiobuttonSegmentCapillary,'value',0);
                set(handles.radiobuttonSegmentVein,'value',1);
            else
                set(handles.radiobuttonSegmentArterial,'value',0);
                set(handles.radiobuttonSegmentCapillary,'value',0);
                set(handles.radiobuttonSegmentVein,'value',0);
            end
        else
            set(handles.radiobuttonSegmentArterial,'value',0);
            set(handles.radiobuttonSegmentCapillary,'value',0);
            set(handles.radiobuttonSegmentVein,'value',0);
        end
    else
        set(handles.textSegmentInfo,'string','');
    end

    updateAxes( handles )

    
elseif get(handles.radiobuttonNewNode,'value')
    %%%%%%%%%%%%%
    % New Node
    
    im.nBflag = 1;
    
    if ~isfield(im,'nodePos')
        if im.ImageView==1
            [foo,zzNew] = max(im.III(yy,xx,zzMin:zzMax),[],3);
        elseif im.ImageView==2
            [foo,zzNew] = max(im.III(zzMin:zzMax,yy,xx),[],1);
        else
            [foo,zzNew] = max(im.III(yy,zzMin:zzMax,xx),[],2);
        end
        zzNew = zzNew + zzMin - 1;
        foo = [xx yy zzNew];
        im.nodePos(1,pIdx) = foo;
        im.nodeDiam(1) = 1;
        im.nodeDiamThetaIdx(1) = 1;
        im.nodeVel(1) = 0;
        im.nodeBC(1) = 0;
        im.nodeBCType(1) = 0;
        im.nodeType(1) = 0;
        im.nodeSegN(1) = 0;
        im.nodeEdges = [];
        im.edgeFlag = [];

    elseif strcmp(get(gcf,'Selectiontype'),'extend')
        if im.ImageView==1
            [foo,zzNew] = max(im.III(yy,xx,zzMin:zzMax),[],3);
        elseif im.ImageView==2
            [foo,zzNew] = max(im.III(zzMin:zzMax,yy,xx),[],1);
        else
            [foo,zzNew] = max(im.III(yy,zzMin:zzMax,xx),[],2);
        end
        zzNew = zzNew + zzMin - 1;
        foo = [xx yy zzNew];
        im.nodePos(end+1,pIdx) = foo;
        im.nodeDiam(end+1) = 1;
        im.nodeDiamThetaIdx(end+1) = 1;
        im.nodeVel(end+1) = 0;
        im.nodeBC(end+1) = 0;
        im.nodeBCType(end+1) = 0;
        im.nodeType(end+1) = 0;
        im.nodeSegN(end+1) = 0;
%        im.nodeEdges = [];
%        im.edgeFlag = [];
        
    else

        nNodes = size(im.nodePos,1);
        foo = ones(nNodes,1)*[xx yy zz];
        d = im.nodePos - foo(:,pIdx2);
%        d = sum(d.^2,2).^0.5;             % this calculates distance in 3D
        d = sum(d(:,pIdx(1:2)).^2,2).^0.5; % this calculates distance in 2D projection
        flag = 1;
        while flag
            [foo,idx] = min(d);
            if im.nodePos(idx,pIdx(3))>=zzMin && im.nodePos(idx,pIdx(3))<=zzMax
                flag = 0;
                d(idx) = max(d);
            else
                d(idx) = max(d);
            end
        end
        flag2 = 0;
        if strcmp(get(gcf,'Selectiontype'),'alt')
            flag2 = 1;
            flag = 1;
            while flag
                [foo,idx2] = min(d);
                if im.nodePos(idx2,pIdx(3))>=zzMin && im.nodePos(idx2,pIdx(3))<=zzMax
                    flag = 0;
                    d(idx2) = max(d);
                else
                    d(idx2) = max(d);
                end
            end
        end
        if im.ImageView==1
            [foo,zzNew] = max(im.III(yy,xx,zzMin:zzMax),[],3);
        elseif im.ImageView==2
            [foo,zzNew] = max(im.III(zzMin:zzMax,yy,xx),[],1);
        else
            [foo,zzNew] = max(im.III(yy,zzMin:zzMax,xx),[],2);
        end
        zzNew = zzNew + zzMin - 1;
        im.nodePos(end+1,pIdx) = [xx yy zzNew];
        im.nodeEdges(end+1,:) = [idx nNodes+1];
        im.nodeDiam(end+1) = 1;
        im.nodeDiamThetaIdx(end+1) = 1;
        im.nodeVel(end+1) = 0;
%        if isfield(im,'nodeBC')
            im.nodeBC(end+1) = 0;
            im.nodeBCType(end+1) = 0;
            im.nodeSegN(end+1) = 0;
%        end
        im.nodeType(end+1) = 0;
        im.edgeFlag(end+1) = 0;
        if flag2
            im.nodeEdges(end+1,:) = [idx2 nNodes+1];
            im.edgeFlag(end+1) = 0;
        end

        pos0 = im.nodePos(idx,:);
        foo = [xx yy zzNew];
        pos1(pIdx) = foo;
        rsep = norm(pos1-pos0);
        if rsep>0
            cxyz = (pos1-pos0) / rsep;
            rstep = 0;
            pos = pos0;
            while rstep<rsep
                im.III(round(pos(2)),round(pos(1)),round(pos(3))) = 251;
                pos = pos + cxyz*0.5;
                rstep = rstep + 0.5;
            end
            im.III(round(pos0(2)),round(pos0(1)),round(pos0(3))) = 255;
            im.III(round(pos1(2)),round(pos1(1)),round(pos1(3))) = 255;
        end
        if flag2
            pos0 = im.nodePos(idx2,:);
            foo = [xx yy zzNew];
            pos1(pIdx) = foo;
            rsep = norm(pos1-pos0);
            if rsep>0
                cxyz = (pos1-pos0) / rsep;
                rstep = 0;
                pos = pos0;
                while rstep<rsep
                    im.III(round(pos(2)),round(pos(1)),round(pos(3))) = 251;
                    pos = pos + cxyz*0.5;
                    rstep = rstep + 0.5;
                end
            end
            im.III(round(pos0(2)),round(pos0(1)),round(pos0(3))) = 255;
            im.III(round(pos1(2)),round(pos1(1)),round(pos1(3))) = 255;
        end
        set(handles.textNumEdges,'string',sprintf('%d edges',size(im.nodeEdges,1)))
%        set(handles.pushbuttonImaris,'enable','on')
    end

    im.nodeSelected = size(im.nodePos,1);
    
    updateAxes( handles );
    set(handles.uipanelGraph,'title',sprintf('Graph (%d nodes)',size(im.nodePos,1)));
    hold on
    if im.ImageView==1
        h=plot(im.nodePos(end,pIdx(1)),im.nodePos(end,pIdx(2)),'cp');
    else
        h=plot(im.nodePos(end,pIdx(2)),im.nodePos(end,pIdx(1)),'cp');
    end
    hold off
    
    set(handles.pushbuttonNodeUndo,'visible','on');

    
elseif get(handles.radiobuttonSeedPoint,'value')
    %%%%%%%%%%%%%
    % New seed

    if im.ImageView==1
        [foo,zIdx] = max(im.III(yy,xx,zzMin:zzMax),[],3);
    elseif im.ImageView==2
        [foo,zIdx] = max(im.III(zzMin:zzMax,yy,xx),[],1);
    else
        [foo,zIdx] = max(im.III(yy,zzMin:zzMax,xx),[],2);
    end
    zz = zzMin+zIdx-1;

    if ~im.seedDirectionFlag
        if foo>0
            nSeed = nSeed + 1;
            im.nSeed = nSeed;
            im.seed(nSeed).pos(pIdx) = [xx yy zz];
            im.seed(nSeed).pos2 = [];
            im.seed(nSeed).bidirectional = 1;
            im.seedDirectionFlag = 1;
            im.seed(nSeed).seedLaunched = 0;
            im.seed(nSeed).seedProcessed = 0;

            im.seed(nSeed).nTrks2Run = str2num(get(handles.editAutoSeedNtrks,'string'));
            im.seed(nSeed).Hxy = str2num(get(handles.editSeedHxy,'string'));
            im.seed(nSeed).Hz = str2num(get(handles.editSeedHz,'string'));
            im.seed(nSeed).funcVal = str2num(get(handles.editSeedFuncVal,'string'));

            im.seed(nSeed).nTrksRan = 0;
            set(handles.editNtrks,'string',sprintf('%d',im.seed(nSeed).nTrks2Run))

            im.curSeed = nSeed;

            set(handles.uipanelSeed,'title',sprintf('Seed %d of %d',nSeed,nSeed))
            set(handles.textSeedInfo,'string','SET DIRECTION');
            set(handles.checkboxBidirectional,'value',im.seed(im.curSeed).bidirectional);
            set(handles.pushbuttonRunSeeds,'Enable', 'on' );
            set(handles.uipanelRunAll,'visible','on')
        end

        updateAxes( handles )
    else
        im.seed(nSeed).pos2(pIdx) = [xx yy zz];
        im.seed(nSeed).cxyz = im.seed(nSeed).pos2 - im.seed(nSeed).pos;
        im.seedDirectionFlag = 0;

        set(handles.textSeedInfo,'string','');
        set(handles.pushbuttonRunSeeds,'Enable', 'on' );
        set(handles.uipanelRunAll,'visible','on')

        updateAxes( handles )
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbuttonSeedButton_Callback(hObject, eventdata, handles)
global im

if im.curSeed==0
    return
end

if im.seedDirectionFlag==1
    im.nSeed = im.nSeed - 1;
    im.curSeed = im.nSeed;
    im.seedDirectionFlag = 0;
elseif eventdata==1 % decrease current seed number
    im.curSeed = max(im.curSeed-1,1);
elseif eventdata==2 % increase current seed number
    im.curSeed = min(im.curSeed+1,im.nSeed);
elseif eventdata==-1 % decrease current seed number
    im.curSeed = max(im.curSeed-10,1);
elseif eventdata==-2 % increase current seed number
    im.curSeed = min(im.curSeed+10,im.nSeed);
elseif eventdata==3 % delete seed
    ii = menu(sprintf('Delete node %d?',im.curSeed),'Yes','No');
    if ii==1
        for ii=im.curSeed:im.nSeed-1
            im.seed(ii) = im.seed(ii+1);
        end
        im.nSeed = im.nSeed - 1;
        % DON'T FORGET TO MOVE THE TRACK FILES
    end
end

set(handles.uipanelSeed,'title',sprintf('Seed %d of %d',im.curSeed,im.nSeed))
set(handles.textSeedInfo,'string','');
set(handles.checkboxBidirectional,'value',im.seed(im.curSeed).bidirectional);
set(handles.editNtrks,'string',sprintf('%d',im.seed(im.curSeed).nTrks2Run))
set(handles.editSeedHxy,'string',sprintf('%d',im.seed(im.curSeed).Hxy))
set(handles.editSeedHz,'string',sprintf('%d',im.seed(im.curSeed).Hz))
set(handles.editSeedFuncVal,'string',sprintf('%.2f ',im.seed(im.curSeed).funcVal))
if im.seed(im.curSeed).seedProcessed
    set(handles.textSeedInfo,'string',sprintf('Processed %d',im.seed(im.curSeed).nTrksRan));
elseif im.seed(im.curSeed).seedLaunched
    set(handles.textSeedInfo,'string','Seed Launched');
end

z = round(get(handles.slider1,'value'));
nS = str2num(get(handles.editNslices,'string'));
if im.ImageView==1
    zSeed = im.seed(im.curSeed).pos(3);
    zNew = min(max(zSeed-floor(nS/2),1),im.nZ-nS);
elseif im.ImageView==2
    zSeed = im.seed(im.curSeed).pos(2);
    zNew = min(max(zSeed-floor(nS/2),1),im.nY-nS);
else
    zSeed = im.seed(im.curSeed).pos(1);
    zNew = min(max(zSeed-floor(nS/2),1),im.nX-nS);
end
set(handles.slider1,'value',zNew);
set(handles.textSlices,'string',sprintf('Slice %d to %d',zNew,zNew+nS-1))

updateAxes( handles );


% --- Executes on button press in checkboxBidirectional.
function checkboxBidirectional_Callback(hObject, eventdata, handles)
global im
if im.curSeed>0 
    im.seed(im.curSeed).bidirectional = get(handles.checkboxBidirectional,'value');
end




% --- Executes on button press in pushbuttonSave.
function pushbuttonSave_Callback(hObject, eventdata, handles)
global im

[fname fpath] = uiputfile( '*.seed','Save seed file');
if fname==0
    return
end
imView3d_SaveData( [], [], fpath, fname)



% --- Executes on button press in pushbuttonLoad.
function pushbuttonLoad_Callback(hObject, eventdata, handles)
global im

[fname fpath] = uigetfile( '*.seed','Load seed file');
if fname==0
    return
end
load([fpath fname],'-mat','im2');
im = im2;
im.gcf = gcf;

updateImageOverlay(0, handles);
UpdateSeeds( handles )


% copied from opening function
if isfield(im,'gui')
    guis = fieldnames(im.gui);
    for ii=1:length(guis)
        foos = guis{ii};
        try
            if strcmpi(foos(1:4),'push')
                eval( sprintf('set(handles.%s,''enable'',''%s'');',foos,getfield(im.gui,foos)) );
            elseif strcmpi(foos(1:2),'ui')
                set(handles.uipanelRunAll,'visible',im.gui.uipanelRunAll);
            elseif strcmpi(foos(1:4),'menu')
                eval( sprintf('set(handles.%s,''checked'',''%s'');',foos,getfield(im.gui,foos)) );
            elseif strcmpi(foos(1:4),'chec')
                eval( sprintf('set(handles.%s,''value'',''%d'');',foos,getfield(im.gui,foos)) );
            else
                eval( sprintf('set(handles.%s,''string'',''%s'');',foos,getfield(im.gui,foos)) );
            end
        catch
        end
    end
end

if strcmpi(get(handles.menu_Flow_useSegments,'checked'),'on')
    im.flagUseSegments = 1;
else
    im.flagUseSegments = 0;
end

set(handles.vesselGraph,'Name',sprintf('vesselGraph %s',im.filenm))

set(handles.slider1,'max',im.nZ)
set(handles.slider1,'value',1)
set(handles.slider1,'min',1)
set(handles.slider1,'sliderstep',[1/im.nZ .1])

axes( handles.axes1 );
xlim([1 size(im.III,2)]);
ylim([1 size(im.III,1)]);

updateAxes( handles );

set( handles.axes1,'ButtonDownFcn','imView3d(''axes1ButtonDown_Callback'',gcbo,[],guidata(gcbo))' );
set(get(handles.axes1,'children'), ...
    'ButtonDownFcn', 'imView3d(''axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))' );

if strcmpi(get(handles.menu_ViewSeedPanels,'checked'),'off')
    set(handles.uipanel_AutoSeed,'visible','off')
    set(handles.uipanelSeed,'visible','off')
    set(handles.uipanel_SeedConnectivity,'visible','off')
    set(handles.radiobuttonSeedPoint,'visible','off')
    set(handles.uipanelRunAll,'visible','off')
else
    set(handles.uipanel_AutoSeed,'visible','on')
    set(handles.uipanelSeed,'visible','on')
    set(handles.uipanel_SeedConnectivity,'visible','on')
    set(handles.radiobuttonSeedPoint,'visible','on')
    set(handles.uipanelRunAll,'visible','on')
end


% not copied from opening function
set(handles.pushbuttonImageGraph,'enable','on')
if isfield(im,'nodeDiam')
    set(handles.radiobuttonSelSegment,'value',1)
end
if isfield(im,'RunPath')
    set(handles.pushbuttonRunPath,'string',sprintf('Path: %s',im.RunPath));
end

pushbuttonImageGraph_Callback([], [], handles)



function UpdateSeeds( handles )
global im

if im.curSeed>0
    set(handles.uipanelSeed,'title',sprintf('Seed %d of %d',im.curSeed,im.nSeed))
    set(handles.textSeedInfo,'string','');
    set(handles.checkboxBidirectional,'value',im.seed(im.curSeed).bidirectional);
    set(handles.editNtrks,'string',sprintf('%d',im.seed(im.curSeed).nTrks2Run))
    set(handles.editSeedHxy,'string',sprintf('%d',im.seed(im.curSeed).Hxy))
    set(handles.editSeedHz,'string',sprintf('%d',im.seed(im.curSeed).Hz))
    set(handles.editSeedFuncVal,'string',sprintf('%.2f ',im.seed(im.curSeed).funcVal))
    if im.seed(im.curSeed).seedProcessed
        set(handles.textSeedInfo,'string',sprintf('Processed %d',im.seed(im.curSeed).nTrksRan));
    elseif im.seed(im.curSeed).seedLaunched
        set(handles.textSeedInfo,'string','Seed Launched');
    end

    set(handles.checkboxBidirectional,'value',im.seed(im.curSeed).bidirectional);
end

nS = str2num(get(handles.editNslices,'string'));
if ~isfield(im,'ImageView')
    im.ImageView = 1;
end
if im.curSeed>0
    if im.ImageView==1
        zSeed = im.seed(im.curSeed).pos(3);
        zNew = min(max(zSeed-floor(nS/2),1),im.nZ-nS);
        nZ = im.nZ;
    elseif im.ImageView==2
        zSeed = im.seed(im.curSeed).pos(2);
        zNew = min(max(zSeed-floor(nS/2),1),im.nY-nS);
        nZ = im.nY;
    else
        zSeed = im.seed(im.curSeed).pos(1);
        zNew = min(max(zSeed-floor(nS/2),1),im.nX-nS);
        nZ = im.nX;
    end
else
    zNew = 1;
    if im.ImageView==1
        nZ = im.nZ;
    elseif im.ImageView==2
        nZ = im.nY;
    else
        nZ = im.nX;
    end
end
set(handles.slider1,'max',nZ)
set(handles.slider1,'min',1)
set(handles.slider1,'value',zNew);
set(handles.slider1,'sliderstep',[1/nZ .1])
set(handles.textSlices,'string',sprintf('Slice %d to %d',zNew,zNew+nS-1))

set(handles.vesselGraph,'Name',sprintf('vesselGraph %s',im.filenm))

if isfield(im,'nodePos')
    set(handles.uipanelGraph,'title',sprintf('Graph (%d nodes)',size(im.nodePos,1)));
end

if isfield(im,'nodeEdges')
    set(handles.textNumEdges,'string',sprintf('%d edges',size(im.nodeEdges,1)))
%    set(handles.pushbuttonPruneLoops,'enable','on')
%    set(handles.pushbuttonPruneLoops4,'enable','on');
end

axes( handles.axes1 );
if im.ImageView==1
    xlim([1 size(im.III,2)]);
    ylim([1 size(im.III,1)]);
elseif im.ImageView==2
    xlim([1 size(im.III,3)]);
    ylim([1 size(im.III,1)]);
else
    xlim([1 size(im.III,3)]);
    ylim([1 size(im.III,2)]);
end

updateAxes( handles )




% --- Executes on button press in pushbuttonRunSeeds.
function pushbuttonRunSeeds_Callback(hObject, eventdata, handles)
global im

wd = cd;

%%%%%%%%%%%%%%%%%%
if 0 % run locally
    cd('Seeds')
    
    flag = 0;
    fid = fopen('launchScript.bat','w');
    for ii=1:im.nSeed
        if ~im.seed(ii).seedLaunched & im.seed(ii).nTrksRan==0
            fprintf( fid, sprintf( 'del %s_seed%d.trk\n', im.filenm, ii ) );
        end
    end
    for ii=1:10:im.nSeed
        fprintf( fid, sprintf( 'DTItrack'));
        for jj=ii:min(ii+9,im.nSeed)
            if ~im.seed(jj).seedLaunched
                fprintf( fid, sprintf( ' %s_seed%d', im.filenm, jj ) );

                fid2 = fopen( sprintf('%s_seed%d.inp',im.filenm,jj), 'w' );
                fprintf( fid2, '%d %d %d\n', size(im.III,2), size(im.III,1), size(im.III,3) );
                fprintf( fid2, '%s.img\n', im.filenm );
                fprintf( fid2, '%d %d\n%d\n', im.seed(jj).nTrks2Run, im.seed(jj).nTrksRan,...
                    -im.seed(jj).bidirectional );
                fprintf( fid2, '-563980456\n' );  % RANDOM NUMBER SEED, GENERATE FROM CLOCK!
                fprintf( fid2, '%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\n', ...
                    im.seed(jj).pos(1)-1, im.seed(jj).pos(2)-1, im.seed(jj).pos(3)-1,...
                    im.seed(jj).cxyz(1), im.seed(jj).cxyz(2), im.seed(jj).cxyz(3) );
                fprintf( fid2, '%d %d\n', im.seed(jj).Hxy, im.seed(jj).Hz );
                fprintf( fid2, '%.2f %.2f\n', im.seed(jj).funcVal(1), im.seed(jj).funcVal(2));
                fclose( fid2 );

                im.seed(jj).seedLaunched = 1;
                im.seed(jj).seedProcessed = 0;
                im.seed(jj).seedGraphed = 0;
                flag = 1;
            end
        end
        fprintf( fid, sprintf( '\n') );
    end
    fprintf( fid, 'exit\n' );
    fclose(fid);

    if flag
        !launchScript.bat &

        set(handles.pushbuttonRunSeeds,'Enable', 'off' );
        set(handles.uipanelRunAll,'visible','off')
        set(handles.pushbuttonProcSeeds,'Enable','On')
    end

%%%%%%%%%%%%%%%%%%
else % run on cluster
    % make sure launcher.sh is running on seychelles
    linuxWD = im.RunPath;
    cd(linuxWD);

    flag = 0;
    fid = fopen('delScript.bat','w');
    for ii=1:im.nSeed
        if ispc
            if ~im.seed(ii).seedLaunched & im.seed(ii).nTrksRan==0
                fprintf( fid, sprintf( 'del %s_seed%d.trk\n', im.filenm, ii ) );
            end
        else
            if ~im.seed(ii).seedLaunched & im.seed(ii).nTrksRan==0
                fprintf( fid, sprintf( 'rm %s_seed%d.trk\n', im.filenm, ii ) );
            end
        end
    end
    fclose(fid);
    if isunix
        !chmod +x delScript.bat
    end
    
    fid = fopen('launchScript2.bat','w');
    for ii=1:10:im.nSeed
        fprintf( fid, sprintf( 'DTItrack'));
        for jj=ii:min(ii+9,im.nSeed)
            if ~im.seed(jj).seedLaunched
                fprintf( fid, sprintf( ' %s_seed%d', im.filenm, jj ) );

                fid2 = fopen( sprintf('%s_seed%d.inp',im.filenm,jj), 'w' );
                fprintf( fid2, '%d %d %d\n', size(im.III,2), size(im.III,1), size(im.III,3) );
                fprintf( fid2, '%s.img\n', im.filenm );
                fprintf( fid2, '%d %d\n%d\n', im.seed(jj).nTrks2Run, im.seed(jj).nTrksRan,...
                    -im.seed(jj).bidirectional );
                fprintf( fid2, '-563980456\n' );  % RANDOM NUMBER SEED, GENERATE FROM CLOCK!
                fprintf( fid2, '%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\n', ...
                    im.seed(jj).pos(1)-1, im.seed(jj).pos(2)-1, im.seed(jj).pos(3)-1,...
                    im.seed(jj).cxyz(1), im.seed(jj).cxyz(2), im.seed(jj).cxyz(3) );
                fprintf( fid2, '%d %d\n', im.seed(jj).Hxy, im.seed(jj).Hz );
                fprintf( fid2, '%.2f %.2f\n', im.seed(jj).funcVal(1), im.seed(jj).funcVal(2));
                fclose( fid2 );

                im.seed(jj).seedLaunched = 1;
                im.seed(jj).seedProcessed = 0;
                im.seed(jj).seedGraphed = 0;
                flag = 1;
            end
        end
        fprintf( fid, sprintf( '\n') );
    end
    fclose(fid);

    if flag
        !delScript.bat
        save foo.mat flag
        if ispc
            !copy foo.mat dothings
        else
            !cp foo.mat dothings
        end

        set(handles.pushbuttonRunSeeds,'Enable', 'off' );
        set(handles.uipanelRunAll,'visible','off')
        set(handles.pushbuttonProcSeeds,'Enable','On')
        set(handles.pushbuttonRemoveSeeds,'visible','off')
    end
    
    
end

cd(wd);


% --- Executes on button press in pushbuttonProcSeeds.
function pushbuttonProcSeeds_Callback(hObject, eventdata, handles)
global im

wd = cd;
if 0
    cd('Seeds');
else
    linuxWD = im.RunPath;
    cd(linuxWD);
end

hWait = waitbar(0,'Processing seeds');

flag = 0;
flag2 = 1;
for ii=1:im.nSeed
    waitbar(ii/im.nSeed, hWait);
    if im.seed(ii).seedLaunched & ~im.seed(ii).seedProcessed
        filenm = sprintf('%s_seed%d.trk',im.filenm,ii);
        if exist(filenm,'file')
            fid = fopen(filenm,'rb');
            n =  fread(fid,2,'int');
            trk = fread(fid,n(1)*n(2)*3,'int');
            fclose(fid);
            trk = reshape(trk,[3 n(2) n(1)]);
            im.seed(ii).trkPos = permute(trk,[2 1 3]) + 1;
            im.seed(ii).seedProcessed = 1;
            im.seed(ii).seedGraphed = 0;
            im.seed(ii).nTrksRan = n(1);
            flag = 1;
        else
            flag2 = 0;
        end
    end
end

close(hWait);

cd(wd)

updateAxes( handles )
set(handles.uipanelSeed,'title',sprintf('Seed %d of %d',im.curSeed,im.nSeed))
set(handles.textSeedInfo,'string','');
set(handles.checkboxBidirectional,'value',im.seed(im.curSeed).bidirectional);
set(handles.editNtrks,'string',sprintf('%d',im.seed(im.curSeed).nTrks2Run))
set(handles.editSeedHxy,'string',sprintf('%d',im.seed(im.curSeed).Hxy))
set(handles.editSeedHz,'string',sprintf('%d',im.seed(im.curSeed).Hz))
set(handles.editSeedFuncVal,'string',sprintf('%.2f ',im.seed(im.curSeed).funcVal))
if im.seed(im.curSeed).seedProcessed
    set(handles.textSeedInfo,'string',sprintf('Processed %d',im.seed(im.curSeed).nTrksRan));
elseif im.seed(im.curSeed).seedLaunched
    set(handles.textSeedInfo,'string','Seed Launched');
end

if flag2
    %set(handles.pushbuttonRunSeeds,'Enable','On')
    %set(handles.pushbuttonGraph,'Enable','On')
    set(handles.pushbuttonProcSeeds,'Enable','Off')
end
if flag
    set(handles.pushbuttonSeedConnectivity,'Enable','On')
end

% --- Executes on button press in pushbuttonSeedConnectivity.
function pushbuttonSeedConnectivity_Callback(hObject, eventdata, handles)
global im

h = str2num(get(handles.editSeedConnHxy,'string'));
hz = str2num(get(handles.editSeedConnHz,'string'));
nThresh = str2num(get(handles.editSeedConnThresh,'string'));

for ii=1:im.nSeed
    pos(ii,:) = im.seed(ii).pos;
end

hWait = waitbar( 0, 'Calculating seed connectivity...' );

clear M
M = zeros(im.nSeed);
if isfield(im,'SeedConnectivityMatrix')
    M = im.SeedConnectivityMatrix;
    Mgrps = im.ConnectivityGrps;
end

for ii=1:im.nSeed
    waitbar(ii/im.nSeed, hWait);
    if im.seed(ii).seedProcessed & ~im.seed(ii).seedGraphed

        im.seed(ii).seedGraphed = 0;

        im.seed(ii).trkConn = zeros(1,im.seed(ii).nTrksRan);
        xTrks = squeeze(im.seed(ii).trkPos(:,1,:));
        yTrks = squeeze(im.seed(ii).trkPos(:,2,:));
        zTrks = squeeze(im.seed(ii).trkPos(:,3,:));
        for jj=1:im.nSeed
            if jj~=ii
                [lstR,lstC] = find(xTrks>=pos(jj,1)-h & xTrks<=pos(jj,1)+h & ...
                    yTrks>=pos(jj,2)-h & yTrks<=pos(jj,2)+h & ...
                    zTrks>=pos(jj,3)-hz & zTrks<=pos(jj,3)+hz );
                lstTrks = unique(lstC);
                if ~isempty(lstTrks)
                    M(ii,jj) = length(lstTrks);
                    if length(lstTrks)>=nThresh
                        for kk=1:length(lstTrks)
                            ll = find(lstC==lstTrks(kk));
                            [foo,mm] = max(lstR(ll));  % maybe this should be min
                            if lstR(ll(mm)) > im.seed(ii).trkConn(lstTrks(kk))
                                im.seed(ii).trkConn(lstTrks(kk)) = lstR(ll(mm));
                            end
                        end
                    end
                else
                    M(ii,jj) = 0;
                end
                % I still need to record the last trk pt (lstR) that connects
                % with a seed to cut graphing at that point.
            else
                M(ii,jj) = -1;
                if nThresh==0
                    im.seed(ii).trkConn = size(im.seed(ii).trkPos,1)*ones(1,size(im.seed(ii).trkPos,3));
                end
            end
        end
    end
end
close(hWait)

% find the number of unique groups of connected seeds
Mgrps = zeros(1,im.nSeed);
nGrps = 0;
for ii=1:size(M,1)
    if Mgrps(ii)==0
        lst = unique([find(M(ii,:)>0) (find(M(:,ii)>0))']);
        nGrps = nGrps + 1;
        Mgrps(ii) = nGrps;
        Mgrps(lst) = nGrps;
        
        while ~isempty(lst)
            jj=lst(end);
            lst = lst(1:end-1);
            lst2 = unique([find(M(jj,:)>0) (find(M(:,jj)>0))']);
            lst3 = find(Mgrps>0);
            Mgrps(lst2) = nGrps;
            lst = unique([lst setdiff(lst2,lst3)]);
        end
    end
end

im.SeedConnectivityMatrix = M;
im.ConnectivityGrps = Mgrps;

set(handles.pushbuttonGraph,'Enable','On')
set(handles.pushbuttonSeedConnectivity,'enable','on');

%M
Mgrps




% --- Executes on button press in pushbuttonGraph.
function pushbuttonGraph_Callback(hObject, eventdata, handles)
global im


% Create the graph
if ~isfield(im,'nodePos')
    nNodes = 0;
    nEdges = 0;
    nodePos = zeros(10000,3);%[0 0 0];
    nodeEdges = zeros(50000,2);%[0 0];
else
    nodePos = im.nodePos;
    nNodes = size(nodePos,1);
    nodeEdges = im.nodeEdges;
    nEdges = size(nodeEdges,1);
    nodePos = [nodePos; zeros(nNodes,3)];
    nodeEdges = [nodeEdges; zeros(nEdges,2)];
end



hWait = waitbar( 0, 'Graphing tracks...' );

hxy = str2num(get(handles.editGraphHxy,'string'));
hz = str2num(get(handles.editGraphHz,'string'));
lstSeedNodes = [];
nEdgesTossed = 0;
for ii=1:im.nSeed
    waitbar(ii/im.nSeed,hWait);
    if im.seed(ii).seedProcessed & ~im.seed(ii).seedGraphed
        [ii-1 nNodes nEdges nEdgesTossed]
        nTrkPosMax = size(im.seed(ii).trkPos,1);
        nTrks = size(im.seed(ii).trkPos,3);
        for jj=1:nTrks
            kk = 1;
            lastNode = 0;
            lastPtSum = 0;
            lastPos = 0;
            while im.seed(ii).trkPos(kk,1,jj)>0 & kk<=im.seed(ii).trkConn(jj) & kk<nTrkPosMax
                pos = squeeze(im.seed(ii).trkPos(kk,:,jj));
                lst = find(pos(1)>=(nodePos(1:nNodes,1)-hxy) & pos(1)<=(nodePos(1:nNodes,1)+hxy) & ...
                    pos(2)>=(nodePos(1:nNodes,2)-hxy) & pos(2)<=(nodePos(1:nNodes,2)+hxy) & ...
                    pos(3)>=(nodePos(1:nNodes,3)-hz) & pos(3)<=(nodePos(1:nNodes,3)+hz) );
                if isempty(lst)
                    nNodes = nNodes+1;
                    nodePos(nNodes,:) = pos;
                    if jj==1 & kk==1
                        lstSeedNodes(end+1) = nNodes;
                    end
                    if norm(pos-lastPos)<2*hxy & lastNode>0
%                        if lastNode>0
%                            nodeConn{lastNode}.next(end+1) = nNodes;
%                        end
                        nEdges = nEdges + 1;
                        nodeEdges(nEdges,:) = [lastNode nNodes];
%                        nodeConn{nNodes}.prev(1) = lastNode;
%                        nodeConn{nNodes}.next = [];
                    elseif lastNode>0
                        %need to find closest node and check if it follows back
                        %to last node. IF NOT THEN USE LAST NODE!!!
                        d=sum((ones(nNodes-1,1)*pos-nodePos(1:nNodes-1,:)).^2,2).^0.5;
                        [foo,idx]=min(d);
                        if ~isempty(idx)
                            lastNode = idx;
                            nEdges = nEdges + 1;
                            nodeEdges(nEdges,:) = [lastNode nNodes];
%                            nodeConn{lastNode}.next(end+1) = nNodes;
%                            nodeConn{nNodes}.prev(1) = lastNode;
%                            nodeConn{nNodes}.next = [];
                        else
                            disp('HELP!!!')
                            keyboard
                        end
                    end
                    lastNode = nNodes;    
                    lastPos = pos;
                    lastPt = kk;
                    lastPtSum = 0;
                else
                    % This code is behaving better now that I added
                    % nodesConnected() and I am not adding redundant nodes
                    if length(lst)>1
                        clear d
                        for iLst=1:length(lst)
                            d(iLst) = norm(pos-nodePos(lst(iLst),:));
                        end
                        [foo, closestNode] = min(d);
                    else
                        closestNode = 1;
                    end
                    newNode = lst(closestNode);
                    if newNode~=lastNode & lastNode~=0 
                        if isempty(find((nodeEdges(1:nEdges,1)==lastNode & nodeEdges(1:nEdges,2)==newNode) | ...
                                        (nodeEdges(1:nEdges,2)==lastNode & nodeEdges(1:nEdges,1)==newNode) ))
                            flag = nodesConnected( nodeEdges(1:nEdges,:), newNode, lastNode, 5 );
%                            flag = 0;
                            if ~flag
                                nEdges = nEdges + 1;
                                if lastNode<newNode
                                    nodeEdges(nEdges,:) = [lastNode newNode];
                                else
                                    nodeEdges(nEdges,:) = [newNode lastNode];
                                end
                            else
                                % I could iterate this until nB(ii)~=1
                                [ir,ic] = find(nodeEdges(1:nEdges,:)==lastNode);
                                if length(ir)==1 & isempty(find(lastNode==lstSeedNodes))
                                    foo = ones(nEdges,1);
                                    foo(ir) = 0;
                                    nEdges = nEdges - 1;
                                    nodeEdges(1:nEdges,:) = nodeEdges(find(foo==1),:);

                                    foo = ones(nNodes,1);
                                    foo(lastNode) = 0;
                                    nNodes = nNodes - 1;
                                    nodePos(1:nNodes,:) = nodePos(find(foo==1),:);

                                    lst = find(nodeEdges>lastNode);
                                    nodeEdges(lst) = nodeEdges(lst) - 1;

                                    lst = find(lstSeedNodes>lastNode);
                                    if ~isempty(lst)
                                        lstSeedNodes(lst) = lstSeedNodes(lst) - 1;
                                    end
                                    
                                    if newNode>lastNode
                                        newNode = newNode - 1;
                                    end

                                    nEdgesTossed = nEdgesTossed + 1;
                                end
                            end
                        else
                            % if the lastNode has only 1 edge and that edge 
                            % is between newNode and lastnode then I can
                            % delete it. 
                            % What if edge is not between newNode and
                            % lastNode? This needs to go above
                            [ir,ic] = find(nodeEdges(1:nEdges,:)==lastNode);
                            if length(ir)==1 & isempty(find(lastNode==lstSeedNodes))
                                foo = ones(nEdges,1);
                                foo(ir) = 0;
                                nEdges = nEdges - 1;
                                nodeEdges(1:nEdges,:) = nodeEdges(find(foo==1),:);
                                
                                foo = ones(nNodes,1);
                                foo(lastNode) = 0;
                                nNodes = nNodes - 1;
                                nodePos(1:nNodes,:) = nodePos(find(foo==1),:);
                                
                                lst = find(nodeEdges>lastNode);
                                nodeEdges(lst) = nodeEdges(lst) - 1;
                                
                                lst = find(lstSeedNodes>lastNode);
                                if ~isempty(lst)
                                    lstSeedNodes(lst) = lstSeedNodes(lst) - 1;
                                end
                                
                                if newNode>lastNode
                                    newNode = newNode - 1;
                                end

                                nEdgesTossed = nEdgesTossed + 1;
                            end
                        end
%                        nodeConn{lastNode}.next(end+1) = lst(closestNode);
%                        nodeConn{lst(closestNode)}.next(end+1) = lastNode;

%                        nodeConn{lastNode}.next = setdiff(unique(nodeConn{lastNode}.next),[nodeConn{lastNode}.prev lastNode]);
%                        nodeConn{lst(closestNode)}.next = setdiff(unique(nodeConn{lst(closestNode)}.next),[nodeConn{lst(closestNode)}.prev lst(closestNode)]);
                    end
                    lastNode = newNode;
                    lastPos = nodePos(lastNode,:);
                end
                kk = kk + 1;
            end
        end
        
        im.seed(ii).seedGraphed = 1;
    end
end

nodePos = nodePos(1:nNodes,:);
nodeEdges = nodeEdges(1:nEdges,:);

%%%%%%%%%%%%%%
% prune edges - still need to handle small loops
%
% point edges
nodeEdges = nodeEdges(find(nodeEdges(:,1)~=nodeEdges(:,2)),:);
% redundant edges
sE = cell(size(nodeEdges,1),1);
for ii=1:length(nodeEdges)
    sE{ii} = sprintf('%05d%05d',nodeEdges(ii,1),nodeEdges(ii,2));
end
[b,i,j]=unique(sE);
nodeEdges = nodeEdges(sort(i),:);


close(hWait);

im.nodePos = nodePos;
im.nodeEdges = nodeEdges;
im.edgeFlag = zeros(1,size(nodeEdges,1));
im.nodeBC = zeros(1,size(nodePos,1));
im.nodeBCType = zeros(1,size(nodePos,1));
im.nodeSegN = zeros(1,size(nodePos,1));
im.nodeType = zeros(1,size(nodePos,1));
im.nodeVel = zeros(1,size(nodePos,1));

set(handles.pushbuttonGraph,'Enable','Off')
set(handles.pushbuttonImageGraph,'Enable','On')
%set(handles.pushbuttonImaris,'Enable','On')
set(handles.uipanelGraph,'title',sprintf('Graph (%d nodes)',nNodes));
set(handles.textNumEdges,'string',sprintf('%d edges',size(im.nodeEdges,1)))
%set(handles.pushbuttonPruneLoops,'enable','on')
%set(handles.pushbuttonPruneLoops4,'enable','on');




% --- Executes on button press in pushbuttonImageGraph.
function pushbuttonImageGraph_Callback(hObject, eventdata, handles)
global im

updateImageOverlay(0, handles);

hWait = waitbar( 0, 'Imaging graph...' );

nodePos = im.nodePos;
nodeEdges = im.nodeEdges;
if isfield(im,'edgeFlag')
    edgeFlag = im.edgeFlag;
else
    edgeFlag = zeros(size(nodeEdges,1),1);
end
nNodes = size(nodePos,1);
nEdges = size(nodeEdges,1);

if get(handles.checkboxDisplayGrps,'value') & isfield(im,'nodeGrp')
    if ~get(handles.checkboxHighlightGrp,'value')
        grp = im.nodeGrp;
    else
        grp = 2*ones(nNodes,1);
        grp(find(im.nodeGrp==str2num(get(handles.editHighlightGrp,'string')))) = 1;
    end
else
    grp = ones(nNodes,1);
end

for ii=1:nEdges
    waitbar(ii/nEdges,hWait);
    pos0 = max(nodePos(nodeEdges(ii,1),:),1);
    pos1 = max(nodePos(nodeEdges(ii,2),:),1);
    rsep = norm(pos1-pos0);
    if rsep>0
        cxyz = (pos1-pos0) / rsep;
        rstep = 0;
        pos = pos0;
        while rstep<rsep
            im.III(round(pos(2)),round(pos(1)),max(round(pos(3)),1)) = min(250 - edgeFlag(ii) + grp(nodeEdges(ii,1)),254);
            pos = pos + cxyz*0.5;
            if pos(1)<2 & pos(2)<2
                keyboard
            end
            rstep = rstep + 0.5;
        end
    end
    im.III(round(pos0(2)),round(pos0(1)),round(pos0(3))) = 255;
    im.III(round(pos1(2)),round(pos1(1)),round(pos1(3))) = 255;
end

close(hWait);
updateAxes( handles );

set(handles.pushbuttonImageGraph,'enable','off')
im.Stats.E1nodeIdx = [];
set(handles.textE1info,'string','');
set(handles.pushbuttonSelBranch,'enable','off')
set(handles.pushbuttonPlotTree,'enable','off')
set(handles.textStatInfo,'string','');

if get(handles.checkboxDisplayGrps,'value') & isfield(im,'nodeGrp') & get(handles.checkboxHighlightGrp,'value')
    grpN = str2num(get(handles.editHighlightGrp,'string'));
    ch = menu(sprintf('Move to plane containing group %d?',grpN),'Yes','No');
    if ch==1
        if im.ImageView==1
            pIdx = [1 2 3];
            xlim([1 size(im.III,2)]);
            ylim([1 size(im.III,1)]);
            nIs = size(im.III,3);
        elseif im.ImageView==2
            pIdx = [1 3 2];
            xlim([1 size(im.III,2)]);
            ylim([1 size(im.III,3)]);
            nIs = size(im.III,2);
        else
            pIdx = [2 3 1];
            xlim([1 size(im.III,1)]);
            ylim([1 size(im.III,3)]);
            nIs = size(im.III,1);
        end
        nS = str2num(get(handles.editNslices,'string'));

        nIdx = find(im.nodeGrp==grpN);
        if ~isempty(nIdx)
            zz = min(max(mean(im.nodePos(nIdx,pIdx(3)),1)-nS/2,1),nIs-nS);
            zz = min(zz,get(handles.slider1,'Max'));
            set(handles.slider1,'value',zz);
            set(handles.textSlices,'string',sprintf('Slice %d to %d',round(zz),round(zz+nS-1)))

            updateAxes(handles);
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateImageOverlay( withTracks, handles )
global im

hWait = waitbar( 0,'Loading Original Image' );

if ~isfield(im,'I')
    if strcmp(im.filetype,'img')
        filenm = sprintf( '%s.img',im.filenm );
        fid = fopen( filenm, 'rb' );
        n = fread( fid, 3, 'int');
        nx = n(1); ny = n(2); nz = n(3);
        I = fread( fid, nx*ny*nz, 'float' );
        I = reshape(I, [ny nx nz] );
        fclose(fid);
    elseif strcmp(im.filetype,'tif') | strcmp(im.filetype,'tiff')
        finfo = imfinfo(sprintf('%s.%s',im.filenm,im.filetype));
        ny = finfo(1).Height;
        nx = finfo(1).Width;
        nz = length(finfo);
        I = zeros(ny,nx,nz);
        for ii=1:nz
            I(:,:,ii) = imread(sprintf('%s.%s',im.filenm,im.filetype),ii);
        end
        I = I / max(I(:));
    end
    
    I = I-0.0;
    I = 32*I;
    I = I/(1-0.0);
    
    im.I = uint8(I);
    
    im.nX = nx;
    im.nY = ny;
    im.nZ = nz;
else
    I = uint8(im.I);
end

if ~isfield(im,'III')
    im.III = uint8(zeros(size(I)));
end
im.III = max(uint8(I),1);
%if withTracks
%    lst = find(im.II>0);
%    im.III(lst) = 32+32*min(im.II(lst),4)/4;
%end

% overlay Vessel Mask or Velocity
if get(handles.radiobuttonOverlayVesselType,'value');
    if isfield(im,'Vm')
        im.III = im.III + (im.Vm-1)*32;
    end
elseif get(handles.radiobuttonOverlayVel,'value');
    if isfield(im,'Vvel');
        lst = find(im.Vvel>0);
        im.III(lst) = 160 + im.Vvel(lst);
    end
elseif get(handles.radiobuttonOverlayFlow,'value');
    if isfield(im,'Vflow');
        lst = find(im.Vflow>0);
        im.III(lst) = 160 + im.Vflow(lst);
    end
elseif get(handles.radiobuttonOverlayPressure,'value');
    if isfield(im,'Vpres');
        lst = find(im.Vpres>0);
        im.III(lst) = 160 + im.Vpres(lst);
    end
end

close(hWait);







% --- Executes on button press in pushbuttonImaris.
function pushbuttonImaris_Callback(hObject, eventdata, handles)
global im

% start Imaris
if ~isfield(im,'vImarisApplication')
    ch=1;
else
    ch = menu('Start Imaris?','Yes','No');
end
if ch==1
    vImarisApplication = actxserver('Imaris.Application');
    vImarisApplication.mVisible = true;
else
    vImarisApplication = im.vImarisApplication;
end
im.vImarisApplication = vImarisApplication;

% Create the data container
ch = menu('Make sure volume is loaded and selected and then hit okay','okay','cancel');
if ch==2
    return;
end
vVolume = vImarisApplication.mFactory.ToVolume(vImarisApplication.mSurpassSelection);
vParent = vImarisApplication.mFactory.CreateDataContainer;
vParent.mName = 'Graph';
vVolume.GetParent.AddChild(vParent);


% get image dimensions
nX=vImarisApplication.mDataSet.get.mSizeX;
nY=vImarisApplication.mDataSet.get.mSizeY;
nZ=vImarisApplication.mDataSet.get.mSizeZ;

hx=(vImarisApplication.mDataSet.get.mExtendMaxX - vImarisApplication.mDataSet.get.mExtendMinX)/nX;
hy=(vImarisApplication.mDataSet.get.mExtendMaxY - vImarisApplication.mDataSet.get.mExtendMinY)/nY;
hz=(vImarisApplication.mDataSet.get.mExtendMaxZ - vImarisApplication.mDataSet.get.mExtendMinZ)/nZ;

% create new filament
vNewFilament = vImarisApplication.mFactory.CreateFilament;
vNewFilament.mIndexT = vImarisApplication.mVisibleIndexT;
vNewFilament.mName = 'foo';
vNewFilament.SetColor(1,0,0, 0.0);
vParent.AddChild(vNewFilament);

nodePos = im.nodePos;
%nodeConn = im.nodeConn;
%nNodes = size(nodePos,1);
scl = [hx hy hz];
imgMin = [vImarisApplication.mDataSet.get.mExtendMinX vImarisApplication.mDataSet.get.mExtendMinY vImarisApplication.mDataSet.get.mExtendMinZ];
pos = (nodePos-1) .* (ones(size(nodePos,1),1)*scl) + ones(size(nodePos,1),1)*imgMin;
rad = ones(size(pos,1),1);
if isfield(im,'nodeDiam')
    rad = (im.nodeDiam/2);% * scl;
end
edges = im.nodeEdges - 1;
%eIdx = 0;
%for ii=1:nNodes
%	for jj=1:length(nodeConn{ii}.next)
%        eIdx = eIdx + 1;
%        edges(eIdx,:) = [ii-1 nodeConn{ii}.next(jj)-1];
%    end
%end

vNewFilament.Set( pos, rad, edges )

% Manually process graph and then get
ch = menu('Manually Center filament and get Diameters. Press Okay when done.','Okay','Cancel');
if ch==2
    return;
end
[posNew, radNew, edgesNew] = vNewFilament.Get;
im.nodePos = (posNew-ones(size(nodePos,1),1)*imgMin)./(ones(size(nodePos,1),1)*scl) + 1;
%im.nodeConn = nodeConn;  % This may be in need of modifying to be like edges
im.nodeDiam = radNew*2;% / scl;
im.nodeEdges = edgesNew + 1;
im.Hvox = [hx hy hz];

set(handles.pushbuttonRegraphNodes,'Enable','On')
set(handles.pushbuttonImageGraph,'Enable','On')
%set(handles.pushbuttonPruneDiam,'Enable','On')
%set(handles.pushbuttonImaris,'Enable','Off')
set(handles.pushbuttonCalcFlowCircuitEq,'enable','on')






% --- Executes on button press in pushbuttonRegraphNodes.
function pushbuttonRegraphNodes_Callback(hObject, eventdata, handles)
global im

% Re-graph after Imaris tuning
nodePos = im.nodePos;
nodeEdges = im.nodeEdges;
nNodes = size(nodePos,1);
nEdges = size(nodeEdges,1);
if isfield(im,'nodeDiam')
    nodeDiam = im.nodeDiam;
else
    nodeDiam = zeros(nNodes,1);
end

hxy = str2num(get(handles.editGraphHxy,'string'));
hz = str2num(get(handles.editGraphHz,'string'));

nNodesUnique = 1;
nodeMap = zeros(nNodes,1);
nodeUnique = zeros(nNodes,1);

nodeMap(1) = 1;
nodePosNew = nodePos(1,:);
nodeUnique(1) = 1;

hwait = waitbar(0,'Regraphing nodes...');
for ii=2:nNodes
    waitbar(ii/nNodes,hwait);
    pos = nodePos(ii,:);
    lst = find(pos(1)>=(nodePosNew(:,1)-hxy) & pos(1)<=(nodePosNew(:,1)+hxy) & ...
        pos(2)>=(nodePosNew(:,2)-hxy) & pos(2)<=(nodePosNew(:,2)+hxy) & ...
        pos(3)>=(nodePosNew(:,3)-hz) & pos(3)<=(nodePosNew(:,3)+hz) );
    if isempty(lst)
        nNodesUnique = nNodesUnique+1;

        nodeMap(ii) = nNodesUnique;
        nodeUnique(ii) = 1;
        
        nodePosNew(nNodesUnique,:) = pos;
        nodeDiamNew(nNodesUnique) = nodeDiam(ii);
    else
        if length(lst)>1
            clear d
            for iLst=1:length(lst)
                d(iLst) = norm(pos-nodePosNew(lst(iLst),:));
            end
            [foo, closestNode] = min(d);
        else
            closestNode = 1;
        end
        nodeMap(ii) = lst(closestNode);
    end
end
close(hwait);

nodeEdgesNew = nodeMap(nodeEdges);
nodeEdges = nodeEdgesNew;

%%%%%%%%%%%%%%
% prune edges - still need to handle small loops
% point edges
nodeEdges = nodeEdges(find(nodeEdges(:,1)~=nodeEdges(:,2)),:);
% redundant edges
sE = cell(size(nodeEdges,1),1);
for ii=1:length(nodeEdges)
    if nodeEdges(ii,1)<nodeEdges(ii,2)
        sE{ii} = sprintf('%05d%05d',nodeEdges(ii,1),nodeEdges(ii,2));
    else
        sE{ii} = sprintf('%05d%05d',nodeEdges(ii,2),nodeEdges(ii,1));
    end
end
[b,i,j]=unique(sE);
nodeEdges = nodeEdges(sort(i),:);

im.nodePos = nodePosNew;
im.nodeEdges = nodeEdges;
im.nodeDiam = nodeDiamNew;
im.edgeFlag = zeros(size(nodeEdges,1),1);

disp(sprintf('Regraph reduced %d nodes to %d, and %d loops to %d',nNodes,size(im.nodePos,1),...
    nEdges-nNodes+1,size(im.nodeEdges,1)-size(im.nodePos,1)+1))

set(handles.pushbuttonRegraphNodes,'Enable','Off')
set(handles.pushbuttonImageGraph,'Enable','On')
set(handles.uipanelGraph,'title',sprintf('Graph (%d nodes)',size(im.nodePos,1)));
set(handles.textNumEdges,'string',sprintf('%d edges',size(im.nodeEdges,1)))
%set(handles.pushbuttonPruneLoops,'enable','on')
%set(handles.pushbuttonPruneLoops4,'enable','on');







function editNtrks_Callback(hObject, eventdata, handles)
global im

if ~isfield(im,'seed')
    return;
end

n2Run = str2double(get(hObject,'String'));

if n2Run > im.seed(im.curSeed).nTrksRan
    im.seed(im.curSeed).nTrks2Run = n2Run;
    im.seed(im.curSeed).seedLaunched = 0;
    set(handles.pushbuttonRunSeeds,'enable','on')
    set(handles.uipanelRunAll,'visible','on')
else
    set(handles.editNtrks,'string',num2str(im.seed(im.curSeed).nTrks2Run))
end









function editSeedConn_Callback(hObject, eventdata, handles)
global im

set(handles.pushbuttonSeedConnectivity,'enable','on');



% --- Executes on button press in pushbuttonAutoSeed.
function pushbuttonAutoSeed_Callback(hObject, eventdata, handles)
global im

if im.ImageView~=1
    menu('Auto Seed only works for XY image view!','Okay');
    return;
end


% filenm = sprintf('%s.img',im.filenm);
% fid = fopen( filenm, 'rb' );
% n = fread( fid, 3, 'int');
% nx = n(1); ny = n(2); nz = n(3);
% I = fread( fid, nx*ny*nz, 'float' );
% I = reshape(I, [ny nx nz] );
% fclose(fid);
% 

hWait = waitbar(0,'Loading original image');

if strcmp(im.filetype,'img')
    filenm = sprintf( '%s.img',im.filenm );
    fid = fopen( filenm, 'rb' );
    n = fread( fid, 3, 'int');
    nx = n(1); ny = n(2); nz = n(3);
    I = fread( fid, nx*ny*nz, 'float' );
    I = reshape(I, [ny nx nz] );
    fclose(fid);
elseif strcmp(im.filetype,'tif') | strcmp(im.filetype,'tiff')
    finfo = imfinfo(sprintf('%s.%s',im.filenm,im.filetype));
    ny = finfo(1).Height;
    nx = finfo(1).Width;
    nz = length(finfo);
    I = zeros(ny,nx,nz);
    for ii=1:nz
        I(:,:,ii) = imread(sprintf('%s.%s',im.filenm,im.filetype),ii);
    end
    I = I / max(I(:));
elseif strcmp(im.filetype,'mat')
    load(im.filenm,'I');
    
    [ny,nx,nz] = size(I);
    I = I / max(I(:));
end

close(hWait);

if get(handles.radiobuttonAutoSeedVisible,'value')
    zzMin = round(get(handles.slider1,'value'));
    nS = str2num(get(handles.editNslices,'string'));
    zzMax = zzMin+nS-1;
    I = I(:,:,zzMin:zzMax);
    nz = size(I,3);
else
    zzMin = 1;
end

h = str2num(get(handles.editAutoSeedHxy,'string'));
hz = str2num(get(handles.editAutoSeedHz,'string'));
Ithresh = str2num(get(handles.editAutoSeedThresh,'string'));
nTrks = str2num(get(handles.editAutoSeedNtrks,'string'));
seedHxy = str2num(get(handles.editSeedHxy,'string'));
seedHz = str2num(get(handles.editSeedHz,'string'));
funcVal = str2num(get(handles.editSeedFuncVal,'string'));

hWait = waitbar( 0, 'Auto Seeding...');

% blank out current seeds
for ii=1:im.nSeed
    if im.seed(ii).seedLaunched
        px = im.seed(ii).pos(1);
        py = im.seed(ii).pos(2);
        pz = im.seed(ii).pos(3)-zzMin+1;
        I( max(py-h,1):min(py+h,ny), max(px-h,1):min(px+h,nx), max(pz-hz,1):min(pz+hz,nz) ) = ...
            0;
    end
end

lst = find(I>Ithresh);
nVox = length(lst);
nMag = floor(log10(nVox));
nSeed = 0;
while ~isempty(lst)
    
    waitbar( min(-log10(length(lst)/nVox),nMag)/nMag, hWait );
    
    [foo,idx] = max(I(lst));
    
    pz = ceil(lst(idx)/(nx*ny));
    px = ceil((lst(idx) - (pz-1)*(nx*ny))/ny);
    py = lst(idx) - (pz-1)*(nx*ny) - (px-1)*ny;
    
    nSeed = nSeed + 1;
    pos(nSeed,:) = [px py pz+zzMin-1];
    
    I( max(py-h,1):min(py+h,ny), max(px-h,1):min(px+h,nx), max(pz-hz,1):min(pz+hz,nz) ) = ...
        0;
    
    lst = find(I>Ithresh);
    [nSeed length(lst)]
end

close( hWait );

[p2, iOrder] = sort(pos(:,3));
ii = 1;
jj = 1;
while ii<=length(iOrder)  % write over seeds not yet launched
    if jj>im.nSeed
        im.seed(jj).seedLaunched = 0;
    end
    if ~im.seed(jj).seedLaunched
        im.seed(jj).pos = pos(iOrder(ii),:);
        im.seed(jj).pos2 = [];
        im.seed(jj).cxyz=[0 0 0];
        im.seed(jj).bidirectional=1;
        im.seed(jj).seedProcessed=0;
        im.seed(jj).seedLaunched=0;
        im.seed(jj).trkPos=[];
        im.seed(jj).nTrks2Run=nTrks;
        im.seed(jj).nTrksRan=0;
        im.seed(jj).Hxy = seedHxy;
        im.seed(jj).Hz = seedHz;
        im.seed(jj).funcVal = funcVal;
        
        ii = ii + 1;
    end
    jj = jj + 1;
end


if ~isfield(im,'curSeed')
    im.curSeed = 1;
elseif im.curSeed == 0;
    im.curSeed = 1;
end
if jj-1>im.nSeed
    im.nSeed = jj-1;
else  % actually it is probably highly unlikely that a seed was launched 
      % out here but seeds before weren't
    ii = jj;
    while ii<=im.nSeed
        if im.seed(ii).seedLaunched
            im.seed(jj) = im.seed(ii);
            jj = jj + 1;
        end
        ii = ii + 1;
    end
    im.nSeed = jj-1;
end

if im.curSeed > im.nSeed
    im.curSeed = 1;
end

set(handles.pushbuttonRunSeeds,'enable','on')
set(handles.uipanelRunAll,'visible','on')
set(handles.pushbuttonRemoveSeeds,'visible','on')

UpdateSeeds( handles )










% --- Executes on button press in pushbuttonSegmentDel.
function pushbuttonSegmentAction_Callback(hObject, eventdata, handles)
global im

if eventdata==1
    % DELETE SEGMENTS
     lst = find(im.edgeFlag==0);
     im.edgeFlag = im.edgeFlag(lst);
     im.nodeEdges = im.nodeEdges(lst,:);
 
     nNodes = size(im.nodePos,1);
     [nB,im] = nBupdate( im );
%      nB=zeros(1,nNodes);
%      for ii=1:nNodes
%          nB(ii)=length(find(im.nodeEdges(:,1)==ii | im.nodeEdges(:,2)==ii));
%      end
     node0 = ones(nNodes,1);
     node0(find(nB==0)) = 0;

     % remove abandoned nodes
    [im.nodePos,im.nodeDiam,im.nodeDiamThetaIdx,im.nodeBC,im.nodeBCType,im.nodeType,im.nodeSegN,im.nodeEdges,im.edgeFlag] = removeNodes( node0, im.nodePos, im.nodeDiam, im.nodeDiamThetaIdx, im.nodeBC, im.nodeBCType, im.nodeType, im.nodeSegN, im.nodeEdges );
    im.nBflag = 1;

    set(handles.pushbuttonUpdateVesselMask,'enable','on')
        
% 
%     nodePos = im.nodePos;
%     nodeEdges = im.nodeEdges;
%     nNodes = size(nodePos,1);
%     nEdges = size(nodeEdges,1);
%     if isfield(im,'nodeDiam')
%         nodeDiam = im.nodeDiam;
%     else
%         nodeDiam = zeros(nNodes,1);
%     end
% 
%     nNodesUnique = 1;
%     nodeMap = zeros(nNodes,1);
%     nodeUnique = zeros(nNodes,1);
% 
%     nodeMap(1) = 1;
%     nodePosNew = nodePos(1,:);  % SHOULDN'T ASSUME nodePos(1,:) is not abandoned
%     nodeUnique(1) = 1;
% 
%     for ii=2:nNodes
%         if ~node0(ii)
%             nNodesUnique = nNodesUnique+1;
% 
%             nodeMap(ii) = nNodesUnique;
%             nodeUnique(ii) = 1;
% 
%             nodePosNew(nNodesUnique,:) = nodePos(ii,:);
%             nodeDiamNew(nNodesUnique) = nodeDiam(ii);
%             nodeBCNew(nNodesUnique) = nodeBC(ii);
%             nodeBCTypeNew(nNodesUnique) = nodeBCType(ii);
%             nodeTypeNew(nNodesUnique) = nodeType(ii);
%         end
%     end
% 
%     nodeEdgesNew = nodeMap(nodeEdges);
%     im.nodeEdges = nodeEdgesNew;
%     im.nodePos = nodePosNew;
%     im.nodeDiam = nodeDiamNew;
% 
    pushbuttonImageGraph_Callback(hObject, eventdata, handles);
    set(handles.textNumEdges,'string',sprintf('%d edges',size(im.nodeEdges,1)))
    set(handles.uipanelGraph,'title',sprintf('Graph (%d nodes)',size(im.nodePos,1)));

elseif eventdata==2
    % UNFLAG SEGMENTS
    lstEdgesSel = find(im.edgeFlag==1)';
    im.edgeFlag(lstEdgesSel) = 0;
    for ii2=1:length(lstEdgesSel)
        ii = lstEdgesSel(ii2);
        pos0 = im.nodePos(im.nodeEdges(ii,1),:);
        pos1 = im.nodePos(im.nodeEdges(ii,2),:);
        rsep = norm(pos1-pos0);
        if rsep>0
            cxyz = (pos1-pos0) / rsep;
            rstep = 0;
            pos = pos0;
            while rstep<rsep
                im.III(max(round(pos(2)),1),max(round(pos(1)),1),max(round(pos(3)),1)) = 251 - im.edgeFlag(ii);
                pos = pos + cxyz*0.5;
                rstep = rstep + 0.5;
            end
            im.III(round(pos0(2)),round(pos0(1)),round(pos0(3))) = 255;
            im.III(round(pos1(2)),round(pos1(1)),round(pos1(3))) = 255;
        end
    end

    updateAxes( handles );

end

set(handles.uipanelSegments,'visible','off');



% --- Executes on button press in pushbuttonViewSegment.
function pushbuttonViewSegment_Callback(hObject, eventdata, handles)
global im

I = im.III;
lst = find(I>=250);
I(lst) = 0;

thresh = ceil(str2num(get(handles.editImageThresh,'string'))*32);
im.TreeThresh = thresh;

nEdges = str2num(get(handles.editViewSegmentN,'string'));
im.TreeNedges = nEdges;
[pMin,pMax,im.Tree] = plotTreeSeg(im.nodeEdges,im.nodePos,im.nodeSelected,nEdges,[1 thresh],I);













function editSeedHxy_Callback(hObject, eventdata, handles)
global im
if im.curSeed>0 
    im.seed(im.curSeed).Hxy = str2num(get(handles.editSeedHxy,'string'));
end


function editSeedHz_Callback(hObject, eventdata, handles)
global im
if im.curSeed>0 
    im.seed(im.curSeed).Hz = str2num(get(handles.editSeedHz,'string'));
end


function editSeedFuncVal_Callback(hObject, eventdata, handles)
global im
if im.curSeed>0 
    im.seed(im.curSeed).funcVal = str2num(get(handles.editSeedFuncVal,'string'));
end



function editImageThresh_Callback(hObject, eventdata, handles)
global im

axes( handles.axes1 );


cm = gray(32);
cm(33:64,:) = hsv2rgb([ones(32,1)*.1 ones(32,1) [1:32]'/32]);
cm(65:255,:) = 0;
cm(65:96,:) = hsv2rgb([ones(32,1)*0 ones(32,1) [1:32]'/32]);
cm(97:128,:) = hsv2rgb([ones(32,1)*0.35 ones(32,1) [1:32]'/32]);
cm(129:160,:) = hsv2rgb([ones(32,1)*0.68 ones(32,1) [1:32]'/32]);
cm(161:192,:) = jet(32);
cm(255,:) = [1 0 0];
cm(251,:) = [1 1 0];
cm(250,:) = [1 0 0];
if get(handles.checkboxDisplayGrps,'value') && isfield(im,'nodeGrp')
    if ~get(handles.checkboxHighlightGrp,'value')
        nGrps = length(unique(im.nodeGrp));
        cm(250+[1:nGrps],:) = jet(nGrps);
    else
        nGrps = 2;
        cm(251,:) = [1 1 0];
        cm(252,:) = [0 1 1];
    end
else
    nGrps = 1;
end
thresh = ceil(str2num(get(handles.editImageThresh,'string'))*32);
cm(1:thresh,2:3) = 0;
colormap(cm)


% cm = gray(32);
% cm(33:64,:) = hsv2rgb([ones(32,1)*0 ones(32,1) [1:32]'/32]);
% cm(65:255,:) = 0;
% cm(65:96,:) = hsv2rgb([ones(32,1)*0 ones(32,1) [1:32]'/32]);
% cm(97:128,:) = hsv2rgb([ones(32,1)*0.35 ones(32,1) [1:32]'/32]);
% cm(129:160,:) = hsv2rgb([ones(32,1)*0.68 ones(32,1) [1:32]'/32]);
% cm(161:192,:) = jet(32);
% cm(255,:) = [1 0 0];
% cm(251,:) = [1 1 0];
% cm(250,:) = [1 0 0];
% if get(handles.checkboxDisplayGrps,'value') && isfield(im,'grpNodes')
%     nGrps = length(unique(im.grpNodes));
%     cm(250+[1:nGrps],:) = jet(nGrps);
% else
%     nGrps = 1;
% end
% thresh = ceil(str2num(get(handles.editImageThresh,'string'))*32);
% cm(1:thresh,2:3) = 0;
% colormap(cm)

% cm = gray(32);
% cm(33:64,:) = hsv2rgb([ones(32,1)*0 ones(32,1) [1:32]'/32]);
% cm(66,:) = [1 1 0];
% cm(65,:) = [1 0 0];
% 
% thresh = ceil(str2num(get(handles.editImageThresh,'string'))*32);
% cm(1:thresh,2:3) = 0;
% 
% colormap(cm)



function radiobuttonAutoSeed_Callback(hObject, eventdata, handles)
if eventdata==1
    set(handles.radiobuttonAutoSeedAll,'value',1);
    set(handles.radiobuttonAutoSeedVisible,'value',0);
else
    set(handles.radiobuttonAutoSeedAll,'value',0);
    set(handles.radiobuttonAutoSeedVisible,'value',1);
end    




function editGraphH_Callback(hObject, eventdata, handles)
global im

if isfield(im,'nodePos')
    set(handles.pushbuttonRegraphNodes,'enable','on')
end


function editNtrksAll_Callback(hObject, eventdata, handles)
global im

n2Run = str2double(get(hObject,'String'));
for ii=1:im.nSeed
    if n2Run > im.seed(ii).nTrksRan
        im.seed(ii).nTrks2Run = n2Run;
        im.seed(ii).seedLaunched = 0;
        set(handles.pushbuttonRunSeeds,'enable','on')
        set(handles.uipanelRunAll,'visible','on')
    end
end
set(handles.editNtrks,'string',sprintf('%d',im.seed(im.curSeed).nTrks2Run))

    


% --- Executes on button press in pushbuttonUpdateStats.
function pushbuttonUpdateStats_Callback(hObject, eventdata, handles)
global im

hwait = waitbar( 0, 'Updating Stats');
nNodes = size(im.nodePos,1);
[nB,im] = nBupdate( im );
% nB=zeros(1,nNodes);
% for ii=1:nNodes
%     nB(ii)=length(find(im.nodeEdges(:,1)==ii | im.nodeEdges(:,2)==ii)); 
% end
set(handles.textEdge1,'string',sprintf('E1: %d',length(find(nB==1))));
set(handles.textEdge2,'string',sprintf('E2: %d',length(find(nB==2))));
set(handles.textEdge3,'string',sprintf('E3: %d',length(find(nB==3))));
set(handles.textEdge4,'string',sprintf('E4+: %d',length(find(nB>=4))));
close(hwait)

ch = menu( sprintf('There are %d E1 nodes. Shall we find those with segment length of 1?',length(find(nB==1))),'Yes','No');
if ch==1
    lstE1 = find(nB==1);
    nodeFlag = ones(nNodes,1);
    edgeFlag2 = ones(size(im.nodeEdges,1),1);
    for ii = 1:length(lstE1)
        nE2B = 0;
        nIdx = lstE1(ii);
        nIdxLst = [];
        eIdx = find(im.nodeEdges(:,1)==nIdx | im.nodeEdges(:,2)==nIdx );
        while length(eIdx)<3 & nE2B<2
            eIdxo = eIdx(1);
            nIdxLst = [nIdxLst nIdx];
            nLst = im.nodeEdges(eIdx,:);
            nIdx = setdiff( nLst(:), nIdxLst );
            if length(nIdx)>1
                warning( 'Problem...' )
                keyboard
            end
            nE2B = nE2B + 1;
            if ~isempty(nIdx)
                eIdx = find(im.nodeEdges(:,1)==nIdx | im.nodeEdges(:,2)==nIdx );
            else
                eIdx = 1:100;
            end
        end
        if nE2B==1
            nodeFlag(lstE1(ii)) = 0;
            edgeFlag2(eIdxo) = 0;
        end
    end

    ch = menu( sprintf('Remove %d E1 nodes with segment length 1?',length(find(nodeFlag==0))), 'Yes','No');
    
    if ch==1
        % remove edges connected to 1 node
        lst = find(edgeFlag2==1);
        im.nodeEdges = im.nodeEdges(lst,:);

        % remove nodes
        [im.nodePos,im.nodeDiam,im.nodeDiamThetaIdx,im.nodeBC,im.nodeBCType,im.nodeType,im.nodeSegN,im.nodeEdges,im.edgeFlag] = removeNodes( nodeFlag, im.nodePos, im.nodeDiam, im.nodeDiamThetaIdx, im.nodeBC, im.nodeBCType, im.nodeType, im.nodeSegN, im.nodeEdges );
        set(handles.pushbuttonUpdateVesselMask,'enable','on')
        im.nBflag = 1;

        % update info in GUI
        set(handles.textNumEdges,'string',sprintf('%d edges',size(im.nodeEdges,1)))
        set(handles.uipanelGraph,'title',sprintf('Graph (%d nodes)',size(im.nodePos,1)));
        set(handles.pushbuttonImageGraph,'enable','on')
        
        nNodes = size(im.nodePos,1);
        [nB,im] = nBupdate( im );
%         nB=zeros(1,nNodes);
%         for ii=1:nNodes
%             nB(ii)=length(find(im.nodeEdges(:,1)==ii | im.nodeEdges(:,2)==ii));
%         end
        set(handles.textEdge1,'string',sprintf('E1: %d',length(find(nB==1))));
        set(handles.textEdge2,'string',sprintf('E2: %d',length(find(nB==2))));
        set(handles.textEdge3,'string',sprintf('E3: %d',length(find(nB==3))));
        set(handles.textEdge4,'string',sprintf('E4+: %d',length(find(nB>=4))));
    end
end


% --- Executes on button press in pushbuttonE1Dec.
function pushbuttonE1_Callback(hObject, eventdata, handles)
global im

if isfield(im,'Stats')
    if ~isfield(im.Stats,'E1Idx')
        im.Stats.E1Idx = 1;
    end
end
E1Idx = im.Stats.E1Idx;

nNodes = size(im.nodePos,1);
[nB,im] = nBupdate( im );
%nB=zeros(1,nNodes);
%for ii=1:nNodes
%    nB(ii)=length(find(im.nodeEdges(:,1)==ii | im.nodeEdges(:,2)==ii)); 
%end
lst = find(nB==1);

if eventdata==1
    E1Idx = E1Idx - 1;
    if E1Idx==0
        E1Idx = length(lst);
    end
else
    E1Idx = E1Idx + 1;
    if E1Idx>length(lst);
        E1Idx = 1;
    end
end
if E1Idx>length(lst);
    E1Idx = length(lst);
end
im.Stats.E1Idx = E1Idx;
im.Stats.E1nodeIdx = lst(E1Idx);
im.nodeSelected = im.Stats.E1nodeIdx;
im.segmentSelectedNodes = im.nodeSelected;

set(handles.textE1info,'string',sprintf('%d\n(%d)',lst(E1Idx),E1Idx))
if ~isfield(im,'nodeBC')
    im.nodeBC = zeros(size(im.nodePos,1),1);
end
set(handles.editE1BC,'string',num2str(im.nodeBC(im.Stats.E1nodeIdx)))

BCtype = im.nodeBCType(im.Stats.E1nodeIdx);
if BCtype==1 || BCtype==3
    set(handles.radiobuttonBCPres,'value',1);
    set(handles.radiobuttonBCVel,'value',0);
elseif BCtype==2 || BCtype==4
    set(handles.radiobuttonBCPres,'value',0);
    set(handles.radiobuttonBCVel,'value',1);
else
    set(handles.radiobuttonBCPres,'value',0);
    set(handles.radiobuttonBCVel,'value',0);
end
if BCtype==3 || BCtype==4
    set(handles.checkboxBClitVal,'value',1);
else
    set(handles.checkboxBClitVal,'value',0);
end

nE2B = 0;
nIdx = lst(E1Idx);
nIdxLst = [];
eIdx = find(im.nodeEdges(:,1)==nIdx | im.nodeEdges(:,2)==nIdx );
while length(eIdx)<3
    nIdxLst = [nIdxLst nIdx];
    nLst = im.nodeEdges(eIdx,:);
    nIdx = setdiff( nLst(:), nIdxLst );
    if length(nIdx)>1
        warning( 'Problem...' )
        keyboard
    end
    nE2B = nE2B + 1;
    if ~isempty(nIdx)
        eIdx = find(im.nodeEdges(:,1)==nIdx | im.nodeEdges(:,2)==nIdx );
    else
        eIdx = 1:100;
    end
end
set(handles.textStatInfo,'string',sprintf('nE2B = %d',nE2B))

nS = str2num(get(handles.editNslices,'string'));

if im.ImageView==1
    pIdx = [1 2 3];
    xlim([1 size(im.III,2)]);
    ylim([1 size(im.III,1)]);
elseif im.ImageView==2
    pIdx = [1 3 2];
    xlim([1 size(im.III,2)]);
    ylim([1 size(im.III,3)]);
else
    pIdx = [2 3 1];
    xlim([1 size(im.III,1)]);
    ylim([1 size(im.III,3)]);
end

ii = lst(E1Idx);
zz = max(im.nodePos(ii,pIdx(3))-nS/2,1);
zz = min(zz,get(handles.slider1,'Max'));
set(handles.slider1,'value',zz);
set(handles.textSlices,'string',sprintf('Slice %d to %d',round(zz),round(zz+nS-1)))

if length(unique(im.nodeType(im.segmentSelectedNodes)))==1
    if unique(im.nodeType(im.segmentSelectedNodes))==1
        set(handles.radiobuttonSegmentCapillary,'value',0);
        set(handles.radiobuttonSegmentVein,'value',0);
        set(handles.radiobuttonSegmentArterial,'value',1);
    elseif unique(im.nodeType(im.segmentSelectedNodes))==2
        set(handles.radiobuttonSegmentArterial,'value',0);
        set(handles.radiobuttonSegmentVein,'value',0);
        set(handles.radiobuttonSegmentCapillary,'value',1);
    elseif unique(im.nodeType(im.segmentSelectedNodes))==3
        set(handles.radiobuttonSegmentArterial,'value',0);
        set(handles.radiobuttonSegmentCapillary,'value',0);
        set(handles.radiobuttonSegmentVein,'value',1);
    else
        set(handles.radiobuttonSegmentArterial,'value',0);
        set(handles.radiobuttonSegmentCapillary,'value',0);
        set(handles.radiobuttonSegmentVein,'value',0);
    end
else
    set(handles.radiobuttonSegmentArterial,'value',0);
    set(handles.radiobuttonSegmentCapillary,'value',0);
    set(handles.radiobuttonSegmentVein,'value',0);
end
set(handles.uipanelSegments,'visible','on');

vSel = im.nodeSelected;
foos = sprintf( 'Diam (min avg max), nE\n%.1f  %.1f  %.1f, %d', ...
    min(im.nodeDiam(vSel)), ...
    mean(im.nodeDiam(vSel)), ...
    max(im.nodeDiam(vSel)), 1 );
if isfield( im,'nodePressure' )
    foos = sprintf('%s\nVel %.1f, Pres %.1f', foos,...
        mean(im.nodeVel(vSel))/1e3, mean(im.nodePressure(vSel)) );
end
if isfield( im, 'segDiam' )
    if length(vSel)>2
        iSeg = median(im.nodeSegN(vSel));
    else
        iSeg = im.nodeSegN(vSel(1));
    end
    foos = sprintf( '%s\nSeg %d: d%.1f', foos, iSeg, im.segDiam(iSeg) );
end

set(handles.textSegmentInfo,'string',foos);

updateAxes(handles);
hold on
h=plot(im.nodePos(ii,pIdx(1)),im.nodePos(ii,pIdx(2)),'cp');
set(h,'markersize',16);
set(h,'markerfacecolor','c')
for ii=1:6
    set(h,'marker','none')
    pause(0.1)
    set(h,'marker','p')
    pause(0.1)
end
hold off

set(handles.pushbuttonSelBranch,'enable','on')
set(handles.pushbuttonPlotTree,'enable','on')


% --- Executes on button press in pushbuttonGroupStats.
function pushbuttonGroupStats_Callback(hObject, eventdata, handles)
global im

im = nodeGrps( im );
set(handles.textGrpStats,'string',sprintf('%d Groups, %d Segments',max(im.nodeGrp),length(im.segNedges)) )
if get(handles.checkboxDisplayGrps,'value');
    pushbuttonImageGraph_Callback(hObject, eventdata, handles);
end



% --- Executes on button press in checkboxDisplayGrps.
function checkboxDisplayGrps_Callback(hObject, eventdata, handles)
set(handles.pushbuttonImageGraph,'enable','on');




% --- Executes on button press in checkboxHighlightGrp.
function checkboxHighlightGrp_Callback(hObject, eventdata, handles)
set(handles.pushbuttonImageGraph,'enable','on');



function editHighlightGrp_Callback(hObject, eventdata, handles)
set(handles.pushbuttonImageGraph,'enable','on');



function editE1BC_Callback(hObject, eventdata, handles)
global im

if ~isfield(im,'nodeBC')
    im.nodeBC = zeros(size(im.nodePos,1),1);
end
im.nodeBC(im.Stats.E1nodeIdx) = str2num(get(handles.editE1BC,'string'));
set(handles.pushbuttonCalcFlowCircuitEq,'enable','on')




% REMOVED 6/5/08
% % --- Executes on button press in pushbuttonSelBranch.
% function pushbuttonSelBranch_Callback(hObject, eventdata, handles)
% global im
% 
% nIdx = im.Stats.E1nodeIdx;
% 
% nE2B = 0;
% nIdxLst = [];
% eIdxLst = [];
% eIdx = find(im.nodeEdges(:,1)==nIdx | im.nodeEdges(:,2)==nIdx );
% while length(eIdx)<3
%     nIdxLst = [nIdxLst nIdx];
%     eIdxLst = [eIdxLst eIdx'];
%     nLst = im.nodeEdges(eIdx,:);
%     nIdx = setdiff( nLst(:), nIdxLst );
%     if length(nIdx)>1
%         warning( 'Problem...' )
%         kayboard
%     end
%     nE2B = nE2B + 1;
%     if ~isempty(nIdx)
%         eIdx = find(im.nodeEdges(:,1)==nIdx | im.nodeEdges(:,2)==nIdx );
%     else
%         eIdx = 1:100;
%     end
% end
% 
% lstEdgesSel = eIdxLst;
% im.edgeFlag(lstEdgesSel) = 1;%xor(im.edgeFlag(lstEdgesSel),1);
% 
% for ii=lstEdgesSel
%     pos0 = im.nodePos(im.nodeEdges(ii,1),:);
%     pos1 = im.nodePos(im.nodeEdges(ii,2),:);
%     rsep = norm(pos1-pos0);
%     if rsep>0
%         cxyz = (pos1-pos0) / rsep;
%         rstep = 0;
%         pos = pos0;
%         while rstep<rsep
%             im.III(round(pos(2)),round(pos(1)),round(pos(3))) = 251 - im.edgeFlag(ii);
%             pos = pos + cxyz*0.5;
%             rstep = rstep + 0.5;
%         end
%     end
% end
% 
% updateAxes( handles );
% 
% set(handles.uipanelSegments,'visible','on');
% 
% vSel = unique(im.nodeEdges(lstEdgesSel,:));
% if ~isempty(vSel)
%     foos = sprintf( 'diameter (min avg max), nE\n%.1f  %.1f  %.1f, %d', ...
%         min(im.nodeDiam(vSel)), ...
%         mean(im.nodeDiam(vSel)), ...
%         max(im.nodeDiam(vSel)), length(lstEdgesSel) );
%     set(handles.textSegmentInfo,'string',foos);
% else
%     set(handles.textSegmentInfo,'string','');
% end


% REMOVED 6/5/08
% % --- Executes on button press in pushbuttonPlotTree.
% function pushbuttonPlotTree_Callback(hObject, eventdata, handles)
% global im
% 
% if ~isempty(im.Stats.E1nodeIdx)
%     plotTreeSeg(im.nodeEdges,im.nodePos,im.Stats.E1nodeIdx,8,1,double(im.I));
% end



% --- Executes on button press in pushbuttonImageTracks.
function pushbuttonImageTracks_Callback(hObject, eventdata, handles)
global im

updateImageOverlay(0,handles);

hWait = waitbar( 0, 'Imaging tracks...' );

[ny,nx,nz] = size(im.III);
I = zeros(ny,nx,nz);

for ii=1:length(im.seed)
    waitbar(ii/length(im.seed),hWait);
    for jj=1:length(im.seed(ii).trkConn)
        nPts = im.seed(ii).trkConn(jj);
%        nPts = size(im.seed(ii).trkPos,1);
        if nPts>0
            for kk=1:nPts
                xx = min(max(round(im.seed(ii).trkPos(kk,1,jj)),1),nx);
                yy = min(max(round(im.seed(ii).trkPos(kk,2,jj)),1),ny);
                zz = min(max(round(im.seed(ii).trkPos(kk,3,jj)),1),nz);
                I(yy,xx,zz) = I(yy,xx,zz) + 1;
            end
        end
    end
end
close(hWait)

if get(handles.checkboxSaveTracks,'value')
    [filenm,filepath] = uiputfile('*.tif','Save Tracks as a TIFF file');
    wd = cd;
    cd(filepath);
    imwrite(I(:,:,1),filenm,'tiff','compression','none');
    for ii=2:nz
        imwrite(I(:,:,ii),filenm,'tiff','compression','none','writemode','append');
    end
    cd(wd);
end

lst=find(I>0);
im.III(lst) = uint8(min(I(lst)*3,32)+32);

updateAxes( handles )


set(handles.pushbuttonImageGraph,'enable','on');





% --- Executes on button press in pushbuttonRunPath.
function pushbuttonRunPath_Callback(hObject, eventdata, handles)
global im

if ~isfield(im,'RunPath')  
    foo=uigetdir('.','Set Path to directory for running seed tracking');
else
    foo=uigetdir(im.RunPath,'Set Path to directory for running seed tracking');    
end
if foo~=0
    im.RunPath = foo;
    set(handles.pushbuttonRunPath,'string',sprintf('Path: %s',im.RunPath));
    
    needs = 'This directory needs';
    flag = 0;
    if exist(sprintf('%s\\qqsub.pl',im.RunPath),'file')~=2 & exist(sprintf('%s/qqsub.pl',im.RunPath),'file')~=2
        needs=sprintf('%s qqsub.pl',needs);
        flag = 1;
    end
    if exist(sprintf('%s\\launcher.sh',im.RunPath),'file')~=2 & exist(sprintf('%s/launcher.sh',im.RunPath),'file')~=2
        needs=sprintf('%s launcher.sh',needs);
        flag = 1;
    end
    if exist(sprintf('%s\\log',im.RunPath),'file')~=7 & exist(sprintf('%s/log',im.RunPath),'file')~=7
        needs=sprintf('%s \\log',needs);
        flag = 1;
    end
    if exist(sprintf('%s\\DTItrack',im.RunPath),'file')~=2 & exist(sprintf('%s/DTItrack',im.RunPath),'file')~=2
        needs=sprintf('%s DTItrack',needs);
        flag = 1;
    end
    if flag==1
        ch = menu(needs,'Okay');
    end
    
    if exist(sprintf('%s\\%s.img',im.RunPath,im.filenm),'file')~=2 & exist(sprintf('%s/%s.img',im.RunPath,im.filenm),'file')~=2
        ch = menu('This directory needs the .img file for DTItrack.\nShall I create it?','Yes','No');
        if ch==1
            hWait = waitbar( 0,'Copying Image' );

            if strcmp(im.filetype,'img')
                filenm = sprintf( '%s.img',im.filenm );
                fid = fopen( filenm, 'rb' );
                n = fread( fid, 3, 'int');
                nx = n(1); ny = n(2); nz = n(3);
                I = fread( fid, nx*ny*nz, 'float' );
                I = reshape(I, [ny nx nz] );
                fclose(fid);
            elseif strcmp(im.filetype,'tif') | strcmp(im.filetype,'tiff')
                finfo = imfinfo(sprintf('%s.%s',im.filenm,im.filetype));
                ny = finfo(1).Height;
                nx = finfo(1).Width;
                nz = length(finfo);
                I = zeros(ny,nx,nz);
                for ii=1:nz
                    I(:,:,ii) = imread(sprintf('%s.%s',im.filenm,im.filetype),ii);
                end
                I = I / max(I(:));
            elseif strcmp(im.filetype,'mat')
                load(im.filenm,'I');

                [ny,nx,nz] = size(I);
                I = I / max(I(:));
            end



            wd = cd;
            cd(im.RunPath)
            filenm = sprintf('%s.img',im.filenm);
            fid = fopen(filenm,'wb');
            fwrite(fid,[nx ny nz],'int');
            fwrite(fid,I,'float');
            fclose(fid)
            cd(wd);
            close(hWait)
        end
    end
    
end





% --- Executes on button press in pushbuttonRunAll.
function pushbuttonRunAll_Callback(hObject, eventdata, handles)
global im

hWait = waitbar(0,'Running All');

pushbuttonRunSeeds_Callback(hObject, eventdata, handles);
wd = cd;
linuxWD = im.RunPath;
cd(linuxWD);
if ispc
    !del done
else
    !rm done
end

set(handles.pushbuttonRunSeeds,'Enable', 'off' );
set(handles.uipanelRunAll,'visible','on')

% process seeds until no more to process
while ~exist('done','file') %strcmpi(get(handles.pushbuttonProcSeeds,'enable'),'on')
    pause(10);
end
cd(wd);
pushbuttonProcSeeds_Callback(hObject, eventdata, handles);

if get(handles.checkboxRunAllSeedConn,'value')
    pushbuttonSeedConnectivity_Callback(hObject, eventdata, handles);
end

if get(handles.checkboxRunAllGraphTracks,'value')
    pushbuttonGraph_Callback(hObject, eventdata, handles);
end

if get(handles.checkboxRunAllImageTracks,'value')
    pushbuttonImageTracks_Callback(hObject, eventdata, handles);
end

set(handles.uipanelRunAll,'visible','off')
close(hWait)







function editViewSegmentN_Callback(hObject, eventdata, handles)
% hObject    handle to editViewSegmentN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editViewSegmentN as text
%        str2double(get(hObject,'String')) returns contents of editViewSegmentN as a double


% --- Executes during object creation, after setting all properties.
function editViewSegmentN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editViewSegmentN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in pushbuttonSegmentCenter.
function pushbuttonSegmentCenter_Callback(hObject, eventdata, handles)
global im

E = im.nodeEdges;
V = im.nodePos;
n0 = im.nodeSelected;

hWait = waitbar(0,'Centering segment...');

I = loadImage( im );


figure(12);
clf
nsp = length(im.segmentSelectedNodes);
nr = 1;
nc = 4;
if nsp>4
    nr = ceil(nsp/4);
    nc = 4;
end

for idxNode = 1:length(im.segmentSelectedNodes)
    n0 = im.segmentSelectedNodes(idxNode);
    
    lst = find( E(:,1)==n0 | E(:,2)==n0 );
    if length(lst)==2
        n1 = setdiff(E(lst(1),:),n0);
        n2 = setdiff(E(lst(2),:),n0);

        x = repmat([-10:10],[21 1]);
        y = repmat([-10:10]',[1 21]);
        z = zeros(21,21);
        xx = x(:);
        yy = y(:);
        zz = z(:);

        r1 = V(n1,:);
        r2 = V(n2,:);
        dr = r1-r2;
        dx = r1(1) - r2(1);
        dz = r1(3) - r2(3);
        r = norm(dr);
        rho = norm(dr(1:2));
        beta = acos(dz/r);
        gamma = -acos(dx/rho);
        Ry = [cos(beta) 0 -sin(beta); 0 1 0; sin(beta) 0 cos(beta)];
        Rz = [cos(gamma) sin(gamma) 0; -sin(gamma) cos(gamma) 0; 0 0 1];
        foo = Ry*[xx';yy';zz'];
        foo = Rz*foo;
        xx = foo(1,:)';
        yy = foo(2,:)';
        zz = foo(3,:)';
        [ny nx nz] = size(I);
        %    lst = round((yy+V(pt,2)) + (xx+V(pt,1))*ny + (zz+V(pt,3))*nx*ny);
        xo = V(n0,1);
        yo = V(n0,2);
        zo = V(n0,3);
        for ii=1:size(x(:))
            A(ii) = I(min(max(round(yy(ii)+yo),1),ny),min(max(round((xx(ii)+xo)),1),nx),min(max(round((zz(ii)+zo)),1),nz));
        end

        [foo,xyi] = max(A(:));

        im.nodePos(n0,:) = round([xo+xx(xyi) yo+yy(xyi) zo+zz(xyi)]);

        subplot(nr,nc,idxNode)
        imagesc(reshape(A,[21 21]))
        axis image
        hold on
        plot(11,11,'r*')
        xi = floor(xyi/size(x,1))+1;
        yi = rem(xyi,size(x,1));
        plot(xi,yi,'g*')
        hold off
        colormap(gray)
        title( sprintf('b=%.2f   g=%.2f',beta*180/pi,gamma*180/pi) )
    end
end
close(hWait)

set(handles.pushbuttonImageGraph,'enable','on')





% --- Executes on button press in pushbuttonNodeUndo.
function pushbuttonNodeUndo_Callback(hObject, eventdata, handles)
global im
    
nodeFlag = ones(size(im.nodePos,1),1);
nodeFlag(im.nodeSelected) = 0;

[im.nodePos,im.nodeDiam,im.nodeDiamThetaIdx,im.nodeBC,im.nodeBCType,im.nodeType,im.nodeSegN,im.nodeEdges,im.edgeFlag] = removeNodes( nodeFlag, im.nodePos, im.nodeDiam, im.nodeDiamThetaIdx, im.nodeBC, im.nodeBCType, im.nodeType, im.nodeSegN,im.nodeEdges );
set(handles.pushbuttonUpdateVesselMask,'enable','on')
im.nBflag = 1;

set(handles.pushbuttonNodeUndo,'visible','off');

pushbuttonImageGraph_Callback(hObject, eventdata, handles);




function radiobuttonImageView_Callback(hObject, eventdata, handles)
global im

lst = [3 2 1];
im.ImageView = eventdata;
if isfield(im,'nodeSelected')
    targetS = im.nodePos(im.nodeSelected,lst(im.ImageView));
else
    targetS = 1;
end

axes( handles.axes1 );

iS = get(handles.slider1,'value');
nS = str2num(get(handles.editNslices,'string'));
if eventdata==1
    xlim([1 size(im.III,2)]);
    ylim([1 size(im.III,1)]);
    set(handles.slider1,'max',im.nZ-nS)
    set(handles.slider1,'value',min(im.nZ-nS,targetS))
    set(handles.slider1,'sliderstep',[1/im.nZ .1])
elseif eventdata ==2
    xlim([1 size(im.III,2)]);
    ylim([1 size(im.III,3)]);
    set(handles.slider1,'max',im.nY-nS)
    set(handles.slider1,'value',min(im.nY-nS,targetS))
    set(handles.slider1,'sliderstep',[1/im.nY .1])
else
    xlim([1 size(im.III,1)]);
    ylim([1 size(im.III,3)]);
    set(handles.slider1,'max',im.nX-nS)
    set(handles.slider1,'value',min(im.nX-nS,targetS))
    set(handles.slider1,'sliderstep',[1/im.nX .1])
end
iS = get(handles.slider1,'value');
set(handles.textSlices,'string',sprintf('Slice %d to %d',round(iS),round(iS+nS-1)))

updateAxes( handles );


% --- Executes on button press in radiobuttonSegmentArterial.
function radiobuttonSegmentACV_Callback(hObject, eventdata, handles)
global im

im.nodeType(im.segmentSelectedNodes) = eventdata*get(hObject,'value');
if ~isfield(im,'nodeTypeUpdated')
    im.nodeTypeUpdated = zeros(length(im.nodeType),1);
end
im.nodeTypeUpdated(im.segmentSelectedNodes) = 1;

if eventdata==1
    set(handles.radiobuttonSegmentCapillary,'value',0);
    set(handles.radiobuttonSegmentVein,'value',0);
elseif eventdata==2
    set(handles.radiobuttonSegmentArterial,'value',0);
    set(handles.radiobuttonSegmentVein,'value',0);
elseif eventdata==3
    set(handles.radiobuttonSegmentArterial,'value',0);
    set(handles.radiobuttonSegmentCapillary,'value',0);
end
set(handles.pushbuttonUpdateVesselMask,'enable','on')
set(handles.pushbuttonCalcFlowCircuitEq,'enable','on')



% --- Executes on button press in pushbuttonLabelCapillary.
function pushbuttonLabelCapillary_Callback(hObject, eventdata, handles)
global im

if ~isfield(im,'nodeDiam')
    menu( 'Node Diameters have not been calculated!','Okay');
    return;
end
if ~isfield(im,'segDiam')
    menu('You need to update Group Stats. Press the [Groups] button.','Okay');
    return;
end

diamThresh = str2num(get(handles.editCapillaryMaxDiam,'string'));

hwait = waitbar(0,'Labeling capillaries');
for iSeg = 1:max(im.nodeSegN)
    waitbar(iSeg/max(im.nodeSegN),hwait);
    lst = find(im.nodeSegN==iSeg);
    if median(im.nodeDiam(lst)) <= diamThresh
        im.nodeType(lst) = 2;
    end
end
close(hwait)

% for iNode = 1:size(im.nodePos,1)
%     if im.nodeDiam(iNode) <= diamThresh
%         im.nodeType(iNode) = 2;
%     end
% end

%imView3d_vesselMask( ); % why was this here?

%set(handles.pushbuttonUpdateVesselOverlay,'enable','on');
set(handles.pushbuttonUpdateVesselMask,'enable','on')


% --- Executes on button press in pushbuttonUpdateVesselMask.
function pushbuttonUpdateVesselMask_Callback(hObject, eventdata, handles)
global im

imView3d_vesselMask( );

set(handles.pushbuttonUpdateVesselOverlay,'enable','on');
set(handles.pushbuttonUpdateVesselMask,'enable','off')

    
    


% --- Executes on button press in radiobuttonBCPres.
function radiobuttonBCType_Callback(hObject, eventdata, handles)
global im

if eventdata==1
    set(handles.radiobuttonBCVel,'value',0);
elseif eventdata==2
    set(handles.radiobuttonBCPres,'value',0);
end

litval = get(handles.checkboxBClitVal,'value');
im.nodeBCType(im.Stats.E1nodeIdx) = (eventdata + litval*2);% * get(hObject,'value');
set(handles.pushbuttonCalcFlowCircuitEq,'enable','on')


% --- Executes on button press in checkboxBClitVal.
function checkboxBClitVal_Callback(hObject, eventdata, handles)
global im

BCtype =  im.nodeBCType(im.Stats.E1nodeIdx);
if BCtype > 0
    if get(hObject,'value')
        im.nodeBCType(im.Stats.E1nodeIdx) = mod( BCtype-1,2)+3;
    else
        im.nodeBCType(im.Stats.E1nodeIdx) = mod( BCtype-1,2)+1;
    end
end
set(handles.pushbuttonCalcFlowCircuitEq,'enable','on')



function Exit_CallbackFcn(hObject, eventdata, handles)
global im

stop(im.timer)
delete(im.timer)

    


% --- Executes on button press in pushbuttonUpdateVesselOverlay.
function pushbuttonUpdateVesselOverlay_Callback(hObject, eventdata, handles)
pushbuttonImageGraph_Callback(hObject, eventdata, handles)
set(handles.pushbuttonUpdateVesselOverlay,'enable','off');


% --- Executes on button press in pushbuttonCalcFlowCircuitEq.
function pushbuttonCalcFlowCircuitEq_Callback(hObject, eventdata, handles)
imView3d_flowCircuitEq()
set(handles.pushbuttonUpdateVesselMask,'enable','on')
set(handles.pushbuttonCalcFlowCircuitEq,'enable','off')


% --- Executes on button press in radiobuttonOverlayVel.
function radiobuttonOverlay_Callback(hObject, eventdata, handles)
global im

if eventdata==1
    set(handles.radiobuttonOverlayVel,'value',0);
    set(handles.radiobuttonOverlayFlow,'value',0);
    set(handles.radiobuttonOverlayPressure,'value',0);
elseif eventdata==2
    set(handles.radiobuttonOverlayVesselType,'value',0);
    set(handles.radiobuttonOverlayFlow,'value',0);
    set(handles.radiobuttonOverlayPressure,'value',0);
    set(handles.textOverlayInfo,'string', sprintf('Max Vel = %.1f mm/s',1e-3*max(abs(im.edgeVel))) );
    figure(11)
    imagesc(1e-3*[32:-1:1]'*max(abs(im.edgeVel))/32,1e-3*[1 32]*max(abs(im.edgeVel))/32)
    colormap(jet(32))
    colorbar
elseif eventdata==3
    set(handles.radiobuttonOverlayVesselType,'value',0);
    set(handles.radiobuttonOverlayVel,'value',0);
    set(handles.radiobuttonOverlayPressure,'value',0);
    set(handles.textOverlayInfo,'string', sprintf('Max Flow = %.2e\num^3/s',max(abs(im.edgeFlow))) );
    figure(11)
    imagesc([32:-1:1]'*max(abs(im.edgeFlow))/32,[1 32]*max(abs(im.edgeFlow))/32)
    colormap(jet(32))
    colorbar
elseif eventdata==4
    set(handles.radiobuttonOverlayVesselType,'value',0);
    set(handles.radiobuttonOverlayVel,'value',0);
    set(handles.radiobuttonOverlayFlow,'value',0);
    set(handles.textOverlayInfo,'string', sprintf('Max Pressure = %.0f',max(abs(im.nodePressure))) );
    figure(11)
    imagesc([32:-1:1]'*max(abs(im.nodePressure))/32,[1 32]*max(abs(im.nodePressure))/32)
    colormap(jet(32))
    colorbar
end
set(handles.pushbuttonUpdateVesselOverlay,'enable','on');


% --- Executes on button press in pushbuttonLabelVessels.
function pushbuttonLabelVessels_Callback(hObject, eventdata, handles)
global im

hwait = waitbar(0,'Labeling segments');
[nB,im] = nBupdate( im );
%nB = zeros(size(im.nodePos,1),1);
nN = size(im.nodePos,1);
%for ii=1:nN
%    nB(ii) = length(find(im.nodeEdges(:,1)==ii | im.nodeEdges(:,2)==ii));
%end

nodeType = im.nodeType(im.nodeSelected);
if nodeType~=1 & nodeType~=3
    menu('Selected node must be an arteriole or vein','okay');
    close(hwait);
    return;
end

diamThresh = str2num(get(handles.editCapillaryMaxDiam,'string'));

iN = im.nodeSelected;
iSeg = im.nodeSegN(iN);

lst = find(im.nodeSegN==iSeg);
if median(im.nodeDiam(lst))<=diamThresh
    menu('Selected segment meets the capillary criterion, so nothing done.','Okay');
    close(hwait);
    return;
end
for ii=1:length(lst)
    iEtmp = find(im.nodeEdges(:,1)==lst(ii) | im.nodeEdges(:,2)==lst(ii));
    foo = im.nodeEdges(iEtmp,:);
    lst = [lst; foo(:)];
end
lst = unique(lst);
im.nodeType(lst) = nodeType;
iNB = lst(find(nB(lst)>2));
eLst = [];
for ii=1:length(iNB)
    eLst = [eLst; find(im.nodeEdges(:,1)==iNB(ii) | im.nodeEdges(:,2)==iNB(ii))];
end
lstSeg = setdiff(unique(im.edgeSegN(eLst)),iSeg);
oldSeg = iSeg;


while ~isempty(lstSeg)
    iSeg = lstSeg(end);
    lstSeg = lstSeg(1:end-1);
    oldSeg = [oldSeg; iSeg];

    lst = find(im.nodeSegN==iSeg);
    if median(im.nodeDiam(lst))>diamThresh
        for ii=1:length(lst)
            iEtmp = find(im.nodeEdges(:,1)==lst(ii) | im.nodeEdges(:,2)==lst(ii));
            foo = im.nodeEdges(iEtmp,:);
            lst = [lst; foo(:)];
        end
        lst = unique(lst);
        im.nodeType(lst) = nodeType;
        if ~isfield(im,'nodeTypeUpdated')
            im.nodeTypeUpdated = zeros(length(im.nodeType),1);
        end
        im.nodeTypeUpdated(lst) = 1;
        
        iNB = lst(find(nB(lst)>2));
        eLst = [];
        for ii=1:length(iNB)
            eLst = [eLst; find(im.nodeEdges(:,1)==iNB(ii) | im.nodeEdges(:,2)==iNB(ii))];
        end
        lstSeg = unique([lstSeg; setdiff(unique(im.edgeSegN(eLst)),oldSeg)]);
    end
end
close(hwait)

set(handles.pushbuttonUpdateVesselMask,'enable','on');


% --------------------------------------------------------------------
function menu_ViewSeedPanels_Callback(hObject, eventdata, handles)

if strcmp(get(hObject,'checked'),'off')
    set(handles.uipanel_AutoSeed,'visible','on')
    set(handles.uipanelSeed,'visible','on')
    set(handles.uipanel_SeedConnectivity,'visible','on')
    set(handles.radiobuttonSeedPoint,'visible','on')
    set(handles.uipanelRunAll,'visible','on')
    set(hObject,'checked','on')
else
    set(handles.uipanel_AutoSeed,'visible','off')
    set(handles.uipanelSeed,'visible','off')
    set(handles.uipanel_SeedConnectivity,'visible','off')
    set(handles.radiobuttonSeedPoint,'visible','off')
    set(handles.uipanelRunAll,'visible','off')
    set(hObject,'checked','off')
end



% --------------------------------------------------------------------
function menu_ListCommand_Callback(hObject, eventdata, handles)
disp(sprintf('\n======================================\nCommands' ))
disp(sprintf('======================================' ))
disp(sprintf('  a - scroll up'))
disp(sprintf('  z - scroll down'))
disp(sprintf('\nSELECT SEGMENT'))
disp(sprintf('           left mouse - select 1 edge'))
disp(sprintf('  <shift>  left mouse - select entire segment from branch pt to branch pt'))
disp(sprintf('\nADD NODES'))
disp(sprintf('           left mouse - add node with single edge to nearest node'))
disp(sprintf('          right mouse - add node with 2 edges to 2 nearest nodes'))
disp(sprintf('  <shift> right mouse - add node without an edge'))


% --- Executes on button press in pushbuttonEstDiam.
function pushbuttonEstDiam_Callback(hObject, eventdata, handles)

reEstimated = imView3d_GetDiam();
if reEstimated
    set(handles.pushbuttonCalcFlowCircuitEq,'enable','on')
%    set(handles.pushbuttonPruneDiam,'enable','on')
%    set(handles.pushbuttonRegraphNodes,'enable','on')
%    set(handles.pushbuttonImageGraph,'enable','on')
end

% --- Executes on button press in pushbuttonCenterNodes.
function pushbuttonCenterNodes_Callback(hObject, eventdata, handles)
global im

if strcmpi(get(handles.menu_VisualizeCentering,'checked'),'off')
    flagVisualize = 0;
else
    flagVisualize = 1;
end

%if eventdata==1 % get step size
    centerStep = get(handles.checkboxCenterStep,'value');
    %str2num(get(handles.editCenterStep,'string')) / 100;
%else
%    centerStep = 1;
%end

thresh = ceil(str2num(get(handles.editImageThresh,'string'))*32);

reEstimated = imView3d_CenterNodes( eventdata, centerStep, flagVisualize, 1, thresh );
if reEstimated==1
    set(handles.pushbuttonRegraphNodes,'enable','on')
    set(handles.pushbuttonImageGraph,'enable','on')
elseif reEstimated==2
    set(handles.pushbuttonImageGraph,'enable','on')
    set(handles.textNumEdges,'string',sprintf('%d edges',size(im.nodeEdges,1)))
    set(handles.uipanelGraph,'title',sprintf('Graph (%d nodes)',size(im.nodePos,1)));
end


% --- Executes on button press in pushbuttonRemoveSeeds.
function pushbuttonRemoveSeeds_Callback(hObject, eventdata, handles)
global im

ch = menu('Do you want to remove unlaunched seeds?','Yes','No');

% find seeds not launched and remove
if ch==1
    lst = [];
    for ii=1:im.nSeed
        if im.seed(ii).seedLaunched
            lst(end+1) = ii;
        end
    end
    im.seed = im.seed(lst);
    im.nSeed = length(lst);
    if length(lst)>0
        im.curSeed = 1;
    else
        im.curSeed = 0;
    end
    
    UpdateSeeds( handles )
    set(handles.pushbuttonRemoveSeeds,'visible','off')
end


% --- Executes on button press in pushbuttonFillNodes.
function pushbuttonFillNodes_Callback(hObject, eventdata, handles)
global im

nE = size(im.nodeEdges,1);
hxy = str2num(get(handles.editGraphHxy,'string')) * mean(im.Hvox(1:2));
hz = str2num(get(handles.editGraphHz,'string')) * im.Hvox(3);
h = min(hxy,hz);

nN = size(im.nodePos,1);
nE = size(im.nodeEdges,1);
flag = 0;
hwait = waitbar(0,'Filling in sparse edges');
nEo = nE;
for iE = 1:nEo
    waitbar(iE/nEo,hwait)
    n1 = im.nodeEdges(iE,1);
    n2 = im.nodeEdges(iE,2);
    len = abs( im.nodePos(n1,:) - im.nodePos(n2,:) );
    if sum(len>2*[hxy hxy hz])
        flag = 1;
        nStep = max(floor(len ./ [hxy hxy hz]) - 1);
        rStep = (im.nodePos(n2,:) - im.nodePos(n1,:)) / (nStep+1);
        pos = im.nodePos(n1,:);
        for jj=1:nStep
            pos = pos + rStep;
            im.nodePos(end+1,:) = pos;
            nN = nN + 1;
            im.nodeEdges(end+1,:) = [n1 nN];
            nE = nE + 1;
            n1 = nN;

            im.nodeDiam(nN) = 0;
            im.nodeVel(nN) = 0;
            im.nodeBC(nN) = 0;
            im.nodeBCType(nN) = 0;
            im.nodeType(nN) = 0;
            im.nodeSegN(nN) = 0;
            im.edgeFlag(nE) = 0;
        end
        im.nodeEdges(iE,:) = [nN n2];
    end
end
close(hwait)

if flag
    set(handles.pushbuttonRegraphNodes,'enable','on')
    set(handles.pushbuttonImageGraph,'enable','on')
    set(handles.uipanelGraph,'title',sprintf('Graph (%d nodes)',size(im.nodePos,1)));
    set(handles.textNumEdges,'string',sprintf('%d edges',size(im.nodeEdges,1)))
end


% --------------------------------------------------------------------
function menu_VisualizeCentering_Callback(hObject, eventdata, handles)

if strcmpi(get(hObject,'checked'),'off')
    set(hObject,'checked','on')
else
    set(hObject,'checked','off')
end







% --- Executes on button press in pushbuttonCenterIterate.
function pushbuttonCenterIterate_Callback(hObject, eventdata, handles)
global im

if strcmpi(get(handles.menu_VisualizeCentering,'checked'),'off')
    flagVisualize = 0;
else
    flagVisualize = 1;
end

%if eventdata==1 % get step size
    centerStep = get(handles.checkboxCenterStep,'value');
%else
%    centerStep = 1;
%end

thresh = ceil(str2num(get(handles.editImageThresh,'string'))*32);



ch = menu('Iterate:','S5RI','[S5FR]I','[CZSFR]I','[(CZS)5FR]I on 1-nodes','(CSR)5I','S5R Z S5R E S5 RI','(CZSFR)5I','(CERFR)5I','Cancel');
drawnow


reEstimated = 0;
switch ch
    case 1 % straighten 5 times
        hwait = waitbar(0,'Straightening 5 times then Regraph and Image');
        for ii=1:5
            waitbar(ii/5,hwait)
            foo = imView3d_CenterNodes( 2, centerStep, flagVisualize, 1, thresh );
            reEstimated = reEstimated | foo;
        end
        if reEstimated
            pushbuttonRegraphNodes_Callback(handles.pushbuttonRegraphNodes, [], handles);
            pushbuttonImageGraph_Callback(handles.pushbuttonImageGraph, [], handles)
            reEstimated = 0;
        end
        close(hwait)
        
    case 2 % [S5RFR]I
        nNodes0 = 9e9;
        iter = 0;
        while abs(nNodes0-size(im.nodePos,1))>10
            iter = iter + 1;
            hwait = waitbar(0,sprintf('[S5FR]I iteration %d',iter));
            nNodes0 = size(im.nodePos,1);
            for ii=1:5
                waitbar(ii/5,hwait)
                imView3d_CenterNodes( 2, centerStep, flagVisualize, 1, thresh );
            end
%            pushbuttonRegraphNodes_Callback(handles.pushbuttonRegraphNodes, [], handles);
            pushbuttonFillNodes_Callback(hObject, eventdata, handles);
            pushbuttonRegraphNodes_Callback(handles.pushbuttonRegraphNodes, [], handles);
            close(hwait)
        end
        pushbuttonImageGraph_Callback(handles.pushbuttonImageGraph, [], handles)
        reEstimated = 0;
        
    case 3 % [CZSFR]I
        nNodes0 = 9e9;
        iter = 0;
        while abs(nNodes0-size(im.nodePos,1))>5
            iter = iter + 1;
            hwait = waitbar(0,sprintf('[CZSFR]I iteration %d',iter));
            nNodes0 = size(im.nodePos,1);
            flagZero = 1;
            imView3d_CenterNodes( 1, centerStep, flagVisualize, flagZero, thresh );
            imView3d_CenterNodesXYZ( centerStep, flagVisualize, thresh );
            imView3d_CenterNodes( 2, centerStep, flagVisualize, flagZero, thresh );
            pushbuttonFillNodes_Callback(hObject, eventdata, handles);
            pushbuttonRegraphNodes_Callback(handles.pushbuttonRegraphNodes, [], handles);
            close(hwait)
        end
        pushbuttonImageGraph_Callback(handles.pushbuttonImageGraph, [], handles)
        reEstimated = 0;

    case 4 %  [(CZS)5FR]I    [CFR]I
        nNodes0 = 9e9;
        nDangling0 = 9e9;
        nDangling = 0;
        n1Nodes0 = 9e9;
        n1Nodes = 0;
        iter = 0;
        
        while abs(nNodes0-size(im.nodePos,1))>5 | abs(n1Nodes0-n1Nodes)>0
            if iter>0
                disp( sprintf('%d) #1-nodes = %d,   #dangling = %d     #1-nodes to branch point = %d', iter, length(find(nB==1)), nDangling, length(nLst)) )
                hwait = waitbar(0,sprintf('[(CZS)5FR]I iteration %d on %d 1-nodes',iter,length(nLst)));
                nNodes0 = size(im.nodePos,1);
                flagZero = 1;
                for ii = 1:5
                    waitbar(ii/5,hwait);
                    imView3d_CenterNodes( 1, centerStep, flagVisualize, flagZero, thresh, nLst );
                    imView3d_CenterNodesXYZ( centerStep, flagVisualize, thresh, nLst );
                    imView3d_CenterNodes( 2, centerStep, flagVisualize, flagZero, thresh, nLst );
                end
                pushbuttonFillNodes_Callback(hObject, eventdata, handles);
                pushbuttonRegraphNodes_Callback(handles.pushbuttonRegraphNodes, [], handles);
                close(hwait)
                drawnow
            end
            
            [nB,im] = nBupdate( im );
%            nB = zeros(size(im.nodePos,1),1);
%            nN = size(im.nodePos,1);
%            for ii=1:nN
%                nB(ii) = length(find(im.nodeEdges(:,1)==ii | im.nodeEdges(:,2)==ii));
%            end
            
            nLst = find(nB==1);
            nDangling0 = nDangling;
            nDangling = 0;
            n1Nodes0 = n1Nodes;
            n1Nodes = length(nLst);
            for ii=1:length(nLst)
                n1 = nLst(ii);
                lst = find(im.nodeEdges(:,1)==n1 | im.nodeEdges(:,2)==n1);
                n2 = setdiff(im.nodeEdges(lst,:),n1);
                nLstSub = [];
                nSub = 0;
                while nB(n2)==2 & nSub<10
                    nLstSub(end+1) = n2;
                    n1(end+1) = n2;
                    nSub = nSub + 1;
                    lst = find(im.nodeEdges(:,1)==n1(end) | im.nodeEdges(:,2)==n1(end));
                    foo = im.nodeEdges(lst,:);
                    n2 = setdiff(foo(:),n1);
                end
                if nSub<10
                    nLst(end+[1:nSub]) = nLstSub;
                    nDangling = nDangling + 1;
                end
            end
        end
        disp( sprintf('%d) #1-nodes = %d,   #dangling = %d     #1-nodes to branch point = %d', iter, length(find(nB==1)), nDangling, length(nLst)) )
        pushbuttonImageGraph_Callback(handles.pushbuttonImageGraph, [], handles)
        reEstimated = 0;
        
    case 5 % center and straighten 5 times
        hwait = waitbar(0,'Center, Straighten, and Regraph 5 times then Image');
        for ii=1:5
            waitbar(ii/5,hwait)
            flagZero = 1;
            foo = imView3d_CenterNodes( 1, centerStep, flagVisualize, flagZero, thresh );
            reEstimated = reEstimated | foo;
            foo = imView3d_CenterNodes( 2, centerStep, flagVisualize, 1, thresh );
            reEstimated = reEstimated | foo;
            pushbuttonRegraphNodes_Callback(handles.pushbuttonRegraphNodes, [], handles);
        end
        if reEstimated
            pushbuttonImageGraph_Callback(handles.pushbuttonImageGraph, [], handles)
            reEstimated = 0;
        end
        close(hwait)
    case 6
        flagZero = 1;
        
        hwait = waitbar(0,'S5R Z S5R E S5 then Regraph and Image');
        for ii=1:5
            waitbar(ii/17,hwait)
            imView3d_CenterNodes( 2, centerStep, flagVisualize, 1, thresh );
        end
        pushbuttonRegraphNodes_Callback(handles.pushbuttonRegraphNodes, [], handles);
%        pushbuttonImageGraph_Callback(handles.pushbuttonImageGraph, [], handles)
        
        waitbar(6/17,hwait)
        imView3d_CenterNodesXYZ( centerStep, flagVisualize, thresh );
        
        for ii=1:5
            waitbar((6+ii)/17,hwait)
            imView3d_CenterNodes( 2, centerStep, flagVisualize, 1, thresh );
        end
        pushbuttonRegraphNodes_Callback(handles.pushbuttonRegraphNodes, [], handles);
%        pushbuttonImageGraph_Callback(handles.pushbuttonImageGraph, [], handles)
        
        waitbar(12/17,hwait)
        imView3d_erodeCenterLine( thresh );
        
        for ii=1:5
            waitbar((12+ii)/17,hwait)
            imView3d_CenterNodes( 2, centerStep, flagVisualize, 1, thresh );
        end
        pushbuttonRegraphNodes_Callback(handles.pushbuttonRegraphNodes, [], handles);
        pushbuttonImageGraph_Callback(handles.pushbuttonImageGraph, [], handles)
        
        close(hwait)
        reEstimated = 0;
        
    case 7 % (CZSFR)5I
        hwait = waitbar(0,'(CZSFR)5I');
        for ii=1:5
            waitbar(ii/5,hwait)
            flagZero = 1;
            imView3d_CenterNodes( 1, centerStep, flagVisualize, flagZero, thresh );
            imView3d_CenterNodesXYZ( centerStep, flagVisualize, thresh );
            imView3d_CenterNodes( 2, centerStep, flagVisualize, flagZero, thresh );
            pushbuttonFillNodes_Callback(hObject, eventdata, handles);
            pushbuttonRegraphNodes_Callback(handles.pushbuttonRegraphNodes, [], handles);
        end
        pushbuttonImageGraph_Callback(handles.pushbuttonImageGraph, [], handles)
        close(hwait)

    case 8
        hwait = waitbar(0,'(CERFR)5I');
        for ii = 1:5
            waitbar(ii/5,hwait)
            flagZero = 1;
            imView3d_CenterNodes( 1, centerStep, flagVisualize, flagZero, thresh );
            imView3d_erodeCenterLine( thresh );
            pushbuttonRegraphNodes_Callback(handles.pushbuttonRegraphNodes, [], handles);
            pushbuttonFillNodes_Callback(hObject, eventdata, handles);
            pushbuttonRegraphNodes_Callback(handles.pushbuttonRegraphNodes, [], handles);
        end
        pushbuttonImageGraph_Callback(handles.pushbuttonImageGraph, [], handles)
        close(hwait)
        
end


set(handles.textNumEdges,'string',sprintf('%d edges',size(im.nodeEdges,1)))
set(handles.uipanelGraph,'title',sprintf('Graph (%d nodes)',size(im.nodePos,1)));

if reEstimated==1
    set(handles.pushbuttonRegraphNodes,'enable','on')
    set(handles.pushbuttonImageGraph,'enable','on')
elseif reEstimated==2
    set(handles.pushbuttonImageGraph,'enable','on')
end


% --------------------------------------------------------------------
function menu_ViewSeeds_Callback(hObject, eventdata, handles)

if strcmp(get(hObject,'checked'),'off')
    set(hObject,'checked','on')
    updateAxes( handles )
else
    set(hObject,'checked','off')
    updateAxes( handles )
end


% --- Executes on button press in pushbuttonErodeCenterLine.
function pushbuttonErodeCenterLine_Callback(hObject, eventdata, handles)

thresh = ceil(str2num(get(handles.editImageThresh,'string'))*32);

imView3d_erodeCenterLine( thresh );
set(handles.pushbuttonRegraphNodes,'enable','on')
set(handles.pushbuttonImageGraph,'enable','on')


% --------------------------------------------------------------------
function menu_viewSegments_Callback(hObject, eventdata, handles)
global im

menu_viewDanglingSegmentsSub(hObject, 2, handles );
return

if ~isfield(im,'segNedges')
    menu('You must get segment and group information before using this tool. Push the [groups] button.','Okay');
    return;
end

lenSeg = menu('View segments of length (# edges)',...
    sprintf('1 (# = %d)',length(find(im.segNedges==1))),...
    sprintf('2 (# = %d)',length(find(im.segNedges==2))),...
    sprintf('3 (# = %d)',length(find(im.segNedges==3))),...
    sprintf('4 (# = %d)',length(find(im.segNedges==4))),...
    sprintf('5 (# = %d)',length(find(im.segNedges==5))),...
    'Cancel');
if lenSeg==6
    return;
end
lst = find(im.segNedges==lenSeg);


[nB,im] = nBupdate( im );
%nB = zeros(size(im.nodePos,1),1);
%nN = size(im.nodePos,1);
%for ii=1:nN
%    nB(ii) = length(find(im.nodeEdges(:,1)==ii | im.nodeEdges(:,2)==ii));
%end

% info update axes view
if im.ImageView==1
    pIdx = [1 2 3];
    xlim([1 size(im.III,2)]);
    ylim([1 size(im.III,1)]);
    nIs = size(im.III,3);
elseif im.ImageView==2
    pIdx = [1 3 2];
    xlim([1 size(im.III,2)]);
    ylim([1 size(im.III,3)]);
    nIs = size(im.III,2);
else
    pIdx = [2 3 1];
    xlim([1 size(im.III,1)]);
    ylim([1 size(im.III,3)]);
    nIs = size(im.III,1);
end
nS = str2num(get(handles.editNslices,'string'));




edgeSegFlag = ones(size(im.edgeSegN));
nodeFlag = ones(size(im.nodePos));

for ii=1:length(lst)
    eLst = find(im.edgeSegN==lst(ii));
    eIdx = eLst(1);
    nIdx = im.nodeEdges(eIdx,1);

    % update tree seg view
    lst2 = find(im.edgeSegN==lst(ii));
    edgeSegFlag(lst2) = 0;
    plotTreeSeg( im.nodeEdges, im.nodePos, nIdx, 10, 1, im.I, [], [], [], edgeSegFlag );
    pos = round(im.nodePos(nIdx,:));
    title( sprintf('Segment %d of %d     pos = (%d,%d,%d)', ii, length(lst),pos(1),pos(2),pos(3)) )
    rotate3d on
    
    ch = 2;
    while ch==2
        ch = menu('Option?','Okay','View 2D image','Delete','Cancel');
        if ch==2
            % update axes view
            zz = min(max(im.nodePos(nIdx,pIdx(3))-nS/2,1),nIs-nS);
            zz = min(zz,get(handles.slider1,'Max'));
            set(handles.slider1,'value',zz);
            set(handles.textSlices,'string',sprintf('Slice %d to %d',round(zz),round(zz+nS-1)))

            updateAxes(handles);
            hold on
            h=plot(im.nodePos(nIdx,pIdx(1)),im.nodePos(nIdx,pIdx(2)),'cp');
            set(h,'markersize',16);
            set(h,'markerfacecolor','c')
            for jj=1:6
                set(h,'marker','none')
                pause(0.1)
                set(h,'marker','p')
                pause(0.1)
            end
            hold off
        end
    end
    edgeSegFlag(lst2) = 1;
    if ch==3
        edgeSegFlag(lst2) = 2;
        nLst = unique(im.nodeEdges(lst2,:));
        nLst = nLst(find(nB(nLst)<=2));
        nodeFlag(nLst) = 0;
    end
    if ch==4
        break
    end
end

% remove edges
hwait = waitbar(0,'Removing edges and nodes');
lst = find(edgeSegFlag==2);
im.nodeEdges(lst,:) = [];

% remove nodes
[im.nodePos,im.nodeDiam,im.nodeDiamThetaIdx,im.nodeBC,im.nodeBCType,im.nodeType,im.nodeSegN,im.nodeEdges,im.edgeFlag] = removeNodes( nodeFlag, im.nodePos, im.nodeDiam, im.nodeDiamThetaIdx, im.nodeBC, im.nodeBCType, im.nodeType, im.nodeSegN,im.nodeEdges );
set(handles.pushbuttonUpdateVesselMask,'enable','on')
im.nBflag = 1;

set(handles.textNumEdges,'string',sprintf('%d edges',size(im.nodeEdges,1)))
set(handles.uipanelGraph,'title',sprintf('Graph (%d nodes)',size(im.nodePos,1)));
close(hwait)

pushbuttonGroupStats_Callback(handles.pushbuttonGroupStats, eventdata, handles)

set(handles.pushbuttonImageGraph,'enable','on')


% --------------------------------------------------------------------
function menu_pruneLoops_Callback(hObject, eventdata, handles)
global im

[nB,im] = nBupdate( im );
%nB = zeros(size(im.nodePos,1),1);
%nN = size(im.nodePos,1);
%for ii=1:nN
%    nB(ii) = length(find(im.nodeEdges(:,1)==ii | im.nodeEdges(:,2)==ii));
%end

% I am looping through all branching nodes
% identify the loops
% delete 1 edge if it meets the condition
lstN = 1:nN;
lstNbr = find(nB>2);
nNbr = length(lstNbr);
nEdel = 0;
hwait = waitbar(0,'Pruning loops');
n1 = 0; % keep track of number of nodes for which no loops removed
n2 = 0; % keep track of number of loops not pruned for skipped nodes
while ~isempty(lstNbr)
    iN = lstNbr(end);
    lstNbr = lstNbr(1:end-1);
    waitbar((nNbr-length(lstNbr))/nNbr,hwait);

    % GET SUB GRAPH
    % find nodes within a certain distance of node iN
    pos = im.nodePos(iN,:);
    r = sum( [(im.nodePos(:,1)-ones(nN,1)*pos(1)) ...
        (im.nodePos(:,2)-ones(nN,1)*pos(2))...
        (im.nodePos(:,3)-ones(nN,1)*pos(3)) ].^2, 2 ).^0.5;
    lstNsub = find(r<20);
    
    % identify edges
    lstE = [];
    for jj=1:length(lstNsub)
        [ir,ic] = find(im.nodeEdges==lstNsub(jj));
        lstE = [lstE;ir];
    end
    lstE = unique(lstE);
    
    % sub graph
    e = im.nodeEdges(lstE,:);
    for jj=1:length(lstNsub)
        e(find(e==lstNsub(jj))) = jj;
    end
    % find dangling nodes
    lst = find(e>length(lstNsub));
    for jj = 1:length(lst)
        kk=find(e(lst(jj))==lstNsub);
        if isempty(kk)
            lstNsub(end+1) = e(lst(jj));
            e(lst(jj)) = length(lstNsub);
        else
            e(lst(jj)) = kk;
        end
    end
    
    % find edges connected to nodes with 3+ edges
    E3p = nB(lstNsub(e(:,1)))>2 & nB(lstNsub(e(:,2)))>2;

    %
    % identify unique loops if there are any
    nLoops = size(e,1) - length(lstNsub) + 1;
    if nLoops>0
        c = grCycleBasis( e );
        
        % find loops with a unique 3+ edge
        lst = find( (sum(c,2)'.*E3p')==1 );
        % remove no more than one of these edges
        % from each loop and no nodes will be abandoned
        %
        % COULD THIS LEAVE A NODE WITH ONE EDGE IF I INADVERTANTLY REMOVE
        % EDGES FRO A GIVEN NODE TO REDUCE nB TO 1. I SHOULD CHECK THAT
        % nB>2 BEFORE REMOVING AN EDGE... EASY CHECK
        %
        if ~isempty(lst)
            lst3 = [];
            for iLoop = 1:nLoops
                lst2 = find(c(lst,iLoop).*E3p(lst)==1);
                if ~isempty(lst2)
                    iE = lst(lst2(1));
                    if nB(lstNsub(e(iE,1)))>2 & nB(lstNsub(e(iE,2)))>2
                        nB(lstNsub(e(iE,1))) = nB(lstNsub(e(iE,1))) - 1;
                        nB(lstNsub(e(iE,2))) = nB(lstNsub(e(iE,2))) - 1;
                        lst3(end+1) = iE;
                        nEdel = nEdel + 1;
                    end
                end
            end
            im.nodeEdges(lstE(lst3),:) = [];
        else
            n1 = n1 + 1;
            n2 = n2 + nLoops;
        end
                
    end
end
close(hwait)
disp( sprintf('No loops removed for %d nodes for a total of %d loops skipped',n1,n2) )

set(handles.textNumEdges,'string',sprintf('%d edges',size(im.nodeEdges,1)))
set(handles.uipanelGraph,'title',sprintf('Graph (%d nodes)',size(im.nodePos,1)));

set(handles.pushbuttonImageGraph,'enable','on')


    
% --------------------------------------------------------------------
function menu_resetNodeVesselType_Callback(hObject, eventdata, handles)
global im

ch = menu('Are you sure you want to reset the vessel type for all nodes?','Yes','No');

if ch==1
    im.nodeType = zeros(size(im.nodePos,1),1);
    set(handles.pushbuttonUpdateVesselMask,'enable','on');
end


% --------------------------------------------------------------------
function menu_resetGraph_Callback(hObject, eventdata, handles)
global im

ch = menu('Do you want to reset the graph?','Yes','No');
if ch==2
    return
end

for ii=1:length(im.seed)
    im.seed(ii).seedGraphed = 0;
end
if isfield(im,'nodePos')
    im = rmfield(im,'nodePos');
end
if isfield(im,'nodeDiam')
    im = rmfield(im,'nodeDiam');
end
if isfield(im,'nodeEdges')
    im = rmfield(im,'nodeEdges');
end
if isfield(im,'edgeFlag')
    im = rmfield(im,'edgeFlag');
end
if isfield(im,'SeedConnectivityMatrix')
    im = rmfield(im,'SeedConnectivityMatrix');
end
if isfield(im,'ConnectivityGrps')
    im = rmfield(im,'ConnectivityGrps');
end

set(handles.pushbuttonGraph,'enable','on')
set(handles.pushbuttonSeedConnectivity,'enable','on')
set(handles.pushbuttonImageGraph,'enable','off')
set(handles.pushbuttonRegraphNodes,'enable','off')
%set(handles.pushbuttonImaris,'enable','off')
set(handles.uipanelGraph,'title',sprintf('Graph (0 nodes)'));
set(handles.textNumEdges,'string','0 edges');
%set(handles.pushbuttonPruneLoops,'enable','off');
%set(handles.pushbuttonPruneLoops4,'enable','off');
%set(handles.pushbuttonPruneDiam,'enable','off');


% --------------------------------------------------------------------
function menu_viewE4pNodes_Callback(hObject, eventdata, handles)
global im


[nB,im] = nBupdate( im );
%nB = zeros(size(im.nodePos,1),1);
%nN = size(im.nodePos,1);
%for ii=1:nN
%    nB(ii) = length(find(im.nodeEdges(:,1)==ii | im.nodeEdges(:,2)==ii));
%end

% info update axes view
if im.ImageView==1
    pIdx = [1 2 3];
    xlim([1 size(im.III,2)]);
    ylim([1 size(im.III,1)]);
    nIs = size(im.III,3);
elseif im.ImageView==2
    pIdx = [1 3 2];
    xlim([1 size(im.III,2)]);
    ylim([1 size(im.III,3)]);
    nIs = size(im.III,2);
else
    pIdx = [2 3 1];
    xlim([1 size(im.III,1)]);
    ylim([1 size(im.III,3)]);
    nIs = size(im.III,1);
end
nS = str2num(get(handles.editNslices,'string'));



lst = find(nB>=4);

for ii=1:length(lst);
    nIdx = lst(ii);
    plotTreeSeg( im.nodeEdges, im.nodePos, nIdx, 10, 1, im.I, [], [], [], ones(size(im.nodeEdges)) );
    title( sprintf('Node %d of %d', ii, length(lst)) )
    rotate3d on
    
    ch = 2;
    while ch==2
        ch = menu('Option?','Okay','View 2D image','Cancel');
        if ch==2
            % update axes view
            zz = min(max(im.nodePos(nIdx,pIdx(3))-nS/2,1),nIs-nS);
            zz = min(zz,get(handles.slider1,'Max'));
            set(handles.slider1,'value',zz);
            set(handles.textSlices,'string',sprintf('Slice %d to %d',round(zz),round(zz+nS-1)))

            updateAxes(handles);
            hold on
            h=plot(im.nodePos(nIdx,pIdx(1)),im.nodePos(nIdx,pIdx(2)),'cp');
            set(h,'markersize',16);
            set(h,'markerfacecolor','c')
            for jj=1:6
                set(h,'marker','none')
                pause(0.1)
                set(h,'marker','p')
                pause(0.1)
            end
            hold off
        end
    end
    if ch==3
        break
    end
end

    


% --------------------------------------------------------------------
function menu_labelSegments_Callback(hObject, eventdata, handles)
global im

if ~isfield(im,'segNedges')
    menu('You must get segment and group information before using this tool. Push the [groups] button.','Okay');
    return;
end


[nB,im] = nBupdate( im );
%nB = zeros(1,size(im.nodePos,1));
%nN = size(im.nodePos,1);
%for ii=1:nN
%    nB(ii) = length(find(im.nodeEdges(:,1)==ii | im.nodeEdges(:,2)==ii));
%end

% info update axes view
if im.ImageView==1
    pIdx = [1 2 3];
    xlim([1 size(im.III,2)]);
    ylim([1 size(im.III,1)]);
    nIs = size(im.III,3);
elseif im.ImageView==2
    pIdx = [1 3 2];
    xlim([1 size(im.III,2)]);
    ylim([1 size(im.III,3)]);
    nIs = size(im.III,2);
else
    pIdx = [2 3 1];
    xlim([1 size(im.III,1)]);
    ylim([1 size(im.III,3)]);
    nIs = size(im.III,1);
end
nS = str2num(get(handles.editNslices,'string'));



nSeg = length(im.segNedges);
nUnlabeled = 0;
for iSeg = 1:nSeg
    lstN = find(im.nodeSegN==iSeg);
    lstN0 = find(im.nodeType(lstN)==0 & nB(lstN)<=2);
    if ~isempty(lstN0)
        nUnlabeled = nUnlabeled + 1;
    end
end

nFound = 0;
for iSeg = 1:nSeg
    lstN = find(im.nodeSegN==iSeg);
    lstN0 = find(im.nodeType(lstN)==0 & nB(lstN)<=2);
    if ~isempty(lstN0)
        nFound = nFound + 1;
        
        % update axes view
        zz = min(max(mean(im.nodePos(lstN,pIdx(3)),1)-nS/2,1),nIs-nS);
        zz = min(zz,get(handles.slider1,'Max'));
        set(handles.slider1,'value',zz);
        set(handles.textSlices,'string',sprintf('Slice %d to %d',round(zz),round(zz+nS-1)))

        updateAxes(handles);
        hold on
        h=plot(im.nodePos(lstN,pIdx(1)),im.nodePos(lstN,pIdx(2)),'cp');
        set(h,'markersize',16);
        set(h,'markerfacecolor','c')
        for jj=1:6
            set(h,'marker','none')
            pause(0.1)
            set(h,'marker','p')
            pause(0.1)
        end
        hold off
        
        ch = menu(sprintf('Label Segment %d of %d?\n Diam (min median max) = %.1f %.1f %.1f',nFound,nUnlabeled,...
            min(im.nodeDiam(lstN)),...
            median(im.nodeDiam(lstN)),max(im.nodeDiam(lstN))),...
            'Arteriole','Capillary','Vein','Skip','Cancel');
       if ch==5
            return;
        end
        if ch<4
            im.nodeType(lstN)=ch;
            set(handles.pushbuttonUpdateVesselMask,'enable','on');
        end
    end
end


% --------------------------------------------------------------------
function menu_pruneGroups_Callback(hObject, eventdata, handles)
global im

ch = menu(sprintf('Are you sure you want to prune all but the largest group of nodes?\nDid you first update Group Stats?'),'Yes','No');
if ch==2
    return;
end

nGrps = max(im.nodeGrp);
for ii=1:nGrps
    nNg(ii) = length(find(im.nodeGrp==ii));
end
[foo,idxG] = max(nNg);

nodeFlag = zeros(size(im.nodePos,1),1);
lst = find(im.nodeGrp==idxG);
nodeFlag(lst) = 1;
% remove nodes
[im.nodePos,im.nodeDiam,im.nodeDiamThetaIdx,im.nodeBC,im.nodeBCType,im.nodeType,im.nodeSegN,im.nodeEdges,im.edgeFlag] = removeNodes( nodeFlag, im.nodePos, im.nodeDiam, im.nodeDiamThetaIdx, im.nodeBC, im.nodeBCType, im.nodeType, im.nodeSegN, im.nodeEdges );
im.nBflag = 1;

set(handles.pushbuttonImageGraph,'enable','on')
set(handles.pushbuttonUpdateVesselMask,'enable','on');
set(handles.textNumEdges,'string',sprintf('%d edges',size(im.nodeEdges,1)))
set(handles.uipanelGraph,'title',sprintf('Graph (%d nodes)',size(im.nodePos,1)));


% --------------------------------------------------------------------
function menu_Flow_LitPressure_Callback(hObject, eventdata, handles)
global im

ch = menu(sprintf('Do you want to set the boundary condition\nfor E1 nodes to literature value for pressure?'),...
    'yes','no');

if ch==2
    return;
end

[nB,im] = nBupdate( im );
%nB = zeros(1,size(im.nodePos,1));
%nN = size(im.nodePos,1);
%for ii=1:nN
%    nB(ii) = length(find(im.nodeEdges(:,1)==ii | im.nodeEdges(:,2)==ii));
%end

lstNe1 = find(nB==1);
for ii = 1:length(lstNe1)
    iN = lstNe1(ii);
    
    im.nodeBCType(iN) = 3;
end

set(handles.pushbuttonCalcFlowCircuitEq,'enable','on')


% --------------------------------------------------------------------
function menu_Flow_useSegments_Callback(hObject, eventdata, handles)
global im

if strcmpi(get(handles.menu_Flow_useSegments,'checked'),'on')
    set(handles.menu_Flow_useSegments,'checked','off')
    im.flagUseSegments = 0;
else
    set(handles.menu_Flow_useSegments,'checked','on')
    im.flagUseSegments = 1;
end


% --- Executes on button press in pushbuttonCenterNodesXYZ.
function pushbuttonCenterNodesXYZ_Callback(hObject, eventdata, handles)

centerStep1vox = get(handles.checkboxCenterStep,'value');

if strcmpi(get(handles.menu_VisualizeCentering,'checked'),'off')
    flagVisualize = 0;
else
    flagVisualize = 1;
end

thresh = ceil(str2num(get(handles.editImageThresh,'string'))*32);

imView3d_CenterNodesXYZ( centerStep1vox, flagVisualize, thresh );

set(handles.pushbuttonRegraphNodes,'enable','on')
set(handles.pushbuttonImageGraph,'enable','on')


% --------------------------------------------------------------------
function menu_nNodesBelowThresh_Callback(hObject, eventdata, handles)
global im

Ithresh = ceil(str2num(get(handles.editImageThresh,'string'))*32);

n0 = 0;
nN = size(im.nodePos,1);

hwait = waitbar(0,'Counting nodes below image threshold');
for ii=1:nN
    waitbar(ii/nN,hwait);
    pos = round(im.nodePos(ii,:));
    if im.I(pos(2),pos(1),pos(3))<Ithresh
        n0 = n0 + 1;
    end
end
close(hwait)

% [nIx,nIy,nIz] = size(im.I);
% lstIdx = min(round(im.nodePos(:,2)),nIy) + (min(round(im.nodePos(:,1)),nIx)-1)*nIy + (min(round(im.nodePos(:,3)),nIz)-1)*nIx*nIy;
% n0 = find( im.I( lstIdx )==0 );

menu(sprintf('There are %d nodes below Image Thresh',n0),'Okay');
disp(sprintf('There are %d nodes below Image Thresh',n0));


% --------------------------------------------------------------------
function menu_viewDanglingSegments_Callback(hObject, eventdata, handles)
% view end point nodes
% with segLen < segDiam
% and query to manually remove
global im

menu_viewDanglingSegmentsSub(hObject, 1, handles );



function menu_viewDanglingSegmentsSub(hObject, eventdata, handles );
global im

[nB,im] = nBupdate( im );
%hwait = waitbar( 0, 'Preparing to view segments');
%nB = zeros(1,size(im.nodePos,1));
%nN = size(im.nodePos,1);
%for ii=1:nN
%    nB(ii) = length(find(im.nodeEdges(:,1)==ii | im.nodeEdges(:,2)==ii));
%end
%close(hwait)

if eventdata==1 % view dangling segments with len < diam
    lst = find(nB==1);
    lstSeg = im.nodeSegN(lst);
    lst2 = find(im.segLen(lstSeg)<im.segDiam(lstSeg));
    lst = lst(lst2);
    
elseif eventdata==2
    lenSeg = menu('View segments of length (# edges)',...
        sprintf('1 (# = %d)',length(find(im.segNedges==1))),...
        sprintf('2 (# = %d)',length(find(im.segNedges==2))),...
        sprintf('3 (# = %d)',length(find(im.segNedges==3))),...
        sprintf('4 (# = %d)',length(find(im.segNedges==4))),...
        sprintf('5 (# = %d)',length(find(im.segNedges==5))),...
        'Cancel');
    if lenSeg==6
        return;
    end
    lstSeg = find(im.segNedges==lenSeg);
    lst = im.segEndNodes(lstSeg,1);
end


nodeFlag = ones(size(im.nodePos,1),1);
eFlag = ones(size(im.nodeEdges,1),1);
for ii=1:length(lst)
    if eventdata==1
        iSeg = im.nodeSegN(lst(ii));
    elseif eventdata==2
        iSeg = lstSeg(ii);
    end
    lst3 = find(im.edgeSegN==iSeg);

    ch = 3;
    nEdges = 15;
    while ch>2
        eFlag(lst3) = 2;
        [pMin,pMax,Tree]=plotTreeSeg(im.nodeEdges,im.nodePos,lst(ii),nEdges,[1 1],im.I,[],[],[],eFlag);
        iSeg = im.nodeSegN(lst(ii));
        title( sprintf('(%d of %d) Node%d  Seg%d  diam=%.0f  len=%.0f',ii,length(lst),lst(ii),iSeg,im.segDiam(iSeg),im.segLen(iSeg)) )

        figure(1)
        z1 = floor(min(im.nodePos(im.segEndNodes(iSeg,:),3)));
        z2 = ceil(max(im.nodePos(im.segEndNodes(iSeg,:),3)));
        imagesc( max(im.I(:,:,z1:z2),[],3) )
        hold on
        for jj=1:length(Tree)
            n1 = im.nodeEdges(Tree(jj),1);
            n2 = im.nodeEdges(Tree(jj),2);
            hl=plot([im.nodePos(n1,1) im.nodePos(n2,1)],[im.nodePos(n1,2) im.nodePos(n2,2)],'r-');
            set(hl,'linewidth',2)
            if eFlag(Tree(jj))~=2
                set(hl,'color','y')
            end
        end
        plot([im.nodePos(lst(ii),1)],[im.nodePos(lst(ii),2)],'r*')
        hold off
        title( sprintf('Z = %.1f',im.segPos(iSeg,3)) )
        colormap gray

        eFlag(lst3) = 1;
        
        ch = menu('What to do?','Keep','Delete','20 edges','25 edges','+5 edges','Quit');
        if ch==3
            nEdges = 20;
        elseif ch==4
            nEdges = 25;
        elseif ch==5
            nEdges = nEdges + 5;
        elseif ch==6
            break
        elseif ch==2
            eFlag(lst3) = 0;
            nLst = im.nodeEdges(lst3,:);
            nLst = unique(nLst(:));
            foo = find(nB(nLst)==1 | nB(nLst)==2);
            nodeFlag(nLst(foo)) = 0;
        end
    end
    if ch==6
        break
    end
end

ch = menu('Remove deleted nodes and recalculate Segments','Yes','No');
drawnow
if ch==1
    [im.nodePos,im.nodeDiam,im.nodeDiamThetaIdx,im.nodeBC,im.nodeBCType,im.nodeType,...
        im.nodeSegN,im.nodeEdges,im.edgeFlag] = removeNodes( nodeFlag, im.nodePos, ...
        im.nodeDiam, im.nodeDiamThetaIdx, im.nodeBC, im.nodeBCType, im.nodeType, im.nodeSegN, im.nodeEdges );
    im.nBflag = 1;
    im = nodeGrps(im);
    set(handles.textGrpStats,'string',sprintf('%d Groups, %d Segments',max(im.nodeGrp),length(im.segNedges)) )
end


