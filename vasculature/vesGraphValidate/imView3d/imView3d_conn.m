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

% Last Modified by GUIDE v2.5 08-Aug-2007 10:10:42

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

% Choose default command line output for imView3d
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

if length(varargin)==1
    im.filenm = varargin{1};
else
    [fname fpath] = uigetfile( '*.img','Load file');
    if fname==0
        return
    end
    cd(fpath)
    im.filenm = fname(1:end-4);
end


filenm = sprintf( '%s.img',im.filenm );
fid = fopen( filenm, 'rb' );
n = fread( fid, 3, 'int');
nx = n(1); ny = n(2); nz = n(3);
im.I = fread( fid, nx*ny*nz, 'float' );
im.I = reshape(im.I, [ny nx nz] );
fclose(fid);
im.III = zeros(size(im.I));
im.III = 32*(im.I-0.1)/(1-0.1);
im.nZ = size(im.I,3);

im.nSeed = 0;
im.seedDirectionFlag = 0;
im.curSeed = 0;

set(handles.figure1,'Name',sprintf('vesselGraph %s',im.filenm))

set(handles.slider1,'max',im.nZ)
set(handles.slider1,'value',1)
set(handles.slider1,'min',1)

axes( handles.axes1 );
xlim([1 size(im.I,2)]);
ylim([1 size(im.I,1)]);

updateAxes( handles );

set( handles.axes1,'ButtonDownFcn','imView3d(''axes1ButtonDown_Callback'',gcbo,[],guidata(gcbo))' );
set(get(handles.axes1,'children'), ...
    'ButtonDownFcn', 'imView3d(''axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))' );


% UIWAIT makes imView3d wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = imView3d_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
global im

z = round(get(handles.slider1,'value'));
nS = str2num(get(handles.editNslices,'string'));
if (z+nS-1)>size(im.I,3)
    z = size(im.I,3)-nS+1;
    set(handles.slider1,'value',z)
end
set(handles.textSlices,'string',sprintf('Slice %d to %d',z,z+nS-1))

updateAxes( handles );





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateAxes( handles )
global im

axes( handles.axes1 );
xrng = xlim();
yrng = ylim();


zzMin = round(get(handles.slider1,'value'));
nS = str2num(get(handles.editNslices,'string'));
zzMax = zzMin+nS-1;

imagesc( max(im.III(:,:,zzMin:zzMax),[],3), [1 66] )
xlim(xrng)
ylim(yrng)

cm = gray(32);
% purple
%cm(33:64,:) = hsv2rgb([ones(32,1)*0.8333 ones(32,1) [1:32]'/32]);
% yellow
%cm(33:64,:) = hsv2rgb([ones(32,1)*0.1667 ones(32,1) [1:32]'/32]);
% red
cm(33:64,:) = hsv2rgb([ones(32,1)*0 ones(32,1) [1:32]'/32]);
cm(65,:) = [1 1 0];
cm(66,:) = [1 0 0];

colormap(cm)

set( handles.axes1,'ButtonDownFcn','imView3d(''axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))' );
set(get(handles.axes1,'children'), ...
    'ButtonDownFcn', 'imView3d(''axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))' );

for ii=1:im.nSeed
    if im.seed(ii).pos(3)>=zzMin & im.seed(ii).pos(3)<=zzMax
        ht=text(im.seed(ii).pos(1),im.seed(ii).pos(2),num2str(ii));
        set(ht,'color',[0.1667 1 1])
        if ~isempty(im.seed(ii).pos2)
            hl = line([im.seed(ii).pos(1) im.seed(ii).pos2(1)],[im.seed(ii).pos(2) im.seed(ii).pos2(2)]);
            set(hl,'color',[0.1667 1 1])
        end
    end
end
    





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function editNslices_Callback(hObject, eventdata, handles)
global im

z = get(handles.slider1,'value');
nS = str2num(get(handles.editNslices,'string'));

set(handles.slider1,'value',min(z,im.nZ-nS+1));
set(handles.slider1,'max',im.nZ-nS+1);

z=round(z);
set(handles.textSlices,'string',sprintf('Slice %d to %d',z,z+nS-1))

updateAxes( handles );





% --- Executes on button press in radiobuttonZoom.
function radiobuttonZoomPan_Callback(hObject, eventdata, handles)
switch eventdata
    case 1
        pan off
        zoom on
        set(handles.radiobuttonZoom,'value',1);
        set(handles.radiobuttonPan,'value',0);
        set(handles.radiobuttonSeedPoint,'value',0);
    case 2
        zoom off
        pan on
        set(handles.radiobuttonZoom,'value',0);
        set(handles.radiobuttonPan,'value',1);
        set(handles.radiobuttonSeedPoint,'value',0);
    case 3
        zoom off
        pan off
        set(handles.radiobuttonZoom,'value',0);
        set(handles.radiobuttonPan,'value',0);
        set(handles.radiobuttonSeedPoint,'value',1);
        set( handles.axes1,'ButtonDownFcn','imView3d(''axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))' );
        set(get(handles.axes1,'children'), ...
            'ButtonDownFcn', 'imView3d(''axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))' );
end




% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
global im

if strcmp(get(gcf,'Selectiontype'),'normal')
    pos = get(gca,'CurrentPoint');
else
    return
end

nSeed = im.nSeed;

xx = round(pos(1,1));
yy = round(pos(1,2));

z = round(get(handles.slider1,'value'));
nS = str2num(get(handles.editNslices,'string'));
zzMin = z;
zzMax = z+nS-1;

% select segment
h=20;
lst = find(xx>im.nodePos(:,1)-h & xx<im.nodePos(:,1)+h & ...
     yy>im.nodePos(:,2)-h & yy<im.nodePos(:,2)+h & ...
     im.nodePos(:,3)>=zzMin & im.nodePos(:,3)<=zzMax);
if isempty(lst)
    return;
end
for ii=1:length(lst)
    rsep(ii) = norm([im.nodePos(lst(ii),1)-xx  im.nodePos(lst(ii),2)-yy]);
end
[foo,nIdx] = min(rsep);
nIdx = lst(nIdx);

nodeConn = im.nodeConn;
lstNodes = [];
nLst = 0;
while length(nodeConn{nIdx}.prev)==1 & length(nodeConn{nIdx}.next)==1
    nLst = nLst + 1;
    lstNodes(nLst) = nIdx;
    nIdx = nodeConn{nIdx}.prev;
    im.nodeFlag(nIdx) = xor(im.nodeFlag(nIdx),1);
end
if nLst>0
    nIdx = lstNodes(1);
    lstNodes = lstNodes(end:-1:1);
    im.nodeFlag(nIdx) = xor(im.nodeFlag(nIdx),1);
end
while length(nodeConn{nIdx}.prev)==1 & length(nodeConn{nIdx}.next)==1
    nLst = nLst + 1;
    lstNodes(nLst) = nIdx;
    nIdx = nodeConn{nIdx}.next;
    if nLst>1
        im.nodeFlag(nIdx) = xor(im.nodeFlag(nIdx),1);
    end
end
pushbuttonImageGraph_Callback(hObject, eventdata, handles);

return





[foo,zIdx] = max(im.III(yy,xx,zzMin:zzMax),[],3);
zz = zzMin+zIdx-1;

if ~im.seedDirectionFlag
    if foo>0
        nSeed = nSeed + 1;
        im.nSeed = nSeed;
        im.seed(nSeed).pos = [xx yy zz];
        im.seed(nSeed).pos2 = [];
        im.seed(nSeed).bidirectional = 1;
        im.seedDirectionFlag = 1;
        im.seed(nSeed).seedLaunched = 0;
        im.seed(nSeed).seedProcessed = 0;

        im.seed(nSeed).nTrks2Run = str2num(get(handles.editAutoSeedNtrks,'string'));

        im.seed(nSeed).nTrksRan = 0;
        set(handles.editNtrks,'string',sprintf('%d',im.seed(nSeed).nTrks2Run))

        im.curSeed = nSeed;
        
        set(handles.uipanelSeed,'title',sprintf('Seed %d of %d',nSeed,nSeed))
        set(handles.textSeedInfo,'string','SET DIRECTION');
        set(handles.checkboxBidirectional,'value',im.seed(im.curSeed).bidirectional);
    end
else
    im.seed(nSeed).pos2 = [xx yy zz];
    im.seed(nSeed).cxyz = im.seed(nSeed).pos2 - im.seed(nSeed).pos;
    im.seedDirectionFlag = 0;

    set(handles.textSeedInfo,'string','');
    set(handles.pushbuttonRunSeeds,'Enable', 'on' );
end

updateAxes( handles )


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
if im.seed(im.curSeed).seedProcessed
    set(handles.textSeedInfo,'string',sprintf('Processed %d',im.seed(im.curSeed).nTrksRan));
elseif im.seed(im.curSeed).seedLaunched
    set(handles.textSeedInfo,'string','Seed Launched');
end

z = round(get(handles.slider1,'value'));
nS = str2num(get(handles.editNslices,'string'));
zSeed = im.seed(im.curSeed).pos(3);
zNew = min(max(zSeed-floor(nS/2),1),im.nZ-nS);
set(handles.slider1,'value',zNew);
set(handles.textSlices,'string',sprintf('Slice %d to %d',zNew,zNew+nS-1))

updateAxes( handles );


% --- Executes on button press in checkboxBidirectional.
function checkboxBidirectional_Callback(hObject, eventdata, handles)
global im

im.seed(im.curSeed).bidirectional = get(handles.checkboxBidirectional,'value');




% --- Executes on button press in pushbuttonSave.
function pushbuttonSave_Callback(hObject, eventdata, handles)
global im

[fname fpath] = uiputfile( '*.seed','Save seed file');
if fname==0
    return
end
save([fpath fname],'im');



% --- Executes on button press in pushbuttonLoad.
function pushbuttonLoad_Callback(hObject, eventdata, handles)
global im

[fname fpath] = uigetfile( '*.seed','Load seed file');
if fname==0
    return
end
load([fpath fname],'-mat','im');

UpdateSeeds( handles )



function UpdateSeeds( handles )
global im

set(handles.uipanelSeed,'title',sprintf('Seed %d of %d',im.curSeed,im.nSeed))
set(handles.textSeedInfo,'string','');
set(handles.checkboxBidirectional,'value',im.seed(im.curSeed).bidirectional);
set(handles.editNtrks,'string',sprintf('%d',im.seed(im.curSeed).nTrks2Run))
if im.seed(im.curSeed).seedProcessed
    set(handles.textSeedInfo,'string',sprintf('Processed %d',im.seed(im.curSeed).nTrksRan));
elseif im.seed(im.curSeed).seedLaunched
    set(handles.textSeedInfo,'string','Seed Launched');
end

set(handles.checkboxBidirectional,'value',im.seed(im.curSeed).bidirectional);

nS = str2num(get(handles.editNslices,'string'));
zSeed = im.seed(im.curSeed).pos(3);
zNew = min(max(zSeed-floor(nS/2),1),im.nZ-nS);
set(handles.slider1,'max',im.nZ)
set(handles.slider1,'min',1)
set(handles.slider1,'value',zNew);
set(handles.textSlices,'string',sprintf('Slice %d to %d',zNew,zNew+nS-1))

set(handles.figure1,'Name',sprintf('vesselGraph %s',im.filenm))

if isfield(im,'nodePos')
    set(handles.uipanelGraph,'title',sprintf('Graph (%d nodes)',size(im.nodePos,1)));
end


axes( handles.axes1 );
xlim([1 size(im.I,2)]);
ylim([1 size(im.I,1)]);

updateAxes( handles )




% --- Executes on button press in pushbuttonRunSeeds.
function pushbuttonRunSeeds_Callback(hObject, eventdata, handles)
global im

flag = 0;
fid = fopen('launchScript.bat','w');
for ii=1:im.nSeed
    if ~im.seed(ii).seedLaunched & im.seed(ii).nTrksRan==0
        fprintf( fid, sprintf( 'del %s_seed%d.trk\n', im.filenm, ii ) );
    end
end
for ii=1:im.nSeed
    if ~im.seed(ii).seedLaunched
        fprintf( fid, sprintf( 'DTItrack %s_seed%d\n', im.filenm, ii ) );
        
        fid2 = fopen( sprintf('%s_seed%d.inp',im.filenm,ii), 'w' );
        fprintf( fid2, '%d %d %d\n', size(im.I,2), size(im.I,1), size(im.I,3) );
        fprintf( fid2, '%s.img\n', im.filenm );
        fprintf( fid2, '%d %d\n%d\n', im.seed(ii).nTrks2Run, im.seed(ii).nTrksRan,...
            -im.seed(ii).bidirectional );
        fprintf( fid2, '-563980456\n' );  % RANDOM NUMBER SEED, GENERATE FROM CLOCK!
        fprintf( fid2, '%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f', ...
            im.seed(ii).pos(1)-1, im.seed(ii).pos(2)-1, im.seed(ii).pos(3)-1,...
            im.seed(ii).cxyz(1), im.seed(ii).cxyz(2), im.seed(ii).cxyz(3) );
        fclose( fid2 );
        
        im.seed(ii).seedLaunched = 1;
        im.seed(ii).seedProcessed = 0;
        im.seed(ii).seedGraphed = 0;
        flag = 1;
    end
end
fprintf( fid, 'exit\n' );
fclose(fid);

if flag
    !launchScript.bat &

    set(handles.pushbuttonRunSeeds,'Enable', 'off' );
    set(handles.pushbuttonProcSeeds,'Enable','On')
end



% --- Executes on button press in pushbuttonProcSeeds.
function pushbuttonProcSeeds_Callback(hObject, eventdata, handles)
global im
% 
% IIsum = im.II;
% for ii=1:im.nSeed
%     if ~im.seed(ii).seedLaunched
%         disp(sprintf('Running track %d...',ii) );
%         [II,info,trkPos] = DTItrack( im.I, ...
%             [im.seed(ii).pos(1) im.seed(ii).pos(2) im.seed(ii).pos(3)], ...
%             [im.seed(ii).cxyz(1) im.seed(ii).cxyz(2) im.seed(ii).cxyz(3)], ...
%             10*(im.seed(ii).bidirectional+1), -im.seed(ii).bidirectional);
%         eval(sprintf('save track%d.mat II info trkPos',ii));        
%         
%         IIsum = IIsum + II;
% 
%         im.seed(ii).seedLaunched = 1;
%         im.seed(ii).trkPos = trkPos;
%     end
% end
% im.II = IIsum;
% 
% updateImageOverlay(1);

flag = 0;
flag2 = 1;
for ii=1:im.nSeed
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

updateAxes( handles )
set(handles.uipanelSeed,'title',sprintf('Seed %d of %d',im.curSeed,im.nSeed))
set(handles.textSeedInfo,'string','');
set(handles.checkboxBidirectional,'value',im.seed(im.curSeed).bidirectional);
set(handles.editNtrks,'string',sprintf('%d',im.seed(im.curSeed).nTrks2Run))
if im.seed(im.curSeed).seedProcessed
    set(handles.textSeedInfo,'string',sprintf('Processed %d',im.seed(im.curSeed).nTrksRan));
elseif im.seed(im.curSeed).seedLaunched
    set(handles.textSeedInfo,'string','Seed Launched');
end

if flag2
    %set(handles.pushbuttonRunSeeds,'Enable','On')
    %set(handles.pushbuttonGraph,'Enable','On')
    set(handles.pushbuttonProcSeeds,'Enable','Off')
    im.pushButtonGraphEnable = 1;  %WHAT IS THIS USED FOR? MAYBE NOT NEEDED.
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
                            [foo,mm] = max(lstR(ll));
                            im.seed(ii).trkConn(lstTrks(kk)) = lstR(ll(mm));
                        end
                    end
                else
                    M(ii,jj) = 0;
                end
                % I still need to record the last trk pt (lstR) that connects
                % with a seed to cut graphing at that point.
            else
                M(ii,jj) = -1;
            end
        end
    end
end
close(hWait)

% find the number of unique groups of connected seeds
Mgrps = zeros(1,im.nSeed);
nGrps = 0;
for ii=1:im.nSeed
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
    nodePos = [0 0 0];
else
    nodePos = im.nodePos;
    nNodes = size(nodePos,1);
    nodeConn = im.nodeConn;
end

hWait = waitbar( 0, 'Graphing tracks...' );

hxy = str2num(get(handles.editGraphHxy,'string'));
hz = str2num(get(handles.editGraphHz,'string'));
for ii=1:im.nSeed
    waitbar(ii/im.nSeed,hWait);
    if im.seed(ii).seedProcessed & ~im.seed(ii).seedGraphed
        [ii nNodes]
        nTrkPosMax = size(im.seed(ii).trkPos,1);
        nTrks = size(im.seed(ii).trkPos,3);
        for jj=1:nTrks
            kk = 1;
            lastNode = 0;
            lastPtSum = 0;
            lastPos = 0;
            while im.seed(ii).trkPos(kk,1,jj)>0 & kk<=im.seed(ii).trkConn(jj) & kk<nTrkPosMax
                pos = squeeze(im.seed(ii).trkPos(kk,:,jj));
                lst = find(pos(1)>=(nodePos(:,1)-hxy) & pos(1)<=(nodePos(:,1)+hxy) & ...
                    pos(2)>=(nodePos(:,2)-hxy) & pos(2)<=(nodePos(:,2)+hxy) & ...
                    pos(3)>=(nodePos(:,3)-hz) & pos(3)<=(nodePos(:,3)+hz) );
                if isempty(lst)
                    nNodes = nNodes+1;
                    nodePos(nNodes,:) = pos;
                    if norm(pos-lastPos)<2*hxy | lastNode==0
                        if lastNode>0
                            nodeConn{lastNode}.next(end+1) = nNodes;
                        end
                        nodeConn{nNodes}.prev(1) = lastNode;
                        nodeConn{nNodes}.next = [];
                    else
                        %need to find closest node and check if it follows back
                        %to last node. IF NOT THEN USE LAST NODE!!!
                        d=sum((ones(nNodes-1,1)*pos-nodePos(1:end-1,:)).^2,2).^0.5;
                        [foo,idx]=min(d);
                        if ~isempty(idx)
                            lastNode = idx;
                            nodeConn{lastNode}.next(end+1) = nNodes;
                            nodeConn{nNodes}.prev(1) = lastNode;
                            nodeConn{nNodes}.next = [];
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
                    % this is assuming length(lst)==1, and other assumptions
                    % that should be okay, but I should worry about it
                    if length(lst>1)
                        clear d
                        for iLst=1:length(lst)
                            d(iLst) = norm(pos-nodePos(lst(iLst),:));
                        end
                        [foo, closestNode] = min(d);
                    else
                        closestNode = 1;
                    end
                    if lst(closestNode)~=lastNode & lastNode~=0
                        nodeConn{lastNode}.next(end+1) = lst(closestNode);
                        nodeConn{lst(closestNode)}.next(end+1) = lastNode;

                        nodeConn{lastNode}.next = setdiff(unique(nodeConn{lastNode}.next),[nodeConn{lastNode}.prev lastNode]);
                        nodeConn{lst(closestNode)}.next = setdiff(unique(nodeConn{lst(closestNode)}.next),[nodeConn{lst(closestNode)}.prev lst(closestNode)]);
                    end
                    lastNode = lst(closestNode);
                    lastPos = nodePos(lastNode,:);
                end
                kk = kk + 1;
            end
        end
        
        im.seed(ii).seedGraphed = 1;
    end
end

close(hWait);

im.nodePos = nodePos;
im.nodeConn = nodeConn;
im.nodeFlag = zeros(1,nNodes);

set(handles.pushbuttonGraph,'Enable','Off')
set(handles.pushbuttonImageGraph,'Enable','On')
set(handles.pushbuttonImaris,'Enable','On')
im.pushButtonGraphEnable = 0;
set(handles.uipanelGraph,'title',sprintf('Graph (%d nodes)',nNodes));




% --- Executes on button press in pushbuttonImageGraph.
function pushbuttonImageGraph_Callback(hObject, eventdata, handles)
global im

updateImageOverlay(0);%get(handles.checkboxWithTracks,'value'));

%if get(handles.checkboxWithGraph,'value')
    nodePos = im.nodePos;
    nodeConn = im.nodeConn;
    if isfield(im,'nodeFlag')
        nodeFlag = im.nodeFlag;
    else
        nodeFlag = zeros(size(nodePos,1),1);
    end
    nNodes = size(nodePos,1);

    for ii=1:nNodes
        for jj=1:length(nodeConn{ii}.next)
            pos0 = nodePos(ii,:);
            pos1 = nodePos(nodeConn{ii}.next(jj),:);
            rsep = norm(pos1-pos0);
            if rsep>0
                cxyz = (pos1-pos0) / rsep;
                rstep = 0;
                pos = pos0;
                while rstep<rsep
                    im.III(round(pos(2)),round(pos(1)),round(pos(3))) = 65 + nodeFlag(ii);
                    pos = pos + cxyz*0.5;
                    rstep = rstep + 0.5;
                end
            end
        end
    end
%end

updateAxes( handles );

set(handles.pushbuttonImageGraph,'enable','off')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateImageOverlay( withTracks )
global im

im.III = zeros(size(im.I));
im.III = 32*(im.I-0.1)/(1-0.1);
if withTracks
    lst = find(im.II>0);
    im.III(lst) = 32+32*min(im.II(lst),4)/4;
end









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

% create new filament
vNewFilament = vImarisApplication.mFactory.CreateFilament;
vNewFilament.mIndexT = vImarisApplication.mVisibleIndexT;
vNewFilament.mName = 'foo';
vNewFilament.SetColor(1,0,0, 0.0);
vParent.AddChild(vNewFilament);

nodePos = im.nodePos;
nodeConn = im.nodeConn;
nNodes = size(nodePos,1);
pos = (nodePos-1) * 352.778;
rad = ones(size(pos,1),1);
eIdx = 0;
for ii=1:nNodes
	for jj=1:length(nodeConn{ii}.next)
        eIdx = eIdx + 1;
        edges(eIdx,:) = [ii-1 nodeConn{ii}.next(jj)-1];
    end
end

vNewFilament.Set( pos, rad, edges )

% Manually process graph and then get
ch = menu('Manually Center filament and get Diameters. Press Okay when done.','Okay','Cancel');
if ch==2
    return;
end
[posNew, radNew, edgesNew] = vNewFilament.Get;
im.nodePos = posNew/352.778 + 1;
im.nodeConn = nodeConn;  % This may be in need of modifying to be like edges
im.nodeDiam = radNew*2;

set(handles.pushbuttonRegraphNodes,'Enable','On')
set(handles.pushbuttonImaris,'Enable','Off')






% --- Executes on button press in pushbuttonRegraphNodes.
function pushbuttonRegraphNodes_Callback(hObject, eventdata, handles)
global im

% OLD ROUTINE, DECIDED TO COPY THE FIRST GRAPH ROUTINE
%
% % Re-graph after Imaris tuning
% nodePos = im.nodePos;
% nodeConn = im.nodeConn;
% nNodes = size(nodePos,1);
% nodeDiam = im.nodeDiam;
% 
% hxy = str2num(get(handles.editGraphHxy,'string'));
% hz = str2num(get(handles.editGraphHz,'string'));
% 
% nNodesUnique = 1;
% nodeMap = zeros(nNodes,1);
% nodeUnique = zeros(nNodes,1);
% 
% nodeMap(1) = 1;
% nodePosNew = nodePos(1,:);
% nodeUnique(1) = 1;
% 
% for ii=2:nNodes
%     pos = nodePos(ii,:);
%     lst = find(pos(1)>=(nodePosNew(:,1)-hxy) & pos(1)<=(nodePosNew(:,1)+hxy) & ...
%         pos(2)>=(nodePosNew(:,2)-hxy) & pos(2)<=(nodePosNew(:,2)+hxy) & ...
%         pos(3)>=(nodePosNew(:,3)-hz) & pos(3)<=(nodePosNew(:,3)+hz) );
%     if isempty(lst)
%         nNodesUnique = nNodesUnique+1;
% 
%         nodeMap(ii) = nNodesUnique;
%         nodeUnique(ii) = 1;
%         
%         nodePosNew(nNodesUnique,:) = pos;
%         nodeDiamNew(nNodesUnique) = nodeDiam(ii);
%     else
%         if length(lst>1)
%             clear d
%             for iLst=1:length(lst)
%                 d(iLst) = norm(pos-nodePosNew(lst(iLst),:));
%             end
%             [foo, closestNode] = min(d);
%         else
%             closestNode = 1;
%         end
%         nodeMap(ii) = lst(closestNode);
%     end
% end
% 
% 
% nodeConnNew = [];
% for ii=1:nNodes
%     if nodeUnique(ii)
%         if nodeConn{ii}.prev~=0
%             nodeConnNew{nodeMap(ii)}.prev = nodeMap(nodeConn{ii}.prev);
%         else
%             nodeConnNew{nodeMap(ii)}.prev = 0;
%         end
%         nodeConnNew{nodeMap(ii)}.next = [];
%     end
%     nodeConnNew{nodeMap(ii)}.next = setdiff(unique([nodeConnNew{nodeMap(ii)}.next nodeMap(nodeConn{ii}.next)']),[nodeConnNew{nodeMap(ii)}.prev nodeMap(ii)]);
% end
% 
% im.nodePos = nodePosNew;
% im.nodeConn = nodeConnNew;
% im.nodeDiam = nodeDiamNew;
% im.nodeFlag = zeros(1,size(nodePosNew,1));
% 
% set(handles.pushbuttonRegraphNodes,'Enable','Off')
% set(handles.pushbuttonImageGraph,'Enable','On')
% set(handles.uipanelGraph,'title',sprintf('Graph (%d nodes)',size(im.nodePos,1)));





% Re-graph after Imaris tuning
nodePos = im.nodePos;
nodeConn = im.nodeConn;
nNodes = size(nodePos,1);
nodeDiam = im.nodeDiam;

hxy = str2num(get(handles.editGraphHxy,'string'));
hz = str2num(get(handles.editGraphHz,'string'));

nNodesNew = 1;
nodeMap = zeros(nNodes,1);
%nodeUnique = zeros(nNodes,1);

nodeMap(1) = 1;
nodePosNew = nodePos(1,:);
nodeConnNew{1}.prev = 0;
nodeConnNew{1}.next = [];
%nodeUnique(1) = 1;


hWait = waitbar( 0, 'Re-Graphing Nodes...' );

lastPos = nodePos(1,:);
lastNode =1;
for ii=2:nNodes
    waitbar(ii/nNodes,hWait);
    pos = nodePos(ii,:);
    lst = find(pos(1)>=(nodePosNew(:,1)-hxy) & pos(1)<=(nodePosNew(:,1)+hxy) & ...
        pos(2)>=(nodePosNew(:,2)-hxy) & pos(2)<=(nodePosNew(:,2)+hxy) & ...
        pos(3)>=(nodePosNew(:,3)-hz) & pos(3)<=(nodePosNew(:,3)+hz) );
    if isempty(lst)
        nNodesNew = nNodesNew+1;
        nodePosNew(nNodesNew,:) = pos;
        nodeDiamNew(nNodesNew) = nodeDiam(ii);
        
        nodeMap(ii) = nNodesNew;
        
        nodeConnNew{nNodesNew}.prev = nodeConn{ii}.prev;
        nodeConnNew{nNodesNew}.next = nodeConn{ii}.next;
    else
        clear d
        for iLst=1:length(lst)
            d(iLst) = norm(pos-nodePosNew(lst(iLst),:));
        end
        [foo, closestNode] = min(d);
        nodeMap(ii) = lst(closestNode);

        lst2 = nodeConn{ii}.prev;
        for jj=1:length(lst2)
            if lst2(jj)<ii & lst2(jj)~=0
                kk = nodeMap(lst2(jj));
                lst3 = find(nodeConnNew{kk}.next~=ii);
                nodeConnNew{kk}.next = [nodeConnNew{kk}.next(lst3) nodeConn{ii}.next];
            elseif lst2(jj)~=0
                kk = lst2(jj);
                lst3 = find(nodeConn{kk}.next~=ii);
                nodeConn{kk}.next = [nodeConn{kk}.next(lst3) nodeConn{ii}.next];
            end
        end
        lst2 = nodeConn{ii}.next;
        for jj=1:length(lst2)
            if lst2(jj)<ii
                kk = nodeMap(lst2(jj));
                lst3 = find(nodeConnNew{kk}.prev~=ii);
                nodeConnNew{kk}.prev = [nodeConnNew{kk}.prev(lst3) nodeConn{ii}.prev];
            else
                kk = lst2(jj);
                lst3 = find(nodeConn{kk}.prev~=ii);
                nodeConn{kk}.prev = [nodeConn{kk}.prev(lst3) nodeConn{ii}.prev];
            end
        end
    end
end

for ii=1:nNodesNew
    nodeConnNew{ii}.prev = nodeMap(nodeConnNew{ii}.prev);
    nodeConnNew{ii}.next = nodeMap(nodeConnNew{ii}.next);
end

close(hWait);

im.nodePos = nodePosNew;
im.nodeConn = nodeConnNew;
im.nodeDiam = nodeDiamNew;
im.nodeFlag = zeros(1,nNodes);

set(handles.pushbuttonRegraphNodes,'Enable','Off')
set(handles.pushbuttonImageGraph,'Enable','On')
set(handles.uipanelGraph,'title',sprintf('Graph (%d nodes)',size(im.nodePos,1)));





function editNtrks_Callback(hObject, eventdata, handles)
global im

n2Run = str2double(get(hObject,'String'));

if n2Run > im.seed(im.curSeed).nTrksRan
    im.seed(im.curSeed).nTrks2Run = n2Run;
    im.seed(im.curSeed).seedLaunched = 0;
    set(handles.pushbuttonRunSeeds,'enable','on')
else
    set(handles.editNtrks,'string',num2str(im.seed(im.curSeed).nTrks2Run))
end









function editSeedConn_Callback(hObject, eventdata, handles)

set(handles.pushbuttonSeedConnectivity,'enable','on');




% --- Executes on button press in pushbuttonAutoSeed.
function pushbuttonAutoSeed_Callback(hObject, eventdata, handles)
global im

filenm = sprintf('%s.img',im.filenm);
fid = fopen( filenm, 'rb' );
n = fread( fid, 3, 'int');
nx = n(1); ny = n(2); nz = n(3);
I = fread( fid, nx*ny*nz, 'float' );
I = reshape(im.I, [ny nx nz] );
fclose(fid);

h = str2num(get(handles.editAutoSeedHxy,'string'));
hz = str2num(get(handles.editAutoSeedHz,'string'));
Ithresh = str2num(get(handles.editAutoSeedThresh,'string'));
nTrks = str2num(get(handles.editAutoSeedNtrks,'string'));

hWait = waitbar( 0, 'Auto Seeding...');

% blank out current seeds
for ii=1:im.nSeed
    if im.seed(ii).seedLaunched
        px = im.seed(ii).pos(1);
        py = im.seed(ii).pos(2);
        pz = im.seed(ii).pos(3);
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
    pos(nSeed,:) = [px py pz];
    
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
        
        ii = ii + 1;
    end
    jj = jj + 1;
end

% 
% for ii=1:length(iOrder);
%     im.seed(ii).pos = pos(iOrder(ii),:);
%     im.seed(ii).pos2 = [];
%     im.seed(ii).cxyz=[0 0 0];
%     im.seed(ii).bidirectional=1;
%     im.seed(ii).seedProcessed=0;
%     im.seed(ii).seedLaunched=0;
%     im.seed(ii).trkPos=[];
%     im.seed(ii).nTrks2Run=10;
%     im.seed(ii).nTrksRan=0;
% end

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

UpdateSeeds( handles )








% --- Executes on button press in pushbuttonGraphReset.
function pushbuttonGraphReset_Callback(hObject, eventdata, handles)
global im

for ii=1:length(im.seed)
    im.seed(ii).seedGraphed = 0;
end
im = rmfield(im,'nodePos');

set(handles.pushbuttonGraph,'enable','on')
set(handles.pushbuttonSeedConnectivity,'enable','on')
set(handles.pushbuttonImageGraph,'enable','off')
set(handles.pushbuttonRegraphNodes,'enable','off')
set(handles.pushbuttonImaris,'enable','off')

