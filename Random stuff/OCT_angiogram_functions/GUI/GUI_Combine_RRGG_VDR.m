function varargout = GUI_Combine_RRGG_VDR(varargin)
% GUI_COMBINE_RRGG_VDR MATLAB code for GUI_Combine_RRGG_VDR.fig
%      GUI_COMBINE_RRGG_VDR, by itself, creates a new GUI_COMBINE_RRGG_VDR or raises the existing
%      singleton*.
%
%      H = GUI_COMBINE_RRGG_VDR returns the handle to a new GUI_COMBINE_RRGG_VDR or the handle to
%      the existing singleton*.
%
%      GUI_COMBINE_RRGG_VDR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_COMBINE_RRGG_VDR.M with the given input arguments.
%
%      GUI_COMBINE_RRGG_VDR('Property','Value',...) creates a new GUI_COMBINE_RRGG_VDR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_Combine_RRGG_VDR_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_Combine_RRGG_VDR_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_Combine_RRGG_VDR

% Last Modified by GUIDE v2.5 09-Aug-2019 13:50:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Combine_RRGG_VDR_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Combine_RRGG_VDR_OutputFcn, ...
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


% --- Executes just before GUI_Combine_RRGG_VDR is made visible.
function GUI_Combine_RRGG_VDR_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_Combine_RRGG_VDR (see VARARGIN)

% Choose default command line output for GUI_Combine_RRGG_VDR
handles.output = hObject;

handles.defaultpathVDR='/projectnb/npboctiv/ns/';
handles.defaultpathGG='/projectnb/npboctiv/ns/';
handles.defaultpathVz='/projectnb/npboctiv/ns/';

%% add MATLAB functions' path
% addpath('D:\OCT Imaging\Data Process CODE\CODE-BU\Functions') % Path on JTOPTICS
addpath('/projectnb/npboctiv/ns/Jianbo/BU-SCC/Functions') % subFunction Path on SCC server
addpath('/projectnb/npboctiv/ns/Jianbo/BU-SCC/GUI') % GUI path on SCC server
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_Combine_RRGG_VDR wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_Combine_RRGG_VDR_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btnRRGG.
function btnRRGG_Callback(hObject, eventdata, handles)
% hObject    handle to btnRRGG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% add MATLAB functions' path
% addpath('D:\OCT Imaging\Data Process CODE\CODE-BU\Functions') % Path on JTOPTICS
addpath('/projectnb/npboctiv/ns/Jianbo/BU-SCC/Functions') % subFunction Path on SCC server
addpath('/projectnb/npboctiv/ns/Jianbo/BU-SCC/GUI') % GUI path on SCC server

%%%%% load file, select any one of the same segmentation files %%%%%%%%%%
[filename0,datapath0]=uigetfile(handles.defaultpathGG);
handles.defaultpathGG=datapath0;
guidata(hObject, handles);
pathInfo=strsplit(datapath0,'/');
fileInfo=pathInfo{end-1};
fileParameter=strsplit(fileInfo,'-');
fIndex=str2num(fileParameter{end});

%% Number of files to process %%%%%%%%%%%%%%%%%%
prompt={'Start file index','Number of files to Load and generate jobs'};
inputNum=inputdlg(prompt,'', 1,{num2str(fIndex),'1'});
Start_file=str2num(inputNum{1});  % number of segments
N_file=str2num(inputNum{2});  % number of segments

%%
Combine_RRGG(filename0,datapath0, Start_file,N_file);
msgbox('RRGG combined!')


% --- Executes on button press in btnVDR.
function btnVDR_Callback(hObject, eventdata, handles)
% hObject    handle to btnVDR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% add MATLAB functions' path
% addpath('D:\OCT Imaging\Data Process CODE\CODE-BU\Functions') % Path on JTOPTICS
addpath('/projectnb/npboctiv/ns/Jianbo/BU-SCC/Functions') % subFunction Path on SCC server
addpath('/projectnb/npboctiv/ns/Jianbo/BU-SCC/GUI') % GUI path on SCC server

[filename,datapath0]=uigetfile(handles.defaultpathVDR);
handles.defaultpathVDR=datapath0;
guidata(hObject, handles);
pathparts=strsplit(datapath0,'/');
folder1=pathparts{end-3};
if folder1(end-1) =='-'
    fileindex=str2num(folder1(end));
elseif folder1(end-2) =='-'
    fileindex=str2num(folder1(end-1:end));
end

%% Number of files to process %%%%%%%%%%%%%%%%%%
prompt={'Start file index','Number of files to Load and generate jobs'};
inputNum=inputdlg(prompt,'', 1,{num2str(fileindex),'1'});
Start_file=str2num(inputNum{1});  % number of segments
N_file=str2num(inputNum{2});  % number of segments

Combine_VDR(filename,datapath0, Start_file,N_file);
msgbox('VDR combined!')

% --- Executes on button press in AVG_AG.
function AVG_AG_Callback(hObject, eventdata, handles)
% hObject    handle to AVG_AG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% add MATLAB functions' path
% addpath('D:\OCT Imaging\Data Process CODE\CODE-BU\Functions') % Path on JTOPTICS
addpath('/projectnb/npboctiv/ns/Jianbo/BU-SCC/Functions') % subFunction Path on SCC server
addpath('/projectnb/npboctiv/ns/Jianbo/BU-SCC/GUI') % GUI path on SCC server

[filename,datapath0]=uigetfile(handles.defaultpathVDR);
handles.defaultpathVDR=datapath0;
guidata(hObject, handles);
fileparts=strsplit(filename,'-');
fileinfo=fileparts{end};
fileindex=fileinfo(1:end-4);

%% Number of files to process %%%%%%%%%%%%%%%%%%
prompt={'Start file index','Number of files to Load and generate jobs'};
inputNum=inputdlg(prompt,'', 1,{fileindex,'10'});
Start_file=str2num(inputNum{1});  % number of segments
N_file=str2num(inputNum{2});  % number of segments

for ifile=Start_file:Start_file+N_file-1
    if exist([datapath0,filename(1:end-5),num2str(ifile),'.mat'],'file')==2
        if ifile==Start_file
            AG=abs(LoadMAT(datapath0,[filename(1:end-5),num2str(ifile)]));
        else
            AG=AG+abs(LoadMAT(datapath0,[filename(1:end-5),num2str(ifile)]));
        end
        display(num2str(ifile));
    else
        display([num2str(ifile),' skiped!']);
    end
end
AG=AG/N_file;
save([datapath0,'/',filename(1:end-6),'-AVG',num2str(N_file),'.mat'],'-v7.3','AG');
figure,imshow(log(squeeze(max(AG(:,:,:),[],1))),[])
msgbox('Angiograms are averaged and saved!')



% --- Executes on button press in CombVz.
function CombVz_Callback(hObject, eventdata, handles)
% hObject    handle to CombVz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% add MATLAB functions' path
% addpath('D:\OCT Imaging\Data Process CODE\CODE-BU\Functions') % Path on JTOPTICS
addpath('/projectnb/npboctiv/ns/Jianbo/BU-SCC/Functions') % subFunction Path on SCC server
addpath('/projectnb/npboctiv/ns/Jianbo/BU-SCC/GUI') % GUI path on SCC server

%%%%% load file, select any one of the same segmentation files %%%%%%%%%%
[filename0,datapath0]=uigetfile(handles.defaultpathVz);
handles.defaultpathVz=datapath0;
guidata(hObject, handles);
if datapath0(end-3) =='-'
    fileindex=str2num(datapath0(end-2:end-1));
else
    fileindex=str2num(datapath0(end-1));
end

%% Number of files to process %%%%%%%%%%%%%%%%%%
prompt={'Start file index','Number of files to Load and generate jobs'};
inputNum=inputdlg(prompt,'', 1,{num2str(fileindex),'1'});
Start_file=str2num(inputNum{1});  % number of segments
N_file=str2num(inputNum{2});  % number of segments
%%
% Combine_PHVz(filename0,datapath0, Start_file,N_file);
Combine_prVzNg1OCTA(filename0,datapath0, Start_file,N_file);
msgbox('Vz combined!')

% --- Executes on button press in btnReset.
function btnReset_Callback(hObject, eventdata, handles)
% hObject    handle to btnReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
addpath('/autofs/space/blinky_003/users/Jianbo/CODE/GUI') % Path on server
handles.defaultpathVDR='/projectnb/npboctiv/ns/';
handles.defaultpathGG='/projectnb/npboctiv/ns/';
handles.defaultpathVz='/projectnb/npboctiv/ns/';

guidata(hObject, handles);
