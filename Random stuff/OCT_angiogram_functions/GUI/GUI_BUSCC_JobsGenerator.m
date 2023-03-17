function varargout = GUI_BUSCC_JobsGenerator(varargin)
% GUI_BUSCC_JOBSGENERATOR MATLAB code for GUI_BUSCC_JobsGenerator.fig
%      GUI_BUSCC_JOBSGENERATOR, by itself, creates a new GUI_BUSCC_JOBSGENERATOR or raises the existing
%      singleton*.10w GUI_BUSCC_JOBSGENERATOR or raises the existing
%      singleton*.10
%w GUI_BUSCC_JOBSGENERATOR or raises the existing
%      singleton*.10
%
%      H = GUI_BUSCC_JOBSGENERATOR returns the handle to a new GUI_BUSCC_JOBSGENERATOR or the handle to
%      the existing singleton*.qs
%
%      GUI_BUSCC_JOBSGENERATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_BUSCC_JOBSGENERATOR.M with the given input arguments.
%
%      GUI_BUSCC_JOBSGENERATOR('Property','Value',...) creates a new GUI_BUSCC_JOBSGENERATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_BUSCC_JobsGenerator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_BUSCC_JobsGenerator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singlet1on)".
%
%See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_BUSCC_JobsGenerator

% Last Modified by GUIDE v2.5 09-Aug-2019 11:17:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_BUSCC_JobsGenerator_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_BUSCC_JobsGenerator_OutputFcn, ...
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


% --- Executes just before GUI_BUSCC_JobsGenerator is made visible.
function GUI_BUSCC_JobsGenerator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_BUSCC_JobsGenerator (see VARARGIN)

% Choose default command line output for GUI_BUSCC_JobsGenerator
handles.output = hObject;

handles.defaultpathDat='/projectnb/npboctiv/ns/';
handles.defaultpathGG='/projectnb/npboctiv/ns/';
handles.defaultpathk='/projectnb/npboctiv/ns/';
%% add MATLAB functions' path
% addpath('D:\OCT Imaging\Data Process CODE\CODE-BU\Functions') % Path on JTOPTICS
addpath('/projectnb/npboctiv/ns/Jianbo/BU-SCC/Functions') % subFunction Path on SCC server
addpath('/projectnb/npboctiv/ns/Jianbo/BU-SCC/GUI') % GUI path on SCC server
handles.jobfolder=uigetdir('/projectnb/npboctiv/ns/Jianbo/BU-SCC/');
handles.jobfolder=[handles.jobfolder, '/'];
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_BUSCC_JobsGenerator wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_BUSCC_JobsGenerator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btnDat2GG.
function btnDat2GG_Callback(hObject, eventdata, handles)
% hObject    handle to btnDat2GG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% add MATLAB functions' path
% addpath('D:\PJ DLS-OCT\CODE\Functions') % Path on JTOPTICS
addpath('/projectnb/npboctiv/ns/Jianbo/BU-SCC/Functions') % subFunction Path on SCC server
addpath('/projectnb/npboctiv/ns/Jianbo/BU-SCC/GUI') % GUI path on SCC server

[filename0,datapath0]=uigetfile(fullfile(handles.defaultpathDat,'*.dat'),'select file');
handles.defaultpathDat=datapath0;
guidata(hObject, handles);

[dim, fNameBase,fIndex]=GetNameInfoRaw(filename0);
%% Number of files to process %%%%%%%%%%%%%%%%%%
prompt={'Start file index','Number of files to Load and generate jobs'};
inputNum=inputdlg(prompt,'', 1,{num2str(fIndex),'1'});
Start_file=str2num(inputNum{1});  % number of segments
N_file=str2num(inputNum{2});  % number of segments

JOBS_Dat2GG(filename0,datapath0, Start_file,N_file,handles.jobfolder);
msgbox('Dat2GG job generated and saved!')

% --- Executes on button press in btnGG2VDR.
function btnGG2VDR_Callback(hObject, eventdata, handles)
% hObject    handle to btnGG2VDR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Generate job file for launchpad calculating GG2VDR
%% add MATLAB functions' path
% addpath('D:\PJ DLS-OCT\CODE\Functions') % Path on JTOPTICS
addpath('/projectnb/npboctiv/ns/Jianbo/BU-SCC/Functions') % subFunction Path on SCC server
addpath('/projectnb/npboctiv/ns/Jianbo/BU-SCC/GUI') % GUI path on SCC server

%%%%% load file, select any one of the same segmentation files %%%%%%%%%%
[filename0,datapath0]=uigetfile(handles.defaultpathGG);
handles.defaultpathGG=datapath0;
guidata(hObject, handles);

pathInfo=strsplit(datapath0,'/');
fileInfo=pathInfo{end-2};
fileParameter=strsplit(fileInfo,'-');
fIndex=str2num(fileParameter{end});

%% Number of files to process %%%%%%%%%%%%%%%%%%
prompt={'Start file index','Number of files to Load and generate jobs'};
inputNum=inputdlg(prompt,'', 1,{num2str(fIndex),'1'});
Start_file=str2num(inputNum{1});  % number of segments
N_file=str2num(inputNum{2});  % number of segments

JOBS_GG2VDR(filename0,datapath0, Start_file,N_file,handles.jobfolder);
msgbox('GG2VDR job generated and saved!')


% --- Executes on button press in Dat2AG.
function Dat2AG_Callback(hObject, eventdata, handles)
% hObject    handle to Dat2AG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% add MATLAB functions' path
% addpath('D:\PJ DLS-OCT\CODE\Functions') % Path on JTOPTICS
addpath('/projectnb/npboctiv/ns/Jianbo/BU-SCC/Functions') % subFunction Path on SCC server
addpath('/projectnb/npboctiv/ns/Jianbo/BU-SCC/GUI') % GUI path on SCC server

[filename0,datapath0]=uigetfile(fullfile(handles.defaultpathDat,'*.dat'),'select file');
handles.defaultpathDat=datapath0;
guidata(hObject, handles);

[dim, fNameBase,fIndex]=GetNameInfoRaw(filename0);

%% Number of files to process %%%%%%%%%%%%%%%%%%
prompt={'Start file index','Number of files to Load and generate jobs','intDk(-0.19 for BU Anes, -0.235 for BU Awake)'};
inputNum=inputdlg(prompt,'', 1,{num2str(fIndex),'10','-0.19'});
Start_file=str2num(inputNum{1});  % number of segments
N_file=str2num(inputNum{2});  % number of segments
intDk=str2num(inputNum{3});

JOBS_Dat2AG(filename0,datapath0, Start_file,N_file,intDk, handles.jobfolder);
msgbox('Dat2AG job generated and saved!')   


% --- Executes on button press in btnReset.
function btnReset_Callback(hObject, eventdata, handles)
% hObject    handle to btnReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.defaultpathDat='/projectnb/npboctiv/ns/';
handles.defaultpathGG='/projectnb/npboctiv/ns/';
handles.defaultpathk='/projectnb/npboctiv/ns/';
% Update handles structure
guidata(hObject, handles);
