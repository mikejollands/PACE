function varargout = PACE_mainscreen(varargin)
% PACE_MAINSCREEN MATLAB code for PACE_mainscreen.fig
%      PACE_MAINSCREEN, by itself, creates a new PACE_MAINSCREEN or raises the existing
%      singleton*.
%
%      H = PACE_MAINSCREEN returns the handle to a new PACE_MAINSCREEN or the handle to
%      the existing singleton*.
%
%      PACE_MAINSCREEN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PACE_MAINSCREEN.M with the given input arguments.
%
%      PACE_MAINSCREEN('Property','Value',...) creates a new PACE_MAINSCREEN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PACE_mainscreen_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PACE_mainscreen_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PACE_mainscreen

% Last Modified by GUIDE v2.5 10-Jan-2020 21:27:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PACE_mainscreen_OpeningFcn, ...
                   'gui_OutputFcn',  @PACE_mainscreen_OutputFcn, ...
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
end

% --- Executes just before PACE_mainscreen is made visible.
function PACE_mainscreen_OpeningFcn(hObject, eventdata, handles, varargin)
set( handles.figure1, 'Units', 'pixels' );
screenSize = get(0, 'ScreenSize');
position = get( handles.figure1, 'Position' );
position(3) = 380;
position(4) = 206;
position(1) = (screenSize(3)-position(3))/2; % put the GUI in the middle of the screen
position(2) = (screenSize(4)-position(4))/2; % put the GUI in the middle of the screen
set( handles.figure1,'Position', position );


% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PACE_mainscreen (see VARARGIN)

% Choose default command line output for PACE_mainscreen
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
end
% UIWAIT makes PACE_mainscreen wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PACE_mainscreen_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

% PACE
function pushbutton1_Callback(hObject, eventdata, handles)
PACE_1_5
end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
PACE_GD_1_2
end
