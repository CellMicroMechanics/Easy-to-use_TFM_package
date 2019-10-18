% Copyright (C) 2010 - 2019, Sabass Lab
%
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option) 
% any later version. This program is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
% Public License for more details. You should have received a copy of the 
% GNU General Public License along with this program.
% If not, see <http://www.gnu.org/licenses/>.

function varargout = calculation_done(varargin)
% CALCULATION_DONE MATLAB code for calculation_done.fig
%      CALCULATION_DONE, by itself, creates a new CALCULATION_DONE or raises the existing
%      singleton*.
%
%      H = CALCULATION_DONE returns the handle to a new CALCULATION_DONE or the handle to
%      the existing singleton*.
%
%      CALCULATION_DONE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALCULATION_DONE.M with the given input arguments.
%
%      CALCULATION_DONE('Property','Value',...) creates a new CALCULATION_DONE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before calculation_done_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to calculation_done_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help calculation_done

% Last Modified by GUIDE v2.5 18-Oct-2019 11:12:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @calculation_done_OpeningFcn, ...
                   'gui_OutputFcn',  @calculation_done_OutputFcn, ...
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


% --- Executes just before calculation_done is made visible.
function calculation_done_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to calculation_done (see VARARGIN)

% Choose default command line output for calculation_done
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes calculation_done wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = calculation_done_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in OK.
function OK_Callback(hObject, eventdata, handles)
delete (handles.figure1);
