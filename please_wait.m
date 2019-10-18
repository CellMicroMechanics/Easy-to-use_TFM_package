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

function varargout = please_wait(varargin)
% PLEASE_WAIT MATLAB code for please_wait.fig
%      PLEASE_WAIT, by itself, creates a new PLEASE_WAIT or raises the existing
%      singleton*.
%
%      H = PLEASE_WAIT returns the handle to a new PLEASE_WAIT or the handle to
%      the existing singleton*.
%
%      PLEASE_WAIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLEASE_WAIT.M with the given input arguments.
%
%      PLEASE_WAIT('Property','Value',...) creates a new PLEASE_WAIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before please_wait_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to please_wait_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help please_wait

% Last Modified by GUIDE v2.5 08-Jan-2019 11:00:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @please_wait_OpeningFcn, ...
                   'gui_OutputFcn',  @please_wait_OutputFcn, ...
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


% --- Executes just before please_wait is made visible.
function please_wait_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to please_wait (see VARARGIN)

% Choose default command line output for please_wait
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes please_wait wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = please_wait_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
