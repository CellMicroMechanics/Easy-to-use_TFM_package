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

function varargout = Select_noise(varargin)
% SELECT_NOISE MATLAB code for Select_noise.fig
%      SELECT_NOISE, by itself, creates a new SELECT_NOISE or raises the existing
%      singleton*.
%
%      H = SELECT_NOISE returns the handle to a new SELECT_NOISE or the handle to
%      the existing singleton*.
%
%      SELECT_NOISE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECT_NOISE.M with the given input arguments.
%
%      SELECT_NOISE('Property','Value',...) creates a new SELECT_NOISE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Select_noise_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Select_noise_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Select_noise

% Last Modified by GUIDE v2.5 23-Dec-2018 15:57:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Select_noise_OpeningFcn, ...
                   'gui_OutputFcn',  @Select_noise_OutputFcn, ...
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


% --- Executes just before Select_noise is made visible.
function Select_noise_OpeningFcn(hObject, eventdata, handles, varargin)
handles.data = varargin{1};

frame = getappdata(0,'frame_1');
bild_datei_index = getappdata(0,'bild_datei_index_1');
bild_dateien = getappdata(0,'bild_dateien_1');

%If we do not get a cell array, we need to make one (eg if only one entry)
if ~isempty(bild_dateien) && ~iscellstr(bild_dateien)
    bild_dateien =cellstr(bild_dateien);
end
no_image = false;
if ~isempty(bild_dateien) && (length(bild_dateien) >= bild_datei_index) && ~strcmp(bild_dateien{bild_datei_index},'No image available')
    bild_datei = bild_dateien{bild_datei_index};
    try
        handles.data.bild = imread(fullfile(handles.data.imagedir_name, bild_datei));
    catch
        no_image = true;
        handles.data.bild = [];
    end
else
    no_image = true;
end

axes(handles.axes1);
cla; axis equal, colormap(gca,gray); hold on;
if ~no_image
    imagesc(handles.data.bild);
end
quiver(handles.data.displacement(frame).pos(:,1),handles.data.displacement(frame).pos(:,2),...
    handles.data.displacement(frame).vec(:,1),handles.data.displacement(frame).vec(:,2),'r');
set(gca, 'DataAspectRatio', [1,1,50],'YDir','reverse','XTick',[],'YTick',[])
hold off;

handles.output = hObject;

guidata(hObject, handles);

% UIWAIT makes Select_noise wait for user response (see UIRESUME)
% uiwait(handles.select_noise);


% --- Outputs from this function are returned to the command line.
function varargout = Select_noise_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


% --- Executes on selection change in preview_image.
function preview_image_Callback(hObject, eventdata, handles)


% Hints: contents = cellstr(get(hObject,'String')) returns preview_image contents as cell array
%        contents{get(hObject,'Value')} returns selected item from preview_image


% --- Executes during object creation, after setting all properties.
function preview_image_CreateFcn(hObject, eventdata, handles)
% hObject    handle to preview_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
  
 
