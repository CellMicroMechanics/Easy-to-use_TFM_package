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

function varargout = get_data(varargin)
%GET_DATA MATLAB code file for get_data.fig
%      GET_DATA, by itself, creates a new GET_DATA or raises the existing
%      singleton*.
%
%      H = GET_DATA returns the handle to a new GET_DATA or the handle to
%      the existing singleton*.
%
%      GET_DATA('Property','Value',...) creates a new GET_DATA using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to get_data_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      GET_DATA('CALLBACK') and GET_DATA('CALLBACK',hObject,...) call the
%      local function named CALLBACK in GET_DATA.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help get_data

% Last Modified by GUIDE v2.5 14-Jan-2019 09:50:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @get_data_OpeningFcn, ...
                   'gui_OutputFcn',  @get_data_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before get_data is made visible.
function get_data_OpeningFcn(hObject, eventdata, handles, varargin)

handles.data.strain_noise_file_name = 0;
handles.data.imagedir_name = 0;

handles.data.young = 10000;
handles.data.poisson = 0.5;
handles.data.pix_durch_my = 0.067;
handles.data.zdepth = 0.0;
handles.data.displacement = []; 
handles.data.noise = []; 

set(handles.figure1, 'units', 'normalized');

% Choose default command line output for get_data
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes get_data wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = get_data_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Continue.
function Continue_Callback(hObject, eventdata, handles)

blah = get(handles.selectmethod, 'SelectedObject');
method =get(blah,'String');

if isempty(handles.data.strain_noise_file_name) || isequal(handles.data.strain_noise_file_name,0) || ...
         isempty(handles.data.imagedir_name) || isequal(handles.data.imagedir_name,0)
        errordlg('Some of the file/directory names have not been specified properly.','Error');
        return;
end

hilf  = load('-mat', handles.data.strain_noise_file_name);
hilf_name = fieldnames(hilf);

if isempty(hilf) || all(strcmp(hilf_name, 'input_data') == 0) 
        errordlg('Specified file with displacement data does not contain data or is empty.','Error');
        return;
else
    hilf_name = fieldnames(hilf.input_data.displacement);
    if all(strcmp(hilf_name, 'pos') == 0) || all(strcmp(hilf_name, 'vec') == 0)
        errordlg('Specified file with displacement data does not contain the fields .pos/.vec.','Error');
        return;
    else
        handles.data.displacement = hilf.input_data.displacement;
    end
end
 
 hilf_n  = load('-mat', handles.data.strain_noise_file_name);
 hilf_name_n = fieldnames(hilf_n);

if isempty(hilf_n) || all(strcmp(hilf_name_n, 'input_data') == 0) 
        errordlg('Specified file with noise data does not contain data or is empty.','Error');
        return;
else
 
     if isfield(hilf.input_data,'noise') == 0
        handles.data.noise  = 1; 
     else
 
     hilf_name_n = fieldnames(hilf.input_data.noise);
       if all(strcmp(hilf_name_n, 'pos') == 0) || all(strcmp(hilf_name_n, 'vec') == 0)
         errordlg('Specified file with noise data does not contain the fields .pos/.vec.','Error');
       return;
       else
       handles.data.noise = hilf.input_data.noise;
       end
     end  
end
 
switch method
    case 'Regularization'
          Regularization(handles.data);
    case 'Bayesian regularization'
          Bayesian_regularization(handles.data);   
end

delete(handles.figure1);

return;


% --- Executes on button press in abort.
function abort_Callback(hObject, eventdata, handles)
delete(handles.figure1);
return;


% --- Executes when selected object is changed in selectmethod.
function selectmethod_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in selectmethod 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in radiobutton9.
function radiobutton9_Callback(hObject, eventdata, handles)
% hObject    handle to the selected object in selectmethod 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function young_win_Callback(hObject, eventdata, handles)
written = str2double(get(hObject,'String'));
if ~isnan(written) && written > 0
    handles.data.young = written*1000;
else
    errordlg('The young modulus must be given as a positive number.','Error');
    set(handles.young_win,'String', num2str(handles.data.young/1000));
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function young_win_CreateFcn(hObject, eventdata, handles)
% hObject    handle to young_win (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function poisson_win_Callback(hObject, eventdata, handles)
written = str2double(get(hObject,'String'));
if ~isnan(written) && written >= 0 && written <= 0.5
    handles.data.poisson = written;
else
    errordlg('The poisson ratio must be between 0 and 0.5.','Error');
    set(handles.poisson_win,'String', num2str(handles.data.poisson));
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function poisson_win_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poisson_win (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pix_durch_my_win_Callback(hObject, eventdata, handles)
written = str2double(get(hObject,'String'));
if ~isnan(written) && written > 0
    handles.data.pix_durch_my = written;
else
    errordlg('Please enter a positive number here.','Error');
    set(handles.pix_durch_my_win,'String', num2str(handles.data.pix_durch_my));
end
guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of pix_durch_my_win as text
%        str2double(get(hObject,'String')) returns contents of pix_durch_my_win as a double


% --- Executes during object creation, after setting all properties.
function pix_durch_my_win_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
 set(hObject,'BackgroundColor','white');
end



function depth_in_substrate_win_Callback(hObject, eventdata, handles)
written = str2double(get(hObject,'String'));
if ~isnan(written) && written > 0
    handles.data.zdepth = written;
else
    errordlg('Please enter a positive number or zero here.','Error');
    set(handles.depth_in_substrate_win,'String', num2str(handles.data.zdepth));
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function depth_in_substrate_win_CreateFcn(hObject, eventdata, handles)
% hObject    handle to depth_in_substrate_win (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in strain_noise_file_browse.
function strain_noise_file_browse_Callback(hObject, eventdata, handles)

[filename, pathname] = uigetfile('*.mat', 'Select a .mat file');
if ~isequal(filename,0)
    handles.data.strain_noise_file_name = fullfile(pathname, filename);
    if size(handles.data.strain_noise_file_name,2) > 55
        disp_string = ['... ',handles.data.strain_noise_file_name(end-55:end)];
    else
        disp_string = handles.data.strain_noise_file_name;
    end
        set(handles.strainfile_win,'string',disp_string);
end


guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function strainfile_win_CreateFcn(hObject, eventdata, handles)
% hObject    handle to strainfile_win (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function noisefile_win_CreateFcn(hObject, eventdata, handles)
% hObject    handle to strainfile_win (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on button press in imagedir_browse.
function imagedir_browse_Callback(hObject, eventdata, handles)

handles.data.imagedir_name = uigetdir('','Location of cell images');

dir_struct = vertcat(dir(fullfile(handles.data.imagedir_name,'*.tif*')),dir(fullfile(handles.data.imagedir_name,'*.jpg*')));
  if isempty(dir_struct) || all(strcmp(dir_struct.name))
        errordlg('Specified file with cell images do not contain the fields .tif/.jpg.','Error');
        return;
  end 
     
if ~isequal(handles.data.imagedir_name,0)
    if size(handles.data.imagedir_name,2) > 55
        disp_string = ['... ',handles.data.imagedir_name(end-55:end)];
    else
        disp_string = handles.data.imagedir_name;
    end
     set(handles.imagedir_win,'String',disp_string);
end
guidata(hObject, handles);

% --- Executes on button press in radiobutton10.
function radiobutton10_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
