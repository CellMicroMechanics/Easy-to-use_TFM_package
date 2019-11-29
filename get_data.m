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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%% Setting initial parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
handles.data.strain_noise_file_name = 0;
handles.data.imagedir_name = 0;
handles.data.young = 10000;   %% Young's modulus
handles.data.poisson = 0.5;   %%  Poisson's ratio
handles.data.pix_durch_my = 0.067;  %% size of one pixel in micrometer
handles.data.zdepth = 0.0;    %% depth of focal plane below the surface
handles.data.displacement = []; 
handles.data.noise = []; 
handles.data.data_nonempty = [];

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


% --- Executes on pressing the button ''Continue''.
function Continue_Callback(hObject, eventdata, handles)

%%check if all parameters are provided correctly
if isnan(handles.data.young) || (handles.data.young< 0)
        errordlg(' ''Young modulus'' must be a positive number.','Error');
        return;
end     
if isnan(handles.data.pix_durch_my) || (handles.data.pix_durch_my<0)
    errordlg(' ''Micrometer per pixel'' must be a positive number.','Error');
    return;
end
if isnan(handles.data.poisson) || (handles.data.poisson<0) || (handles.data.poisson>0.5)
    errordlg(' ''Poisson ratio'' must be in the range (0-0.5).','Error');
    return;
end
if isnan(handles.data.zdepth) || (handles.data.zdepth<0)
    errordlg(' ''Depth of z-plane'' should be a positive number (the distance between gel surface and measurement plane).','Error');
    return;
end

%%check the existence for inpput_dat file and images folder 
if isempty(handles.data.strain_noise_file_name) || isequal(handles.data.strain_noise_file_name,0) 
        errordlg('Det data file has not been specified properly.','Error');
        return;
end

%%load file with displacements
try
    hilf  = load('-mat', handles.data.strain_noise_file_name);
catch
    errordlg('The displacement data must be provided as a Matlab file (.mat)','Error');
    return;    
end

%%check for the correct file structure and transfer data 
hilf_name = fieldnames(hilf);
if isempty(hilf) || all(strcmp(hilf_name, 'input_data') == 0) 
        errordlg('The input file must contain a structure with the name ''input_data''.','Error');
        return;
else
    hilf_name = fieldnames(hilf.input_data);
    if isempty(hilf) || all(strcmp(hilf_name, 'displacement') == 0) || isempty(hilf.input_data.displacement) || ~isstruct(hilf.input_data.displacement)
            errordlg('The input file must contain a non-empty structure named ''input_data.displacement()''.','Error');
            return;
    end    
    hilf_name = fieldnames(hilf.input_data.displacement);
    if all(strcmp(hilf_name, 'pos') == 0) || all(strcmp(hilf_name, 'vec') == 0)
        errordlg('The provided data structure does not contain the fields ''.pos()'' and ''.vec()''.','Error');
        return;
    else
        handles.data.displacement = hilf.input_data.displacement;
    end

    if isfield(hilf.input_data,'noise') == 0 || ~isstruct(hilf.input_data.noise) || isempty(hilf.input_data.noise) 
        handles.data.noise  = []; 
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

%identify the data sets that contain displacement data
%This information is stored in a variable called 'data_nonempty'
framenumber = length(handles.data.displacement);
data_nonempty = false(framenumber,1);
for i= 1: framenumber
    size_pos = size(handles.data.displacement(i).pos);
    size_vec = size(handles.data.displacement(i).vec);
    
    if size_pos(1) > 1 && size_pos(2) >= 2 && size_vec(1)>1  && size_vec(2) >= 2 
       data_nonempty(i,1) = true;
    else 
       data_nonempty(i,1) = false ;  
    end 
end
handles.data.data_nonempty = [data_nonempty  (1: framenumber)']; 

%check all displacement datasets to see if there is any data
data_nonempty =  handles.data.data_nonempty;
first_index = find(data_nonempty(:,1),1);
if isempty(first_index)
    errordlg('The dataset only contains only empty or too small dispacement fields.','Error');
    return;
end


%select the next menu that will open now
blah = get(handles.selectmethod, 'SelectedObject');
method =get(blah,'String');
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
handles.data.young = written*1000;
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
handles.data.poisson = written;
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
handles.data.pix_durch_my = written;
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
handles.data.zdepth = written;
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
%read the name of the file containing input data
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
%Read the file names for the images
function imagedir_browse_Callback(hObject, eventdata, handles)

handles.data.imagedir_name = uigetdir('','Location of cell images');

%%%%%%%%%%%%%%%%%%% change
% dir_struct = vertcat(dir(fullfile(handles.data.imagedir_name,'*.tif*')),dir(fullfile(handles.data.imagedir_name,'*.jpg*')));
%%%%%%%%%%%%%%%%%%%

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
