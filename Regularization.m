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

function varargout = Regularization(varargin)
%REGULARIZATION MATLAB code file for Regularization.fig
%      REGULARIZATION, by itself, creates a new REGULARIZATION or raises the existing
%      singleton*.
%
%      H = REGULARIZATION returns the handle to a new REGULARIZATION or the handle to
%      the existing singleton*.
%
%      REGULARIZATION('Property','Value',...) creates a new REGULARIZATION using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to Regularization_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      REGULARIZATION('CALLBACK') and REGULARIZATION('CALLBACK',hObject,...) call the
%      local function named CALLBACK in REGULARIZATION.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Regularization

% Last Modified by GUIDE v2.5 16-Aug-2019 08:25:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Regularization_OpeningFcn, ...
                   'gui_OutputFcn',  @Regularization_OutputFcn, ...
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


% --- Executes just before Regularization is made visible.
function Regularization_OpeningFcn(hObject, eventdata, handles, varargin)
handles.data = varargin{1};

dir_struct = vertcat(dir(fullfile(handles.data.imagedir_name,'*.tif*')),dir(fullfile(handles.data.imagedir_name,'*.jpg*')));
[sorted_names,sorted_index] = sortrows({dir_struct.name}');
set(handles.preview_image,'String',sorted_names,'Value',1)
for i = 1:length(handles.data.displacement)
    hilf_cell{i,1} = num2str(i);
end
set(handles.preview_frame,'String',hilf_cell,'Value',1);
set(handles.figure1, 'units', 'normalized');

max_eck(1:2) = [max(handles.data.displacement(1).pos(:,1)), max(handles.data.displacement(1).pos(:,2))];
min_eck(1:2) = [min(handles.data.displacement(1).pos(:,1)), min(handles.data.displacement(1).pos(:,2))];
  
handles.data.meshsize = round(sqrt((max_eck(1)-min_eck(1))*(max_eck(2)-min_eck(2))/size(handles.data.displacement(1).pos,1)));
set(handles.meshsize_win,'String',num2str(handles.data.meshsize));

handles.data.regparam = 20;
set(handles.regparam_win,'String',num2str(handles.data.regparam,'%f'));

handles.data.pos = 0;
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Regularization wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Regularization_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on mouse press over axes background.
function axes2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in abort.
function abort_Callback(hObject, eventdata, handles)
delete(handles.figure1);
return;



% --- Executes on button press in analyze.
function analyze_Callback(hObject, eventdata, handles)

grid_mat = [];

haha = waitbar(0,'Please wait while data is being assembled..','WindowStyle','modal');
framenumber = length(handles.data.displacement);

%%% scale the Young modulus to 1
ne = handles.data.young;
young= 1;

for frame = 1:framenumber
    if isempty(handles.data.displacement(frame).pos) || isempty(handles.data.displacement(frame).vec)||...
       length(handles.data.displacement(frame).pos) < 2
         disp(['There dataset of frame ', num2str(frame),' was not valid. Skipping this frame!']);
        continue;
    end
   waitbar(frame/framenumber);
   [grid_mat,u, i_max,j_max] = interp_vec2grid(handles.data.displacement(frame).pos+handles.data.displacement(frame).vec,...
                                handles.data.displacement(frame).vec,handles.data.meshsize, grid_mat);

   [TFM_results(frame).pos,TFM_results(frame).traction,TFM_results(frame).traction_magnitude,f] = ...
                         reg_fourier_TFM(grid_mat,u,young,handles.data.poisson, handles.data.pix_durch_my,...
                         handles.data.zdepth, handles.data.meshsize, i_max,j_max, handles.data.regparam);
   TFM_results(frame).traction=ne*TFM_results(frame).traction;
   TFM_results(frame).traction_magnitude=ne*TFM_results(frame).traction_magnitude;
   
   bnd = 6;
   f = f.*ne;
   TFM_results(frame).energy = 1/2*sum(sum(u(bnd:end-bnd+1,bnd:end-bnd+1,1).*f(bnd:end-bnd+1,bnd:end-bnd+1,1) +...
         u(bnd:end-bnd+1,bnd:end-bnd+1,2).*f(bnd:end-bnd+1,bnd:end-bnd+1,2)))*...
         (handles.data.meshsize)^2*handles.data.pix_durch_my^3/10^6; 
  
    pu(:,1) = reshape(u(:,:,1),i_max*j_max,1);
    pu(:,2) = reshape(u(:,:,2),i_max*j_max,1);
    TFM_results(frame).displacement= pu;  
end

TFM_settings.poisson = handles.data.poisson;
TFM_settings.young = handles.data.young;
TFM_settings.micrometer_per_pix = handles.data.pix_durch_my;
TFM_settings.regularization_parameter = handles.data.regparam;
TFM_settings.meshsize = handles.data.meshsize;
TFM_settings.zdepth = handles.data.zdepth;
TFM_settings.i_max = i_max;
TFM_settings.j_max = j_max;
close(haha);


 [filepath,name,ext] = fileparts(handles.data.strain_noise_file_name); 
 handles.data.targetdir_name = filepath;

savefile_name = fullfile(handles.data.targetdir_name,['Reg-FTTC_results_',datestr(now, 'dd-mm-yy'),'.mat']);
if exist(savefile_name)
    button = questdlg('The file exists already. Overwrite?','Error','Yes');
    if strcmpi(button,'No') || strcmpi(button,'')
            return;
    end
end
save(savefile_name,'TFM_results','TFM_settings','-mat');

calculation_done;

return;


% --- Executes on selection change in preview_image.
function preview_image_Callback(hObject, eventdata, handles)
% hObject    handle to preview_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function preview_image_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in preview_frame.
function preview_frame_Callback(hObject, eventdata, handles)
% hObject    handle to preview_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns preview_frame contents as cell array
%        contents{get(hObject,'Value')} returns selected item from preview_frame


% --- Executes during object creation, after setting all properties.
function preview_frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to preview_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in preview.
function preview_Callback(hObject, eventdata, handles)

frame = get(handles.preview_frame,'value');
if isempty(handles.data.displacement(frame).pos) || isempty(handles.data.displacement(frame).vec)||...
length(handles.data.displacement(frame).pos) < 2
    errordlg('The selected dataset is empty or not valid.','Error');
    return;
end

bild_datei_index = get(handles.preview_image,'value');
bild_dateien = get(handles.preview_image,'string');
bild_datei = bild_dateien{bild_datei_index};

%%% scale the Young modulus to 1
ne = handles.data.young;
young= 1;

[grid_mat,u, i_max,j_max] = interp_vec2grid(handles.data.displacement(frame).pos+handles.data.displacement(frame).vec,...
                            handles.data.displacement(frame).vec,handles.data.meshsize);

[handles.data.pos,traction,traction_magnitude,f] = ...
                         reg_fourier_TFM(grid_mat,u,young,handles.data.poisson, handles.data.pix_durch_my,...
                         handles.data.zdepth, handles.data.meshsize, i_max,j_max, handles.data.regparam);
 handles.data.bild = imread(fullfile(handles.data.imagedir_name, bild_datei));
handles.data.force=ne.*traction;
handles.data.f_mat=ne.*f;

axes(handles.axes1);
cla; axis equal, hold on; colormap gray, imagesc(handles.data.bild);

hilf = get(handles.show_vectors,'Value');
if hilf
    quiver(handles.data.pos(:,1),handles.data.pos(:,2),handles.data.force(:,1),handles.data.force(:,2),2,'r');
end
set(gca, 'DataAspectRatio', [1,1,50],'YDir','reverse','XTick',[],'YTick',[])
hold off;

axes(handles.axes2);
cla; hold on; colormap jet;
fnorm = (handles.data.f_mat(:,:,2).^2 + handles.data.f_mat(:,:,1).^2).^0.5;
surf(grid_mat(:,:,1), grid_mat(:,:,2),fnorm),view(0,90),shading interp, axis equal;
set(gca, 'DataAspectRatio', [1,1,50],'YDir','reverse','XTick',[],'YTick',[]),hold off;
guidata(hObject, handles);



function meshsize_win_Callback(hObject, eventdata, handles)
written = str2double(get(hObject,'String'));
if ~isnan(written) && written > 0
    handles.data.meshsize = written;
else
    errordlg('The mesh size must be given as a positive number.','Error');
    set(handles.meshsize_win,'String', num2str(handles.data.meshsize));
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function meshsize_win_CreateFcn(hObject, eventdata, handles)
% hObject    handle to meshsize_win (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function regparam_win_Callback(hObject, eventdata, handles)

written = str2double(get(hObject,'String'));
if ~isnan(written) && written >= 0
    handles.data.regparam = written;
else
    errordlg('The smoothing parameter must be larger or equal zero.','Error');
    set(handles.regparam_win,'String', num2str(handles.data.regparam));
end
guidata(hObject, handles);



% Hints: get(hObject,'String') returns contents of regparam_win as text
%        str2double(get(hObject,'String')) returns contents of regparam_win as a double


% --- Executes during object creation, after setting all properties.
function regparam_win_CreateFcn(hObject, eventdata, handles)
% hObject    handle to regparam_win (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in show_vectors.
function show_vectors_Callback(hObject, eventdata, handles)

if handles.data.pos ~= 0
    axes(handles.axes1);
    cla; axis equal, hold on;  imagesc(handles.data.bild);
    if (get(hObject,'Value') == get(hObject,'Max'))
        quiver(handles.data.pos(:,1),handles.data.pos(:,2),handles.data.force(:,1),handles.data.force(:,2),2,'r');
    end
    hold off;
end

% Hint: get(hObject,'Value') returns toggle state of show_vectors
