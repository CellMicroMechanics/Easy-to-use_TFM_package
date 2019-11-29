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

% Last Modified by GUIDE v2.5 21-Nov-2019 19:47:44

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%% Setting initial parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
handles.data = varargin{1};

dir_struct = vertcat(dir(fullfile(handles.data.imagedir_name,'*.tif*')),dir(fullfile(handles.data.imagedir_name,'*.jpg*')));
[sorted_names,sorted_index] = sortrows({dir_struct.name}');
if ~isempty(sorted_names)
    set(handles.preview_image,'String',sorted_names,'Value',1)
else
    set(handles.preview_image,'String','No image available','Value',1)
end

for i = 1:length(handles.data.displacement)
    hilf_cell{i,1} = num2str(i);
end
set(handles.preview_frame,'String',hilf_cell,'Value',1);

%set length units for figure
set(handles.figure1, 'units', 'normalized');

% use the first nonempty dataset to calculate a grid spacing that is
% suggested in the menu
data_nonempty =  handles.data.data_nonempty;
first_index = find(data_nonempty(:,1),1);
if ~isempty(first_index)
    first_useful_frame = data_nonempty(first_index,2);
else
    errordlg('The dataset only contains empty or too small dispacement fields.','Error');
end
max_eck(1:2) = [max(handles.data.displacement(first_useful_frame).pos(:,1)), max(handles.data.displacement(first_useful_frame).pos(:,2))];
min_eck(1:2) = [min(handles.data.displacement(first_useful_frame).pos(:,1)), min(handles.data.displacement(first_useful_frame).pos(:,2))];
handles.data.meshsize = round(sqrt((max_eck(1)-min_eck(1))*(max_eck(2)-min_eck(2))/size(handles.data.displacement(first_useful_frame).pos,1)));
set(handles.meshsize_win,'String',num2str(handles.data.meshsize));
set(handles.meshsize_win,'String',num2str(handles.data.meshsize));

%provide arbitrary inital suggestion for regularization parameter
handles.data.regparam = 100;
set(handles.regparam_win,'String',num2str(handles.data.regparam,'%f'));

%initialize variables
handles.data.pos = [];
handles.data.bild = [];

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Performs Regularized TFM calculations for the data sequence and save result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%IMPORTANT: for strain energy calculations, we crop the grid by
%Crop_rim_perc percent from all sides to avoid boundary effects. Change if desired
Crop_rim_perc = 15; %value in the rang 0-50 [%] 

%initialize waitbar
haha = waitbar(0,'Please wait while data is being assembled..','WindowStyle','modal');

%init
framenumber = length(handles.data.displacement);
data_nonempty =  handles.data.data_nonempty;
grid_mat = [];
 
for frame = 1:framenumber    
     if  ~data_nonempty(frame,1)
         continue;
     end

       waitbar(frame/framenumber);
       
       Orig_pos = handles.data.displacement(frame).pos(:,1:2);  % position of displacement
       vec = handles.data.displacement(frame).vec(:,1:2);       % value of displacement
       New_pos = Orig_pos +vec;                            % shifted positions
       
       [grid_mat,u, i_max,j_max] = interp_vec2grid(New_pos,vec,handles.data.meshsize,grid_mat);
       rescale = handles.data.young;
       E = 1;                       %%% scale the Young modulus to 1            
       Ftu(:,:,1) = fft2(u(:,:,1));
       Ftu(:,:,2) = fft2(u(:,:,2));

      %calculate traction forces
      [TFM_results(frame).pos,TFM_results(frame).traction,TFM_results(frame).traction_magnitude,f_n_m,~,~] = ...
                          reg_fourier_TFM(Ftu(:,:,1), Ftu(:,:,2), handles.data.regparam, E, handles.data.poisson,...
                          handles.data.meshsize, i_max, j_max, grid_mat, handles.data.pix_durch_my, handles.data.zdepth);                  
  
       %scale tractions with correct Young's modulus
       TFM_results(frame).traction = rescale*TFM_results(frame).traction;
       TFM_results(frame).traction_magnitude = rescale*TFM_results(frame).traction_magnitude; 
       f_mat = f_n_m.*rescale;
       
       bnd = floor(min(size(u(:,:,1)))*Crop_rim_perc/100+1);
       %calculate strain energy. Note the cropping of the field by bnd.
       TFM_results(frame).energy = 1/2*sum(sum(u(bnd:end-bnd+1,bnd:end-bnd+1,1).*f_mat(bnd:end-bnd+1,bnd:end-bnd+1,1) +...
                                   u(bnd:end-bnd+1,bnd:end-bnd+1,2).*f_mat(bnd:end-bnd+1,bnd:end-bnd+1,2)))*...
                                   (handles.data.meshsize)^2*handles.data.pix_durch_my^3/10^6; 

        u_reshape(:,1) = reshape(u(:,:,1),i_max*j_max,1);
        u_reshape(:,2) = reshape(u(:,:,2),i_max*j_max,1);
        TFM_results(frame).displacement = u_reshape;      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%% save TFM_settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

 msgbox('Calculation completed.');

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Preview the Regularized TFM results for individual frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if we have enough data in the dataset number
frame = get(handles.preview_frame,'value');
if ~handles.data.data_nonempty(frame,1)
    errordlg('The dataset number refers to an empty or too small displacement field.','Error');
    return;
end

rescale = handles.data.young;
E = 1;    %%% scale the Young modulus to 1

%assigning an empty grid_mat variable will trigger automatic determination of
%a new grid size and location
grid_mat = [];

Orig_pos = handles.data.displacement(frame).pos(:,1:2);  % position of displacement
vec = handles.data.displacement(frame).vec(:,1:2);       % value of displacement
New_pos = Orig_pos +vec;                            % shifted positions

% interpolate displacement data onto a grid
[grid_mat,u, i_max,j_max] = interp_vec2grid(New_pos,vec,handles.data.meshsize,grid_mat);

%Fourier-transform the displacement field
Ftu(:,:,1) = fft2(u(:,:,1));
Ftu(:,:,2) = fft2(u(:,:,2));

%calculate traction forces
[handles.data.pos,f_nm_2,traction_magnitude,f_n_m,~,~] = ...
    reg_fourier_TFM(Ftu(:,:,1), Ftu(:,:,2), handles.data.regparam, E, handles.data.poisson, handles.data.meshsize,...
    i_max,j_max, grid_mat,handles.data.pix_durch_my, handles.data.zdepth);

%scale traction forces with correct Young's modulus
handles.data.force = rescale.*f_nm_2;
handles.data.f_mat = rescale.* f_n_m;

%Next we show the preview results in the GUI
%get file name of image that we display with the TFM preview
bild_datei_index = get(handles.preview_image,'value');
bild_dateien = get(handles.preview_image,'string');
%If we do not get a cell array, we need to make one (eg if only one entry)
if ~isempty(bild_dateien) && ~iscellstr(bild_dateien)
   bild_dateien =cellstr(bild_dateien);
end
%check if we have a proper image and load it
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
%display the cell image and plot vectors
axes(handles.axes1);
cla; axis equal, hold on; colormap(gca,gray);
if ~no_image
    imagesc(handles.data.bild);
end
hilf = get(handles.show_vectors,'Value');
if hilf
    quiver(handles.data.pos(:,1),handles.data.pos(:,2),handles.data.force(:,1),handles.data.force(:,2),2,'r');
end
set(gca, 'DataAspectRatio', [1,1,50],'YDir','reverse','XTick',[],'YTick',[],'YColor','w','XColor','w')
hold off;

%display the heatmap showing traction magnitude
axes(handles.axes2);
colorbar off;
cla; hold on; colormap(gca,jet);
fnorm = (handles.data.f_mat(:,:,2).^2 + handles.data.f_mat(:,:,1).^2).^0.5;
surf(grid_mat(:,:,1), grid_mat(:,:,2),fnorm),view(0,90),shading interp, axis equal;
set(gca, 'DataAspectRatio', [1,1,50],'YDir','reverse','XTick',[],'YTick',[],'YColor','w','XColor','w');
hilf = get(handles.showColorbar,'Value');
if hilf
    colorbar('location','East','YColor','w','XColor','w');
end
hold off;

 %pass variables back (if needed)
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
    errordlg('The smoothing regularization parameter must be bigger or equal to zero.','Error');
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

if ~isempty(handles.data.pos)
    axes(handles.axes1);
    cla; axis equal, hold on; 
    if ~isempty(handles.data.bild)
        imagesc(handles.data.bild);
    end
    if (get(hObject,'Value') == get(hObject,'Max'))
        quiver(handles.data.pos(:,1),handles.data.pos(:,2),handles.data.force(:,1),handles.data.force(:,2),2,'r');
    end
    hold off;
end

% Hint: get(hObject,'Value') returns toggle state of show_vectors


% --- Executes on button press in showColorbar.
function showColorbar_Callback(hObject, eventdata, handles)
% hObject    handle to showColorbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showColorbar

axes(handles.axes2);
if (get(hObject,'Value') == get(hObject,'Max'))
    hold on;
    colorbar('location','East','YColor','w','XColor','w');
    hold off;
else
    colorbar off;
end
