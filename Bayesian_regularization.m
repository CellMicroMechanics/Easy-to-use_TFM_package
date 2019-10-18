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

function varargout = Bayesian_regularization(varargin)
%BAYESIAN_REGULARIZATION MATLAB code file for Bayesian_regularization.fig
%      BAYESIAN_REGULARIZATION, by itself, creates a new BAYESIAN_REGULARIZATION or raises the existing
%      singleton*.
%
%      H = BAYESIAN_REGULARIZATION returns the handle to a new BAYESIAN_REGULARIZATION or the handle to
%      the existing singleton*.
%
%      BAYESIAN_REGULARIZATION('Property','Value',...) creates a new BAYESIAN_REGULARIZATION using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to Bayesian_regularization_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      BAYESIAN_REGULARIZATION('CALLBACK') and BAYESIAN_REGULARIZATION('CALLBACK',hObject,...) call the
%      local function named CALLBACK in BAYESIAN_REGULARIZATION.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Bayesian_regularization

% Last Modified by GUIDE v2.5 16-Aug-2019 08:22:39

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Bayesian_regularization_OpeningFcn, ...
                   'gui_OutputFcn',  @Bayesian_regularization_OutputFcn, ...
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


% --- Executes just before Bayesian_regularization is made visible.
function Bayesian_regularization_OpeningFcn(hObject, eventdata, handles, varargin)
    
handles.data = varargin{1};

dir_struct = vertcat(dir(fullfile(handles.data.imagedir_name,'*.tif*')),dir(fullfile(handles.data.imagedir_name,'*.jpg*')));
[sorted_names,sorted_index] = sortrows({dir_struct.name}');
set(handles.preview_image,'String',sorted_names,'Value',1)
for i = 1:length(handles.data.displacement)
    hilf_cell{i,1} = num2str(i);
end
set(handles.preview_frame,'String',hilf_cell,'Value',1);

max_eck(1:2) = [max(handles.data.displacement(1).pos(:,1)), max(handles.data.displacement(1).pos(:,2))];
min_eck(1:2) = [min(handles.data.displacement(1).pos(:,1)), min(handles.data.displacement(1).pos(:,2))];

for i = 1:length(handles.data.noise)
    hilf_cell_noise{i,1} = num2str(i);
end

handles.data.meshsize = round(sqrt((max_eck(1)-min_eck(1))*(max_eck(2)-min_eck(2))/size(handles.data.displacement(1).pos,1)));
set(handles.meshsize_win,'String',num2str(handles.data.meshsize));

set(handles.figure1, 'units', 'normalized');

handles.data.regparam = 0.2/handles.data.young;

handles.data.pos = 0;
handles.output = hObject;

handles.data.leftalpha= 1;
handles.data.rightalpha = 1e7;
handles.data.stepalpha = 1;
handles.data.iterationnumber = 200;

% Update handles structure
guidata(hObject, handles);




% UIWAIT makes Bayesian_regularization wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Bayesian_regularization_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in abort.
function abort_Callback(hObject, eventdata, handles)
delete(handles.figure1);
return;


% --- Executes on button press in analyze.
function analyze_Callback(hObject, eventdata, handles)
handles.data.zdepth = 0.0; 
grid_mat = [];
framenumber = length(handles.data.displacement);

%%% scale the Young modulus to 1
ne = handles.data.young;
young= 1;

reg_parameter=[];
for frame = 1:framenumber
    if isempty(handles.data.displacement(frame).pos) || isempty(handles.data.displacement(frame).vec)||...
       length(handles.data.displacement(frame).pos) < 2
         disp(['There dataset of frame ', num2str(frame),' was not valid. Skipping this frame!']);
         continue;
    end
    
  [grid_mat,u, i_max,j_max] = interp_vec2grid(handles.data.displacement(frame).pos+...
              handles.data.displacement(frame).vec, handles.data.displacement(frame).vec,handles.data.meshsize);
  kx_vec = 2*pi/i_max/handles.data.meshsize.*[0:(i_max/2-1) (-i_max/2:-1)];
  ky_vec = 2*pi/j_max/handles.data.meshsize.*[0:(j_max/2-1) (-j_max/2:-1)];
  kx = repmat(kx_vec',1,j_max);
  ky = repmat(ky_vec,i_max,1);
    
  kx(1,1) = 1;
  ky(1,1) = 1;
  k = sqrt(kx.^2+ky.^2);
  
  conf = 2.*(1+handles.data.poisson)./(young.*k.^3);
  Ginv_xx = conf .* ((1-handles.data.poisson).*k.^2+handles.data.poisson.*ky.^2);
  Ginv_xy = conf .* (-handles.data.poisson.*kx.*ky);
  Ginv_yy = conf .* ((1-handles.data.poisson).*k.^2+handles.data.poisson.*kx.^2);

  Ginv_xx(1,1) = 0;
  Ginv_yy(1,1) = 0;
  Ginv_xy(1,1) = 0;
  
  Ginv_xy(i_max/2+1,:) = 0;
  Ginv_xy(:,j_max/2+1) = 0;

  G1 = reshape(Ginv_xx,[1,i_max*j_max]);
  G2 = reshape(Ginv_yy,[1,i_max*j_max]);
  X1 = reshape([G1; G2], [], 1)';
  
  G3 = reshape(Ginv_xy,[1,i_max*j_max]);
  G4 = zeros(1, i_max*j_max); 
  X2 = reshape([G4; G3], [], 1)';
  X3 = X2(1,2:end);
   
  X4 = diag(X1);
  X5 = diag(X3,1);
  X6 = diag(X3,-1);
  X=X4+X5+X6;
 
  Ftu(:,:,1) = fft2(u(:,:,1));
  Ftu(:,:,2) = fft2(u(:,:,2));

  fux1=reshape(Ftu(:,:,1),i_max*j_max,1);
  fuy1=reshape(Ftu(:,:,2),i_max*j_max,1);

  fuu(1:2:size(fux1)*2,1) = fux1;
  fuu(2:2:size(fuy1)*2,1) = fuy1;
  
  blah = get(handles.noiseway, 'SelectedObject');
  noise_way=get(blah,'String');
  
  switch noise_way
    case 'Noise from input file'

        size_noise = size(handles.data.noise);
        if  size_noise(2) ==1                      
           errordlg('No noise from input file. Using manual noise selection','Error');
           return;
        end
  
        if  frame==1
         using_noise = handles.data.noise(frame);          
         haha = waitbar(0,'Please wait while data is being assembled..','WindowStyle','modal');
        else 
         using_noise = handles.data.noise(frame);
        end
    case 'Manual noise selection' 
        
       if  frame==1

        bild_dateien_1 = get(handles.preview_image,'string');
        frame1 = get(handles.preview_frame,'value');
        setappdata(0,'frame1', frame1 );
        setappdata(0,'bild_dateien_1', bild_dateien_1 );
        Select_noise(handles.data)
        [xx,yy]=getline('closed');
        
       if isempty(xx) || isempty(yy)
        errordlg('The selected noise is empty.','Error');
        return;
       end     
        delete(Select_noise(handles.data));
     
        indata = inpolygon(handles.data.displacement(frame).pos(:,1),handles.data.displacement(frame).pos(:,2),xx,yy);
        using_noise.pos = handles.data.displacement(frame).pos(indata,:);
        using_noise.vec = handles.data.displacement(frame).vec(indata,:);
     
     haha = waitbar(0,'Please wait while data is being assembled..','WindowStyle','modal');
     else
         
     using_noise.pos = handles.data.displacement(frame).pos(indata,:);
     using_noise.vec = handles.data.displacement(frame).vec(indata,:);
     
     end
      
  end
  
   waitbar(frame/framenumber);

  noise_u(1:2:size(using_noise.vec,1)*2,1) = using_noise.vec(:,1);
  noise_u(2:2:size(using_noise.vec,1)*2,1) = using_noise.vec(:,2);
  beta = 1/var(noise_u);
  
   
  [L]= optimal_lambda_sequence(beta,fuu,Ftu(:,:,1),Ftu(:,:,2),...
                          young,handles.data.poisson,handles.data.meshsize,i_max, j_max,X);

  handles.data.regparam=[];
  handles.data.regparam = L;
   
   [TFM_results(frame).pos,TFM_results(frame).traction,TFM_results(frame).traction_magnitude,f] = ...
                         reg_fourier_TFM(grid_mat,u,young,handles.data.poisson, handles.data.pix_durch_my,...
                         handles.data.zdepth, handles.data.meshsize, i_max,j_max, L);

   TFM_results(frame).traction=ne*TFM_results(frame).traction;
   TFM_results(frame).traction_magnitude=ne*TFM_results(frame).traction_magnitude;
   
   bnd = 6;
   f = f.*ne;
   TFM_results(frame).energy = 1/2*sum(sum(u(bnd:end-bnd+1,bnd:end-bnd+1,1).*f(bnd:end-bnd+1,bnd:end-bnd+1,1) +...
         u(bnd:end-bnd+1,bnd:end-bnd+1,2).*f(bnd:end-bnd+1,bnd:end-bnd+1,2)))*(handles.data.meshsize)^2*...
         handles.data.pix_durch_my^3/10^6; 
     
       
    pu(:,1) = reshape(u(:,:,1),i_max*j_max,1);
    pu(:,2) = reshape(u(:,:,2),i_max*j_max,1);
    TFM_results(frame).displacement= pu; 
    reg_parameter(end+1,1)=L;
end

TFM_settings.poisson = handles.data.poisson;
TFM_settings.young = handles.data.young;
TFM_settings.micrometer_per_pix = handles.data.pix_durch_my;
TFM_settings.regularization_parameter = reg_parameter;
TFM_settings.meshsize = handles.data.meshsize;
TFM_settings.zdepth = handles.data.zdepth;

TFM_settings.i_max = i_max;
TFM_settings.j_max = j_max;

reg_parameter=[];

 close(haha);

[filepath,name,ext] = fileparts(handles.data.strain_noise_file_name); 
handles.data.targetdir_name = filepath;                               

savefile_name = fullfile(handles.data.targetdir_name,['Bay-FTTC_results_',datestr(now, 'dd-mm-yy'),'.mat']);
if exist(savefile_name)
    button = questdlg('The file exists already. Overwrite?','Error','Yes');
    if strcmpi(button,'No') || strcmpi(button,'')
            return;
    end
end

save(savefile_name,'TFM_results','TFM_settings','-mat');

clear XX;
clear yy;

return;


% --- Executes on button press in show_vectors.
function show_vectors_Callback(hObject, eventdata, handles)
if handles.data.pos ~= 0
    axes(handles.axes2);
    cla; axis equal, hold on;  imagesc(handles.data.bild);
    if (get(hObject,'Value') == get(hObject,'Max'))
        quiver(handles.data.pos(:,1),handles.data.pos(:,2),handles.data.force(:,1),handles.data.force(:,2),2,'r');
    end
    hold off;
end

% Hint: get(hObject,'Value') returns toggle state of show_vectors


% --- Executes on selection change in preview_image.
function preview_image_Callback(hObject, eventdata, handles)
% hObject    handle to preview_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns preview_image contents as cell array
%        contents{get(hObject,'Value')} returns selected item from preview_image


% --- Executes during object creation, after setting all properties.
function preview_image_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in preview_frame.
function preview_frame_Callback(hObject, eventdata, handles)
% hObject    handle to preview_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function preview_frame_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in preview.
function preview_Callback(hObject, eventdata, handles)
    
handles.data.zdepth = 0.0;    

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

[grid_mat,u, i_max,j_max] = interp_vec2grid(handles.data.displacement(frame).pos+...
              handles.data.displacement(frame).vec, handles.data.displacement(frame).vec,handles.data.meshsize);

  kx_vec = 2*pi/i_max/handles.data.meshsize.*[0:(i_max/2-1) (-i_max/2:-1)];
  ky_vec = 2*pi/j_max/handles.data.meshsize.*[0:(j_max/2-1) (-j_max/2:-1)];
  kx = repmat(kx_vec',1,j_max);
  ky = repmat(ky_vec,i_max,1);
    
  kx(1,1) = 1;
  ky(1,1) = 1;
  k = sqrt(kx.^2+ky.^2);
  
  conf = 2.*(1+handles.data.poisson)./(young.*k.^3);
  Ginv_xx = conf .* ((1-handles.data.poisson).*k.^2+handles.data.poisson.*ky.^2);
  Ginv_xy = conf .* (-handles.data.poisson.*kx.*ky);
  Ginv_yy = conf .* ((1-handles.data.poisson).*k.^2+handles.data.poisson.*kx.^2);

  Ginv_xx(1,1) = 0;
  Ginv_yy(1,1) = 0;
  Ginv_xy(1,1) = 0;
  
  Ginv_xy(i_max/2+1,:) = 0;
  Ginv_xy(:,j_max/2+1) = 0;
  
  G1 = reshape(Ginv_xx,[1,i_max*j_max]);
  G2 = reshape(Ginv_yy,[1,i_max*j_max]);
  X1 = reshape([G1; G2], [], 1)';
  
  G3 = reshape(Ginv_xy,[1,i_max*j_max]);
  G4 = zeros(1, i_max*j_max); 
  X2 = reshape([G4; G3], [], 1)';
  X3 = X2(1,2:end);
   
  X4 = diag(X1);
  X5 = diag(X3,1);
  X6 = diag(X3,-1);
  X=X4+X5+X6;
  
  Ftu(:,:,1) = fft2(u(:,:,1));
  Ftu(:,:,2) = fft2(u(:,:,2));

  fux1=reshape(Ftu(:,:,1),i_max*j_max,1);
  fuy1=reshape(Ftu(:,:,2),i_max*j_max,1);

  fuu(1:2:size(fux1)*2,1) = fux1;
  fuu(2:2:size(fuy1)*2,1) = fuy1;
  
  blah = get(handles.noiseway, 'SelectedObject');
  noise_way=get(blah,'String');
  
  switch noise_way
    case 'Noise from input file'
        size_noise = size(handles.data.noise);
        if  size_noise(2) ==1               
           errordlg('No noise from input file. Using manual noise selection','Error');
           return;
        end
          using_noise = handles.data.noise(frame);
      
    case 'Manual noise selection' 
        
     bild_dateien_1 = get(handles.preview_image,'string');
     frame1 = get(handles.preview_frame,'value');
     setappdata(0,'frame1', frame1 );
     setappdata(0,'bild_dateien_1', bild_dateien_1 );
     Select_noise(handles.data)
  
     [xx,yy]=getline('closed');
     
     if isempty(xx) || isempty(yy)
       errordlg('The selected noise is empty.','Error');
       return;
     end

     delete(Select_noise(handles.data));
     
     nxx=size(xx);
     if isempty(xx) || nxx(1)==2
        errordlg('Noise does not select','Error');
        return;
     end
      
     indata = inpolygon(handles.data.displacement(frame).pos(:,1),handles.data.displacement(frame).pos(:,2),xx,yy);
     using_noise.pos = handles.data.displacement(frame).pos(indata,:);
     using_noise.vec = handles.data.displacement(frame).vec(indata,:);
     
  end
  
   clear XX;
   clear yy;
  
  
  noise_u(1:2:size(using_noise.vec,1)*2,1) = using_noise.vec(:,1);
  noise_u(2:2:size(using_noise.vec,1)*2,1) = using_noise.vec(:,2);
  beta = 1/var(noise_u);

  please_wait;
 
 [L evidencep evidence_one]= optimal_lambda(beta,fuu,Ftu(:,:,1),Ftu(:,:,2),...
                             young,handles.data.poisson,handles.data.meshsize,i_max, j_max,X);
 
 handles.data.evidencep=[];
 handles.data.evidence_one=[];
 handles.data.evidencep = evidencep;
 handles.data.evidence_one = evidence_one;

 handles.data.regparam=[];
 handles.data.regparam = L;

 [handles.data.pos,traction,traction_magnitude,f] = ...
                           reg_fourier_TFM(grid_mat,u,young,handles.data.poisson, handles.data.pix_durch_my,...
                           handles.data.zdepth, handles.data.meshsize, i_max,j_max, L);

 handles.data.bild = imread(fullfile(handles.data.imagedir_name, bild_datei));

 handles.data.force=ne.*traction;
 handles.data.f_mat=ne.*f;

 axes(handles.axes1)
 cla;hold on;
 plot(handles.data.evidencep(1,:), handles.data.evidencep(2,:),'o');
 xlabel('\lambda_2E^2');
 ylabel('Log(Evidence)');
 plot(handles.data.regparam, handles.data.evidence_one,'r*','MarkerSize',10);
 ps_t=[handles.data.regparam handles.data.evidence_one];
 strValues = num2str(ps_t(1),'%.4f');
 text(ps_t(1), 1*ps_t(2),['\lambda_2E^2=', strValues],'FontSize',10,'VerticalAlignment','bottom');
 hold off;

 axes(handles.axes2);
 cla; axis equal, hold on; colormap gray, imagesc(handles.data.bild);
 hilf = get(handles.show_vectors,'Value');
 if hilf
    quiver(handles.data.pos(:,1),handles.data.pos(:,2),handles.data.force(:,1),handles.data.force(:,2),2,'r');
 end
 set(gca, 'DataAspectRatio', [1,1,50],'YDir','reverse','XTick',[],'YTick',[]),hold off;

 axes(handles.axes3);
 cla; hold on; colormap jet;
 fnorm = (handles.data.f_mat(:,:,2).^2 + handles.data.f_mat(:,:,1).^2).^0.5;
 surf(grid_mat(:,:,1), grid_mat(:,:,2),fnorm),view(0,90),shading interp, axis equal;
 set(gca, 'DataAspectRatio', [1,1,50],'YDir','reverse','XTick',[],'YTick',[]),hold off;

 delete(please_wait);

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
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in Take_noise.
 function Take_noise_Callback(hObject, eventdata, handles)

 Select_noise(handles.data);   


% --- Executes on button press in noise_draw.
function noise_draw_Callback(hObject, eventdata, handles)
% hObject    handle to noise_draw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in noise_frame.
function noise_frame_Callback(hObject, eventdata, handles)
% hObject    handle to noise_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function noiseway_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noiseway (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


 
