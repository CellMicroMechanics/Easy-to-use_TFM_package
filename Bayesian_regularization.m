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

% Last Modified by GUIDE v2.5 26-Nov-2019 22:02:00

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%% Setting initial parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
  
handles.data = varargin{1};

%set figure units
set(handles.figure1, 'units', 'normalized');

%get image file names and display them in the drop down menu
dir_struct = vertcat(dir(fullfile(handles.data.imagedir_name,'*.tif*')),dir(fullfile(handles.data.imagedir_name,'*.jpg*')));
[sorted_names,sorted_index] = sortrows({dir_struct.name}');
if ~isempty(sorted_names)
    set(handles.preview_image,'String',sorted_names,'Value',1)
else
    set(handles.preview_image,'String','No image available','Value',1)
end
%get dataset numbers
for i = 1:length(handles.data.displacement)
    hilf_cell{i,1} = num2str(i);
end
set(handles.preview_frame,'String',hilf_cell,'Value',1);

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

%initialize variables
handles.data.pos = [];
handles.data.bild = [];
handles.data.selected_region_with_noise_xy = [];
handles.data.selected_region_with_noise_frame = [];
handles.output = hObject;

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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Performs Bayesian TFM calculations for the data sequence and save result
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% IMPORTANT: for strain energy calculations, we crop the grid by
        %Crop_rim_perc percent from all sides to avoid boundary effects. Change if desired
        Crop_rim_perc = 15; %value between 0 [%] and 50 [%]  
        
        %With Bayesian TFM, we assume for simplicity z=0. This can be 
        %changed, but functions need to be adapted
        handles.data.zdepth = 0.0;
        
        %get all the data sets that contain nonempty displacement fields
        framenumber = length(handles.data.displacement);
        data_nonempty =  handles.data.data_nonempty;
        indices = find(data_nonempty(:,1));
        
        %init
        reg_parameter = [];          %regularization parameter
        grid_mat = [];               %reguar grid
        beta = zeros(framenumber,1); %inverse noise variance
        
        % Select source of noise data
        blah = get(handles.noiseway, 'SelectedObject');
        noise_way = get(blah,'Tag');
        switch noise_way
            case 'inputfile'
                noise_from_inputfile = 1;
                noise_framenumber = length(handles.data.noise);
            otherwise
                noise_from_inputfile = 0;
        end
        
        %calculate the noise variance in the whole data sequence
        for frame = 1:framenumber
            if  ~data_nonempty(frame,1)  %%% skip if no displacement data available
                continue;
            end
            
            if frame == indices(1)   %see if we are dealing with the first frame
                first_frame = 1;
            else
                first_frame = 0;
            end
                  
            
            if noise_from_inputfile %noise comes from input data
               if ~isfield(handles.data,'noise') || isempty(handles.data.noise)
                    %display error message if 'noise from input file'
                    %was selected, but input does not contain noise
                    %sample
                    errordlg('The input file contains no noise data. Select noise sample manually.','Error');
                    return;
               else
                   %check each frame length of noise
                    size_noise_pos = size(handles.data.noise(frame).pos);
                    size_noise_vec = size(handles.data.noise(frame).vec);            
               end
               if  first_frame %is this the first dataset?
                    %search for first dataset containing noise data
                    first_noise_frame_index = 0;
                    frame2 = frame;
                    while (first_noise_frame_index == 0) && (frame2 <= noise_framenumber)
                        size_noise_vec1 = size(handles.data.noise(frame2).vec);
                        %does the variable really contain enough data?
                        if (size_noise_vec1(1) >= 1) && (size_noise_vec1(2) >= 2) && length(handles.data.noise(frame2).vec) > 2
                            %This index is used to calculate noise variance
                            %for frames that have no noise data
                            first_noise_frame_index = frame2;
                        end
                        frame2 = frame2 +1;
                    end
                    if first_noise_frame_index==0
                        %display error message if 'noise from input file'
                        %was selected, but input does not contain noise
                        %sample
                        errordlg('The input file contains no noise data. Select noise sample manually.','Error');
                        return;
                    end
                end
                
                if (frame > noise_framenumber) || size_noise_pos(1) < 1 || size_noise_pos(2) < 2 || size_noise_vec(1) < 1 || size_noise_vec(2) < 2
                    disp(['NOTICE: No noise sample in dataset number ', num2str(frame),'. Using noise in dataset number ', num2str(first_noise_frame_index)]);
                    current_noise_frame_index = first_noise_frame_index;
                else
                    current_noise_frame_index = frame;
                end
                
                used_noise_sample = handles.data.noise(current_noise_frame_index).vec(:,1:2);
                beta(frame) = 1/var(used_noise_sample(:));
                
            else  %noise comes from manually selected region
                  if ~isempty(handles.data.selected_region_with_noise_xy) %Region of interest for noise has been specified
                    
                     %select the noise sample in the polygonal roi that was
                     %manually chosen
                     pxy     = handles.data.selected_region_with_noise_xy;
                     indata  = inpolygon(handles.data.displacement(frame).pos(:,1),handles.data.displacement(frame).pos(:,2),pxy(:,1),pxy(:,2));
                     if isempty(indata) ||(nnz(indata) <2)
                        errordlg('In dataset no. ', num2str(frame),': The ROI for noise selection contains to few displacements.','Error');
                        return;
                     end

                     used_noise_sample = handles.data.displacement(frame).vec(indata,1:2);
                     beta(frame) = 1/var(used_noise_sample(:));
                    
                else %Region of interest for noise is not yet determined. Must select Roi manually now.
                    
                        noise_frame = get(handles.preview_frame,'value');   %%frame for noise selection is the current one in the GUI.                       
                        message = sprintf(['Choose a ROI containing noise in the current dataset no. ', num2str(noise_frame)]);
                        uiwait(msgbox(message, 'Notice'));
                        
                        % check if data is provided
                        if ~data_nonempty(noise_frame,1)                            
                            errordlg('Can not select noise sample. The currently selected data set has too few displacements.','Error');
                            return;
                        end
                        

                        bild_dateien_1 = get(handles.preview_image,'string');
                        bild_datei_index_1 = get(handles.preview_image,'value');

                        % Transfer information to Select_noise.mat
                        setappdata(0,'frame_1', noise_frame);
                        setappdata(0,'bild_datei_index_1', bild_datei_index_1 );
                        setappdata(0,'bild_dateien_1', bild_dateien_1 );
                                                
                        %choose the ROI with the mouse
                        Select_noise(handles.data)
                        try
                            [xx,yy] = getline('closed');
                        catch
                            errordlg('The selected noise data is empty.','Error');
                            return;
                        end
                        
                        %check if ROI contains data
                        if isempty(xx) || isempty(yy)
                            errordlg('The selected noise data is empty.','Error');
                            delete(Select_noise(handles.data));
                            return;
                        end
                        indata = inpolygon(handles.data.displacement(noise_frame).pos(:,1),handles.data.displacement(noise_frame).pos(:,2),xx,yy);                        
                        if (length(xx)<3) ||(nnz(indata) <2)
                            errordlg('The selected area contains to few displacements.','Error');
                            delete(Select_noise(handles.data));
                            return;
                        end
                        
                        delete(Select_noise(handles.data));
                        
                        %save the polygon to select noise in other frames
                        handles.data.selected_region_with_noise_xy = [xx yy];
                        handles.data.selected_region_with_noise_frame= noise_frame;
                        
                        %get the noise sample and calculate the variance
                        used_noise_sample = handles.data.displacement(noise_frame).vec(indata,1:2);
                        beta(frame) = 1/var(used_noise_sample(:));
                        
                   end
            end
        end
               
        E = 1;      %set Young modulus to 1 in calculations
        E_rescale = handles.data.young; %Young modulus for final rescaling
        s  =  handles.data.poisson; %Poisson's modulus
        meshsize = handles.data.meshsize;
        
        haha = waitbar(0,'Please wait while data is being assembled..','WindowStyle','modal');
        for frame = 1:framenumber
            waitbar(frame/framenumber);
            if  ~data_nonempty(frame,1)  %%% displacement data available?
                continue;
            end
            
            %calculate Green's function and fourier transformed
            %displacement field
            orig_pos = handles.data.displacement(frame).pos(:,1:2);
            vec = handles.data.displacement(frame).vec(:,1:2);
            [grid_mat,i_max,j_max, X,fuu,Ftux,Ftuy,u] = fourier_X_u(orig_pos,vec, meshsize, E, s,grid_mat);
            
            %%%Calculate the optimal regularization parameter
            [L,~,~] = optimal_lambda(beta(frame),fuu,Ftux,Ftuy,...
                E,s,handles.data.meshsize,i_max, j_max,X,1);
            
            %calculate traction forces
            [TFM_results(frame).pos,TFM_results(frame).traction,TFM_results(frame).traction_magnitude,f_n_m,~,~] = ...
                reg_fourier_TFM(Ftux,Ftuy,L,E,s, handles.data.meshsize,...
                i_max,j_max, grid_mat,handles.data.pix_durch_my, handles.data.zdepth);
            
            %%% save TFM_results
            TFM_results(frame).traction = E_rescale*TFM_results(frame).traction;
            TFM_results(frame).traction_magnitude = E_rescale*TFM_results(frame).traction_magnitude;
            
            %calculate strain energy. Note that we crop the outer rim of the field by 'bnd' grid
            %units to avoid edge effects
            f_mat=E_rescale.*f_n_m;
            bnd = floor(min(size(u(:,:,1)))*Crop_rim_perc/100+1);
            TFM_results(frame).energy = 1/2*sum(sum(u(bnd:end-bnd+1,bnd:end-bnd+1,1).*f_mat(bnd:end-bnd+1,bnd:end-bnd+1,1) +...
                u(bnd:end-bnd+1,bnd:end-bnd+1,2).*f_mat(bnd:end-bnd+1,bnd:end-bnd+1,2)))*(handles.data.meshsize)^2*...
                handles.data.pix_durch_my^3/10^6;
            
            u_reshape(:,1) = reshape(u(:,:,1),i_max*j_max,1);
            u_reshape(:,2) = reshape(u(:,:,2),i_max*j_max,1);
            TFM_results(frame).displacement = u_reshape;
            reg_parameter(frame,1) = L;
        end
        
        %%% save TFM_settings
        TFM_settings.poisson = handles.data.poisson;
        TFM_settings.young = handles.data.young;
        TFM_settings.micrometer_per_pix = handles.data.pix_durch_my;
        TFM_settings.regularization_parameter = reg_parameter;
        TFM_settings.meshsize = handles.data.meshsize;
        TFM_settings.zdepth = handles.data.zdepth;
        TFM_settings.i_max = i_max;
        TFM_settings.j_max = j_max;
        
        %document the noise source in the saved settings
        if noise_from_inputfile
            TFM_settings.type_noise = ('Noise from input_file');  %% save type of noise
        else
            TFM_settings.type_noise = ['Region of interest for noise selected manually'];  %% save type of noise
        end
        
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
        
        msgbox('Calculation completed.');
        %optional, if we want to pass data back after calculations
        %guidata(hObject, handles);
return;
        

% --- Executes on button press in show_vectors.
function show_vectors_Callback(hObject, eventdata, handles)
if ~isempty(handles.data.pos)
    axes(handles.axes2);
    cla; axis equal; hold on;  
    if ~isempty(handles.data.bild)
        imagesc(handles.data.bild);
    end
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Preview the Bayesian TFM results for individual frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %for the Bayesian method, we only assume z=0 for simplicity
    handles.data.zdepth = 0.0;
    
    %check if we have enough data in the dataset number
    frame = get(handles.preview_frame,'value');
    if ~handles.data.data_nonempty(frame,1)
        errordlg('The dataset number refers to an empty or too small displacement field.','Error');
        return;
    end

    E_rescale = handles.data.young; % initial Young modulus
    E = 1;    %%% scale the Young modulus to 1
    Orig_pos = handles.data.displacement(frame).pos(:,1:2);  %% position of displacement
    vec = handles.data.displacement(frame).vec(:,1:2);  %% value of displacement
    meshsize = handles.data.meshsize;
    s  =  handles.data.poisson;
    grid_mat = [];

    %construct Green's function and Fourier transform the displacements
    [grid_mat,i_max,j_max, X,fuu,Ftux,Ftuy,~] = fourier_X_u(Orig_pos,vec, meshsize,E, s, grid_mat);

    %select source for noise data
    blah = get(handles.noiseway, 'SelectedObject');
    noise_way = get(blah,'Tag');
    switch noise_way
        case 'inputfile'
              if ~isfield(handles.data,'noise') || isempty(handles.data.noise) || (frame > length(handles.data.noise)) ||  isempty(handles.data.noise(frame).vec) || size(handles.data.noise(frame).vec,2) < 2 || size(handles.data.noise(frame).vec,2) < 2
                  errordlg(['Current dataset no. ', num2str(frame), ' contains no noise data. Use manual noise selection.'],'Error');
                  return;
              end
              used_noise_sample = handles.data.noise(frame).vec(:,1:2);

        case 'noisedraw'

            bild_dateien_1 = get(handles.preview_image,'string');
            bild_datei_index_1 = get(handles.preview_image,'value');

            % Transfer information to Select_noise.mat
            setappdata(0,'frame_1', frame);
            setappdata(0,'bild_datei_index_1', bild_datei_index_1 );
            setappdata(0,'bild_dateien_1', bild_dateien_1 );

            Select_noise(handles.data)

            try
                [xx,yy] = getline('closed');
            catch
                errordlg('The selected noise data is empty.','Error');
                delete(Select_noise(handles.data));
                return;
            end
            delete(Select_noise(handles.data)); %% delete the Select_noise GUI
            
            %check if ROI contains displacement noise information
            if isempty(xx) || isempty(yy)
                errordlg('The selected noise data is empty.','Error');
                return;
            end
            indata = inpolygon(handles.data.displacement(frame).pos(:,1),handles.data.displacement(frame).pos(:,2),xx,yy);
            if (length(xx)<3) || (nnz(indata) <2)
                errordlg('The selected area contains to few displacements.','Error');
                return;
            end

            %get the displacement data in the ROI
            used_noise_sample = handles.data.displacement(frame).vec(indata,1:2);
            
            %%% save the polygonal area containing the noise data
            handles.data.selected_region_with_noise_xy = [xx yy];
            handles.data.selected_region_with_noise_frame = frame;
    end

    beta = 1/var(used_noise_sample(:));


    % Main calculation starts here, show waitbar
    wait_f = waitbar(0.5,'Please wait...');
    %Calculate optimal regularization parameter L
    [L evidencep evidence_one] = optimal_lambda(beta,fuu,Ftux,Ftuy,...
        E,s,meshsize,i_max, j_max,X);
    %Calculate tractions
    [handles.data.pos,f_nm_2,traction_magnitude,f_n_m,~,~] = ...
        reg_fourier_TFM(Ftux,Ftuy,L,E,s,meshsize,...
        i_max,j_max, grid_mat,handles.data.pix_durch_my, handles.data.zdepth);
    handles.data.force = E_rescale.*f_nm_2;
    f_mat = E_rescale.*f_n_m;


    % Show results on the GUI
    %plot evidence curve
    cla(handles.axes1,'reset');
    axes(handles.axes1)
    hold on;
    plot(evidencep(1,:), evidencep(2,:),'o');
    set(gca,'FontSize',9);
    xlabel('\lambda_2E^2','FontSize',9);
    ylabel('Log(Evidence)','FontSize',9);
    plot(L, evidence_one,'r*','MarkerSize',10);
    ps_t=[L evidence_one];
    strValues = num2str(ps_t(1),'%.4f');
    lamtexth =text(ps_t(1), 1*ps_t(2),['\lambda_2E^2=', strValues],'FontSize',9,'VerticalAlignment','bottom');
    Txpos = get(lamtexth,'Extent');
    ayLims = get(gca,'YLim');
    if Txpos(2)+Txpos(4) >= ayLims(2)
        set(gca,'YLim', [ayLims(1), Txpos(2)+Txpos(4)*1.5]);
    end
    hold off;

    %Plot cel image and vectors
    %first gather information from the GUI
    bild_datei_index = get(handles.preview_image,'value');
    bild_dateien = get(handles.preview_image,'string');
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

    %plot images and traction vectors
    axes(handles.axes2);
    cla; axis equal, colormap(gca,gray); hold on; 
    if ~no_image 
        imagesc(handles.data.bild);
    end
    hilf = get(handles.show_vectors,'Value');
    if hilf
        quiver(handles.data.pos(:,1),handles.data.pos(:,2),handles.data.force(:,1),handles.data.force(:,2),2,'r');
    end
    set(gca, 'DataAspectRatio', [1,1,50],'YDir','reverse','XTick',[],'YTick',[],'YColor','w','XColor','w'),hold off;

    %plot the heatmap for visualizing traction magnitude
    axes(handles.axes3);
    colorbar off;
    cla; hold on; colormap(gca,jet);
    fnorm = (f_mat(:,:,2).^2 + f_mat(:,:,1).^2).^0.5;
    surf(grid_mat(:,:,1), grid_mat(:,:,2),fnorm),view(0,90),shading interp, axis equal;
    set(gca, 'DataAspectRatio', [1,1,50],'YDir','reverse','XTick',[],'YTick',[],'YColor','w','XColor','w');
    hilf = get(handles.showColorbar,'Value');
    if hilf
        colorbar('location','East','YColor','w','XColor','w');
    end
    hold off;

    %close waitbar
    close(wait_f);

    %pass variables back
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



% --- Executes on button press in noisedraw.
function noisedraw_Callback(hObject, eventdata, handles)


% --- Executes on button press in inputfile.
function inputfile_Callback(hObject, eventdata, handles)
% hObject    handle to inputfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in showColorbar.
function showColorbar_Callback(hObject, eventdata, handles)
% hObject    handle to showColorbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showColorbar
axes(handles.axes3);
if (get(hObject,'Value') == get(hObject,'Max'))
    hold on;
    colorbar('location','East','YColor','w','XColor','w');
    hold off;
else
    colorbar off;
end
