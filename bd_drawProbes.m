function bd_drawProbes(tv, av, st, registeredImage, outputDir, screenToUse)
% based on:
% AP_get_probe_histology(tv,av,st,slice_im_path)
%
% Get probe trajectory in histology and convert to ccf
% Andy Peters (peters.andrew.j@gmail.com)

if ~exist('im_type', 'var')
    im_type = '';
end
if ~exist('fullImg', 'var')
    fullImg = [];
end
% Initialize guidata
gui_data = struct;
gui_data.tv = tv;
gui_data.av = av;
gui_data.st = st;

% Query number of probes from user
gui_data.n_probes = str2num(cell2mat(inputdlg('How many probes?')));

% Load in slice images

gui_data.slice_im_path = outputDir;
    
for curr_slice = 1:size(registeredImage,3)
    gui_data.slice_im{curr_slice} = registeredImage(:,:,curr_slice);
end
    

% Load corresponding CCF slices
ccf_slice_fn = [outputDir filesep '/manual/histology_ccf.mat'];
load(ccf_slice_fn);
gui_data.histology_ccf = histology_ccf;

% Load histology/CCF alignment
ccf_alignment_fn = [outputDir filesep '/manual/atlas2histology_tform.mat'];
load(ccf_alignment_fn);
gui_data.histology_ccf_alignment = atlas2histology_tform;

% Warp area labels by histology alignment
gui_data.histology_aligned_av_slices = cell(length(gui_data.slice_im),1);
for curr_slice = 1:length(gui_data.histology_ccf)
    curr_av_slice = gui_data.histology_ccf(curr_slice).av_slices;
    curr_av_slice(isnan(curr_av_slice)) = 1;
    curr_slice_im = gui_data.slice_im{curr_slice};
    
    tform = affine2d;
    tform.T = gui_data.histology_ccf_alignment{curr_slice};
    tform_size = imref2d([size(curr_slice_im,1),size(curr_slice_im,2)]);
    gui_data.histology_aligned_av_slices{curr_slice} = ...
        imwarp(curr_av_slice,tform,'nearest','OutputView',tform_size);
end

% Create figure, set button functions
gui_fig = figure('KeyPressFcn',@keypress, 'Color', 'k');

SCRSZ = screensize(screenToUse); %Get user's screen size
set(gui_fig, 'Position', SCRSZ);
screenPortrait = SCRSZ(4) > SCRSZ(3);
if screenPortrait
    gui_data.histology_ax = axes('Position', [0.05, 0.5, 0.9, 0.45], 'YDir','reverse');
    gui_button_position = [SCRSZ(3)*0.1, SCRSZ(4)*0.45, 100, 40];
else
    gui_data.histology_ax = axes('Position', [0.1, 0.1, 0.5, 0.5], 'YDir','reverse');
    gui_button_position = [SCRSZ(3)*0.1, SCRSZ(4)*0.45, 100, 40]; %QQ change

end

% 'add probe' button 
gui_data.add_probe_btn = uicontrol('Style', 'pushbutton',...
'String', 'Add a probe',...
'Position',gui_button_position,...
'BackgroundColor', rgb('SpringGreen'),...
'CallBack', @(varargin) addProbeButtonPushed(gui_fig));

% 'toggle visibility all' button 
gui_data.toggle_probe_btn = uicontrol('Style', 'pushbutton',...
'String', 'Hide all other probes',...
'Position',gui_button_position + [200, 0, 50, 0],...
'BackgroundColor', rgb('LightCyan'),...
'CallBack', @(varargin) toggleAllProbeButtonPushed(gui_fig));


gui_data.curr_slice = 1;

% Set up axis for histology image
%gui_data.histology_ax = axes();
hold on; colormap(gray); axis image off;
gui_data.slice_im{1}(gui_data.slice_im{1}>1200) = 0;
img = imadjust(gui_data.slice_im{1}, [0.1, 0.8]);
img(img == 1) = 0;
gui_data.histology_im_h = imagesc(img, ...
    'Parent',gui_data.histology_ax);

%caxis([min(min(gui_data.histology_im_h.CData )), max(max(gui_data.histology_im_h.CData ))])
colormap(gray)
% Create title to write area in
gui_data.histology_ax_title = title(gui_data.histology_ax,'','FontSize',14,'Color','white');

% Initialize probe points
gui_data.probe_color = lines(gui_data.n_probes);
gui_data.probe_points_histology = cell(length(gui_data.slice_im),gui_data.n_probes);
gui_data.probe_lines = gobjects(gui_data.n_probes,1);


% initialize probe buttons
nProbes_fit = floor((SCRSZ(4) - gui_button_position(2)) / 50);
nCols = ceil(gui_data.n_probes/ nProbes_fit);
colSpacing = (SCRSZ(3)-100) / nCols;
for iProbe = 1:gui_data.n_probes
    %if iProbe <= nProbes_fit
    nextCol = ceil((iProbe)/(nProbes_fit));
    gui_data.select_probe_btns(iProbe) = uicontrol('Style', 'pushbutton',...
    'String', ['Probe', num2str(iProbe)],...
    'Position', gui_button_position - [60 - (nextCol -1) * colSpacing , 40*(iProbe - ((nProbes_fit) * (nextCol-1))), 0,10],...
    'BackgroundColor', gui_data.probe_color(iProbe,:),...
    'CallBack', @(varargin) selectProbeButtonPushed(gui_fig));

    gui_data.resest_probe_btns(iProbe) = uicontrol('Style', 'pushbutton',...
    'String', 'reset',...
    'Position', gui_button_position - [60 - (nextCol -1) * colSpacing - 100, 40*(iProbe - ((nProbes_fit) * (nextCol-1))), 50,10],...
    'BackgroundColor', gui_data.probe_color(iProbe,:),...
    'CallBack', @(varargin) resetProbeButtonPushed(gui_fig));

    gui_data.del_probe_btns(iProbe) = uicontrol('Style', 'pushbutton',...
    'String', 'hide',...
    'Position', gui_button_position - [60 - (nextCol -1) * colSpacing - 150, 40*(iProbe - ((nProbes_fit) * (nextCol-1))), 60,10],...
    'BackgroundColor', gui_data.probe_color(iProbe,:),...
    'CallBack',  @(varargin) toggleVisiblityProbeButtonPushed(gui_fig));
    
     gui_data.fourShank_probe_btns(iProbe) = uicontrol('Style', 'pushbutton',...
    'String', '4shank with n+3',...
    'Position',gui_button_position - [60 - (nextCol -1) * colSpacing - 190, 40*(iProbe - ((nProbes_fit) * (nextCol -1))), -20,10],...
    'BackgroundColor', gui_data.probe_color(iProbe,:),...
    'CallBack',  @(varargin) fourShankProbeButtonPushed(gui_fig));

end


% Upload gui data
guidata(gui_fig,gui_data);

% Update the slice
update_slice(gui_fig);

end

function keypress(gui_fig,eventdata)

% Get guidata
gui_data = guidata(gui_fig);

switch eventdata.Key
    
    % left/right: move slice
    case 'leftarrow'
        gui_data.curr_slice = max(gui_data.curr_slice - 1,1);
        guidata(gui_fig,gui_data);
        update_slice(gui_fig);
        
    case 'rightarrow'
        gui_data.curr_slice = ...
            min(gui_data.curr_slice + 1,length(gui_data.slice_im));
        guidata(gui_fig,gui_data);
        update_slice(gui_fig);
        
    % Number: add coordinates for the numbered probe
    case [cellfun(@num2str,num2cell(0:9),'uni',false),cellfun(@(x) ['numpad' num2str(x)],num2cell(0:9),'uni',false)]
        keyData = eventdata;
        if isempty(keyData.Modifier)
             curr_probe = str2num(eventdata.Key(end));
        elseif strcmp(keyData.Modifier{:}, 'shift') == 1
             curr_probe = str2num(eventdata.Key(end)) + 10;
        elseif strcmp(keyData.Modifier{:}, 'alt') == 1
            curr_probe = str2num(eventdata.Key(end)) + 20;
        elseif strcmp(keyData.Modifier{:}, 'control') == 1
            curr_probe = str2num(eventdata.Key(end)) + 30;
        end
            
        
        if curr_probe > gui_data.n_probes
           disp(['Probe ' eventdata.Key ' selected, only ' num2str(gui_data.n_probes) ' available']);
           return
        end

        if curr_probe == 0 %quirk in my keyboard
            curr_probe = 4;
        end
        
        update_curr_probe(gui_fig, curr_probe)
        
        % Upload gui data
        guidata(gui_fig,gui_data);
    case 'insert'% if probe # > 9, need an input dialog
        
        probeN = str2num(cell2mat(inputdlg('Probe #: ')));
        curr_probe = probeN;
        gui_data.curr_probe = curr_probe;
        update_curr_probe(gui_fig, curr_probe)
       
        
    case 'escape'
        opts.Default = 'Yes';
        opts.Interpreter = 'tex';
        user_confirm = questdlg('\fontsize{15} Save and quit?','Confirm exit',opts);
        if strcmp(user_confirm,'Yes')
            
            % Initialize structure to save
            probe_ccf = struct( ...
                'points',cell(gui_data.n_probes,1), ...
                'trajectory_coords',cell(gui_data.n_probes,1), ....
                'trajectory_areas',cell(gui_data.n_probes,1));
            
            % Convert probe points to CCF points by alignment and save                     
            for curr_probe = 1:gui_data.n_probes             
                for curr_slice = find(~cellfun(@isempty,gui_data.probe_points_histology(:,curr_probe)'))
                    
                    % Transform histology to atlas slice
                    tform = affine2d;
                    tform.T = gui_data.histology_ccf_alignment{curr_slice};
                    % (transform is CCF -> histology, invert for other direction)
                    tform = invert(tform);
                   
                    % Transform and round to nearest index
                    [probe_points_atlas_x,probe_points_atlas_y] = ...
                        transformPointsForward(tform, ...
                        gui_data.probe_points_histology{curr_slice,curr_probe}(:,1), ...
                        gui_data.probe_points_histology{curr_slice,curr_probe}(:,2));
                   
                    probe_points_atlas_x = round(probe_points_atlas_x);
                    probe_points_atlas_y = round(probe_points_atlas_y);

                    % Get CCF coordinates corresponding to atlas slice points
                    % (CCF coordinates are in [AP,DV,ML])
                    use_points = find(~isnan(probe_points_atlas_x) & ~isnan(probe_points_atlas_y));
                    for curr_point = 1:length(use_points)
                        ccf_ap = gui_data.histology_ccf(curr_slice). ...
                            plane_ap( probe_points_atlas_y(curr_point), ...
                            probe_points_atlas_x(curr_point));
                        ccf_ml = gui_data.histology_ccf(curr_slice). ...
                            plane_dv(probe_points_atlas_y(curr_point), ...
                            probe_points_atlas_x(curr_point));
                        ccf_dv = gui_data.histology_ccf(curr_slice). ...
                            plane_ml(probe_points_atlas_y(curr_point), ...
                            probe_points_atlas_x(curr_point));
                        probe_ccf(curr_probe).points = ...
                            vertcat(probe_ccf(curr_probe).points,[ccf_ap,ccf_dv,ccf_ml]);
                    end              
                end
                
                % Sort probe points by DV (probe always top->bottom)
                [~,dv_sort_idx] = sort(probe_ccf(curr_probe).points(:,2));
                probe_ccf(curr_probe).points = probe_ccf(curr_probe).points(dv_sort_idx,:);
                
            end
            
            % Get areas along probe trajectory            
            for curr_probe = 1:gui_data.n_probes
                
                % Get best fit line through points as probe trajectory
                r0 = mean(probe_ccf(curr_probe).points,1);
                xyz = bsxfun(@minus,probe_ccf(curr_probe).points,r0);
                [~,~,V] = svd(xyz,0);
                histology_probe_direction = V(:,1);
                % (make sure the direction goes down in DV - flip if it's going up)
                if histology_probe_direction(2) < 0
                    histology_probe_direction = -histology_probe_direction;
                end
                
                line_eval = [-1000,1000];
                probe_fit_line = bsxfun(@plus,bsxfun(@times,line_eval',histology_probe_direction'),r0)';
                
                % Get the positions of the probe trajectory
                trajectory_n_coords = max(abs(diff(probe_fit_line,[],2)));
                [trajectory_ap_ccf,trajectory_dv_ccf,trajectory_ml_ccf] = deal( ...
                    round(linspace(probe_fit_line(1,1),probe_fit_line(1,2),trajectory_n_coords)), ...
                    round(linspace(probe_fit_line(2,1),probe_fit_line(2,2),trajectory_n_coords)), ...
                    round(linspace(probe_fit_line(3,1),probe_fit_line(3,2),trajectory_n_coords)));
                
                trajectory_coords_outofbounds = ...
                    any([trajectory_ap_ccf;trajectory_dv_ccf;trajectory_ml_ccf] < 1,1) | ...
                    any([trajectory_ap_ccf;trajectory_dv_ccf;trajectory_ml_ccf] > size(gui_data.av)',1);
                              
                trajectory_coords = ...
                    [trajectory_ap_ccf(~trajectory_coords_outofbounds)' ...
                    trajectory_dv_ccf(~trajectory_coords_outofbounds)', ...
                    trajectory_ml_ccf(~trajectory_coords_outofbounds)'];
                
                trajectory_coords_idx = sub2ind(size(gui_data.av), ...
                    trajectory_coords(:,1),trajectory_coords(:,2),trajectory_coords(:,3));
                
                trajectory_areas_uncut = gui_data.av(trajectory_coords_idx)';
                
                % Get rid of NaN's and start/end 1's (non-parsed)
                trajectory_areas_parsed = find(trajectory_areas_uncut > 1);
                use_trajectory_areas = trajectory_areas_parsed(1): ...
                    trajectory_areas_parsed(end);
                trajectory_areas = reshape(trajectory_areas_uncut(use_trajectory_areas),[],1);
                
                probe_ccf(curr_probe).trajectory_coords = double(trajectory_coords(use_trajectory_areas,:));
                probe_ccf(curr_probe).trajectory_areas = double(trajectory_areas);
                
            end
            
            % Save probe CCF points
            save_fn = [gui_data.slice_im_path filesep 'probe_ccf.mat'];
            save(save_fn,'probe_ccf');
            disp(['Saved probe locations in ' save_fn])
            
            % Close GUI
            %close(gui_fig)
            
            % Plot probe trajectories
            plot_probe(gui_data,probe_ccf);
            
        end
end

end

function update_curr_probe(gui_fig, curr_probe)
% Get guidata
gui_data = guidata(gui_fig);

        set(gui_data.histology_ax_title,'String',['Draw probe ' num2str(curr_probe)]);
        curr_line = imline;
        % If the line is just a click, don't include
        curr_line_length = sqrt(sum(abs(diff(curr_line.getPosition,[],1)).^2));
        if curr_line_length == 0
            return
        end
        gui_data.probe_points_histology{gui_data.curr_slice,curr_probe} = ...
            curr_line.getPosition;
        set(gui_data.histology_ax_title,'String', ...
            ['Arrows to move, Number to draw probe [' num2str(1:gui_data.n_probes) '], Esc to save/quit']);
        
        % Delete movable line, draw line object
        curr_line.delete;
        if size(gui_data.probe_points_histology{gui_data.curr_slice,curr_probe},2)>1 % delete any previous line
            gui_data.probe_points_histology{gui_data.curr_slice,curr_probe}(1,:) = [];
        end
        gui_data.probe_lines(curr_probe) = ...
            line(gui_data.probe_points_histology{gui_data.curr_slice,curr_probe}(:,1), ...
            gui_data.probe_points_histology{gui_data.curr_slice,curr_probe}(:,2), ...
            'linewidth',3,'color',gui_data.probe_color(curr_probe,:));
        gui_data.curr_probe = curr_probe;
        % Upload gui data
        guidata(gui_fig,gui_data);
end

function addProbeButtonPushed(gui_fig)
% Get guidata
gui_data = guidata(gui_fig);

% Add a probe
gui_data.n_probes = gui_data.n_probes + 1;
gui_data.probe_color = lines(gui_data.n_probes);
gui_data.probe_points_histology(:,gui_data.n_probes) = cell(length(gui_data.slice_im),1);
gui_data.probe_lines(gui_data.n_probes) = gobjects(1,1);

% Update all buttons 

nProbes_fit = floor((SCRSZ(4) - gui_button_position(2)) / 50);
nCols = ceil(gui_data.n_probes/ nProbes_fit);
colSpacing = (SCRSZ(3)-100) / nCols;
for iProbe = 1:gui_data.n_probes
    %if iProbe <= nProbes_fit
    nextCol = ceil((iProbe)/(nProbes_fit));
    gui_data.select_probe_btns(iProbe) = uicontrol('Style', 'pushbutton',...
    'String', ['Probe', num2str(iProbe)],...
    'Position', gui_button_position - [60 - (nextCol -1) * colSpacing , 40*(iProbe - ((nProbes_fit) * (nextCol-1))), 0,10],...
    'BackgroundColor', gui_data.probe_color(iProbe,:),...
    'CallBack', @(varargin) selectProbeButtonPushed(gui_fig));

    gui_data.resest_probe_btns(iProbe) = uicontrol('Style', 'pushbutton',...
    'String', 'reset',...
    'Position', gui_button_position - [60 - (nextCol -1) * colSpacing - 100, 40*(iProbe - ((nProbes_fit) * (nextCol-1))), 50,10],...
    'BackgroundColor', gui_data.probe_color(iProbe,:),...
    'CallBack', {@resetProbeButtonPushed, gui_fig});

    gui_data.del_probe_btns(iProbe) = uicontrol('Style', 'pushbutton',...
    'String', 'del',...
    'Position', gui_button_position - [60 - (nextCol -1) * colSpacing - 150, 40*(iProbe - ((nProbes_fit) * (nextCol-1))), 60,10],...
    'BackgroundColor', gui_data.probe_color(iProbe,:),...
    'CallBack', {@delProbeButtonPushed, gui_fig});
    
     gui_data.fourShank_probe_btns(iProbe) = uicontrol('Style', 'pushbutton',...
    'String', '4shank with n+3',...
    'Position',gui_button_position - [60 - (nextCol -1) * colSpacing - 190, 40*(iProbe - ((nProbes_fit) * (nextCol -1))), -20,10],...
    'BackgroundColor', gui_data.probe_color(iProbe,:),...
    'CallBack', {@fourShankProbeButtonPushed, gui_fig});

end

% Upload gui data
guidata(gui_fig,gui_data);

end

function toggleAllProbeButtonPushed(gui_fig)
% Get guidata
gui_data = guidata(gui_fig);

% hide/show 
if contains(gui_data.toggle_probe_btn.String, 'Hide')
    gui_data.toggle_probe_btn.String = 'Show all other probes';
    other_probes = logical(ones(gui_data.n_probes,1));
    other_probes(gui_data.curr_probe) = 0;
    other_probes = find(other_probes); 
    gui_data.visibility = 0;
    for iProbe = 1:size(other_probes,1)
        gui_data.probe_lines(iProbe).Color(4) = 0;
    end
else
    gui_data.toggle_probe_btn.String = 'Hide all other probes';
    other_probes = logical(ones(gui_data.n_probes,1));
    other_probes(gui_data.curr_probe) = 0;
    other_probes = find(other_probes); 
    gui_data.visibility = 1;
    for iProbe = 1:size(other_probes,1)
        gui_data.probe_lines(iProbe).Color(4) = 1;
    end
end
% Upload gui data
guidata(gui_fig,gui_data);
end

function selectProbeButtonPushed(gui_fig)
% Get guidata
gui_data = guidata(gui_fig);

% Get current probe 
curr_probe = find([gui_data.select_probe_btns(:).Value]);
gui_data.curr_probe = curr_probe;
% Change curr probe 
update_curr_probe(gui_fig, curr_probe)

% Upload gui data
guidata(gui_fig,gui_data);

end

function fourShankProbeButtonPushed(gui_fig)
% Get guidata
gui_data = guidata(gui_fig);

% Get current probe 
curr_probe = find([gui_data.select_probe_btns(:).Value]);

% Change curr probe 
update_curr_probe(gui_fig, curr_probe)

% Upload gui data
guidata(gui_fig,gui_data);
end

function toggleVisiblityProbeButtonPushed(gui_fig)
% Get guidata
gui_data = guidata(gui_fig);

% hide/show 
if contains(gui_data.toggle_probe_btn.String, 'Hide')
    gui_data.toggle_probe_btn.String = 'Show all other probes';
    other_probes = logical(ones(gui_data.n_probes,1));
    other_probes(gui_data.curr_probe) = 0;
    other_probes = find(other_probes); 
    gui_data.visibility = 0;
    for iProbe = 1:size(other_probes,1)
        gui_data.probe_lines(iProbe).Color(4) = 0;
    end
else
    gui_data.toggle_probe_btn.String = 'Hide all other probes';
    other_probes = logical(ones(gui_data.n_probes,1));
    other_probes(gui_data.curr_probe) = 0;
    other_probes = find(other_probes); 
    gui_data.visibility = 1;
    for iProbe = 1:size(other_probes,1)
        gui_data.probe_lines(iProbe).Color(4) = 1;
    end
end
% Upload gui data
guidata(gui_fig,gui_data);

end

function update_slice(gui_fig)
% Draw histology and CCF slice

% Get guidata
gui_data = guidata(gui_fig);

% Set next histology slice
set(gui_data.histology_im_h,'CData',gui_data.slice_im{gui_data.curr_slice})

% Clear any current lines, draw probe lines
gui_data.probe_lines.delete;
for curr_probe = find(~cellfun(@isempty,gui_data.probe_points_histology(gui_data.curr_slice,:)))
    if gui_data.visibility == 0 && gui_data.curr_probe ~= curr_probe
        gui_data.probe_lines(curr_probe) = ...
            line(gui_data.probe_points_histology{gui_data.curr_slice,curr_probe}(:,1), ...
            gui_data.probe_points_histology{gui_data.curr_slice,curr_probe}(:,2), ...
            'linewidth',3,'color',[gui_data.probe_color(curr_probe,:),0]);
    else

        gui_data.probe_lines(curr_probe) = ...
            line(gui_data.probe_points_histology{gui_data.curr_slice,curr_probe}(:,1), ...
            gui_data.probe_points_histology{gui_data.curr_slice,curr_probe}(:,2), ...
            'linewidth',3,'color',[gui_data.probe_color(curr_probe,:),1]);
    end
end

set(gui_data.histology_ax_title,'String', ...
            ['Arrows to move, Number to draw probe [' num2str(1:gui_data.n_probes) '], Esc to save/quit']);
        
% Upload gui data
guidata(gui_fig, gui_data);

end

function plot_probe(gui_data,probe_ccf)

% Plot probe trajectories
figure('Name','Probe trajectories');
axes_atlas = axes;
[~, brain_outline] = plotBrainGrid([],axes_atlas);
set(axes_atlas,'ZDir','reverse');
hold(axes_atlas,'on');
axis vis3d equal off manual
view([-30,25]);
caxis([0 300]);
[ap_max,dv_max,ml_max] = size(gui_data.tv);
% xlim([-10,ml_max+10])
% ylim([-10,dv_max+10])
% zlim([-10,ap_max+10])
h = rotate3d(gca);
h.Enable = 'on';

for curr_probe = 1:length(probe_ccf)
    % Plot points and line of best fit
    thesePoints = probe_ccf(curr_probe).points*2.5; % QQ 2.5 to correct for atlas 25 um where we draw probes and 10 um in plotBrainGrid 
    r0 = mean(thesePoints,1);
    xyz = bsxfun(@minus,thesePoints,r0);
    [~,~,V] = svd(xyz,0);
    %V= permute(V, [3, 2, 1]);
    histology_probe_direction = V(:,1);
    % (make sure the direction goes down in DV - flip if it's going up)
    if histology_probe_direction(2) < 0
        histology_probe_direction = -histology_probe_direction;
    end
    
    line_eval = [-1000,1000];
    probe_fit_line = bsxfun(@plus,bsxfun(@times,line_eval',histology_probe_direction'),r0);
    plot3(thesePoints(:,1), ...
        thesePoints(:,2), ...
        thesePoints(:,3), ...
        '.','color',gui_data.probe_color(curr_probe,:),'MarkerSize',20);
    line(probe_fit_line(:,1),probe_fit_line(:,2),probe_fit_line(:,3), ...
        'color',gui_data.probe_color(curr_probe,:),'linewidth',2)
end

% Plot probe areas
figure('Name','Trajectory areas');
% (load the colormap - located in the repository, find by associated fcn)
allenCCF_path = fileparts(which('allenCCFbregma'));
cmap_filename = [allenCCF_path filesep 'allen_ccf_colormap_2017.mat'];
load(cmap_filename);
for curr_probe = 1:length(probe_ccf)
    curr_axes = subplot(1,gui_data.n_probes,curr_probe);
    
    trajectory_area_boundaries = ...
        [1;find(diff(probe_ccf(curr_probe).trajectory_areas) ~= 0);length(probe_ccf(curr_probe).trajectory_areas)];    
    trajectory_area_centers = trajectory_area_boundaries(1:end-1) + diff(trajectory_area_boundaries)/2;
   
%     trajectory_area_labels = gui_data.st.acronym(...
%         ismember(gui_data.st.id, probe_ccf(curr_probe).trajectory_areas(round(trajectory_area_centers))));
%      
    for iArea = 1:size(trajectory_area_centers,1)
        trajectory_area_labels(iArea) = gui_data.st.acronym(gui_data.st.id == ...
            probe_ccf(curr_probe).trajectory_areas(round(trajectory_area_centers(iArea))));
    end
    image(probe_ccf(curr_probe).trajectory_areas);
    colormap(curr_axes,cmap);
    caxis([1,size(cmap,1)])
    set(curr_axes,'YTick',trajectory_area_centers,'YTickLabels',trajectory_area_labels);
    set(curr_axes,'XTick',[]);
    title(['Probe ' num2str(curr_probe)]);
    
end

end









