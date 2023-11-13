function ya_checkAndCorrectAtlasAlignment(tv, av, st, registeredIm, saveDir, screenPortrait)
% based on:
% AP_manual_align_histology_ccf(tv,av,st,slice_im_path)
%
% Align histology slices and matched CCF slices
% Andy Peters (peters.andrew.j@gmail.com)
% things to add: 
% - smoothing
% - brigthness
% - slider
% - raw original image (and trnasformed atlas)
% - structure "selection" and border moving 
% - non-affine transform? 

% Initialize guidata
gui_data = struct;
gui_data.tv = tv;
gui_data.av = av;
gui_data.st = st;

if nargin < 6 || isempty(screenPortrait)
    screenPortrait = 1;
end
% Load in slice images
gui_data.slice_im_path = saveDir;

SCRSZ = get(0,'screensize');
%screenPortrait = SCRSZ(4)>SCRSZ(3);
   
gui_data.slice_im = cell(size(registeredIm,3),1);
for curr_slice = 1:size(registeredIm,3)
    gui_data.slice_im{curr_slice} = registeredIm(:,:,curr_slice);
end

% Load corresponding CCF slices
ccf_slice_fn = [saveDir filesep '/manual/histology_ccf.mat'];
load(ccf_slice_fn);
gui_data.histology_ccf = histology_ccf;

% Load automated alignment
auto_ccf_alignment_fn = [saveDir filesep '/manual/atlas2histology_tform.mat'];
if exist(auto_ccf_alignment_fn,'file')
    load(auto_ccf_alignment_fn);
    gui_data.histology_ccf_auto_alignment = atlas2histology_tform;
end

% Create figure, set button functions
gui_fig = figure('KeyPressFcn',@keypress, 'Color', 'k');
gui_data.curr_slice = 120;

% Set up axis for histology image

if screenPortrait
    gui_data.histology_ax = subplot(2,1,2,'YDir','reverse'); 
    set(gui_data.histology_ax,'Position',[0,0,1.0,0.5]);
else %assume portrait mode 
    gui_data.histology_ax = subplot(1,2,1,'YDir','reverse'); 
    set(gui_data.histology_ax,'Position',[0,0,0.5,0.9]);
end
    
hold on; colormap(gray); axis image off;
gui_data.histology_im_h = imagesc(gui_data.slice_im{gui_data.curr_slice}, ...
    'Parent',gui_data.histology_ax,'ButtonDownFcn',@mouseclick_histology);
colormap(gray)

% Set up histology-aligned atlas overlay
% (and make it invisible to mouse clicks)
histology_aligned_atlas_boundaries_init = ...
    zeros(size(gui_data.slice_im{1},1),size(gui_data.slice_im{gui_data.curr_slice},2));
gui_data.histology_aligned_atlas_boundaries = ...
    imagesc(histology_aligned_atlas_boundaries_init,'Parent',gui_data.histology_ax, ...
    'AlphaData',histology_aligned_atlas_boundaries_init,'PickableParts','none');

structure_aligned_atlas_boundaries_init = ...
    zeros(size(gui_data.slice_im{1},1),size(gui_data.slice_im{gui_data.curr_slice},2));
gui_data.structure_aligned_atlas_boundaries = ...
    imagesc(structure_aligned_atlas_boundaries_init,'Parent',gui_data.histology_ax, ...
    'AlphaData',structure_aligned_atlas_boundaries_init,'PickableParts','none');

% Set up axis for atlas slice
if screenPortrait
    gui_data.atlas_ax = subplot(2,1,1,'YDir','reverse'); 
    set(gui_data.atlas_ax,'Position',[0,0.5,1.0,0.5]);
else %assume portrait mode 
    gui_data.atlas_ax = subplot(1,2,2,'YDir','reverse'); 
    set(gui_data.atlas_ax,'Position',[0.5,0,0.5,0.9]);
end

hold on; axis image off; colormap(gray); caxis([0,400]);
gui_data.atlas_im_h = imagesc(gui_data.slice_im{1}, ...
    'Parent',gui_data.atlas_ax,'ButtonDownFcn',@mouseclick_atlas);

% Set up histology-aligned atlas overlay
% (and make it invisible to mouse clicks)
atlas_aligned_atlas_boundaries_init = ...
    zeros(size(gui_data.histology_ccf(1).tv_slices,1),size(gui_data.histology_ccf(1).tv_slices,2));
gui_data.atlas_aligned_atlas_boundaries = ...
    imagesc(atlas_aligned_atlas_boundaries_init,'Parent',gui_data.atlas_ax, ...
    'AlphaData',atlas_aligned_atlas_boundaries_init,'PickableParts','none');

% Initialize alignment control points and tform matricies
gui_data.histology_control_points = repmat({zeros(0,2)},length(gui_data.slice_im),1);
gui_data.atlas_control_points = repmat({zeros(0,2)},length(gui_data.slice_im),1);

gui_data.histology_control_points_plot = plot(gui_data.histology_ax,nan,nan,'.w','MarkerSize',20);
gui_data.atlas_control_points_plot = plot(gui_data.atlas_ax,nan,nan,'.r','MarkerSize',20);

gui_data.histology_ccf_manual_alignment = gui_data.histology_ccf_auto_alignment;

% Structure selection 
gui_data.structure = 0;
gui_data.structureText = uicontrol('style','text',...                              %Textbox "size" to set filter size. 
    'unit','pix',...      
    'position',[SCRSZ(3)-250 SCRSZ(4)-SCRSZ(2)-250 220 100],...
    'fontsize',16,...
    'string',st.name(gui_data.st.id == gui_data.structure),...
    'BackgroundColor', 'k',...
    'ForegroundColor', [1, 0, 1]);
%set(gui_data.structureText,'Callback',{@selectStructure,gui_data});

% Upload gui data
guidata(gui_fig,gui_data);

% Initialize alignment
align_ccf_to_histology(gui_fig);

% Print controls
CreateStruct.Interpreter = 'tex';
CreateStruct.WindowStyle = 'non-modal';
msgbox( ...
    {'\fontsize{12}' ...
    '\bf Controls: \rm' ...
    '1,2 : switch slice' ...
    'click : set reference points for manual alignment (3 minimum)', ...
    'space : toggle alignment overlay visibility', ...
    'c : clear manually placed points', ...
    's : save', ...
    'Escape: save and close'}, ...
    'Controls',CreateStruct);

end


function keypress(gui_fig,eventdata)

% Get guidata
gui_data = guidata(gui_fig);

switch eventdata.Key
    
    % 1/2: move slice
    case '1'
        gui_data.curr_slice = max(gui_data.curr_slice - 1,1);
        guidata(gui_fig,gui_data);
        update_slice(gui_fig);
        
    case '2'
        gui_data.curr_slice = ...
            min(gui_data.curr_slice + 1,length(gui_data.slice_im));
        guidata(gui_fig,gui_data);
        update_slice(gui_fig);
        
    % O: toggle overlay visibility
    case 'space'
        curr_visibility = ...
            get(gui_data.histology_aligned_atlas_boundaries,'Visible');
        set(gui_data.histology_aligned_atlas_boundaries,'Visible', ...
            cell2mat(setdiff({'on','off'},curr_visibility)))

        curr_visibility_atlas = ...
            get(gui_data.atlas_aligned_atlas_boundaries,'Visible');
        set(gui_data.atlas_aligned_atlas_boundaries,'Visible', ...
            cell2mat(setdiff({'on','off'},curr_visibility_atlas)))
        
    % C: clear current points
    case 'c'
        gui_data.histology_control_points{gui_data.curr_slice} = zeros(0,2);
        gui_data.atlas_control_points{gui_data.curr_slice} = zeros(0,2);
        
        guidata(gui_fig,gui_data);
        update_slice(gui_fig);
        
    % S: save
    case 's'
        atlas2histology_tform = ...
            gui_data.histology_ccf_manual_alignment;
        save_fn = [gui_data.slice_im_path filesep 'atlas2histology_tform.mat'];
        save(save_fn,'atlas2histology_tform');
        disp(['Saved ' save_fn]);
        
    % Escape: save and exit
    case 'escape'
        opts.Default = 'Yes';
        opts.Interpreter = 'tex';
        user_confirm = questdlg('\fontsize{15} Save and quit?','Confirm exit',opts);
        if strcmp(user_confirm,'Yes')            
            atlas2histology_tform = ...
                gui_data.histology_ccf_manual_alignment;
            save_fn = [gui_data.slice_im_path filesep 'atlas2histology_tform.mat'];
            save(save_fn,'atlas2histology_tform');
            disp(['Saved ' save_fn]);
            close(gui_fig);            
        end
        
end

end

% function selectStructure(varargin)
% % Get guidata
% [h,gui_data] = varargin{[1,3]};
% 
% 
% end
function mouseclick_histology(gui_fig,eventdata)
% Draw new point for alignment

% Get guidata
gui_data = guidata(gui_fig);

if eventdata.Button == 3 %right click, toggle nearest structure/global
        selectedStructure = gui_data.histology_ccf(gui_data.curr_slice).av_slices(round(eventdata.IntersectionPoint(2)), ...
            round(eventdata.IntersectionPoint(1)));
        set(gui_data.structureText, 'String', gui_data.st.name(gui_data.st.id ==selectedStructure)) 
        gui_data.st.name(gui_data.st.id == selectedStructure)

        curr_av_slice = gui_data.histology_ccf(gui_data.curr_slice).av_slices;
        curr_av_slice(isnan(curr_av_slice)) = 1;
        
        av_boundaries = curr_av_slice== selectedStructure;
        av_warp_boundaries_orange(:,:,1) = av_boundaries;
        av_warp_boundaries_orange(:,:,2) = zeros(size(av_boundaries));
        av_warp_boundaries_orange(:,:,3) = av_boundaries;
        
        set(gui_data.structure_aligned_atlas_boundaries, ...
            'CData',av_warp_boundaries_orange, ...
            'AlphaData',av_boundaries*0.3);
        
else %eventdata.Button == 1 => left click 
% Add clicked location to control points
gui_data.histology_control_points{gui_data.curr_slice} = ...
    vertcat(gui_data.histology_control_points{gui_data.curr_slice}, ...
    eventdata.IntersectionPoint(1:2));

set(gui_data.histology_control_points_plot, ...
    'XData',gui_data.histology_control_points{gui_data.curr_slice}(:,1), ...
    'YData',gui_data.histology_control_points{gui_data.curr_slice}(:,2));

% Upload gui data
guidata(gui_fig, gui_data);

% If equal number of histology/atlas control points > 3, draw boundaries
if size(gui_data.histology_control_points{gui_data.curr_slice},1) == ...
        size(gui_data.atlas_control_points{gui_data.curr_slice},1) 
    align_ccf_to_histology(gui_fig)
end
end
end


function mouseclick_atlas(gui_fig,eventdata)
% Draw new point for alignment

% Get guidata
gui_data = guidata(gui_fig);

% Add clicked location to control points
gui_data.atlas_control_points{gui_data.curr_slice} = ...
    vertcat(gui_data.atlas_control_points{gui_data.curr_slice}, ...
    eventdata.IntersectionPoint(1:2));

set(gui_data.atlas_control_points_plot, ...
    'XData',gui_data.atlas_control_points{gui_data.curr_slice}(:,1), ...
    'YData',gui_data.atlas_control_points{gui_data.curr_slice}(:,2));

% Upload gui data
guidata(gui_fig, gui_data);

% If equal number of histology/atlas control points > 3, draw boundaries
if size(gui_data.histology_control_points{gui_data.curr_slice},1) == ...
        size(gui_data.atlas_control_points{gui_data.curr_slice},1) || ...
        (size(gui_data.histology_control_points{gui_data.curr_slice},1) > 3 && ...
        size(gui_data.atlas_control_points{gui_data.curr_slice},1) > 3)
    align_ccf_to_histology(gui_fig)
end

end


function align_ccf_to_histology(gui_fig)

% Get guidata
gui_data = guidata(gui_fig);

% Define a weight for blending based on the number of control points
num_points = size(gui_data.histology_control_points{gui_data.curr_slice},1);
weight_manual = min(num_points, 10) / 10; % Weight increases with number of points, max of 1

if num_points >= 3   && numel(gui_data.atlas_control_points{gui_data.curr_slice}) ==...
        numel(gui_data.histology_control_points{gui_data.curr_slice})
    % Use control point alignment for manual transformation
    tform_manual = fitgeotrans(gui_data.atlas_control_points{gui_data.curr_slice}, ...
        gui_data.histology_control_points{gui_data.curr_slice},'affine');
else
    % No manual transformation if less than 3 points
    tform_manual = affine2d(eye(3));
end

% Get automated transformation
if isfield(gui_data,'histology_ccf_auto_alignment')
    tform_auto = affine2d;
    tform_auto.T = gui_data.histology_ccf_auto_alignment{gui_data.curr_slice};
else
    tform_auto = affine2d(eye(3));
end

% Blend the transformations
tform_blended = affine2d(eye(3));
tform_blended.T = (1 - weight_manual) * tform_auto.T + weight_manual * tform_manual.T;

curr_av_slice = squeeze(gui_data.histology_ccf(gui_data.curr_slice).av_slices);
curr_av_slice(isnan(curr_av_slice)) = 1;
curr_slice_im = gui_data.slice_im{gui_data.curr_slice};

tform_size = imref2d([size(curr_slice_im,1),size(curr_slice_im,2)]);
curr_av_slice_warp = imwarp(curr_av_slice, tform_blended, 'OutputView',tform_size);

av_warp_boundaries = round(conv2(curr_av_slice_warp,ones(3)./9,'same')) ~= curr_av_slice_warp;
av_warp_boundaries_red(:,:,1) = av_warp_boundaries;
av_warp_boundaries_red(:,:,2) = zeros(size(av_warp_boundaries));
av_warp_boundaries_red(:,:,3) = zeros(size(av_warp_boundaries));

set(gui_data.histology_aligned_atlas_boundaries, ...
    'CData',av_warp_boundaries_red, ...
    'AlphaData',av_warp_boundaries*0.3);


atlas_boundaries = round(conv2(curr_av_slice,ones(3)./9,'same')) ~= curr_av_slice;
atlas_boundaries_red(:,:,1) = atlas_boundaries;
atlas_boundaries_red(:,:,2) = zeros(size(curr_av_slice));
atlas_boundaries_red(:,:,3) = zeros(size(curr_av_slice));

set(gui_data.atlas_aligned_atlas_boundaries, ...
    'CData',atlas_boundaries_red, ...
    'AlphaData', atlas_boundaries_red(:,:,1)*0.3);

% Update transform matrix
gui_data.histology_ccf_manual_alignment{gui_data.curr_slice} = tform_blended.T;

% Upload gui data
guidata(gui_fig, gui_data);

end


function update_slice(gui_fig)
% Draw histology and CCF slice

% Get guidata
gui_data = guidata(gui_fig);

% Set next histology slice
set(gui_data.histology_im_h,'CData',gui_data.slice_im{gui_data.curr_slice})
maxColVal = arrayfun(@(x) max(max(gui_data.slice_im{x})), 1:size(gui_data.slice_im,1));
thisImage = gui_data.slice_im{gui_data.curr_slice};
thisImage(thisImage>prctile(maxColVal(maxColVal>0), 90)) = prctile(maxColVal(maxColVal>0), 90);
set(gui_data.histology_im_h,'CData',thisImage)

set(gui_data.histology_im_h, 'CData', gui_data.slice_im{gui_data.curr_slice})

% Plot control points for slice
set(gui_data.histology_control_points_plot, ...
    'XData',gui_data.histology_control_points{gui_data.curr_slice}(:,1), ...
    'YData',gui_data.histology_control_points{gui_data.curr_slice}(:,2));
set(gui_data.atlas_control_points_plot, ...
    'XData',gui_data.atlas_control_points{gui_data.curr_slice}(:,1), ...
    'YData',gui_data.atlas_control_points{gui_data.curr_slice}(:,2));

% Reset histology-aligned atlas boundaries if not
histology_aligned_atlas_boundaries_init = ...
    zeros(size(gui_data.slice_im{1},1),size(gui_data.slice_im{1},2));
set(gui_data.histology_aligned_atlas_boundaries, ...
    'CData',histology_aligned_atlas_boundaries_init, ...
    'AlphaData',histology_aligned_atlas_boundaries_init);

structure_aligned_atlas_boundaries_init = ...
    zeros(size(gui_data.histology_ccf(1).tv_slices,1),size(gui_data.histology_ccf(1).tv_slices,2));
set(gui_data.structure_aligned_atlas_boundaries, ...
    'CData',structure_aligned_atlas_boundaries_init, ...
    'AlphaData',structure_aligned_atlas_boundaries_init);

atlas_aligned_atlas_boundaries_init = ...
    zeros(size(gui_data.histology_ccf(1).tv_slices,1),size(gui_data.histology_ccf(1).tv_slices,2));
set(gui_data.atlas_aligned_atlas_boundaries, ...
    'CData',atlas_aligned_atlas_boundaries_init, ...
    'AlphaData',atlas_aligned_atlas_boundaries_init);

% Upload gui data
guidata(gui_fig, gui_data);

% Update atlas boundaries
align_ccf_to_histology(gui_fig)

end



















