
animal = 'JF070';
myPaths;

slice_spacing = 25;
    structure_alpha = 0.2;
% %load probe2epphys eect, plot eqch probe trrqck 
% allen_atlas_path = '/home/julie/Dropbox/Atlas/allenCCF';
% tv = readNPY([allen_atlas_path, filesep, 'template_volume_10um.npy']);
% av = readNPY([allen_atlas_path, filesep, 'annotation_volume_10um_by_index.npy']);
% st = loadStructureTreeJF([allen_atlas_path, filesep, 'structure_tree_safe_2017.csv']);    

load([outputDir '/probe_ccf.mat'])
load([outputDir '/probe2ephys.mat'])
cmap_filename = [allenAtlasPath filesep 'allenCCF/allen_ccf_colormap_2017.mat'];

figure();
load(cmap_filename);
for curr_probe = 1:length(probe_ccf)
    curr_axes = subplot(1,length(probe_ccf),curr_probe);
    
    trajectory_area_boundaries = ...
        [1;find(diff(probe_ccf(curr_probe).trajectory_areas) ~= 0);length(probe_ccf(curr_probe).trajectory_areas)];    
    trajectory_area_centers = trajectory_area_boundaries(1:end-1) + diff(trajectory_area_boundaries)/2;
        for iArea = 1:size(trajectory_area_centers, 1)
        trajectory_area_labels(iArea) = st.acronym(st.id == ...
            probe_ccf(curr_probe).trajectory_areas(round(trajectory_area_centers(iArea))));
        end

    %trajectory_area_labels = st.acronym(probe_ccf(curr_probe).trajectory_areas(round(trajectory_area_centers)));
      
    image(probe_ccf(curr_probe).trajectory_areas);
    colormap(curr_axes,cmap);
    caxis([1,size(cmap,1)])
    set(curr_axes,'YTick',trajectory_area_centers,'YTickLabels',trajectory_area_labels);
    set(curr_axes,'XTick',[]);
%    title([num2str(probe2ephys(curr_probe).day) num2str(probe2ephys(curr_probe).site)]);
    
end

%    probe_ccf(19).points = probe_ccf(16).points;
%    probe_ccf(19).points(:,1) = probe_ccf(16).points(:,1)+60;
%    probe_ccf(19).trajectory_coords = probe_ccf(16).trajectory_coords;
%   probe_ccf(19).trajectory_coords(:,1) = probe_ccf(16).trajectory_coords(:,1)+60;
%    probe_ccf(19).trajectory_areas = probe_ccf(16).trajectory_areas
%  trajectory_coords
[~, brain_outline] = plotBrainGrid([], []);
theseColors = {rgb('DeepSkyBlue');   rgb('DarkOrange'); rgb('Hotpink');rgb('SeaGreen');rgb('Crimson')};
regionsNames = {'CP', 'GPe', 'GPi', 'STN', 'SNr'};
for iRegion = 1:size(regionsNames,2)
    curr_plot_structure = find(strcmp(st.acronym, regionsNames{iRegion}));
 structure_3d = isosurface(permute(av(1:slice_spacing:end, ...
                1:slice_spacing:end, 1:slice_spacing:end) == curr_plot_structure, [3, 1, 2]), 0); 
            hold on;
        axis vis3d equal off manual
        view([-30, 25]);
        caxis([0, 300]);
        [ap_max, dv_max, ml_max] = size(tv);
        xlim([-10, ap_max + 10])
        ylim([-10, ml_max + 10])
        zlim([-10, dv_max + 10])
        structure_patch = patch('Vertices', structure_3d.vertices*slice_spacing, ...
            'Faces', structure_3d.faces, ...
            'FaceColor', theseColors{iRegion}, 'EdgeColor', 'none', 'FaceAlpha', structure_alpha);
       
end
colsThese = lines(length(probe_ccf));
if size(colsThese,1) > 7
%    hashes = qq;
else
    
end

for thisthisProbe = 3%:5%length(probe_ccf)

          r0 = mean(probe_ccf(thisthisProbe).points, 1);
                xyz = bsxfun(@minus, probe_ccf(thisthisProbe).points, r0);
                [~, ~, V] = svd(xyz, 0);
                histology_probe_direction = V(:, 1);

                probe_eval_points = [-1000, 1000];
                probe_line_endpoints = bsxfun(@plus, bsxfun(@times, probe_eval_points', histology_probe_direction'), r0);
                line(probe_line_endpoints(:, 1), probe_line_endpoints(:, 3), probe_line_endpoints(:, 2), 'color', ...
                    colsThese(thisthisProbe,:), 'linewidth', 2)
end



