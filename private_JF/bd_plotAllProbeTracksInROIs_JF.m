patchBrain = 0;

% to do
% add probe depths, types, and mutli color if several regions of interest
% plot only probe in the region of interest?
cl_myPaths;
bregma = [540, 0, 570];
theseAnimals = {'JF058', 'JF059', 'JF088', 'JF089', 'JF101', 'JF108', 'JF109', 'JF110'};%{'JF091', 'JF093', 'JF107', 'JF105'};%{'JF058', 'JF059', 'JF088', 'JF089', 'JF101', 'JF108', 'JF109', 'JF110'};{'JF096', 'JF097', 'JF099', 'JF106'};%
%theseAnimals = {'JF093', 'JF091', 'JF107', 'JF104', 'JF105'};
animalsType = {'Naive'};
regionsNames = {'CP', 'GPe', 'GPi', 'STN', 'SNr'};
regions = {'DMS', 'GPe', 'GPi', 'STN', 'SNr'};
regionPlotLoc = [-1, 1, -1, 1, -1];
%recordingInfo = readtable('/home/julie/Dropbox/Analysis/RecNew.csv');
% all animals probes
cl_myPaths;
allen_atlas_path = '/home/julie/Dropbox/Atlas/allenCCF';
tv = readNPY([allen_atlas_path, filesep, 'template_volume_10um.npy']);
av = readNPY([allen_atlas_path, filesep, 'annotation_volume_10um_by_index.npy']);
st = loadStructureTreeJF([allen_atlas_path, filesep, 'structure_tree_safe_2017.csv']);
    slice_spacing = 10;

[~, ~, st_br, ~] = bd_loadAllenAtlas(atlasBrainRegLocation);
n_str = 0;
n_gpe =0;
n_snr =0;
theseColors =     {[          0 0.7461 1];...
    [0.1797 0.5430 0.3398];...  
    [0.8594 0.0781 0.2344];...
    [     1 0.4102 0.7031];...
    [          1 0.5469 0];...
    [               0 0 0];...
    [0.6445 0.1641 0.1641]};
%add probe types, depths
for iType = 1:size(animalsType, 2)
    %theseTypes = strcmp(recordingInfo.Type, animalsType{iType});
    %theseColors = {rgb('DeepSkyBlue'); rgb('DarkOrange'); rgb('Hotpink'); rgb('SeaGreen'); rgb('Crimson')};

    slice_spacing = 10;
    structure_alpha = 0.2;
    %get colors (overrride allen)
    %figure();
    if patchBrain
        bd_plotBrainSurface(allenAtlasPath)
    else
        [~, brain_outline] = plotBrainGrid([], []);
    end

    %overlay regions
    for iRegion = [1,2,5]
        curr_plot_structure = find(strcmp(st.acronym, regionsNames{iRegion}));
        if regionPlotLoc(iRegion) == -1
            structure_3d = isosurface(permute(av(1:slice_spacing:end, ...
                1:slice_spacing:end, 1:slice_spacing:(1140 / 2)) == curr_plot_structure, [3, 1, 2]), 0);
        else
            structure_3d = isosurface(permute(av(1:slice_spacing:end, ...
                1:slice_spacing:end, (1140 / 2)+1:slice_spacing:end) == curr_plot_structure, [3, 1, 2]), 0);
            structure_3d.vertices(:, 2) = structure_3d.vertices(:, 2) + (1140 / 2 / slice_spacing);
        end
        %midpoint only

        if strcmp(regionsNames{iRegion}, 'STR') == 0
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
                'FaceColor', theseColors{iRegion, :}, 'EdgeColor', 'none', 'FaceAlpha', structure_alpha);
        end
        %plot probe tracks
        %theseProbes = ones((size(recordingInfo.Location, 1)),1);
        %get animal and probe, load track
        %theseAnimals = recordingInfo.Mouse(theseTypes);
        qqqq = unique(theseAnimals);
        for iAnimal = 1:size(theseAnimals, 2)
            %iAnimal = iAnimal + 1;
            probe_ccf_location = AP_cortexlab_filenameJF(theseAnimals{iAnimal}, [], [], 'histo');
            load(probe_ccf_location)

            %thisAnimal = strcmp(recordingInfo.Mouse, qqqq{iAnimal});
            %ttP = (recordingInfo.HistologyProbe( thisAnimal & theseProbes));
            %thisProbe = recordingInfo.HistologyProbe(find(thisAnimal & theseProbes));
            for iProbe = 1:size(probe_ccf, 1)
                


                %iProbe = iProbe +1
                curr_probe = iProbe;
                %animalIdx=find(thisAnimal);
                %ex=recordingInfo.Exclude(animalIdx(iProbe));
                %if  isempty(ex{:})

               struct_curr_br = st_br.id(strcmp(st_br.acronym, regionsNames{iRegion}));
                if (any(probe_ccf(iProbe).trajectory_areas == struct_curr_br)) %|| (iAnimal == 2 && ismember(iProbe,ttP))
                    if iRegion == 1
                    n_str = n_str+1;
                elseif iRegion == 2
                    n_gpe = n_gpe+1;
                else
                    n_snr = n_snr+1;
                end
                    thesePoints = probe_ccf(curr_probe).points * 2.5 ; % QQ 2.5 to correct for atlas 25 um where we draw probes and 10 um in plotBrainGrid
                    if regionPlotLoc(iRegion) == -1
                        thesePoints(thesePoints(:,2)>570,2) =  570 - thesePoints(thesePoints(:,2)>570,2); %528? 320? 456?
                    else
                        thesePoints(thesePoints(:,2)<570,2) =  570 + (570 - thesePoints(thesePoints(:,2)<570,2)); %528? 320? 456?
               
                    end
                     r0 = mean(thesePoints, 1);
                    xyz = bsxfun(@minus, thesePoints, r0);
                    [~, ~, V] = svd(xyz, 0);
                    %V= permute(V, [3, 2, 1]);
                    histology_probe_direction = V(:, 1);
                    % (make sure the direction goes down in DV - flip if it's going up)
                    if histology_probe_direction(2) < 0
                        histology_probe_direction = -histology_probe_direction;
                    end

                    line_eval = [-1000, 1000];
                    probe_fit_line = bsxfun(@plus, bsxfun(@times, line_eval', histology_probe_direction'), r0);
                    plot3(thesePoints(:, 1), ...
                        thesePoints(:, 2), ...
                        thesePoints(:, 3), ...
                        '.', 'color', theseColors{iRegion, :}, 'MarkerSize', 20);
                    line(probe_fit_line(:, 1), probe_fit_line(:, 2), probe_fit_line(:, 3), ...
                        'color', theseColors{iRegion, :}, 'linewidth', 2)
                end
            end
            %                 probe_vector = [probe_ref_vector(:, 1), diff(probe_ref_vector, [], 2) ./ ...
            %                     norm(diff(probe_ref_vector, [], 2)) * probe_length + probe_ref_vector(:, 1)];
            %end
        end

    end
end
view([90, 0])
view([90, 90])
view([0, 0])
%save as .avi rotating vid
set(gcf, 'Color', 'w')
OptionZ.FrameRate = 15;
OptionZ.Duration = 5.5;
OptionZ.Periodic = true;
CaptureFigVid([-20, 10; -110, 10; -190, 80; -290, 10; -380, 10], 'WellMadeVid', OptionZ)


