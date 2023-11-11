function bd_plotAllProbeTracksInROIs(atlasBrainRegLocation, paths, regionNames, regionColors, regionPlotSide, makeMovie, saveMoviePath)

if nargin < 2 || isempty(paths)
    paths = {};
end

if nargin < 3 || isempty(regionNames)
    regionNames = [];
end

if nargin < 4 || isempty(regionColors)
    regionColors = bd_getColors(length(regionNames));
end

if nargin < 5 || isempty(regionPlotSide)
    regionPlotSide = repmat([-1, 1], ceil(length(regionNames)/2));
end

if nargin < 6 || isempty(makeMovie)
    makeMovie = 1;
end

if nargin < 7 || isempty(saveMoviePath)
    saveMoviePath = pwd;
end

allenAtlasPath = fileparts(matlab.desktop.editor.getActiveFilename);

tv = readNPY([allenAtlasPath, filesep, 'template_volume_10um.npy']);
av = readNPY([allenAtlasPath, filesep, 'annotation_volume_10um_by_index.npy']);
st = loadStructureTreeJF([allenAtlasPath, filesep, 'structure_tree_safe_2017.csv']);
bregma = [540, 0, 570];
slice_spacing = 10;
structure_alpha = 0.2;

[~, brain_outline] = plotBrainGrid([], []);

if ~isempty(regionNames)
    [~, ~, st_br, ~] = bd_loadAllenAtlas(atlasBrainRegLocation);
    for iRegion = 1:length(regionNames)
        curr_plot_structure = find(strcmp(st.acronym, regionNames{iRegion}));
        if regionPlotSide(iRegion) == -1
            structure_3d = isosurface(permute(av(1:slice_spacing:end, ...
                1:slice_spacing:end, 1:slice_spacing:(1140 / 2)) == curr_plot_structure, [3, 1, 2]), 0);
        else
            structure_3d = isosurface(permute(av(1:slice_spacing:end, ...
                1:slice_spacing:end, (1140 / 2)+1:slice_spacing:end) == curr_plot_structure, [3, 1, 2]), 0);
            structure_3d.vertices(:, 2) = structure_3d.vertices(:, 2) + (1140 / 2 / slice_spacing);
        end


        if strcmp(regionNames{iRegion}, 'STR') == 0
            hold on;
            axis vis3d equal off manual
            view([-30, 25]);
            clim([0, 300]);
            [ap_max, dv_max, ml_max] = size(tv);
            xlim([-10, ap_max + 10])
            ylim([-10, ml_max + 10])
            zlim([-10, dv_max + 10])
            patch('Vertices', structure_3d.vertices*slice_spacing, ...
                'Faces', structure_3d.faces, ...
                'FaceColor', regionColors(iRegion, :), 'EdgeColor', 'none', 'FaceAlpha', structure_alpha);
        end

        if ~isempty(paths)
            for iAnimal = 1:length(paths, 2)

                load(paths{iAnimal})

                for iProbe = 1:size(probe_ccf, 1)

                    curr_probe = iProbe;
                    thesePoints = probe_ccf(curr_probe).points * 2.5; % QQ 2.5 to correct for atlas 25 um where we draw probes and 10 um in plotBrainGrid
                    struct_curr_br = st_br.id(strcmp(st_br.acronym, regionNames{iRegion}));
                    if (any(probe_ccf(iProbe).trajectory_areas == struct_curr_br))

                        if regionPlotSide(iRegion) == -1
                            thesePoints(thesePoints(:, 2) > 570, 2) = 570 - thesePoints(thesePoints(:, 2) > 570, 2);
                        else
                            thesePoints(thesePoints(:, 2) < 570, 2) = 570 + (570 - thesePoints(thesePoints(:, 2) < 570, 2));

                        end
                        r0 = mean(thesePoints, 1);
                        xyz = bsxfun(@minus, thesePoints, r0);
                        [~, ~, V] = svd(xyz, 0);

                        histology_probe_direction = V(:, 1);

                        if histology_probe_direction(2) < 0
                            histology_probe_direction = -histology_probe_direction;
                        end

                        line_eval = [-1000, 1000];
                        probe_fit_line = bsxfun(@plus, bsxfun(@times, line_eval', histology_probe_direction'), r0);
                        plot3(thesePoints(:, 1), ...
                            thesePoints(:, 2), ...
                            thesePoints(:, 3), ...
                            '.', 'color', regionColors(iRegion, :), 'MarkerSize', 20);
                        line(probe_fit_line(:, 1), probe_fit_line(:, 2), probe_fit_line(:, 3), ...
                            'color', regionColors(iRegion, :), 'linewidth', 2)
                    end
                end
            end
        end
    end
else
    probeColors = bd_getColors(1000); % some large number
    if ~isempty(paths)
        for iAnimal = 1:length(paths, 2)

            load(paths{iAnimal})
            for iProbe = 1:size(probe_ccf, 1)

                curr_probe = iProbe;
                thesePoints = probe_ccf(curr_probe).points * 2.5; % QQ 2.5 to correct for atlas 25 um where we draw probes and 10 um in plotBrainGrid
                r0 = mean(thesePoints, 1);
                xyz = bsxfun(@minus, thesePoints, r0);
                [~, ~, V] = svd(xyz, 0);

                histology_probe_direction = V(:, 1);

                if histology_probe_direction(2) < 0
                    histology_probe_direction = -histology_probe_direction;
                end

                line_eval = [-1000, 1000];
                probe_fit_line = bsxfun(@plus, bsxfun(@times, line_eval', histology_probe_direction'), r0);
                plot3(thesePoints(:, 1), ...
                    thesePoints(:, 2), ...
                    thesePoints(:, 3), ...
                    '.', 'color', probeColors(iProbe, :), 'MarkerSize', 20);
                line(probe_fit_line(:, 1), probe_fit_line(:, 2), probe_fit_line(:, 3), ...
                    'color', probeColors(iProbe, :), 'linewidth', 2)
            end
        end
    end


end


view([90, 0])
view([90, 90])
view([0, 0])

if makeMovie
    %save as .avi rotating vid
    set(gcf, 'Color', 'w')
    OptionZ.FrameRate = 15;
    OptionZ.Duration = 5.5;
    OptionZ.Periodic = true;
    CaptureFigVid([-20, 10; -110, 10; -190, 80; -290, 10; -380, 10], [saveMoviePath, filesep, 'probeTracks'], OptionZ)
end