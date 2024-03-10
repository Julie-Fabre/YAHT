function ya_plotHistoPerMouse_JF(animal, outputDir, st)
cl_myPaths;
slice_spacing = 25;
structure_alpha = 0.2;
scalingFactor = 2.5;
 old_histo = false;
if ~isempty(animal)
    probe2ephys_location = cl_cortexlab_filename(animal, [], [], 'probe2ephys');
    load(probe2ephys_location)
    probe_ccf_location = cl_cortexlab_filename(animal, [], [], 'histo');
    load(probe_ccf_location)
    if isfield(probe2ephys, 'shank')% new histo
        st = readtable([atlasBrainRegLocation, filesep, 'structures.csv']);
        st(size(st,1)+1,1) = {'all'};
        st(size(st,1),2) = {0};
        st(size(st,1),3) = {'all'};
    else
        old_histo = true;
        slice_spacing = 10;
        st = ya_loadStructureTree([allenAtlasPath, filesep, 'structure_tree_safe_2017.csv']);
        scalingFactor = 1;
    end
else
    load([outputDir, '/probe_ccf.mat'])


    if ~isempty(dir([outputDir, '/probe2ephys.mat']))
        load([outputDir, '/probe2ephys.mat'])
    else
        probe2ephys = struct;
    end
    if ~isfield(probe2ephys, 'shank')% new histo
        old_histo = true;
        slice_spacing = 10;
        scalingFactor = 1;
    end
end
cmap_filename = [allenAtlasPath, filesep, '/allen_ccf_colormap_2017.mat'];
%
figure();
load(cmap_filename);
for curr_probe = 1:length(probe_ccf)
    curr_axes = subplot(1, length(probe_ccf), curr_probe);

    trajectory_area_boundaries = ...
        [1; find(diff(probe_ccf(curr_probe).trajectory_areas) ~= 0); length(probe_ccf(curr_probe).trajectory_areas)];
    trajectory_area_centers = trajectory_area_boundaries(1:end-1) + diff(trajectory_area_boundaries) / 2;
    if size(trajectory_area_centers, 1) > 1
        for iArea = 1:size(trajectory_area_centers, 1)
            if old_histo
                trajectory_area_labels(iArea) = st.acronym(...
                    probe_ccf(curr_probe).trajectory_areas(round(trajectory_area_centers(iArea))),:);
            else
                trajectory_area_labels(iArea) = st.acronym(st.id == ...
                    probe_ccf(curr_probe).trajectory_areas(round(trajectory_area_centers(iArea))));
            end
        end

        %trajectory_area_labels = st.acronym(probe_ccf(curr_probe).trajectory_areas(round(trajectory_area_centers)));

        image(probe_ccf(curr_probe).trajectory_areas);
        colormap(curr_axes, cmap);
        caxis([1, size(cmap, 1)])
        set(curr_axes, 'YTick', trajectory_area_centers, 'YTickLabels', trajectory_area_labels);
        set(curr_axes, 'XTick', []);
    end
    title(num2str(curr_probe));
    %    if any( structfun(@isempty, probe2ephys) )
    %       title([num2str(probe2ephys(curr_probe).day) num2str(probe2ephys(curr_probe).site)]);
    %  end

end

[~, brain_outline] = plotBrainGrid([], []);
hold on;
theseColors = {rgb('DeepSkyBlue'); rgb('DarkOrange'); rgb('Hotpink'); rgb('SeaGreen'); rgb('Crimson')};
regionsNames = {'CP', 'GPe', 'GPi', 'STN', 'SNr'};

colsThese = lines(length(probe_ccf));

for curr_probe = 1:length(probe_ccf)


        thesePoints = probe_ccf(curr_probe).points * scalingFactor; % QQ 2.5 to correct for atlas 25 um where we draw probes and 10 um in plotBrainGrid
        
        
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
       if old_histo
             plot3(probe_ccf(curr_probe).points(:,1), ...
        probe_ccf(curr_probe).points(:,3), ...
        probe_ccf(curr_probe).points(:,2), ...
        '.','color',colsThese(curr_probe, :),'MarkerSize',20);
    line(probe_fit_line(:,1),probe_fit_line(:,3),probe_fit_line(:,2), ...
        'color',colsThese(curr_probe, :),'linewidth',2)

       else
             plot3(probe_ccf(curr_probe).points(:,1), ...
        probe_ccf(curr_probe).points(:,2), ...
        probe_ccf(curr_probe).points(:,3), ...
        '.','color',colsThese(curr_probe, :),'MarkerSize',20);
    line(probe_fit_line(:,1),probe_fit_line(:,3),probe_fit_line(:,2), ...
        'color',colsThese(curr_probe, :),'linewidth',2)

       end
        

        v_1 = [probe_fit_line(1, 2), probe_fit_line(1, 3)] - [probe_fit_line(2, 2), probe_fit_line(2, 3)];
        v_2 = [probe_fit_line(1, 2), probe_fit_line(1, 3)] - [probe_fit_line(2, 2), probe_fit_line(2, 2)];
        line_angle(curr_probe, 1) = dot(v_1, v_2) ./ (norm(v_1) * norm(v_2)) * 180 / pi;


        v_1 = [probe_fit_line(1, 2), probe_fit_line(1, 1)] - [probe_fit_line(2, 2), probe_fit_line(2, 1)];
        v_2 = [probe_fit_line(1, 2), probe_fit_line(1, 2)] - [probe_fit_line(2, 2), probe_fit_line(2, 1)];
        line_angle(curr_probe, 2) = dot(v_1, v_2) ./ (norm(v_1) * norm(v_2)) * 180 / pi;


end

% figure();
% subplot(311)
% scatter(1:length(probe_ccf), line_angle(:, 1), 'filled')
% xlabel('probe #')
% ylabel('ap angle')
% 
% subplot(312)
% scatter(1:length(probe_ccf), line_angle(:, 2), 'filled')
% xlabel('probe #')
% ylabel('ml angle')
% 
% subplot(313)
% scatter3(1:length(probe_ccf), line_angle(:, 1), line_angle(:, 2), 'filled')
% xlabel('probe #')
% ylabel('ap angle')
% zlabel('ml angle')

% figure();
% for curr_probe = 1:length(probe_ccf)
% 
%     thesePoints = probe_ccf(curr_probe).points * scalingFactor; % QQ 2.5 to correct for atlas 25 um where we draw probes and 10 um in plotBrainGrid
%     r0 = mean(thesePoints, 1);
%     xyz = bsxfun(@minus, thesePoints, r0);
%     [~, ~, V] = svd(xyz, 0);
%     %V= permute(V, [3, 2, 1]);
%     histology_probe_direction = V(:, 1);
%     % (make sure the direction goes down in DV - flip if it's going up)
%     if histology_probe_direction(2) < 0
%         histology_probe_direction = -histology_probe_direction;
%     end
% 
%     line_eval = [-1000, 1000];
%     probe_fit_line = bsxfun(@plus, bsxfun(@times, line_eval', histology_probe_direction'), r0);
% 
%     subplot(2, length(probe_ccf), curr_probe)
%     %plot(thesePoints(:, 1), ...
%     %   thesePoints(:, 3), ...
%     %  '.', 'color', colsThese(curr_probe, :), 'MarkerSize', 20);
% 
%     line(probe_fit_line(:, 1), probe_fit_line(:, 3), ...
%         'color', colsThese(curr_probe, :), 'linewidth', 2)
%     title(num2str(curr_probe))
%     set(gca, 'xticklabel', {[]})
%     set(gca, 'xtick', [])
%     set(gca, 'yticklabel', {[]})
%     set(gca, 'ytick', [])
%     box off;
% 
%     subplot(2, length(probe_ccf), curr_probe+length(probe_ccf))
%     plot(thesePoints(:, 1), ...
%         thesePoints(:, 2), ...
%         '.', 'color', colsThese(curr_probe, :), 'MarkerSize', 20);
% 
%     line(probe_fit_line(:, 1), probe_fit_line(:, 2), ...
%         'color', colsThese(curr_probe, :), 'linewidth', 2)
%     set(gca, 'xticklabel', {[]})
%     set(gca, 'xtick', [])
%     set(gca, 'yticklabel', {[]})
%     set(gca, 'ytick', [])
%     box off;
% 
% end
% prettify_plot('YLimits', 'all', 'XLimits', 'all')
% 
% 
% lineStyles = {'-', '--', ':', '-.'};
% lineTypes = [repmat(1, 7, 1); repmat(2, 7, 1); repmat(3, 7, 1); repmat(4, 7, 1); repmat(1, 7, 1); repmat(2, 7, 1); repmat(3, 7, 1); repmat(4, 7, 1)];
%
% figure();
% for curr_probe = 1:length(probe_ccf)
%     %try
%
%
%         thesePoints = probe_ccf(curr_probe).points * 2.5; % QQ 2.5 to correct for atlas 25 um where we draw probes and 10 um in plotBrainGrid
%         r0 = mean(thesePoints, 1);
%         xyz = bsxfun(@minus, thesePoints, r0);
%         [~, ~, V] = svd(xyz, 0);
%         %V= permute(V, [3, 2, 1]);
%         histology_probe_direction = V(:, 1);
%         % (make sure the direction goes down in DV - flip if it's going up)
%         if histology_probe_direction(2) < 0
%             histology_probe_direction = -histology_probe_direction;
%         end
%
%         line_eval = [-1000, 1000];
%         probe_fit_line = bsxfun(@plus, bsxfun(@times, line_eval', histology_probe_direction'), r0);
%
%         subplot(3,2,[1,2]);hold on;
%         line(-probe_fit_line(:, 1), -probe_fit_line(:, 3), ...
%             'color', colsThese(curr_probe, :), 'linewidth', 2, 'LineStyle', lineStyles{lineTypes(curr_probe)});
%
%         subplot(3,2,[3,4]);hold on;
%         line(probe_fit_line(:, 2), probe_fit_line(:, 3), ...
%             'color', colsThese(curr_probe, :), 'linewidth', 2, 'LineStyle', lineStyles{lineTypes(curr_probe)});
%
%         subplot(3,2,5);hold on;
%          scatter(probe_fit_line(1, 2), -probe_fit_line(1, 1), 20,...
%              colsThese(curr_probe, :));
%         text(probe_fit_line(1, 2), -probe_fit_line(1, 1),{num2str(curr_probe)},'FontSize',20,'Color', colsThese(curr_probe, :))
%         xlabel('ML')
%         ylabel('AP')
%
%         subplot(3,2,6);hold on;
%         scatter(-probe_points_start_position(curr_probe, 2), -probe_points_start_position(curr_probe, 1), 20,...
%              colsThese(curr_probe, :));
%         text(-probe_points_start_position(curr_probe, 2), -probe_points_start_position(curr_probe, 1),{num2str(curr_probe)},'FontSize',20,'Color', colsThese(curr_probe, :))
%         xlabel('ML')
%         ylabel('AP')
%
%     %catch
%     %end
%     subplot(3,2,[1,2])
%     legend;
% end
%
%
