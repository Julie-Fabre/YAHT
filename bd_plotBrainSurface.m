function [currFig, patchHandle] = bd_plotBrainSurface(allenAtlasPath, figureColor, brainColor, patchAlphaValue)

if nargin < 2 || isempty(figureColor)
    figureColor = 'w';
end

if nargin < 3 || isempty(brainColor)
    brainColor = [0.6, 0.6, 0.6];
end

if nargin < 3 || isempty(patchAlphaValue)
    patchAlphaValue = 0.1;
end

av = readNPY([allenAtlasPath, filesep, 'annotation_volume_10um_by_index.npy']);
st = loadStructureTreeJF([allenAtlasPath, filesep, 'structure_tree_safe_2017.csv']);

slice_spacing = 10;

currFig = figure('Color', figureColor);
brain = find(strcmp(st.acronym, 'root'));
structure_3d = isosurface(permute(av(1:slice_spacing:end, ...
    1:slice_spacing:end, 1:slice_spacing:end) == brain, [3, 1, 2]), 0);

patchHandle = patch('Vertices', structure_3d.vertices*slice_spacing, ...
    'Faces', structure_3d.faces, ...
    'FaceColor', brainColor, 'EdgeColor', 'none', 'FaceAlpha', patchAlphaValue);


set(gca, 'ZDir', 'reverse')
axis(gca, 'equal');
axis(gca, 'vis3d');
axis(gca, 'off');

camlight; 
hold on;

end
