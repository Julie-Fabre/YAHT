function [tv, av, st, bregma] = ya_loadAllenAtlas(atlasLocation)

if contains(atlasLocation, 'brainglobe')
    tv = permute(rot90(loadtiff([atlasLocation, filesep, 'reference.tiff']), 1), [3, 2, 1]);
    av = permute(rot90(loadtiff([atlasLocation, filesep, 'annotation.tiff']), 1), [3, 2, 1]);
    st = readtable([atlasLocation, filesep, 'structures.csv']);
    st(size(st,1)+1,1) = {'all'};
    st(size(st,1),2) = {0};
    st(size(st,1),3) = {'all'};
    bregma = [540,0,570] ./2.5;
else
    tv = readNPY([atlasLocation filesep 'template_volume_10um.npy']);
    av = readNPY([atlasLocation filesep 'annotation_volume_10um_by_index.npy']);
    st = ya_loadStructureTree([atlasLocation filesep 'structure_tree_safe_2017.csv']);
    bregma = [540,0,570];
end



end