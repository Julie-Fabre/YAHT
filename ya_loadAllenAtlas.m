function [tv, av, st, bregma_ap_dv_ml] = ya_loadAllenAtlas(atlasLocation)

if contains(atlasLocation, 'brainglobe') % 25 um * 25um *25 um
    tv = permute(rot90(loadtiff([atlasLocation, filesep, 'reference.tiff']), 1), [3, 2, 1]);
    av = permute(rot90(loadtiff([atlasLocation, filesep, 'annotation.tiff']), 1), [3, 2, 1]);
    st = readtable([atlasLocation, filesep, 'structures.csv']);
    st(size(st,1)+1,1) = {'all'};
    st(size(st,1),2) = {0};
    st(size(st,1),3) = {'all'};
    bregma_ap_dv_ml = ya_getBregma()./2.5; %divide AP by 2.5 because brainreg is 25um/slice in AP and allen atlas is 10. 
elseif contains(atlasLocation, '_v2')
    tv = NaN; % not implemented yet
    av = readNPY([atlasLocation filesep 'annotation_volume_v2_20um_by_index.npy']); % ap x dv x ml
    st = readtable([atlasLocation, filesep, 'UnifiedAtlas_Label_ontology_v2.csv']);
    bregma_ap_dv_ml = ya_getBregma()./2;
else
    tv = readNPY([atlasLocation filesep 'template_volume_10um.npy']);
    av = readNPY([atlasLocation filesep 'annotation_volume_10um_by_index.npy']); % ap x dv x ml
    st = ya_loadStructureTree([atlasLocation filesep 'structure_tree_safe_2017.csv']);
    bregma_ap_dv_ml = ya_getBregma();
end



end