function [tv, av, st, bregma] = bd_loadAllenAtlas(atlasLocation)


tv = loadtiff([atlasLocation, filesep, 'reference.tiff']);
av = loadtiff([atlasLocation, filesep, 'annotation.tiff']);
st = readtable([atlasLocation, filesep, 'structures.csv']);
st(size(st,1)+1,1) = {'all'};
st(size(st,1),2) = {0};
st(size(st,1),3) = {'all'};

bregma = [540,0,570];

end