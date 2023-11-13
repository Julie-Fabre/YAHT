function ya_convertToAPFormat(regIm, tv, av, outputDir)

histology_ccf=struct;
atlas2histology_tform = cell(size(regIm,3),1);
for iSlice = 1:size(regIm,3)
    %histology_ccf(iSlice).plane_ap = 1320 - (1320 - cropAllenLimits(2)) - repmat((iSlice-1)*2, [size(reg,1), size(reg,2)]);%other way? %bregma ?? 
    %histology_ccf(iSlice).plane_ml = repmat(1:2:1140, [size(reg,1),1]);
    %histology_ccf(iSlice).plane_dv = repmat(1:2:800, [ size(reg,2),1])';
    %histology_ccf(iSlice).tv_slices = squeeze(tv((iSlice-1)*2 + cropAllenLimits(1),:,:));
    histology_ccf(iSlice).tv_slices = tv(iSlice, :, :);
    histology_ccf(iSlice).av_slices = av(iSlice, :, :);
    %histology_ccf(iSlice).av_slices = squeeze(av((iSlice-1)*2 + cropAllenLimits(1),:,:));
    atlas2histology_tform{iSlice} = [1 0 0; 0 1 0; 0 0 1]; % no scaling 
end
mkdir([outputDir '/manual'])
save([outputDir '/manual/histology_ccf.mat'], 'histology_ccf','-v7.3')
save([outputDir '/manual/atlas2histology_tform.mat'], 'atlas2histology_tform')
end