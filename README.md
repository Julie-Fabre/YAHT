# allen CCF tools

Some code to work with the Allen Inst Mouse Brain CCF data, specifically the 10µm voxel 2017 version. See the Wiki for detailed instructions.


## Requirements
- A computer mouse with a scroll wheel
- MATLAB (R2017 on a Windows computer was used for all testing)
- This repository (add all folders and subfolders to your MATLAB path)
- [The npy-matlab repository](http://github.com/kwikteam/npy-matlab)
- [The Allen Mouse Brain Atlas volume and annotations](http://data.cortexlab.net/allenCCF/) (download all 4 files from this link)
If you have access, you could also locate the data files at //zserver/Lab/Atlas/allenCCF. Alternatively, see setup_utils to download and preprocess the files yourself. See also https://alleninstitute.github.io/AllenSDK/reference_space.html to access the data directly via the Allen Inst python API.

## Neuropixels trajectory GUI (allen_ccf_npx)
(OLD - new/supported version is here: https://github.com/petersaj/neuropixels_trajectory_explorer)

Plan neuropixels trajectories with an Allen CCF-based GUI

```
>> tv = readNPY('template_volume_10um.npy'); % grey-scale "background signal intensity"
>> av = readNPY('annotation_volume_10um_by_index.npy'); % the number at each pixel labels the area, see note below
>> st = loadStructureTree('structure_tree_safe_2017.csv'); % a table of what all the labels mean

>> allen_ccf_npx(tv,av,st);
```
for 4-shank: 
```
>> allen_ccf_npx_4shank(tv,av,st);
```

## Slice Histology Alignment, Registration, and Probe Track analysis (SHARP-Track)

SHARP-Track is a MATLAB user interface to explore the Allen Mouse Brain Atlas, register asymmetric slice images to the atlas using manual input,
and interactively analyze electrode tracks that span several slices. It can also be used to localize other ROIs such as labelled neurons
or to determine the parameters needed to target particular brain regions with an electrode. All user-oriented scripts can be found in the 'SHARP-Track' folder.

See this repository's wiki for instructions.

*now with scripts to convert  from Allen CCF coordinates to Franklin-Paxinos brain region labels, either through the SHARP-Track pipeline (Convert_Clicked_Points_to_FP_coords.m) or with Allen CCF coordinates in general (Convert_CCF_Coords_to_FP_Regions.m). This uses data from Chon et al. Enhanced and unified anatomical labeling for a common mouse brain atlas (2020).

If you use this tool, please cite [our bioRxiv paper](https://www.biorxiv.org/content/early/2018/10/19/447995).

## To run the Atlas Browser

This is a browser sort of like Paxinos - scroll through coronal slices, see labeled areas, etc.

```
>> tv = readNPY('template_volume_10um.npy'); % grey-scale "background signal intensity"
>> av = readNPY('annotation_volume_10um_by_index.npy'); % the number at each pixel labels the area, see note below
>> st = loadStructureTree('structure_tree_safe_2017.csv'); % a table of what all the labels mean

>> file_save_location = 'C:\Histology\Mouse1'; % where will the probe locations be saved
>> probe_name = 'test'; % name probe to avoid overwriting

>> f = allenAtlasBrowser(tv, av, st, file_save_location, probe_name);
```

## To plot wire mesh of brain
```
>> plotBrainGrid();
```

### Note about annotation volume

The original volume has numbers that correspond to the "id" field in the structure tree, but since I wanted to make a colormap for these, I re-indexed the annotation volume by the row number of the structure tree. So in this version the values correspond to "index"+1. This also allows using uint16 datatype, cutting file size in half. See setup_utils.m.

## Source
© 2015 Allen Institute for Brain Science. Allen Mouse Brain Atlas (2015) with region annotations (2017).
Available from: http://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/

See Allen Mouse Common Coordinate Framework Technical White Paper for details
http://help.brain-map.org/download/attachments/8323525/Mouse_Common_Coordinate_Framework.pdf?version=3&modificationDate=1508178848279&api=v2
