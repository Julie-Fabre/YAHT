# YAHT: Yet Another Histology Tool

YAHT: Yet Another Histology Tool. The goal of this toolbox is to have an all-in-one pipeline to quickly go from histology to highly accurately-tracked probes.

YAHT allows for:
- preprocessing of histology images to get higher quality data,
- automated histology slice registration to a common atlas (using [brainreg](https://github.com/brainglobe/brainreg)),
- manual fine-tuning of the registration to get perfect fits,
- manual drawing of probes in an intuitive and easy to adjust way,  using adjustable bezier curves (we find that this provides the best fit),
- probe to ephys data alignement, and
- generation of high quality 3D probe location summary plots and videos. 

Largely based on [AP_histology](https://github.com/petersaj/AP_histology) and [allenCCF](https://github.com/cortex-lab/allenCCF).

### Installation

YAHT requires: 
- brainreg to be installed locally 
- the files located [here](https://osf.io/fv7ed/#!) to be downloaded locally

Dependencies: 
- 
### Quick start

Adjust the paths in `bd_histologyMain_template` an run through the sections to get started.

### Pipeline steps
#### Step 0: preprocessing
#### Step 1: register your histology using brain reg

#### Step 2: Correct any orientation issues

<img style="float: right;" src="https://github.com/Julie-Fabre/braindraw/blob/master/image/altasOrientationFix.png" width=100% height=100%>

#### Step 3: Fine tune the stretch factors on a slice-by-slice basis 

#### Step 4: Draw probes

![atlasDrawProbes](image/atlasDrawProbes.png)


#### Step 5: Assign probes to ephys days/sites


#### Step 6: Align probes and ephys

#### Step 7: Generate some beautiful plots and videos 




