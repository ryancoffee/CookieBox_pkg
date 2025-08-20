For LCLS-1 data processing, from the root directory (`CookieBox_pkg`)  

# DataSelector
Now there is a DataSelector that composes a training set of up to 512 images that represent a leveled representation of the random-slice-Wasserstein version of nearest neighbor distribution.  


# Bright Images
The runscript is in `./runscript.bash` and it sources the `init.bash` and launches the `./python/Xtcav_preproc.py` with darkpath(input), brightpath(output), and expname and bright run numbers.  
The slurmscript allows batch submitting a list of run numbers as follows...   
```bash
./slurmscript.bash /sdf/data/lcls/ds/amo/amoi0216/scratch /sdf/data/lcls/ds/amo/amoi0216/scratch amoi0216 17 18 19 20 21 22
```

THe alternative for running on directly on the `.xtc` file segments with Xtcav only...   
```bash
for r in 3 7 5 4; do ./slurmscript_xtcOnly.bash /sdf/data/lcls/ds/amo/amoe7615/scratch /sdf/data/lcls/ds/amo/amoe7615/scratch amoe7615 /sdf/data/lcls/ds/amo/amoe7615/xtc/e569-r$(printf "%04d" $r)-s80-c*.xtc;done
for r in 240 241 236 234 229 225 213 212 203 201 182 140 130 125 126 127 128 103 106 97; do ./slurmscript_xtcOnly.bash /sdf/data/lcls/ds/amo/amoh5215/scratch /sdf/data/lcls/ds/amo/amoh5215/scratch amoh5215 /sdf/data/lcls/ds/amo/amoh5215/xtc/e576-r$(printf "%04d" $r)-s80-c*.xtc;done
```


# Dark Images
New way to run directly on the .xtc files of only xtcav.  
```bash
for r in 13 18 26 40 46 58 62 68; do ./slurmdark_xtcOnly.bash amo86815 /sdf/data/lcls/ds/amo/amo86815/xtc/e609-r$(printf "%04d" $r)-s80-c*.xtc;done
for r in 15 19 21 28 41 51 59 70; do ./slurmdark_xtcOnly.bash amoi0216 /sdf/data/lcls/ds/amo/amoi0216/xtc/e826-r$(printf "%04d" $r)-s80-c*.xtc;done
```


Old slurm way to compute the dark images, also storing in single run per darkimage.h5 file.  
```bash
./slurmdark.bash amo86815 13 18 26 40 46 58 62 68
```

This is now used only if running directly from an interactive node.   
```bash
source init.bash
python3 python/Xtcav_store_dark.py amo86815 58 62 68
```


# Training set composition  
Now, for choosing whether a new image is added to the training set or discarded, I'm using ```./python/randomSliceWasserstein.py``` like this...   
```bash
python3 ./python/Xtcav_composeTranin.py /media/coffee/9C33-6BBD/temp_xtcav/amo86815/xtcav_bright_images_*.h5 
```
and this measures random selection from one file with a rendom selection for either the same or a different file.  
The so called "leveling" that I'm doing is to choose whether to keep the sample based on its mean random slice rmse with the comparator.  
In the end, this should be a search through the existing training set and for the minimum rmse and then if that is appropriate, then replace the most well represented sample.

Once the training sets have been gneerated, then we can use them to construct an set of eigenfunctions for use in image decomposition.  
Left to check is that the bright file computation `rolls` the centroid into the middle of the image and then the decomposition and eigenfunction calculation only happens for the 512x512 central crop with a douwnsample of `::4,::4` for both dimensions... so 128x128 images in the end.

```bash
python3 python/Xtcav_trainEigenimages.py /nvme0/temp /nvme0/amo86815/xtcav_bright_images_58_train.h5 /nvme0/amoi0216/xtcav_bright_images_run68_train.h5
```

![Compare Plot](./figures/compare58_58.png)  
Comparing an individual block average image to the run average (mean of block averages) in fun 58.   


![Comparison Plot](./figures/compare68_58.png)   
Comparing a block average in run 68 with the whole run average (mean of block averages) in run 58.   

This indicates that the camera baseline indeed changes over the course of runs at best, and over the course of shots within a run at worst.  

