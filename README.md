For LCLS-1 data processing, from the root directory (`CookieBox_pkg`)  

# DataSelector
Now there is a DataSelector that composes a training set of up to 512 images that represent a leveled representation of the random-slice-Wasserstein version of nearest neighbor distribution.  


# Bright Images
The runscript is in `./runscript.bash` and it sources the `init.bash` and launches the `./python/Xtcav_preproc.py` with darkpath(input), brightpath(output), and expname and bright run numbers.  
The slurmscript allows batch submitting a list of run numbers as follows...   
```bash
./slurmscript.bash /sdf/data/lcls/ds/amo/amoi0216/scratch /sdf/data/lcls/ds/amo/amoi0216/scratch amoi0216 17 18 19 20 21 22
```


# Dark Images
New slurm way to compute the dark images, also storing in single run per darkimage.h5 file.  
```bash
./slurmdark.bash amo86815 13 18 26 40 46 58 62 68
```

This is now used only if running directly from an interactive node.   
```bash
source init.bash
python3 python/Xtcav_store_dark.py amo86815 58 62 68
```
# Training set composition  
Now, for choosing wether a new image is added to the training set or discarded, I'm using ```./python/randomSliceWasserstein.py``` like this...   
```bash
python3 ./python/randomSliceWasserstien.py /media/coffee/9C33-6BBD/temp_xtcav/amo86815/xtcav_bright_images_69.h5 /media/coffee/9C33-6BBD/temp_xtcav/amo86815/xtcav_bright_images_69.h5 
```
and this measures random selection from one file with a rendom selection for either the same or a different file.  
The so called "leveling" that I'm doing is to choose whether to keep the sample based on its mean random slice rmse with the comparator.  
In the end, this should be a search through the existing training set and for the minimum rmse and then if that is appropriate, then replace the most well represented sample.


![Compare Plot](./figures/compare58_58.png)  
Comparing an individual block average image to the run average (mean of block averages) in fun 58.   


![Comparison Plot](./figures/compare68_58.png)   
Comparing a block average in run 68 with the whole run average (mean of block averages) in run 58.   

This indicates that the camera baseline indeed changes over the course of runs at best, and over the course of shots within a run at worst.  

