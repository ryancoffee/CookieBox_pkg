For LCLS-1 data processing, from the root directory (`CookieBox_pkg`)

# Dark Images
```bash
source init.bash
python3 python/Xtcav_store_dark.py amo86815 58 62 68
```

![Compare Plot](./figures/compare58_58.png)  
Comparing an individual block average image to the run average (mean of block averages) in fun 58.  


![Comparison Plot](./figures/compare68_58.png)  
Comparing a block average in run 68 with the whole run average (mean of block averages) in run 58.  

This indicates that the camera baseline indeed changes over the course of runs at best, and over the course of shots within a run at worst.  


# Bright Images
The runscript is in `./runscript.bash` and it sources the `init.bash` and launches the `./python/Xtcav_preproc.py` with dark and expname and runnum.  
The slurmscript allows batch submitting a list of run numbers with the associated dark and expname as follows...   
```bash
./slurmscript.bash /sdf/data/lcls/ds/amo/amo86815/scratch/xtcav_dark_images_68_62_58.h5 amo86815 69 59 61 63 70 71 72 73
```
