For LCLS-1 data processing, from the root directory (`CookieBox_pkg`)

# Dark Images
```bash
source init.bash
python3 python/Xtcav_store_dark.py amo86815 58 62 68
```

![Compare Plot](./figure/compare58_58.png)  
Comparing an individual block average image to the run average (mean of block averages) in fun 58.  


![Comparison Plot](./figures/compare68_58.png)  
Comparing a block average in run 68 with the whole run average (mean of block averages) in run 58.  
