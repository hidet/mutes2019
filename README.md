# Anaysis tools for MUTES project (anaMUTES)

---

# Environment variables
```
export MUTESHOME="$HOME/work/mutes"
export MUTESANADIR="$MUTESHOME/mutes2019"
export MUTESDATADIR="$MUTESHOME/data/TMU_2019G"
```
# Run information and run summary
`RUNINFO` should be updated and set in the `run_mutes_single.py` or `ipy_test.py`
```
RUNINFO="./csv/data_TMU_2019G.csv"
RUNINFO="./csv/data_TMU_2019H.csv"
RUNINFO="./csv/data_TMU_2019I.csv"
```

## --- Updated from 2019 June ---
You can use the options to specify the dataset `TMU_2019X`. You just change  simply as
```
export MUTESDATADIR="$MUTESHOME/data"
```
`RUNINFO` is automatically selected in `run_mutes_single.py` or `ipy_test.py`


# Options
Be careful, the default value is defined as `False` or `None`.
For examples:
 - calibration run (first)
 ```python run_mutes_single.py 76 -R```
 - calibration run (to update)
```python run_mutes_single.py 76 -fsc -R```
 - beam run (first)
```python run_mutes_single.py 58,59 -eg -REG```
 - beam run (to update)
```python run_mutes_single.py 58,59 -fsceg -REG```
 - to delete the hdf5 file
```python run_mutes_single.py 58,59 -egd -REG```
 - to use 'beam off' or 'spill off' categorical cuts for making filter templates and analyzing with more clean pulses 
```python run_mutes_single.py 58,59 -egd -REG --beam=off```
 - to use selections with group trigger data (e.g., sprmc) 
 ```python run_mutes_single.py 58,59 -egd -REG --beam=off --sprmc=on```
 - to specify the dataset `TMU_2019I` without changing the environment variables
```python run_mutes_single.py 1,2 -R --adr=TMU_2019 --cool=I```

 ### options for analysis
```
'-f', '--force',    dest='forceNew',   action='store_true',  help='True to update filter (default=False)'
'-s', '--summary',  dest='summaryNew', action='store_true',  help='True to update summary (default=False)'
'-c', '--calib',    dest='calibNew',   action='store_true',  help='True to update calibration (default=False)'
'-e', '--exttrig',  dest='externTrig', action='store_true',  help='True for calc externTrig (default=False)'
'-g', '--grptrig',  dest='groupTrig',  action='store_true',  help='True for calc groupTrig (default=False)'
'-d', '--delete',   dest='delete',     action='store_true',  help='True to delete hdf5 file (default=False)'
'-R', '--dumproot', dest='dumproot',   action="store_true",  help='dump ROOT except for pulses (default=False)'
'-E', '--rootext',  dest='rootext',    action='store_true',  help='True ROOT with externTrig (default=False)'
'-G', '--rootgrp',  dest='rootgrp',    action='store_true',  help='True ROOT with groupTrig (default=False)'
```

 ### categorical cuts for analyzing average pulse, filter template, drift correction, etc...
```
'--beam',   dest='beam',     action="store",type=str, help='set beam catecut (default=None, on or off)',default="None")
'--sprmc',  dest='sprmc',    action="store",type=str, help='set sprmc catecut (default=None, on or off)',default="None")
'--jbrsc',  dest='jbrsc',    action="store",type=str, help='set jbrsc catecut (default=None, on or off)',default="None")
'--pre',    dest='cut_pre',  action="store",type=int, help='set cut for pre samples',default=0)
'--post',   dest='cut_post', action="store",type=int, help='set cut for post samples',default=0)
```
 ### to specify dataset (from 2019 June)
```
'--adr',    dest='adr',      action="store",type=str, help='set adr tag (default=TMU_2019, TMU_2018, ...)',default="TMU_2019")
'--cool',   dest='cool',     action="store",type=str, help='set cooling tag (default=G, A,B,C...)',default="G")
```
# External trigger timing
The external trigger timing is defined as the time difference between the external triggers and each pulse. The unit is ```row (np.int64 or np.uint64)```, one count is 240 ns in the TMU system (2018-2019). *Watch out overflows when you calculate with* ```row```.

The external trigger data are recorded with an independent client. The file name is usually this like ```runXXXX_extern_trig.hdf5``` . This file is read in ```channel.py``` by ```external_trigger_rowcount()```. [^1]
[^1]: This function `external_trigger_rowcount()` should be run just once by any `ds`. In multiple-run analysis, each numpy array is appended (usually numpy array should not be appended because of slow, but in this case I'd like to keep the dtype to avoid the overflows).


Each pulse has two relative timing to the external triggers,

 1.  ```rows_after_last_external_trigger```
 2.  ```rows_until_next_external_trigger```

The relation of two values is like this;
```
----- ext ----- pulse ------------- ext ---
            ^---after_last   ^---until_next
```
In this case, the pulse is closer to the last external trigger. These are just the integer comparisons to avoid the overflows. The decimal information ```ds.p_rowd``` is obtained by filtering the pulse. To adjust the integer value, there are two options to add ```ds.p_rowp``` or to subtract ```ds.p_rown```.  Unfortunately the sign depends on the version of filter. The added value ```ds.p_rowp``` is created by the newest filter, so you can use it with ```rows_after_last_external_trigger_nrp``` and ```rows_until_next_external_trigger_nrp```.

### New filter or old filter
The new filter parameter is defined in ```mutes_ana.py``` as  ```self.use_new_filters=True```. If you want to try the old filter, just change the value to ```False```. 

Be careful, the timing information derived from old filter has the *opposite sign*. You need to choose the negative-sign values for old-filter analysis in ```mutes_ext.py```, especially you should do manually to calculate ```ds.p_dt``` (time difference with decimal values)
- for new filter
```
pdt[:]   = rows_until_next_external_trigger_nrp - ds.p_rowd
```
- for old filter
```
pdt[:]   = rows_until_next_external_trigger_nrn - ds.p_rowd
```

# Energy calibration with low-energy tail

First of all, the definition of *fitter* for energy calibration was changed from `mass_nov2018`.  It is now in `fluorescence_lines.py`, and the attributes are global.

One of the simplest ways to get a fitter: `fitter = mass.getfitter(linename)`.

- fit parameters
```
NOTE parameters of MultiLorentzianComplexFitter
param_meaning = {
"resolution": 0,
"peak_ph": 1,
"dP_dE": 2,
"amplitude": 3,
"background": 4,
"bg_slope": 5,
"tail_frac": 6,
"tail_length": 7
}
```
Usually, the `dP_dE` is fixed by a guess value at each line. To enable fitting of background and low-energy tail, you need to specify `vary_bg=True` and `vary_tail=True`, respectively. Both `vary_bg` and  `vary_tail` are defined in `mutes_ana.py`, you can change them if you want to switch on/off. [^2] [^3] [^4]

[^2]: Fitting with low-energy tail is necessary to achieve a 0.1-eV class accurate energy calibration. But if it is too slow, you just skip the LE-tail fitting by setting `vary_tail=False`.

[^3]: The optimization procedure is the `MaximumLikelihoodHistogramFitter()`. To fit  with LE-tail, you need to specify boundary conditions for the tail parameters (tail fraction and tail scale length).

[^4]: The original `autocal()` code of mass has no setting of the `vary_tail` parameter. This is why the energy calibration has not been accurate even with choosing the proper X-ray lines.






---


## License

[MIT](http://b4b4r07.mit-license.org)
