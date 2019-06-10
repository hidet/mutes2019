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
You can use the options to specify the dataset `TMU_2019X`, for that you just change as simply
```
export MUTESDATADIR="$MUTESHOME/data"
```
`RUNINFO` is automatically selected in `run_mutes_single.py` or `ipy_test.py`


# Options
Be careful, the default value is `False` or `None`.

For examples:
- calibration run (first)
```python run_mutes_single.py 76 -R```
- calibration run (update)
```python run_mutes_single.py 76 -fsc -R```
- beam run (first)
```python run_mutes_single.py 58,59 -eg -REG```
- beam run (update)
```python run_mutes_single.py 58,59 -fsceg -REG```
- to delete the hdf5 file
```python run_mutes_single.py 58,59 -egd -REG```
- to use beam off data for making filter templates and analyzing with more clean pulses 
```python run_mutes_single.py 58,59 -egd -REG --beam=off```
- to use selections with group trigger data (e.g., sprmc)
- ```python run_mutes_single.py 58,59 -egd -REG --beam=off --sprmc=on```
- to specify the dataset `TMU_2019I` without changing the environment variables
```python run_mutes_single.py 1,2 -R --adr=TMU_2019 --cool=I```

options for analysis
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

categorical cuts for analyzing average pulse, filter template, drift correction, etc...
```
'--beam',   dest='beam',     action="store",type=str, help='set beam catecut (default=None, on or off)',default="None")
'--sprmc',  dest='sprmc',    action="store",type=str, help='set sprmc catecut (default=None, on or off)',default="None")
'--jbrsc',  dest='jbrsc',    action="store",type=str, help='set jbrsc catecut (default=None, on or off)',default="None")
'--pre',    dest='cut_pre',  action="store",type=int, help='set cut for pre samples',default=0)
'--post',   dest='cut_post', action="store",type=int, help='set cut for post samples',default=0)
```
to specify dataset (from 2019 June)
```
'--adr',    dest='adr',      action="store",type=str, help='set adr tag (default=TMU_2019, TMU_2018, ...)',default="TMU_2019")
'--cool',   dest='cool',     action="store",type=str, help='set cooling tag (default=G, A,B,C...)',default="G")
```
# External trigger timing
The external trigger timing is defined as the time difference between the external triggers and each pulse. The unit is ```row (np.int64 or np.uint64)``` count which is 240 ns in the TMU system (2018-2019). <watch out overflows when you calculate with ```row```>

Each pulse has two timing values, (1) ```rows_after_last_external_trigger``` and (2) ```rows_until_next_external_trigger```. The relation of two values is like this;
```
----- ext ----- pulse ------------- ext ---
            ^--after_last .   ^--until_next
```
These are just the integer comparisons to avoid the overflows. The decimal information ```ds.p_rowd``` is obtained by filtering the pulse. To adjust the filtered timing, there are two options, add ```ds.p_rowp``` or subtract ```ds.p_rown```. 





# New filter or old filter
The pulse height and timing information are obtained by filtering each pulse with a template (made from average pulse and noise spectrum). There are two methods in mass code, called new or old filter. The new filter looks better than the old one, but it is complicated to follow.

The new filter is used in ```mutes_ana.py``` with the value ```self.use_new_filters=True```. If you want to try the old filter, just change the value to ```False```. 

Be careful, the timing information of old filter has the  opposite sign. You need to choose the negative-sign values for beam timing analysis in ```mutes_ext.py```, especially you should do manually...
- for new filter
```
pdt[:]   = rows_until_next_external_trigger_nrp - ds.p_rowd
```
- for old filter
```
pdt[:]   = rows_until_next_external_trigger_nrn - ds.p_rowd
```



# Energy calibration with low-energy tail













---





## License

[MIT](http://b4b4r07.mit-license.org)
