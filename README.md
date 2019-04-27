# anaysis tools for MUTES project (anaMUTES)

---

# environment variables
```
export MUTESHOME="$HOME/work/mutes"
export MUTESANADIR="$MUTESHOME/mutes2019"
export MUTESDATADIR="$MUTESHOME/data/TMU_2019G"
```

# options
Be careful, the default value is `False` or `None`.

For examples:
- calibration run (first)
```python run_mutes_single.py 54 -R```
- calibration run (forceNew)
```python run_mutes_single.py 54 -fsc -R```
- beam run (first)
```python run_mutes_single.py 58,59 -eg -REG```
- beam run (forceNew)
```python run_mutes_single.py 58,59 -fsceg -REG```

```
'-f', '--force',    dest='forceNew',   action='store_true',  help='True to update filter (default=False)'
'-s', '--summary',  dest='summaryNew', action='store_true',  help='True to update summary (default=False)'
'-c', '--calib',    dest='calibNew',   action='store_true',  help='True to update calibration (default=False)'
'-e', '--exttrig',  dest='externTrig', action='store_true',  help='True for calc externTrig (default=False)'
'-g', '--grptrig',  dest='groupTrig',  action='store_true',  help='True for calc groupTrig (default=False)'
'-R', '--dumproot', dest='dumproot',   action="store_true",  help='dump ROOT except for pulses (default=False)'
'-E', '--rootext',  dest='rootext',    action='store_true',  help='True ROOT with externTrig (default=False)'
'-G', '--rootgrp',  dest='rootgrp',    action='store_true',  help='True ROOT with groupTrig (default=False)'
```

- categorical cut for averge pulse, template, drift correction, etc...
```
'--beam',   dest='beam',     action="store",type=str, help='set beam catecut (default=None, on or off)',default="None")
'--sprmc',  dest='sprmc',    action="store",type=str, help='set sprmc catecut (default=None, on or off)',default="None")
'--jbrsc',  dest='jbrsc',    action="store",type=str, help='set jbrsc catecut (default=None, on or off)',default="None")
'--pre',    dest='cut_pre',  action="store",type=int, help='set cut for pre samples',default=0)
'--post',   dest='cut_post', action="store",type=int, help='set cut for post samples',default=0)
```


---

## License

[MIT](http://b4b4r07.mit-license.org)
