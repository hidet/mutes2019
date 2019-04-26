# anaysis tools for HEATES project (anaHEATES)


## Description

This is a software packege for E62 HEATES calorimater project. 
It utilizes mass (offline TES analysis package developed by NIST) 
with HEATES-specific modification. 

## Features

- create hdf5 by running mass 
- calculate energy resolution and plot spectra
- merge histogram among runs 
- plot pulse records for a given condition 

## Requirement

- python 2.7 environment 
- mass version 1.0.0 (mass_new) 
-- pip >= version 10 fails with "No module named pip.req". Then downgrade, python -m pip install pip==9.0.3
- pyROOT if you wish to analyze the data in ROOT

## Documantation 

The detail of the software will be documanted in wiki. 

- https://bitbucket.org/heates/ana_2018_jparc/wiki/Home

In this page, only a brief explanation is shown.  

## Usage

### 1. run standard analysis for each run 

    python run_heates_single.py RUNNUM

### 2. merge histgram among runs

    python run_heates_across.py

### 3.plot pulse records

    python run_heates_single.py RUNNUM --pulse_ana --chans=1,7,11

### 4. dump ROOT files

    python run_heates_single.py RUNNUM --dumproot 

It dumps all except for pulse records. 

When you want to dump pulse record, use dumprootpulse option.   

    python run_heates_single.py RUNNUM --dumprootpulse 


### 5. do analysis for each pixel 

    python run_heates_single.py RUNNUM --doiteach

It outputs the spectrum of each pixel with fitting results, and distribution of energy resolution.


### 6. process all production runs of He3 and He4

    ./do_run_heates_single.sh         

The run number of He3 and He4 production are stored in He3.runs and He4.runs. 

### 7. merge hist across runs 

    python run_heates_across.py run20180709
    ....
    ./output/He3_rows_merged_nr8_155_166_run20180709.root is created.

The run number of He3 and He4 production are read from He3.runs and He4.runs. 

    python run_heates_across.py run20180709 -t He4

The option of -t He4 will change He3 to He4 runs. If you add -c, it will merget HDF5 with cut when they are created. 



## Installation

install mass

    $ git clone https://heates@bitbucket.org/heates/mass_new.git

install anaHEATES

    $ git clone https://heates@bitbucket.org/heates/ana_2018_jparc.git

## setup


Define three environmental variables. For example, 

    export HEATESHOME="$HOME/work/ana/HEATES/JPARC201807_E62"
    export HEATESANADIR="$HEATESHOME/ana_2018_jparc"
    export HEATESDATADIR="$HEATESHOME/data/TMU_2018U"

Please change according to your environment. 

Using "lndir" in your analysis directory is recommended. Here is an example in TMU WS. 

    cd $HEATESHOME
    git clone https://heates@bitbucket.org/heates/ana_2018_jparc.git
    mkdir -p data/TMU_2018U                                          
    cd data/TMU_2018U                                          
    lndir /a/heates/data/TMU_2018U .

Set csv files according to run summary. 

- https://docs.google.com/spreadsheets/d/1Ag0jlDUSc5Vn60SxQCBCEc0hCAmHm3_05RbWBxR5cqM/edit#gid=1272501364

If you notice some errors in the log file, please let us know. 


## Things to watch out

- When hdf5 is opened in ipython or somewhere else, mass using the file can not be run.
- When beam is ON, the external trigger sometimes stops. It may hinder beam-off analysis at low photon statistics. 

## The scripts actively maintained

- run_heates_single.py (2018/07/07, v1.2)
- khe_util.py          (2018/07/07, v1.2)
- khe_ana.py           (2018/07/07, v1.2)
- run_heates_across.py (not yet but soon)
- root/heates_root_dumpulse.py (2018/07/07, v1.0)

## History

- 2018.06-07 initial development under E62 run 
- 2018.07.04 compiled from scratch, S.Yamada

## License

[MIT](http://b4b4r07.mit-license.org)
