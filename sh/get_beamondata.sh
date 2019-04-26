#!/bin/bash

name=spill_info.data

scp -p oper@buyu.intra.j-parc.jp:/xvb602/home/oper/scaler_e62_moge/record_scaler_beam_on.data ../csv/
./trim.rb > ../csv/$name
