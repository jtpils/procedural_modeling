#!/bin/bash

logfile="log."

for d in 10 25 50 100 200 500 1000
  do
    :
    python opt.py $d > $logfile$d
  done
