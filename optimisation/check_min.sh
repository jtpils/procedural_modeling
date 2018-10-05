#!/bin/bash
while true;
do awk 'NR ==1 {line=$0; min=$12}
     NR > 1 && $12 < min {line=$0; min=$12}
     END{print line}' func.dat > current_opt.txt;
sleep 60;
done;

  
