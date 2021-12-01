#!/bin/bash


#alp="0.0"
#bet="50e-4"
#path="DIRECT/"
path="env_VERY_SLOW/"
alp="1e-4"
bet="0.0"
gam="1.0"
eps="5.0"
rho="0.1"
N="2000"

for i in $(seq 0 5000);
do
    python simulate.py $path $i $alp $bet $gam $eps $rho $N
done


