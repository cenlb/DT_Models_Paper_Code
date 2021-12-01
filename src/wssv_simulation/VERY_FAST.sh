#!/bin/bash


#alp="0.0"
#bet="50e-4"
#path="DIRECT/"
path="env_VERY_FAST/"
alp="1e-4"
bet="0.0"
gam="1.0"
eps="50000.0"
rho="1000.0"
delta_t="0.00001"
t_max="60.0"
N="2000"

for i in $(seq 0 5000);
do
    echo $i
    ./SIR_P_sim $path $i $alp $bet $gam $eps $rho $delta_t $t_max $N
done


