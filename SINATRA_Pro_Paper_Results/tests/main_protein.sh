#!/bin/python3

python3 ../src/main.py --protA WT --protB R164S \
                --directory "WT_R164S_65_230_2.0" \
                --single \
                --n_sample 10 \
                --offset 0 \
                --struct_file_A "../data/wt.gro" \
                --traj_file_A "../data/wt.xtc" \
                --struct_file_B "../data/r164s.gro" \
                --traj_file_B "../data/r164s.xtc" \
                --selection "protein and resid 65:230" \
                --radius 2.0 \
                --n_cone 1 \
                --n_direction_per_cone 1 \
                --cap_radius 0.80 \
                --ec_type "DECT" \
                --n_filtration 20 \
                --n_mcmc 10000 \
                --parallel --n_core 4 --verbose

