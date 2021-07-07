python3 main.py --protA WT --protB R164S \
                --directory "test" \
                --n_sample 10 \
                --struct_file_A "data/WT/md_0_1.gro" \
                --traj_file_A "data/WT/md_0_1_noPBC.xtc" \
                --struct_file_B "data/R164S/md_0_1.gro" \
                --traj_file_B "data/R164S/md_0_1_noPBC.xtc" \
                --selection "protein and resid 65:230" \
                --radius 2.0 \
                --n_cone 1 \
                --n_direction_per_cone 1 \
                --cap_radius 0.80 \
                --ec_type "DECT" \
                --n_filtration 20 \
                --n_mcmc 10000 \
                --parallel --n_core 4 --verbose
