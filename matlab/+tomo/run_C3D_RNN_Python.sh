#!/bin/bash
module unload python/2.7.9
module unload gcc/4.9.3
module load gcc/7.4.0
module load anaconda3/5.2.0
source activate tensorflow_env
unset MKL_NUM_THREADS
#python3 convert_mat_to_npy.py --mat_dir $1 --npy_dir $2
#echo "Done converting MAT file to NPY."
python3 run_C3D_RNN.py --data_root $3 --c3d_pth $4 --rnn_pth $5 --output_dir $6 --frms $7