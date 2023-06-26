from __future__ import print_function
import os
import os.path as osp
import sys

import numpy as np
from scipy.io import loadmat

sys.path.insert(0, '../')
import config as cfg

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--mat_dir', default='slices_mat_64x64', type=str)
    parser.add_argument('--npy_dir', default='slices_npy_64x64', type=str)
    args = cfg.parse_args(parser)
    args.mat_path = osp.join(args.data_root, args.mat_dir)
    args.npy_path = osp.join(args.data_root, args.npy_dir)

    for date in sorted(os.listdir(args.mat_path)):
        date_mat_path = osp.join(args.mat_path, date)
        date_npy_path = osp.join(args.npy_path, date)

        for session in sorted(os.listdir(date_mat_path)):
            session_mat_path = osp.join(date_mat_path, session)
            session_npy_path = osp.join(date_npy_path, session)
            if not osp.isdir(session_npy_path):
                os.makedirs(session_npy_path)

            for file in sorted(os.listdir(session_mat_path)):
                file_mat_path = osp.join(session_mat_path, file)
                file_npy_path = osp.join(session_npy_path, file.replace('.mat', '.npy'))
                data = loadmat(file_mat_path)
                data = data['fusion'].astype(np.float32)
                np.save(file_npy_path, data)
