from __future__ import print_function
import os
import os.path as osp
import sys
import time
import argparse

import torch
import torch.nn as nn
import torch.utils.data as data
import numpy as np
import numpy.linalg as LA
from scipy.io import loadmat

from rnn import *
from cnn import *
import config as cfg

def default_loader(path, number):
    c3d_data = []
    for index in range(number-2, number+3):
        index = min(max(index, 1), 3332)
        mat = loadmat(osp.join(path, str(index).zfill(5)+'.mat'))
        c3d_data.append(mat['fusion'].astype(np.float32))
    c3d_data = np.array(c3d_data, dtype=np.float32)
    c3d_data = (c3d_data-0.5)/0.5
    c3d_data = c3d_data[np.newaxis, ...]
    rnn_data = c3d_data[0,2,:].copy()
    norm = LA.norm(rnn_data, axis=0)
    rnn_data /= norm[None, :]
    return c3d_data, rnn_data

class DataLayer(data.Dataset):
    def __init__(self, data_root, sessions, loader=default_loader):
        self.data_root = data_root
        self.sessions = sessions
        self.loader = loader

        self.inputs = []
        for session_name in self.sessions:
            session_path = osp.join(self.data_root, 'target', session_name+'.txt')
            session_data = open(session_path, 'r').read().splitlines()
            self.inputs.extend(session_data)

    def __getitem__(self, index):
        data_path, number, air_target, bed_target = self.inputs[index].split()
        c3d_data, rnn_data = self.loader(osp.join(
            self.data_root, 'slices_mat_64x64', data_path), int(number))
        c3d_data, rnn_data = torch.from_numpy(c3d_data), torch.from_numpy(rnn_data)
        air_target = np.array(air_target.split(','), dtype=np.float32)
        air_target = torch.from_numpy(air_target)
        bed_target = np.array(bed_target.split(','), dtype=np.float32)
        bed_target = torch.from_numpy(bed_target)

        return c3d_data, rnn_data, air_target, bed_target, data_path

    def __len__(self):
        return len(self.inputs)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gpu', default='0', type=str)
    parser.add_argument('--batch_size', default=1, type=int)
    parser.add_argument('--num_workers', default=0, type=int)
    parser.add_argument('--c3d_pth', default='pretrained_models/c3d.pth', type=str)
    parser.add_argument('--rnn_pth', default='pretrained_models/rnn.pth', type=str)
    
    # Victor: add command-line arguments for frames to be run
    parser.add_argument('--frms', nargs='+')
    # Victor: add output directory for results file
    parser.add_argument('--output_dir', default='NN_TempOutput.txt', type=str);
    
    args = cfg.parse_args(parser)

    os.environ['CUDA_VISIBLE_DEVICES'] = args.gpu
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    # Victor: open results file to write output to
    file = open(args.output_dir, "w");

    data_set = DataLayer(
        data_root=args.data_root,
        #sessions=args.test_session_set,
        sessions = args.frms,
    )

    data_loader = data.DataLoader(
        data_set,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
    )

    c3d_model = C3D().to(device)
    c3d_model.load_state_dict(torch.load(args.c3d_pth))
    c3d_model.train(False)
    rnn_model = RNN().to(device)
    rnn_model.load_state_dict(torch.load(args.rnn_pth))
    rnn_model.train(False)

    air_criterion = nn.L1Loss().to(device)
    bed_criterion = nn.L1Loss().to(device)
    air_errors = 0.0
    bed_errors = 0.0
    
    frm_init = True;
    frm_name = '';

    start = time.time()
    with torch.set_grad_enabled(False):
        for batch_idx, (c3d_data, rnn_data, air_target, bed_target, data_path) in enumerate(data_loader):
            c3d_data = c3d_data.to(device)
            rnn_data = rnn_data.to(device)
            air_target = air_target.to(device)
            bed_target = bed_target.to(device)

            air_feature, bed_feature = c3d_model.features(c3d_data)

            init = torch.cat((air_feature, bed_feature), 1)
            air_output, bed_output = rnn_model(rnn_data, init)

            # Save these air and bed layer for visualization
            air_layer = (air_output.to('cpu').numpy()+1)*412
            bed_layer = (bed_output.to('cpu').numpy()+1)*412
            
            # Victor: define frame name
            data_path_mod = data_path[0].replace("/", "_");

            # Victor: write results for surface and bottom layers to text file
            file.write(data_path_mod + '\n');
            for idx in range(0, air_layer.shape[1]):
                file.write(str(air_layer[0, idx]) + ' ')
            file.write('\n');
            for idx in range(0, bed_layer.shape[1]):
                file.write(str(bed_layer[0, idx]) + ' ')
            file.write('\n\n\n');

            air_loss = air_criterion(air_output, air_target)
            bed_loss = bed_criterion(bed_output, bed_target)
            
            air_errors += air_loss.item()
            bed_errors += bed_loss.item()
    end = time.time()
    
    # Victor: close file reader
    file.close();

    print('Finish all, errors (air): {:4.2f} (bed): {:4.2f}, '
          'total running time: {:.2f} sec'.format(
              air_errors/len(data_loader.dataset)*412,
              bed_errors/len(data_loader.dataset)*412,
              end-start,
          ))
          
if __name__ == '__main__':
    main()
