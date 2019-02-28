import torch
import torch.nn as nn

__all__ = [
    'RNN',
]

def fc_relu(in_features, out_features, inplace=True):
    return nn.Sequential(
        nn.Linear(in_features, out_features),
        nn.ReLU(inplace=inplace),
        nn.Dropout(p=0.1),
    )

class RNN(nn.Module):
    def __init__(self):
        super(RNN, self).__init__()

        self.hsize = 512

        self.air_rnn_0 = nn.GRUCell(self.hsize, self.hsize)
        self.bed_rnn_0 = nn.GRUCell(self.hsize, self.hsize)
        self.air_rnn_1 = nn.GRUCell(self.hsize, self.hsize)
        self.bed_rnn_1 = nn.GRUCell(self.hsize, self.hsize)

        self.air_fc_in_0 = fc_relu(64, self.hsize)
        self.bed_fc_in_0 = fc_relu(64, self.hsize)
        self.air_fc_in_1 = fc_relu(64, self.hsize)
        self.bed_fc_in_1 = fc_relu(64, self.hsize)
        self.air_fc_out = nn.Linear(self.hsize, 1)
        self.bed_fc_out = nn.Linear(self.hsize, 1)

    def forward(self, data, init):
        air_output = None
        bed_output = None
        air_hidden = [[data.new_zeros((data.shape[0], self.hsize)) for i in range(65)] for j in range(2)]
        bed_hidden = [[data.new_zeros((data.shape[0], self.hsize)) for i in range(65)] for j in range(2)]
        air_hidden[0][0] = init
        bed_hidden[0][0] = init
        air_hidden[1][0] = init
        bed_hidden[1][0] = init

        for i in range(64):
            air_input_0 = self.air_fc_in_0(data[:,:,i])
            bed_input_0 = self.bed_fc_in_0(data[:,:,i])
            air_input_1 = self.air_fc_in_1(data[:,:,63-i])
            bed_input_1 = self.bed_fc_in_1(data[:,:,63-i])

            air_hidden[0][i+1] = self.air_rnn_0(air_input_0, air_hidden[0][i])
            bed_hidden[0][i+1] = self.bed_rnn_0(bed_input_0, bed_hidden[0][i])
            air_hidden[1][i+1] = self.air_rnn_1(air_input_1, air_hidden[1][i])
            bed_hidden[1][i+1] = self.bed_rnn_1(bed_input_1, bed_hidden[1][i])

        for i in range(1, 65):
            air_temp = self.air_fc_out(air_hidden[0][i]+air_hidden[1][65-i])
            bed_temp = self.bed_fc_out(bed_hidden[0][i]+bed_hidden[1][65-i])

            air_output = air_temp if i ==1 else torch.cat((air_output, air_temp), 1)
            bed_output = bed_temp if i ==1 else torch.cat((bed_output, bed_temp), 1)

        return air_output, bed_output
