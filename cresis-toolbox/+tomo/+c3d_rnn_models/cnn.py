import torch
import torch.nn as nn

__all__ = ['C3D']

class Flatten(nn.Module):
    def __init__(self):
        super(Flatten, self).__init__()

    def forward(self, x):
        return x.view(x.shape[0], -1)

class C3D(nn.Module):
    def __init__(self):
        super(C3D, self).__init__()

        self.conv1_s = nn.Sequential(
            nn.Conv3d(1, 16, (3, 5, 3), stride=1, padding=(1, 2, 1)),
            nn.BatchNorm3d(16),
            nn.ReLU(inplace=True),
            nn.MaxPool3d((1, 2, 2)),
        )

        self.conv2_s = nn.Sequential(
            nn.Conv3d(16, 32, (3, 5, 3), stride=1, padding=(1, 2, 1)),
            nn.BatchNorm3d(32),
            nn.ReLU(inplace=True),
            nn.MaxPool3d((1, 2, 2)),
        )

        self.conv3_a = nn.Sequential(
            nn.Conv3d(32, 64, (3, 5, 3), stride=1, padding=(1, 2, 1)),
            nn.BatchNorm3d(64),
            nn.ReLU(inplace=True),
            nn.Conv3d(64, 64, (3, 5, 3), stride=1, padding=(1, 2, 1)),
            nn.BatchNorm3d(64),
            nn.ReLU(inplace=True),
            nn.Conv3d(64, 64, (3, 5, 3), stride=1, padding=(1, 2, 1)),
            nn.BatchNorm3d(64),
            nn.ReLU(inplace=True),
            nn.MaxPool3d((1, 2, 2)),
        )

        self.conv3_b = nn.Sequential(
            nn.Conv3d(32, 64, (3, 5, 3), stride=1, padding=(1, 2, 1)),
            nn.BatchNorm3d(64),
            nn.ReLU(inplace=True),
            nn.Conv3d(64, 64, (3, 5, 3), stride=1, padding=(1, 2, 1)),
            nn.BatchNorm3d(64),
            nn.ReLU(inplace=True),
            nn.Conv3d(64, 64, (3, 5, 3), stride=1, padding=(1, 2, 1)),
            nn.BatchNorm3d(64),
            nn.ReLU(inplace=True),
            nn.MaxPool3d((1, 2, 2)),
        )

        self.conv4_a = nn.Sequential(
            nn.Conv3d(64, 128, (3, 5, 3), stride=1, padding=(1, 2, 1)),
            nn.BatchNorm3d(128),
            nn.ReLU(inplace=True),
            nn.Conv3d(128, 128, (3, 5, 3), stride=1, padding=(1, 2, 1)),
            nn.BatchNorm3d(128),
            nn.ReLU(inplace=True),
            nn.Conv3d(128, 128, (3, 5, 3), stride=1, padding=(1, 2, 1)),
            nn.BatchNorm3d(128),
            nn.ReLU(inplace=True),
            nn.MaxPool3d((5, 2, 2)),
        )

        self.conv4_b = nn.Sequential(
            nn.Conv3d(64, 128, (3, 5, 3), stride=1, padding=(1, 2, 1)),
            nn.BatchNorm3d(128),
            nn.ReLU(inplace=True),
            nn.Conv3d(128, 128, (3, 5, 3), stride=1, padding=(1, 2, 1)),
            nn.BatchNorm3d(128),
            nn.ReLU(inplace=True),
            nn.Conv3d(128, 128, (3, 5, 3), stride=1, padding=(1, 2, 1)),
            nn.BatchNorm3d(128),
            nn.ReLU(inplace=True),
            nn.MaxPool3d((5, 2, 2)),
        )

        self.conv5_a = nn.Sequential(
            nn.Conv3d(128, 256, (1, 4, 4), stride=1, padding=(0, 0, 0)),
            nn.ReLU(inplace=True),
            Flatten(),
        )

        self.conv5_b = nn.Sequential(
            nn.Conv3d(128, 256, (1, 4, 4), stride=1, padding=(0, 0, 0)),
            nn.ReLU(inplace=True),
            Flatten(),
        )

        self.fc6_a = nn.Linear(256, 64)
        self.fc6_b = nn.Linear(256, 64)

    def features(self, img):
        img = self.conv1_s(img)
        img = self.conv2_s(img)
        air = self.conv3_a(img)
        bed = self.conv3_b(img)
        air = self.conv4_a(air)
        bed = self.conv4_b(bed)
        air = self.conv5_a(air)
        bed = self.conv5_b(bed)
        return air, bed

    def forward(self, img):
        air, bed = self.features(img)
        air = self.fc6_a(air)
        bed = self.fc6_b(bed)
        return air, bed
