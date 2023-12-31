import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim

from _negative_binomial import ZeroInflatedNegativeBinomial

from dataclasses import dataclass, field
import numpy as np
import pdb

@dataclass
class ModelConfig:

    """
    Model configuration for crispAI models
        - convolutional layers: conv1d
        - pooling layers: maxpool1d
        - linear layers: linear
        - activation functions: relu
        -
    """

    # input sequence 
    seq_len: int = 23
    seq_encoding_channels: int = 6
    pi_encoding_channels: int = 4
    pi_features: bool = True # whether to use physical features
    channels: int = seq_encoding_channels
    channels += pi_encoding_channels if pi_features else 0
    multitask: bool = False # TODO: add multitask loss

    # conv, use default_factory to avoid mutable default arguments
    conv_seq_configs: list = field(default_factory=lambda: [(128, 3, 1), (32, 1, 1)]) # (out_channels, kernel_size, stride)
    conv_pi_configs: list = field(default_factory=lambda: [(16, 4, 1)]) # (out_channels, kernel_size, stride)   
    conv_dropout: float = 0.1
    conv_batchnorm: bool = True
    #conv_pool: list = field(default_factory=lambda: [(2, 2, 0), (2, 2, 0)]) # (kernel_size, stride, padding)
    conv_pool: list = None
    conv_activation: str = 'relu'

    # lstm
    lstm_bidirectional: bool = True
    lstm_hidden: int = 128
    lstm_dropout: float = 0.1

    # dense
    sequence_dense: list = field(default_factory=lambda: [128, 32])
    pi_dense: list = field(default_factory=lambda: [32])
    output_dense: list = field(default_factory=lambda: [64, 3]) # [concat_dim, ..., (mu, theta, pi)]


class CrispAI_pi(nn.Module):

    def __init__(self, config):
        super().__init__()
        assert config.pi_features, 'pi_features must be True for CrispAI_pi'
        self.pi_features = config.pi_features

        self.config = config
        self.seq_len = config.seq_len
        self.channels_input = config.channels
        self.channels = 0 # temp variable for conv layers
        self.seq_encoding_channels = config.seq_encoding_channels
        self.pi_encoding_channels = config.pi_encoding_channels
        self.multitask = config.multitask
        self.pi_encoding_channels = config.pi_encoding_channels
        self.eps = 1e-6

        # conv_seq layers
        self.conv_seq_layers = nn.ModuleList()
        temp_channels = self.seq_encoding_channels

        for i, (out_channels, kernel_size, stride) in enumerate(config.conv_seq_configs):
            self.conv_seq_layers.append(nn.Conv1d(
                                    in_channels=temp_channels,
                                    out_channels=out_channels,
                                    kernel_size=kernel_size,
                                    stride=stride))
            
            
            if config.conv_batchnorm:
                self.conv_seq_layers.append(nn.BatchNorm1d(out_channels))
              
            if config.conv_pool is not None:
                self.conv_seq_layers.append(nn.MaxPool1d(kernel_size=config.conv_pool[i][0],
                                                        stride=config.conv_pool[i][1],
                                                        padding=config.conv_pool[i][2]))
                                                 
            self.conv_seq_layers.append(nn.Dropout(config.conv_dropout))
            self.conv_seq_layers.append(nn.ReLU())
            temp_channels = out_channels
            conv_seq_out_channels = out_channels

        # conv_pi layers
        self.conv_pi_layers = nn.ModuleList()
        temp_channels = self.pi_encoding_channels

        for out_channels, kernel_size, stride in config.conv_pi_configs:
            self.conv_pi_layers.append(nn.Conv1d(
                                    in_channels=temp_channels,
                                    out_channels=out_channels,
                                    kernel_size=kernel_size,
                                    stride=stride))
            
            if config.conv_batchnorm:
                self.conv_pi_layers.append(nn.BatchNorm1d(out_channels))
            
            if config.conv_pool is not None:
                self.conv_pi_layers.append(nn.MaxPool1d(kernel_size=config.conv_pool[i][0],
                                                            stride=config.conv_pool[i][1],
                                                            padding=config.conv_pool[i][2]))
                                                        
            self.conv_pi_layers.append(nn.ReLU())
            temp_channels = out_channels
            conv_pi_out_channels = out_channels

        # lstm
        self.lstm = nn.LSTM(input_size=conv_seq_out_channels,
                            hidden_size=config.lstm_hidden,
                            batch_first=True,
                            bidirectional=config.lstm_bidirectional,
                            dropout=config.lstm_dropout)
        
        # dense layers
        self.sequence_dense = nn.ModuleList()
        temp_channels = config.lstm_hidden * 2 if config.lstm_bidirectional else config.lstm_hidden
        for out_features in config.sequence_dense:
            self.sequence_dense.append(nn.Linear(in_features=temp_channels,
                                                 out_features=out_features))
            self.sequence_dense.append(nn.ReLU())
            temp_channels = out_features
            seq_dense_out_channels = out_features
        
        self.pi_dense = nn.ModuleList()
        temp_channels = conv_pi_out_channels * (self.seq_len - config.conv_pi_configs[0][1] + 1) 
        for out_features in config.pi_dense:
            self.pi_dense.append(nn.Linear(in_features=temp_channels,
                                                 out_features=out_features))
            self.pi_dense.append(nn.ReLU())
            temp_channels = out_features
            pi_dense_out_channels = out_features
        
        
        self.output_dense = nn.ModuleList()
        temp_channels = seq_dense_out_channels + pi_dense_out_channels
        for i, out_features in enumerate(config.output_dense):
            self.output_dense.append(nn.Linear(in_features=temp_channels,
                                                 out_features=out_features))
            if i != len(config.output_dense) - 1:
                self.output_dense.append(nn.ReLU())
            
            temp_channels = out_features
            output_dense_out_channels = out_features
            

    def forward(self, x):

        '''
        Forward pass of crispAI model
            - x: input sequence (batch_size, seq_len, channels)
            - last self.pi_encoding_channels channels are physical features
            - output: (batch_size, output_dense[-1])
            - sequence features: stacked conv1d + lstm
            - pi features: conv1d
            - fusion features: concat(sequence features, pi features)
        '''

        x_seq = x[:, :, :-self.pi_encoding_channels]
        x_pi = x[:, :, -self.pi_encoding_channels:]

        x_seq = x_seq.permute(0, 2, 1)
        x_pi = x_pi.permute(0, 2, 1)
        
        for layer in self.conv_seq_layers:
            x_seq = layer(x_seq)

        x_seq = x_seq.view(x_seq.shape[0], -1, x_seq.shape[1])
        _, (hn, _) = self.lstm(x_seq)

        hn = hn.permute(1, 0, 2).contiguous().view(hn.shape[1], -1)

        for layer in self.sequence_dense:
            hn = layer(hn)

        for layer in self.conv_pi_layers:
            x_pi = layer(x_pi)

        x_pi = x_pi.view(x_pi.shape[0], -1)

        for layer in self.pi_dense:
            x_pi = layer(x_pi)

        x_cat = torch.cat((hn, x_pi), dim=1)

        for layer in self.output_dense:
            x_cat = layer(x_cat)

        x_mu = torch.exp(x_cat[:, 0]).view(-1, 1)
        x_theta = torch.exp(x_cat[:, 1]).view(-1, 1)
        x_pi = x_cat[:, 2].view(-1, 1)

        x = torch.cat((x_mu, x_theta, x_pi), dim=1)

        return x
    
    def draw_samples(self, x, n_samples = 100):
        '''
        Draw samples from the predicted distribution
        '''
        preds_ = self(x) # forward pass to get lambda
        preds_ = preds_.cpu().detach().numpy()

        mu_ = preds_[:, 0]
        theta_ = preds_[:, 1]
        pi_ = preds_[:, 2]
        mu_ = torch.tensor(mu_, dtype=torch.float32)
        theta_ = torch.tensor(theta_, dtype=torch.float32)
        pi_ = torch.tensor(pi_, dtype=torch.float32)
        zinb = ZeroInflatedNegativeBinomial(mu = mu_, theta = theta_, zi_logits = pi_)
        x_samples = zinb.sample_n(n_samples)
        x_samples = x_samples.cpu().detach().numpy()
        
        return x_samples

    def get_num_params(self):
        return sum(p.numel() for p in self.parameters() if p.requires_grad)