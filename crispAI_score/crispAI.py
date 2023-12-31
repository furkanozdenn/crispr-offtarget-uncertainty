'''
Source code for crispAI module

Use cases:
1:  Off-target cleavage activity for sgRNA-target pair
    - Sample posterior off-target activity distribution N times for a given sgRNA-target pair 
    - Output either point estimate or distribution of predicted off-target activity

2:  Uncertainty aware genome-wide sgRNA specificty prediction 
    - Sample aggregate score distribution N times for a given sgRNA
    - Output either point estimate or distribution of predicted aggregate score

Input for 1:
    - .csv file where columns  sgRNA_sequence,target_sequence,chr,start,end,strand
    - N_samples: number of samples to draw from posterior distribution (default = 0)
    - P: point estimate of off-target activity (default = 1)
    - O: output file name (default = crispAI_output.csv)

Input for 2:
    - .csv file where columns are sgRNAs (required)
    - N_samples: number of samples to draw from posterior distribution (default = 0)
    - P: point estimate of aggregate score (default = 1)
    - O: output file name (default = crispAI_aggregate_output.csv)
    - N_mismatch: for aggregate score - genome-wide search for off-target sites up-to N_mismatch mismatches (default = 4, Max 6)
'''

import argparse
import pandas as pd
import numpy as np
import pickle
import os
import pdb
import subprocess

import datetime

def message(msg):
    #print(f'{datetime.datetime.now()}: {msg}')
    # h-m-s format
    print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: {msg}')


'''
Arguments for CLI:
    - Mode: default 'offt-score', other option is 'agg-score'
    - Input file name: .csv file with sgRNA-target pairs or sgRNAs
    - Number of samples to draw from posterior distribution: default 0
    - if 'agg-score' mode, number of mismatches to search for off-target sites: default 4
    - Point estimate of off-target activity: default 1
    - Output file name: default 'crispAI_output.csv' or 'crispAI_aggregate_output.csv'
'''

parser = argparse.ArgumentParser(description='crispAI: Quantifying uncertainty in off-target activity for CRISPR guide RNAs')
parser.add_argument('--mode', type=str, default='offt-score', help='Mode: offt-score or agg-score')
parser.add_argument('--input_file', type=str, default='input.csv', help='Input file name')
parser.add_argument('--N_samples', type=int, default=1000, help='Number of samples to draw from posterior distribution [100, 2000]')
parser.add_argument('--N_mismatch', type=int, default=4, help='Number of mismatches to search for off-target sites')
parser.add_argument('--O', type=str, default='crispAI_output.csv', help='Output file name')
# optional GPU support, default is CPU if this argument is not provided. Provide cuda number if GPU is available
parser.add_argument('--gpu', type=int, default=-1, help='CUDA device number for GPU support (default: -1 for CPU)')
# plot aggregate score distribution for sgRNAs (optional flag, if provided, plot the distribution)
parser.add_argument('--plot-agg', action='store_true', help='Plot aggregate score distribution for sgRNAs')



args = parser.parse_args()

# check if input file exists
if not os.path.exists(args.input_file):
    raise ValueError('Input file does not exist')

# check if output file exists and warn user
if os.path.exists(args.O):
    message(f'Output file {args.O} already exists. It will be overwritten.')

# check if mode is either offt-score or agg-score
if args.mode not in ['offt-score', 'agg-score']:
    raise ValueError('Mode should be either offt-score or agg-score')

# check if N_samples is valid, min 0 max 2000
if args.N_samples < 100 or args.N_samples > 2000:
    raise ValueError('N_samples should be between 100 and 2000')

# check if N_mismatch is valid, min 0 max 6
if args.N_mismatch < 0 or args.N_mismatch > 6:
    raise ValueError('N_mismatch should be between 0 and 6')


# check if GPU is valid, should be -1 or a valid cuda number
if args.gpu < -1 or args.gpu > 7:
    raise ValueError('GPU should be between -1 and 7')

# if no output file name is provided, use default
# if offt-mode use crispAI_output.csv
if args.mode == 'offt-score' and args.O == 'crispAI_output.csv':
    args.O = 'crispAI_output.csv'
# if agg-mode use crispAI_aggregate_output.csv
elif args.mode == 'agg-score' and args.O == 'crispAI_output.csv':
    args.O = 'crispAI_aggregate_output.csv'
else:
    pass



'''
Obtain off-target sites to run the model on
'''

# read the input file
if args.mode == 'offt-score':
    message(f'Running crispAI in offt-score mode')
    offtscore_input_data = pd.read_csv(args.input_file, header=None)
    offtscore_input_data.columns = ['sgRNA_sequence', 'target_sequence', 'chr', 'start', 'end', 'strand']
    message(f'Input file {args.input_file} read successfully')
    offtscore_input_data['reads'] = None
    offtscore_input_data['mismatch'] = None
    offtarget_sites = offtscore_input_data

else: # elif args.mode == 'agg-score':
    message(f'Running crispAI in agg-score mode')
    aggscore_input_data = pd.read_csv(args.input_file, header=None)
    aggscore_input_data.columns = ['sgRNA']
    message(f'Input file {args.input_file} read successfully')
    CAS_OFFINDER_PATH = './casoffinder'
    UCSC_CHROMS = f'{CAS_OFFINDER_PATH}/ucsc_chroms'
    SEARCH_SEQUENCE = 'N' * 21 + 'GG'
    DEFAULT_DEVICE = args.gpu # defaults to -1 for CPU

    casoffinder_input_path = f'{CAS_OFFINDER_PATH}/casoffinder_input.txt'   
    with open(casoffinder_input_path, 'w') as f:
        f.write(f'{UCSC_CHROMS}\n')
        f.write(f'{SEARCH_SEQUENCE}\n')
        for i, row in aggscore_input_data.iterrows():
            f.write(f'{row.sgRNA}\t{args.N_mismatch}\n')

    num_sgRNAs = len(aggscore_input_data)
    message(f'Running Cas-OFFinder to find off-target sites up-to {args.N_mismatch} mismatches for {num_sgRNAs} sgRNAs. This may take a while.')
    # run Cas-OFFinder
    casoffinder_output_path = f'{CAS_OFFINDER_PATH}/casoffinder_output.txt'

    if args.gpu == -1: # CPU, CasOFFinder does not see CPU's #TODO: fix this
        os.system(f'{CAS_OFFINDER_PATH}/cas-offinder {casoffinder_input_path} C {casoffinder_output_path} > /dev/null 2>&1')
    else: # GPU
        os.system(f'{CAS_OFFINDER_PATH}/cas-offinder {casoffinder_input_path} G[{args.gpu}] {casoffinder_output_path} > /dev/null 2>&1')

    # read the output file
    casoffinder_output_df = pd.read_csv(casoffinder_output_path, sep='\t', header=None)
    casoffinder_output_df.columns = ['sgRNA_sequence', 'chr', 'start', 'target_sequence', 'strand', 'mismatch']
    # off-target_sequence to upper
    casoffinder_output_df['target_sequence'] = casoffinder_output_df['target_sequence'].str.upper()
    # end column is start + 22 (0-based)
    casoffinder_output_df['end'] = casoffinder_output_df['start'] + 22
    # reads column is 'None' (CasOffinder, in-silico) 
    casoffinder_output_df['reads'] = 'None'
    num_obtained_offtargets = len(casoffinder_output_df)
    message(f'Obtained {num_obtained_offtargets} off-target sites')


    offtarget_sites = casoffinder_output_df



'''
Annotate offtarget_sites with physical features
'''
from annotate_pi import annotation_pipeline

annotated_data = annotation_pipeline(offtarget_sites)


'''
Predict the distribution parameters for each sgRNA-target pair
'''
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F

from torch.utils.data import TensorDataset, DataLoader

from utils import preprocess_features
from model import CrispAI, CrispAI_pi, ModelConfig

# load config and model from checkpoint
device = f'cuda:{args.gpu}' if args.gpu >= 0 else 'cpu'
# message device usage
message(f'Using device: {device}')
checkpoint = torch.load('./model_checkpoint/epoch:19-best_valid_loss:0.270.pt')
config = checkpoint['config']
model = CrispAI_pi(config)
model.load_state_dict(checkpoint['model_state_dict'])
model.to(device)
model.eval()

offtarget_data = preprocess_features(df = annotated_data,
                            reads = 'reads',
                            target = 'sgRNA_sequence',
                            offtarget_sequence= 'target_sequence',
                            distance= 'mismatch',
                            read_cutoff= 10,
                            max_reads = 1e4,
                            nupop_occupancy_col= 'NuPoP occupancy',
                            nupop_affinity_col= 'NuPoP affinity',
                            gc_content_col= 'GC flank73',
                            nucleotide_bdm_col= 'nucleotide BDM')

X = np.stack([x.astype(np.float32) for x in offtarget_data['interface_encoding'].values], axis=0)
X_pi = np.stack([x.astype(np.float32) for x in offtarget_data['physical_features'].values], axis=0)

# concat X and X_pi if pi_features
X = np.concatenate([X, X_pi], axis=2)

dataset = TensorDataset(torch.tensor(X))
preds = []
from tqdm import tqdm
pbar = tqdm(DataLoader(dataset, batch_size=128, shuffle=False), desc='Predicting')
with torch.no_grad():
    for x in pbar:
        x = x[0]
        x = x.to(device)
        samples = model.draw_samples(x, n_samples= args.N_samples).T
        preds.append(samples)

preds = np.concatenate(preds, axis=0)
offtarget_data['preds'] = list(preds)

# if offt-score mode, output desired results and exit
if args.mode == 'offt-score':
    print_data = offtarget_data[['target_N', 'chr', 'start', 'end', 'strand', 'target_sequence', 'preds']]
    print_data['mean'] = print_data['preds'].apply(lambda x: np.mean(x))
    print_data['std'] = print_data['preds'].apply(lambda x: np.std(x))

    with open(args.O, 'w') as f:
        f.write(f'target_N\tchr\tstart\tend\tstrand\ttarget_sequence\tmean\tsamples\tstd\n')
        for i, row in print_data.iterrows():
            f.write(f'{row.target_N}\t{row.chr}\t{row.start}\t{row.end}\t{row.strand}\t{row.target_sequence}\t{row["mean"]}\t{",".join(map(lambda x: str(int(x)), row.preds))}\t{row["std"]}\n')


    message(f'Output file {args.O} written successfully')
    exit()

# else if agg-score mode, continue for aggregate score calculation
message(f'Calculating aggregate score for each sgRNA')
# for each unique guide sequence (offtarget_data['target_N']), sum corresponding preds and take log -> aggregate score
dict = {}
for target_ in offtarget_data['target_N'].unique():
    mismatches = offtarget_data.loc[offtarget_data['target_N'] == target_, 'mismatch'].values
    preds_ = offtarget_data.loc[offtarget_data['target_N'] == target_, 'preds'].values # (n,) arrays each of size (n_samples,)
    preds_ = np.stack(preds_, axis=0)
    preds_[preds_ == 0] = 1
    preds_0 = preds_[mismatches == 0]
    if len(preds_[mismatches == 0]) > 1:
        preds_0 = preds_0[0].reshape(1, -1)
    preds_ = preds_ / preds_0
    dict[target_] = np.log(np.sum(preds_, axis=0))


'''
Write the output file using dict, in the format:
columns: sgRNA, aggregate_score_mean, aggregate_score_median, aggregate_score_std [list of n_samples aggregate scores]
'''

with open(args.O, 'w') as f:
    f.write(f'sgRNA\taggregate_score_mean\taggregate_score_median\taggregate_score_std\t{args.N_samples}-samples\n')
    for sgrna in dict.keys():
        f.write(f'{sgrna}\t{round(np.mean(dict[sgrna]), 4)}\t{round(np.median(dict[sgrna]), 4)}\t{round(np.std(dict[sgrna]), 4)}\t{",".join(map(lambda x: str(round(x, 2)), dict[sgrna]))}\n')
message(f'Output file {args.O} written successfully')


# if plot-agg flag is provided, plot the aggregate score distribution for sgRNAs
if args.plot_agg:
    message(f'Plotting aggregate score distributions')
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_style("ticks")

    for sgrna in dict.keys():
        # init fig, ax
        fig, ax = plt.subplots()
        sns.histplot(dict[sgrna], bins=100, color='lightcoral', kde=True, stat='density')
        ax.axvline(np.mean(dict[sgrna]), color='black', linestyle='--')
        ax.set_xlabel('Aggregate score', fontsize=12)
        ax.set_ylabel('Density', fontsize=12)
        ax.tick_params(labelsize=12)
        ax.set_title(f'sgRNA: {sgrna}\nMean aggregate score: {round(np.mean(dict[sgrna]), 4)}', fontsize=12)
        sns.despine()
        plt.tight_layout()
        plt.savefig(f'{sgrna}.png', dpi=500)
        plt.close()

message(f'crispAI finished successfully')