from encoding import lin_or_encoding

import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
import seaborn as sns
import torch
import pdb

def preprocess_features(df,
                    reads:str = None,
                    target:str = None,
                    offtarget_sequence:str = None,
                    distance: str = None,
                    read_cutoff: int = 10,
                    max_reads: int = 1000,
                    nupop_occupancy_col:str = None,
                    nupop_affinity_col:str = None,
                    gc_content_col:str = None,
                    nucleotide_bdm_col:str = None
                    ):


    data_encoded = lin_or_encoding(df[target],df[offtarget_sequence])
    df['interface_encoding'] = data_encoded

    nupop_nan_indexes = df[df[nupop_affinity_col].isna()].index # 7 samples in total
    df = df[~df[nupop_affinity_col].isna()]
    df = df.reset_index(drop=True)
    nupop_occupancy = df[nupop_occupancy_col].apply(lambda x: np.asarray(x.split(',')[73:96], dtype=np.float32))
    gc_flank = df[gc_content_col].apply(lambda x: np.asarray(x.split(',')[73:96], dtype=np.float32))
    nupop_affinity = df[nupop_affinity_col].apply(lambda x: x.split(',')[73:96])
    nucleotide_bdm = df[nucleotide_bdm_col].apply(lambda x: x.split(','))

    # check if any of the nupop_affinity sequences contain 'NA' and replace with 0
    for i in range(len(nupop_affinity)):
        if 'NA' in nupop_affinity[i]:
            nupop_affinity[i] = [0 if x == 'NA' else x for x in nupop_affinity[i]]

    for i in range(len(nucleotide_bdm)):
        if 'NA' in nucleotide_bdm[i]:
            nucleotide_bdm[i] = [0 if x == 'NA' else x for x in nucleotide_bdm[i]]

    nupop_affinity = nupop_affinity.apply(lambda x: np.asarray(x, dtype=np.float32))
    nupop_occupancy = nupop_occupancy.apply(lambda x: np.asarray(x, dtype=np.float32))
    nucleotide_bdm = nucleotide_bdm.apply(lambda x: np.asarray(x, dtype=np.float32))

    df['NuPoP occupancy'] = nupop_occupancy
    df['NuPoP affinity'] = nupop_affinity
    df['GC flank'] = gc_flank
    df['nucleotide BDM'] = nucleotide_bdm

    nupop_lens = df['NuPoP occupancy'].apply(lambda x: len(x))
    nupop_lens = nupop_lens[nupop_lens != 23].index
    df = df[~df['NuPoP occupancy'].apply(lambda x: len(x) != 23)]
    df = df.reset_index(drop=True)

    nupop_affinity = np.stack([x for x in df['NuPoP affinity']], axis=0)  
    nupop_occupancy = np.stack([x for x in df['NuPoP occupancy']], axis=0)
    gc_flank = np.stack([x for x in df['GC flank']], axis=0)
    nucleotide_bdm = np.stack([x for x in df['nucleotide BDM']], axis=0)

    df['NuPoP occupancy'] = (df['NuPoP occupancy'] - nupop_occupancy.min()) / (nupop_occupancy.max() - nupop_occupancy.min())
    df['NuPoP affinity'] = (df['NuPoP affinity'] - nupop_affinity.min()) / (nupop_affinity.max() - nupop_affinity.min())
    df['GC flank'] = (df['GC flank'] - gc_flank.min()) / (gc_flank.max() - gc_flank.min())
    df['nucleotide BDM'] = (df['nucleotide BDM'] - nucleotide_bdm.min()) / (nucleotide_bdm.max() - nucleotide_bdm.min())
    df['physical_features'] = df.apply(lambda row: np.concatenate([row['NuPoP occupancy'], row['NuPoP affinity'], row['GC flank'], row['nucleotide BDM']]).reshape(23, 4), axis=1)

    return df
