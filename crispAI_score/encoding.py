import numpy as np
import pandas as pd

def one_hot_encode(seq, indel = False):
    if indel:
        mapping = dict(zip("ATGC-", range(5)))
        mapping = dict(zip("ATGC-", range(5)))    
        seq_ = [mapping[i] for i in seq]
        return np.eye(5, dtype=int)[seq_]

    else:
        mapping = dict(zip("ATGC", range(4)))    
        seq_ = [mapping[i] for i in seq]
        return np.eye(4, dtype=int)[seq_]


def lin_or_encoding(grna, target):

    def direction_of_bp_loci(grna, target):
        mapping = dict(zip("ATGC", range(4)))
        grna_ = [mapping[i] for i in grna]
        target_ = [mapping[i] for i in target]
        # compare grna_ and target_
        # if the number is same direction is [0,0], if number inreases direction is [1,0] and if number decreases direction is [0,1]
        direction = np.zeros((len(grna_),2))
        for i in range(len(grna_)):
            if grna_[i] == target_[i]:
                direction[i] = [0,0]
            elif grna_[i] < target_[i]:
                direction[i] = [1,0]
            else:
                direction[i] = [0,1]

        return direction       

    if isinstance(grna, pd.Series):
        grna_ohe = grna.apply(one_hot_encode)
        target_ohe = target.apply(one_hot_encode)
        directions = [direction_of_bp_loci(grna_, target_) for grna_, target_ in zip(grna, target)]


        res = np.asarray(list(np.bitwise_or(grna_ohe.to_numpy(), target_ohe.to_numpy())))
        # concatenate directions and res along axis 2
        res = np.concatenate((res, np.asarray(directions)), axis=2)
        return list(res)
        
    elif isinstance(grna, str):
        grna_ohe = one_hot_encode(grna)
        target_ohe = one_hot_encode(target)
        direction = direction_of_bp_loci(grna, target)
        res = np.bitwise_or(grna_ohe, target_ohe)
        # concatenate direction to res along axis 1
        res = np.concatenate((res, direction), axis=1)
        return list(res)

    else:
        return None