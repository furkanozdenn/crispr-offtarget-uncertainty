import genomepy
from pybdm import BDM
import pdb
import os
import numpy as np
import pandas as pd
import warnings

warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")

'''
Annotate off-target dataframe with physical features
    - GC content 
    - Nucleotide BDM
    - NuPoP Affinity and Occupancy scores
- Input: off-target dataframe with columns:
    chrom_col = 'chr'
    start_col = 'start'
    end_col = 'end'
    strand_col = 'strand'
    offtarget_sequence_col = 'target_sequence'
    target_sequence_col = 'sgRNA_sequence'
    mismatch_col = 'mismatch'
'''

def bp_wise_NucleotideBDM(seq, bdm, nuc_dict):
    if 'N' in seq:
        return ','.join(['NA']*23)
    
    if len(seq) != 169:
        return ','.join(['NA']*23)
    
    out = []
    for i in range(23):
        out.append(round(bdm.bdm(np.array([nuc_dict[x] for x in seq[i: i+147]])), 2))
    out = ','.join(map(str, out))
    return out


def gc_content_flank(flank_size, sequence):
    'returns the GC content of the given sequence for each position. Each position is the center of a window of size flank_size bp flank_size on each side' 
    gc_content = []
    sequence = sequence.upper()
    for i in range(flank_size):
        gc_content.append('NA')
    for i in range(flank_size, len(sequence)-flank_size):
        window = sequence[i-flank_size:i+flank_size + 1]
        gc_count = window.count('G') + window.count('C')
        gc_content.append(round(gc_count/(2*flank_size + 1), 3))
    for i in range(flank_size):
        gc_content.append('NA')

    return gc_content

def read_nupop_output(dir, file_id):
    with open(dir + str(file_id) + '.seq_Prediction4.txt') as f:
        df = pd.read_csv(f, sep='\t')
        df = df.to_numpy()
        affinity, occupancy = [], []
        for line in df:
            occup = line[0][19:24]
            aff = line[0][-6:]
            if aff[-2:] == 'NA':
                aff = 'NA'
            elif aff[-6] == ' ':
                aff = aff[-5:]
            affinity.append(aff)
            occupancy.append(occup)
        affinity = ','.join(affinity)
        occupancy = ','.join(occupancy)
        return occupancy, affinity

def get_surrounding_sequence(chrom, start, end, strand, genome, flank=100):
    'if strand is negative, return reverse complement'
    rc = False
    if strand == '-':
        rc = True
    seq = genome.get_seq(name=chrom, start=start-flank+1, end=end+flank, rc=rc)
    return str(seq)

def get_surrounding_sequence_from_row(row, genome, flank, chrom_col, start_col, end_col, strand_col):
    if row[chrom_col].find('_') != -1:
        return 'N/A'    
    else:
        return get_surrounding_sequence(row[chrom_col], row[start_col], row[end_col]+ 1, row[strand_col], genome, flank)
    

def annotation_pipeline(offtarget_data):
    import genomepy # TODO:add other reference versions
    from pybdm import BDM

    g = genomepy.Genome('GRCh38')

    # required columns in offtarget_data
    chrom_col = 'chr'
    start_col = 'start'
    end_col = 'end'
    strand_col = 'strand'
    reads_col = 'reads'
    offtarget_sequence_col = 'target_sequence'
    target_sequence_col = 'sgRNA_sequence'

    offtarget_data['context sequence flank_73'] = offtarget_data.apply(lambda row: get_surrounding_sequence_from_row(row, g, flank=73, chrom_col= chrom_col, start_col=start_col, end_col=end_col, strand_col=strand_col), axis=1)

    print('Flank sequences generated from GRCh38')

    # NuPoP occupancy and affinity annotations (R script runs)
    offtarget_data = offtarget_data[offtarget_data[chrom_col].str.find('_') == -1]
    offtarget_data.reset_index(drop=True, inplace=True)
    context_sequence_flank_73 = offtarget_data['context sequence flank_73']
    chrom = offtarget_data[chrom_col]

    # create directory if it doesn't exist
    if not os.path.exists('./nupop_input_temp/'):
        os.makedirs('./nupop_input_temp/')
    for i in range(len(context_sequence_flank_73)):
        with open('./nupop_input_temp/' + str(i+1) + '.seq', 'w') as f:

            try:
                chrom_ = chrom.iloc[i]
                context_sequence_flank_73_ = context_sequence_flank_73.iloc[i]
                f.write('>' + chrom_ + '\n')
                f.write(context_sequence_flank_73_)

            except:
                pass
    
    nupop_output_dir = './nupop_output_temp/'
    # make /nupop_output_temp/ directory if it doesn't exist
    if not os.path.exists(nupop_output_dir):
        os.makedirs(nupop_output_dir)
    
    # copy nupop.R to nupop_output_temp/
    os.system(f'cp nupop.R {nupop_output_dir}')
    pwd = os.getcwd()
    # cd to nupop_output_temp/ and run nupop.R
    os.chdir(nupop_output_dir)

    # change 2nd line of nupop.R 
    with open('nupop.R', 'r') as f:
        lines = f.readlines()
    lines[1] = "nupop_input_dir <- '../nupop_input_temp/'\n"
    with open('nupop.R', 'w') as f:
        f.writelines(lines)

    # run nupop.R
    os.system('Rscript nupop.R > /dev/null 2>&1')
    os.chdir(pwd)

    # read every file in the directory with '.txt' extension in nupop_output_temp/
    files = [f for f in os.listdir('./nupop_input_temp/') if f.endswith('.seq')]
    len_files = len(files)
    occupancies, affinities = [], []
    for i in range(1, len_files+1):
        try:
            occ, aff = read_nupop_output(nupop_output_dir, i)
            occupancies.append(occ)
            affinities.append(aff)
        except:
            # print(i)
            occupancies.append('NA')
            affinities.append('NA')
    
    offtarget_data['NuPoP occupancy'] = occupancies
    offtarget_data['NuPoP affinity'] = affinities

    os.system('rm -r ./nupop_input_temp/')
    os.system('rm -r ./nupop_output_temp/')
    print('NuPoP occupancy and affinity annotations generated')


    gc_content_flank73 = []
    for i in range(len(context_sequence_flank_73)):
        res_ = gc_content_flank(73, context_sequence_flank_73[i])
        gc_content_flank73.append(','.join(map(str, res_)))
    
    offtarget_data['GC flank73'] = gc_content_flank73

    print('GC content annotation generated')

    bdm = BDM(ndim=1, nsymbols=4)
    nuc_dict = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    offtarget_data['nucleotide BDM'] = context_sequence_flank_73.apply(lambda seq: bp_wise_NucleotideBDM(seq.upper(), bdm, nuc_dict))

    print('Nucleotide BDM annotation generated')

    offtarget_data['target_N'] = offtarget_data[target_sequence_col]
    unique_targets = offtarget_data[target_sequence_col].unique()

    for target_ in unique_targets:
        # if -3rd base is not N then 'Target_Sequence_0mm' is the same as 'target_sequence'
        if target_[-3] != 'N':
            offtarget_data.loc[offtarget_data[target_sequence_col] == target_, 'Target_Sequence_0mm'] = target_
        else:
            try:
                offtarget_data.loc[offtarget_data[target_sequence_col] == target_, 'Target_Sequence_0mm'] = offtarget_data.loc[(offtarget_data[target_sequence_col] == target_) & (offtarget_data[mismatch_col] == 0), offtarget_sequence_col].values[0]
            except:
                # except
                offtarget_data.loc[offtarget_data[target_sequence_col] == target_, 'Target_Sequence_0mm'] = target_.replace('N', 'A')

    offtarget_data[target_sequence_col] = offtarget_data['Target_Sequence_0mm']

    return offtarget_data