import os.path
import pandas as pd

sample_file = "data/test/reference.tsv"
data_dir = "data/test"

sample_dat = pd.read_csv(sample_file, sep = "\t")
# remove empty rows
sample_dat.replace("", float("NaN"), inplace=True)
sample_dat.dropna(subset = ["sample"], inplace=True)
# remove empty columns 
# sample_dat.dropna(axis=1, how='all', inplace=True)
# ensure that paths to files are correct
f1_bnames = [os.path.basename(f) for f in sample_dat['file1'].to_list()]
f1_full_paths = [os.path.join(data_dir, b) for b in f1_bnames]
sample_dat['file1'] = f1_full_paths
if(len(sample_dat.file2.value_counts()) > 0):
    f2_bnames = [os.path.basename(f) for f in sample_dat['file2'].to_list()]
    f2_full_paths = [os.path.join(data_dir, b) for b in f2_bnames]
    sample_dat['file2'] = f2_full_paths
    sample_dat['paired'] = "1"
else:
    sample_dat['paired'] = "0"
sample_dat.to_csv(outf, index=False, sep='\t')