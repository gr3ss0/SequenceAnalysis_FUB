import pandas as pd
import os

# list of input files from Snakemake
files = snakemake.input
sample0=os.path.basename(files[0]).replace(".stats", "")
df0=pd.read_csv(files[0], sep="\t", header=None, names=["ref", "length", sample0, "unmapped"])
df0 = df0.set_index(["ref", "length"])#to make the first two columns as index, so that we can merge the dataframes later
df0 = df0[[sample0]]#keeping only 3rd column



for f in files[1:]:
    sample = os.path.basename(f).replace(".stats", "")#get rid of full path and extension etc

    df = pd.read_csv(
        f,#full path
        sep="\t",
        header=None,
        names=["ref", "length", sample, "unmapped"]#need to add sample name to the already existing column names
    )

    df=df.set_index(["ref", "length"])
    df = df[[sample]] #keeping only 3rd column
    
    # merge with the first dataframe
    df0 = pd.concat([df0, df], axis=1)


df0 = df0.reset_index()

# write output
df0.to_csv(snakemake.output[0], sep="\t", index=False)