from polyscript.brics.BRICS_utils import read_df, apply_BRICS

if __name__ == "__main__":
    smiles = read_df("./source.csv", "SMILES") # replace the csv files path with real content and the column name
    apply_BRICS(smiles)
