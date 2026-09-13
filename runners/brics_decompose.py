from polyscript.brics import BRICSDecompose
import pandas as pd
from rdkit import Chem
import tqdm


if __name__ == "__main__":
    input_path = "/home/kalki/src/polyscript/benchmarking/generator/datasets/merged_smiles_set.csv"
    column_name = "SMILES"
    output_path = "/home/kalki/src/polyscript/benchmarking/generator/datasets/decomp.csv"

    print("[+] Reading input from", input_path)
    input_list = pd.read_csv(input_path)[column_name].tolist()
    print("[+] Starting Decomposition")
    output_list = []
    for smiles in tqdm.tqdm(input_list):
        # print(list(BRICSDecompose(Chem.MolFromSmiles(smiles))))
        output_list.extend(list(BRICSDecompose(Chem.MolFromSmiles(smiles))))

    print("[+] Decomposition complete")
    print("[+] Post processing output ...")

    output_df = pd.DataFrame({"decomp_SMILES": output_list})

    # drop duplicates
    output_df = output_df.drop_duplicates()
    output_df.to_csv(output_path, index=False)
    print("[+] Output saved to", output_path)
