import requests 
import subprocess 
import os
from datetime import datetime
import logging
from alignments.paths import FR3D_PATH


def store_file(data, suffix):
    from random import choice
    now = str(datetime.now().timestamp())
    
    digit = chr(choice(range(ord('a'), ord('z')))) # to make sure the dir is unique
    file_path = os.path.join("/tmp/", f"{digit}_{now}", f'cust.{suffix}')
    # file_path = f'/tmp/cust2_{digit}.{suffix}'
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    # raise Exception(f"File stored: {filepath}")
    with open(file_path, 'wb') as f:
            f.write(data)
            f.close()
            
    return file_path

def download_struct_cif(struct_id, entity_id):
    # url = f"https://rna.bgsu.edu/rna3dhub/rest/getChainSequenceBasePairs?pdb_id={struct_id}&chain={chain}&only_nested={only_nested}"
    # cif_url = f"https://models.rcsb.org/v1/{struct_id}/atoms?label_entity_id={entity_id}&encoding=cif"
    cif_url = f"http://files.rcsb.org/download/{struct_id}.cif"
    
    cif_file_data = requests.get(cif_url).content
    cif_file_path = store_file(cif_file_data, "cif")
    return cif_file_path

def run_fred(filepath, base_name, output, chain_id):
    logger = logging.getLogger("ribovision3-logger")
    
    commands = [
        "/usr/bin/python3", f"{FR3D_PATH}/NA_pairwise_interactions.py", 
        "--input", filepath, base_name, "-o", f"{output}/results/json",
        "-f", "ebi_json", "--chain", chain_id
    ]
    logger.debug(" ".join(commands))
    subprocess.run(commands)
    
def replace_r2dt_base_pairs_with_fred(output, base_name, chain_id, ext):
    rna2d_path = os.path.join(output, "results/json", base_name.replace(ext, f"_{chain_id}_basepair.json"))
    rna2d_bp_path = os.path.join(output, "results/json", "BP_json.json")
    os.rename(rna2d_path, rna2d_bp_path)
    
def get_fred_base_pairs(struct_id, entity_id, chain_id, output, file_path=""):
    logger = logging.getLogger("ribovision3-logger")
    
    if not file_path:
        logger.warning("Structure not provided, Downloading structure.")
        file_path = download_struct_cif(struct_id, entity_id)

    ext = file_path.split(".")[-1]

    dir_path, base_name = os.path.split(file_path)
    logger.debug(f"Exracting basepair from {struct_id} chain {chain_id} using {FR3D_PATH}")
    run_fred(dir_path, base_name, output, chain_id)
    
    replace_r2dt_base_pairs_with_fred(output, base_name, chain_id, f".{ext}")
    logger.info("Base pairs extracted and replaced with R2DT")
    rna2d_path = os.path.join(output, "results/json", base_name.replace(f".{ext}", f"_{chain_id}_basepair.json"))
    return rna2d_path
    
    
if __name__ == "__main__":
    # for 28S
    get_fred_base_pairs("6QZP", "1", "L5", "/home/RiboVision3/R2DT1_4/R2DT1_4/R2DT/R2DT-test20_2024_9_5_17_38_58_540143")

    
    
    
    
    