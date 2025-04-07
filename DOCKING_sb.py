#!/usr/bin/env python3
#docking_sb
import requests
import openbabel as ob
import os
import sys
from time import sleep
from tempfile import NamedTemporaryFile

"""DockingSB: Batch docking script for submitting ligands to docking services."""

DOCKING_SERVICES = [
    {"url": "https://dockthor.lncc.br", "compounds": "/api/compounds", "dock": "/api/dock_batch", "timeout": 60},
    {"url": "https://www.dockingserver.com", "compounds": "/api/ligands", "dock": "/api/dock", "timeout": 60},
]

def fetch_ligands_from_website(service, limit=50):
    """Fetch ligands from a docking service API."""
    try:
        response = requests.get(f"{service['url']}{service['compounds']}?limit={limit}", timeout=10)
        response.raise_for_status()
        data = response.json()
        return data.get('compounds', [])
    except requests.RequestException as e:
        print(f"Failed to fetch ligands: {e}")
        return []

def smiles_to_pdbqt(smiles, output_file):
    """Convert SMILES string to PDBQT format."""
    try:
        conv = ob.OBConversion()
        conv.SetInAndOutFormats("smi", "pdbqt")
        mol = ob.OBMol()
        conv.ReadString(mol, smiles)
        mol.AddHydrogens()
        conv.WriteFile(mol, output_file)
        return True
    except Exception as e:
        print(f"SMILES to PDBQT failed: {e}")
        return False

def fetch_pdb_receptor(pdb_id, output_file="receptor.pdbqt"):
    """Download and convert PDB receptor to PDBQT."""
    try:
        response = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=10)
        response.raise_for_status()
        with NamedTemporaryFile(suffix=".pdb", delete=False) as temp:
            temp.write(response.text.encode())
            temp_path = temp.name
        conv = ob.OBConversion()
        conv.SetInAndOutFormats("pdb", "pdbqt")
        mol = ob.OBMol()
        conv.ReadFile(mol, temp_path)
        mol.AddHydrogens()
        conv.WriteFile(mol, output_file)
        os.remove(temp_path)
        return output_file
    except Exception as e:
        print(f"Failed to fetch receptor: {e}")
        return None

def dock_ligand_batch(service, ligands, receptor_file, use_smiles=True):
    """Submit ligands for batch docking."""
    try:
        if use_smiles:
            data = {"ligands": [ligand.get('smiles') for ligand in ligands if ligand.get('smiles')]}
            with open(receptor_file, "rb") as receptor:
                response = requests.post(f"{service['url']}{service['dock']}", data=data, files={"receptor": receptor}, timeout=service['timeout'])
        else:
            files = []
            with open(receptor_file, "rb") as receptor:
                files.append(("receptor", receptor))
                for ligand in ligands:
                    if 'file' in ligand:
                        files.append(("ligands", open(ligand['file'], "rb")))
                response = requests.post(f"{service['url']}{service['dock']}", files=files, timeout=service['timeout'])
        response.raise_for_status()
        result = response.json()
        if not use_smiles:
            for f in files[1:]:
                f[1].close()
        return result
    except requests.RequestException as e:
        print(f"Docking failed: {e}")
        return {"error": True}

def save_result(mol_name, result, output_file="docking_results.txt"):
    """Save docking results to a file."""
    try:
        with open(output_file, "a") as f:
            f.write(f"{mol_name}: {result}\n")
    except IOError as e:
        print(f"Failed to save result: {e}")

def main(batch_size=50):
    """Main function to run batch docking."""
    receptor_input = input("Enter receptor file path (e.g., receptor.pdbqt) or PDB ID (e.g., 1ABC): ").strip()
    
    if receptor_input.lower().endswith('.pdbqt'):
        receptor_file = receptor_input
        if not os.path.exists(receptor_file):
            print("Receptor file not found!")
            sys.exit(1)
    else:
        receptor_file = fetch_pdb_receptor(receptor_input)
        if not receptor_file:
            print("Failed to fetch receptor from PDB!")
            sys.exit(1)
    
    print("Trying all available batch docking services...")
    
    for service in DOCKING_SERVICES:
        print(f"Using service: {service['url']}")
        start = 0
        use_smiles = True
        
        while True:
            print(f"Fetching ligands {start} to {start + batch_size}")
            ligands = fetch_ligands_from_website(service, batch_size)
            if not ligands:
                print("No ligands or fetch failed. Moving to next service.")
                break
            
            processed_ligands = []
            mol_names = []
            for ligand in ligands:
                smiles = ligand.get('smiles')
                if not smiles:
                    continue
                mol_name = ligand.get('name', ligand.get('id', f"ligand_{start}"))
                if not use_smiles:
                    ligand_file = f"{mol_name}.pdbqt"
                    if smiles_to_pdbqt(smiles, ligand_file):
                        ligand['file'] = ligand_file
                    else:
                        continue
                processed_ligands.append(ligand)
                mol_names.append(mol_name)
            
            if not processed_ligands:
                start += batch_size
                continue
            
            results = dock_ligand_batch(service, processed_ligands, receptor_file, use_smiles)
            if "error" in results and use_smiles:
                use_smiles = False
                for ligand, mol_name in zip(processed_ligands, mol_names):
                    ligand_file = f"{mol_name}.pdbqt"
                    if smiles_to_pdbqt(ligand.get('smiles'), ligand_file):
                        ligand['file'] = ligand_file
                results = dock_ligand_batch(service, processed_ligands, receptor_file, use_smiles)
            
            if "error" not in results:
                scores = results.get("scores", [])
                for i, mol_name in enumerate(mol_names):
                    result = scores[i] if i < len(scores) else "No score"
                    save_result(mol_name, result)
                    print(f"Docked {mol_name}: {result}")
            else:
                print("Docking failed. Trying next service.")
                break
            
            if not use_smiles:
                for ligand in processed_ligands:
                    if 'file' in ligand:
                        try:
                            os.remove(ligand['file'])
                        except OSError as e:
                            print(f"Failed to remove {ligand['file']}: {e}")
            
            sleep(5)
            start += batch_size

if __name__ == "__main__":
    try:
        batch_size = int(input("Enter batch size (default 50): ") or 50)
        if batch_size <= 0:
            raise ValueError
    except ValueError:
        print("Invalid batch size! Using default (50).")
        batch_size = 50
    main(batch_size)


#end