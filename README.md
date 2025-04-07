# DockingSB

DockingSB: A Python script for batch docking ligands against receptors using external services.

## Features
- Fetches ligands from docking service APIs.
- Converts SMILES to PDBQT or downloads receptors from PDB.
- Submits docking jobs in batches and saves results.

## Requirements
- Python 3.6+
- Dependencies: `requests`, `openbabel` (Pybel)

## Installation

1. Clone the repository:

git clone https://github.com/Sangeet01/docking_sb.git
cd docking_sb


2. Install dependencies:
pip install -r requirements.txt



3. Ensure Open Babel is installed (e.g., via `conda install -c conda-forge openbabel`).

## Usage

Run the script and follow prompts:

python dockingsb.py



- Enter a receptor file path (`.pdbqt`) or PDB ID (e.g., `1ABC`).
- Specify batch size (default: 50).

## Notes
- Replace placeholder URLs in `DOCKING_SERVICES` with real docking service endpoints.
- Results are saved to `docking_results.txt`.

## License

MIT License (see `LICENSE`).

## Contributing

Fork, edit, and submit a pull request—keep it simple and functional.

##

PPS: Sangeet’s the name, a daft undergrad splashing through chemistry and code like a toddler—my titrations are a mess, and I’ve used my mouth to pipette.

