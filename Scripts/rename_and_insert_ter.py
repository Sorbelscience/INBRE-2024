import os

def rename_atoms_and_insert_ter(input_file, output_file):
    atom_renaming = {
        'HB1': 'HB2', 'HG1': 'HG2', 'HD1': 'HD2', 'HE1': 'HE2', 'HA1': 'HA2', 'HG11': 'HG13',
        'HD2': 'HD3', 'HD3': 'HD1', 'HG2': 'HG3', 'HG3': 'HG1', 'HA2': 'HA3', 'HB2': 'HB3',
        'HE2': 'HE3', 'HE3': 'HE1', 'HB3': 'HB1', 'HG12': 'HG13', 'HD12': 'HD13', 'HD13': 'HD11'
    }
    
    current_chain = None
    with open(input_file, 'r') as file:
        lines = file.readlines()

    with open(output_file, 'w') as file:
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                parts = line.split()
                atom_name = parts[2]
                chain_id = line[21]

                if current_chain is None:
                    current_chain = chain_id
                elif chain_id != current_chain:
                    file.write("TER\n")
                    current_chain = chain_id

                if atom_name in atom_renaming:
                    line = line.replace(atom_name, atom_renaming[atom_name])

            file.write(line)

        file.write("TER\n")  # Ensure the file ends with a TER

# Path to the input and output files
input_file_path = '/home/storm/INBRE-2024/k1_protonation_states/k1.pdb'
output_file_path = '/home/storm/INBRE-2024/k1_protonation_states/renamed_k1.pdb'

rename_atoms_and_insert_ter(input_file_path, output_file_path)
