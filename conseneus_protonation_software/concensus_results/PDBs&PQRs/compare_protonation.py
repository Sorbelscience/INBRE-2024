from Bio import PDB

###### Extract protonation states from PDB file ######
def extract_protonation_states_pdb(file_path):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', file_path)
    protonation_states = set()
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id()[0] == ' ':  # Exclude heteroatoms and water
                    resname = residue.get_resname()
                    chain_id = chain.get_id()
                    resid = residue.get_id()[1]
                    protonation_states.add((resname, chain_id, resid))
    return protonation_states

###### Extract protonation states from PQR file ######
def extract_protonation_states_pqr(file_path):
    protonation_states = set()
    with open(file_path) as file:
        for line in file:
            if line.startswith('ATOM'):
                resname = line[17:20].strip()
                chain_id = line[21].strip()
                resid = int(line[22:26].strip())
                protonation_states.add((resname, chain_id, resid))
    return protonation_states

###### Compare two sets of protonation states ######
def compare_protonation_states(set1, set2):
    differences = set1.intersection(set2)  # Common elements
    unique_to_set1 = set1 - set2  # Unique to set1
    unique_to_set2 = set2 - set1  # Unique to set2
    return differences, unique_to_set1, unique_to_set2

###### Export comparison results to a file ######
def write_comparison_to_file(file_path, hplusplus_states, propka_states, pypka_states):
    with open(file_path, 'w') as file:
        file.write("Comparison of protonation states between H++, propka, and pyPKa\n")
        file.write("="*60 + "\n\n")

        #### Compare H++ and propka
        file.write("Comparing H++ and propka:\n")
        differences_hp, unique_to_hplusplus, unique_to_propka = compare_protonation_states(hplusplus_states, propka_states)
        file.write(f"Differences: {differences_hp}\n")
        file.write(f"Unique to H++: {unique_to_hplusplus}\n")
        file.write(f"Unique to propka: {unique_to_propka}\n\n")

        #### Compare H++ and pyPKa
        file.write("Comparing H++ and pyPKa:\n")
        differences_hp_py, unique_to_hplusplus_py, unique_to_pypka = compare_protonation_states(hplusplus_states, pypka_states)
        file.write(f"Differences: {differences_hp_py}\n")
        file.write(f"Unique to H++: {unique_to_hplusplus_py}\n")
        file.write(f"Unique to pyPKa: {unique_to_pypka}\n\n")

        #### Compare propka and pyPKa
        file.write("Comparing propka and pyPKa:\n")
        differences_pp, unique_to_propka_py, unique_to_pypka_pp = compare_protonation_states(propka_states, pypka_states)
        file.write(f"Differences: {differences_pp}\n")
        file.write(f"Unique to propka: {unique_to_propka_py}\n")
        file.write(f"Unique to pyPKa: {unique_to_pypka_pp}\n\n")

        #### Detailed interpretation
        file.write("Detailed Interpretation:\n")
        file.write("="*60 + "\n")

        file.write("\nComparing H++ and propka:\n")
        if differences_hp:
            file.write(f"Differences: {differences_hp}\n")
        else:
            file.write("Differences: There are no differences in the protonation states of residues that are present in both H++ and propka.\n")
        file.write(f"Unique to H++: {unique_to_hplusplus}\n")
        file.write(f"Unique to propka: {unique_to_propka}\n\n")

        file.write("Comparing H++ and pyPKa:\n")
        if differences_hp_py:
            file.write(f"Differences: {differences_hp_py}\n")
        else:
            file.write("Differences: There are no differences in the protonation states of residues that are present in both H++ and pyPKa.\n")
        file.write(f"Unique to H++: {unique_to_hplusplus_py}\n")
        file.write(f"Unique to pyPKa: {unique_to_pypka}\n\n")

        file.write("Comparing propka and pyPKa:\n")
        if differences_pp:
            file.write(f"Differences: {differences_pp}\n")
        else:
            file.write("Differences: There are no differences in the protonation states of residues that are present in both propka and pyPKa.\n")
        file.write(f"Unique to propka: {unique_to_propka_py}\n")
        file.write(f"Unique to pyPKa: {unique_to_pypka_pp}\n")

###### Main function to run the comparisons and write the results ######
def main():
    hplusplus_file = 'k1_H++.pdb'
    propka_file = 'k1_propka.pqr'
    pypka_file = 'k1_pypka.pdb'

    #### Extract protonation states from the files
    hplusplus_states = extract_protonation_states_pdb(hplusplus_file)
    propka_states = extract_protonation_states_pqr(propka_file)
    pypka_states = extract_protonation_states_pdb(pypka_file)

    #### Export comparison results to a file
    output_file = 'difference_in_protonation.txt'
    write_comparison_to_file(output_file, hplusplus_states, propka_states, pypka_states)

###### Execute the main function ######
if __name__ == "__main__":
    main()
