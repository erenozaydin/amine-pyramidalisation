import numpy as np
from cov_radii import cov_dictionary

def parse_line(line):
    items = line.split()
    if len(items) >= 4:
        chemicalsymbol = items[0]
        atomcoordinates = np.array(items[1:], dtype=float)
        return [chemicalsymbol, atomcoordinates]
    else:
        return None

def parse_xyz_file(xyzfile):
    """
    This reads and validates .xyz files, extracting atomic information.

    Intake:
        xyzfile (str): path to the .xyz file.

    Output:
        list: parsed list of [chemical symbol, coordinates] for each atom.

    Throws:
        ValueError: If the file is improperly formatted.
        FileNotFoundError: If the file does not exist.
    """
    try:
        with open(xyzfile, 'r') as file:
            lines = file.readlines()
        if len(lines) < 3:
            raise ValueError("This .xyz file is improperly formatted - there are not enough atoms :(")
        try:
            atom_count = int(lines[0].strip())
        except ValueError:
            raise ValueError("This .xyz file is improperly formatted - the number of atoms does not match the number in the first line :(")
        if len(lines) != atom_count + 2:
            raise ValueError("This .xyz file is improperly formatted - the number of lines is wrong :(")
        parsed_lines = [parse_line(line) for line in lines[2:] if parse_line(line) is not None]
        if len(parsed_lines) != atom_count:
            raise ValueError("This .xyz file is improperly formatted - the number of lines is wrong :(")
        return parsed_lines
    except FileNotFoundError:
        raise FileNotFoundError(f"The file {xyzfile} does not exist.")
    except Exception as e:
        raise ValueError("This .xyz file is improperly formatted :(") from e

def find_N(parsexyzfileoutput):
    """
    This searches for the nitrogen atom in the parsed file.

    Intake:
        parsexyzfileoutput (list): list of [symbol, coordinates] from .xyz file.

    Output:
        list: the nitrogen atom [symbol, coordinates], or None if not found.
    """
    for item in parsexyzfileoutput:
        if item[0] == 'N':
            return item
    return None
    
def distance(coord1, coord2):
    return np.linalg.norm(coord1 - coord2)

def bonded(list1, list2):
    """
    Utilises covalent radii from cov_radii to determine which atoms are bonded to a central atom.

    Intake:
        list1 (list): [element, coordinates] of the central atom.
        list2 (list): List of all [element, coordinates] atoms.

    Output:
        list: atoms bonded to list1.
    """
    bondedatoms = []
    firstelementsymbol, firstatomcoordinates = list1
    radius1 = cov_dictionary.get(firstelementsymbol, 0)
    
    for atom in list2:
        if atom == list1:
            continue
        nextelementsymbol, nextatomcoordinates = atom
        radius2 = cov_dictionary.get(nextelementsymbol, 0)
        if distance(firstatomcoordinates, nextatomcoordinates) <= (radius1 + radius2 + 0.4) and distance(firstatomcoordinates, nextatomcoordinates) != 0:
            bondedatoms.append(atom)
    return bondedatoms
    
def angle_between(v1, v2):
    """
    Pretty self explanatory: find the angle between 2 vectors in radians.

    Intake:
        v1 (np.array): first vector.
        v2 (np.array): second vector.

    Output:
        float: angle in radians.
    """
    cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    return np.arccos(np.clip(cos_theta, -1.0, 1.0))

def bond_angles3(list1, list2):
    """
    This computes all unique bond angles between atoms bonded to a central atom.

    Intake:
        list1 (list): central atom [symbol, coordinates].
        list2 (list): list of bonded atoms.

    Output:
        list: list of angles (radians) between each pair of bonded atoms.
    """
    firstchemicalsymbolforbondangles, firstatomcoordinatesforbondangles = list1
    listofvectors = [nextcoordinatesforbondangles - firstatomcoordinatesforbondangles for thenameofthisvariabledoesnotseemtomatter, nextcoordinatesforbondangles in list2]
    
    angles = []
    for i in range(len(listofvectors)):
        for j in range(i + 1, len(listofvectors)):
            angles.append(angle_between(listofvectors[i], listofvectors[j]))
    return angles

def calculatePsigma(angles):
    """
    This calculates P_sigma (a measure of non-planarity in trivalent nitrogen centres).

    Intake:
        angles (list): list of angles in radians.

    Output:
        float: P_sigma value.
    """
    anglesum = sum(angles)
    P_sigma = np.sqrt(2 * np.pi - anglesum)
    return P_sigma
    
#HERE THE MAIN PROGRAM BEGINS
if __name__ == "__main__":

    while True:
        xyzfile = input("Please enter the full path for the .xyz file, including the .xyz suffix: ").strip()
        if not xyzfile.lower().endswith('.xyz'):
            print("The file must have '.xyz' at the end. Ensure all of the letters in the suffix are lowercase and that you have not added any spaces after the suffix, then try again.")
            continue

        try:
            result = parse_xyz_file(xyzfile)
            nitrogen_atom = find_N(result)
            if nitrogen_atom is None:
                print("No nitrogen atom found in the .xyz file.")
                continue

#finding which atoms are bonded to nitrogen
            bonded_atoms = bonded(nitrogen_atom, result)
            if len(bonded_atoms) < 3:
                print("Not enough atoms bonded to nitrogen to calculate three bond angles.")
                continue
            angles = bond_angles3(nitrogen_atom, bonded_atoms)
            if len(angles) < 3:
                print("Could not calculate three bond angles.")
                continue

            P_sigma = calculatePsigma(angles)
            print(f"\nP_Sigma = {P_sigma:.10g}")

        except FileNotFoundError:
            print(f"The file {xyzfile} cannot be located or perhaps does not exist. Make sure you have provided the full path, in a format like this: C:\\Users\\imperialstudentusername\\Downloads\\ethylamine.xyz")
        except Exception as e:
            print(f"An error occurred: {e}.")

#do you want to go again with another .xyz file?
        while True:
            another = input("Would you like to go again with another .xyz file? (yes/no): ").strip().lower()
            if another in ['yes', 'no']:
                break
            print("Please answer with 'yes' or 'no'.")
        
        if another == 'no':
            print("Thank you, I hope you enjoyed this work - Eren Ozaydin :)")
            break
#preventing unsightly abrupt closure of the program
    input("\nPress enter now to exit the program.")

