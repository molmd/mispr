import os
from pymatgen.io.vasp.inputs import Poscar, Potcar

def generate_potcar(poscar_path, potcar_dir, functional):
    """
    Generate POTCAR file from POSCAR and POTCAR directory.

    Parameters
    ----------
    poscar_path : str
        Path to POSCAR file
    potcar_dir : str
        Path to VASP POTCAR directory
    functional : str
        POTCAR functional (e.g., PBE, PBE_54, LDA)
    """
    potcar_paths = []
    poscar = Poscar.from_file(poscar_path)
    elements = [site.specie.symbol for site in poscar.structure]
    elements = sorted(list(set(elements)), key=elements.index)
    #os.environ["PMG_VASP_PSP_DIR"] = potcar_dir
    #potcar = Potcar(symbols=elements, functional=functional)
    #potcar.write_file("POTCAR")
    #TODO: add ability to specify custom POTCAR paths for each element and functional. Add error handling for missing POTCAR files and invalid functionals.
    for el in elements:
        potcar_path = os.path.join(potcar_dir, functional, el, "POTCAR")
        if not os.path.exists(potcar_path):
            raise FileNotFoundError(f"POTCAR not found for {el}: {potcar_path}")
        potcar_paths.append(potcar_path)
        
    with open("POTCAR", "wb") as outfile:
        for path in potcar_paths:
            with open(path, "rb") as infile:
                outfile.write(infile.read())
