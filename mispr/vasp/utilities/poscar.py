from __future__ import annotations

import os
import logging
from xxlimited import Str

from pymatgen.core.surface import SlabGenerator
#from pymatgen.transformations.standard_transformations import SupercellTransformation
from pymatgen.transformations.advanced_transformations import AddAdsorbateTransformation
#from pymatgen.core.adsorption import AdsorbateSiteFinder
#from pymatgen.core.sites import PeriodicSite,Site
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.analysis.interfaces.coherent_interfaces import CoherentInterfaceBuilder

from mispr.gaussian.utilities.mol import process_mol

logger = logging.getLogger(__name__)

def write_poscar(structure, operation_type, **kwargs):
    """
    Write a POSCAR file for a structure.

    Args:
        structure (Structure): The structure to write.
        filename (str, optional): The name of the POSCAR file to write. Default "POSCAR".
        **kwargs: Additional keyword arguments to pass to the Poscar constructor.

    Returns:
        str: The path to the written POSCAR file.
    """
    filename ="POSCAR"
    
    if operation_type == "unitcell":
        poscar = Poscar(structure, **kwargs)
        poscar.write_file(filename)
        return filename
    
    elif operation_type == "supercell":
        multiplicity = kwargs.get("multiplicity") or (2, 2, 1)
        supercell_structure = structure.make_supercell(multiplicity)
        poscar = Poscar(supercell_structure, **kwargs)
        poscar.write_file(filename)
        return filename
    
    elif operation_type == "surface":
        miller_index = kwargs.get("miller_index") or (1, 1, 1)
        min_slab_size = kwargs.get("min_slab_size") or 10
        min_vacuum_size = kwargs.get("min_vacuum_size") or 10
        center_slab = kwargs.get("center_slab") or False
        primitive = kwargs.get("primitive") or False
        max_normal_search = kwargs.get("max_normal_search") or 10
        multiplicity = kwargs.get("multiplicity") or (1, 1, 1)
        slabgen = SlabGenerator(structure, miller_index, min_slab_size=min_slab_size, min_vacuum_size=min_vacuum_size, center_slab=center_slab, primitive=primitive, max_normal_search=max_normal_search,)
        slabs = slabgen.get_slabs()
        struc = slabs[0]
        struc_supercell = struc.make_supercell(multiplicity)
        poscar = Poscar(struc_supercell, comment = f"structure {miller_index[0]} {miller_index[1]} {miller_index[2]} surface")
        poscar.write_file(filename)
        return filename
    
    #TODO: Add Interface strucutre generation and POSCAR writing for a base and a thin slab structure. Add Adsorption structure generation and POSCAR writing for a base and an adsorbate molecule.

    elif operation_type == "interface":
        structure1 = structure
        structure2 = kwargs.get("structure2")
        miller_index_1 = kwargs.get("miller_index_1") or (1, 1, 1)
        miller_index_2 = kwargs.get("miller_index_2") or (1, 1, 1)
        max_lattice_match_length = kwargs.get("max_lattice_match_length") or 20
        max_area_ratio = kwargs.get("max_area_ratio") or 0.5
        max_mismatch = kwargs.get("max_mismatch") or 0.05
        max_angle_diff = kwargs.get("max_angle_diff") or 5
        interface_builder = CoherentInterfaceBuilder(structure1, structure2, miller_index_1, miller_index_2, max_lattice_match_length=max_lattice_match_length, max_area_ratio=max_area_ratio, max_mismatch=max_mismatch, max_angle_diff=max_angle_diff,)
        interfaces = interface_builder.get_interfaces()
        struc = interfaces[0]
        poscar = Poscar(struc, comment=f"interface between {miller_index_1} and {miller_index_2}")
        poscar.write_file(filename)
        return filename
    
    elif operation_type == "adsorption":
        adsorbate_molecule = kwargs.get("adsorbate_mol")
        height = kwargs.get("height") or 2.0
        miller_index = kwargs.get("miller_index") or (1, 1, 1)
        min_slab_size = kwargs.get("min_slab_size") or 10
        min_vacuum_size = kwargs.get("min_vacuum_size") or 10
        center_slab = kwargs.get("center_slab") or False
        primitive = kwargs.get("primitive") or False
        max_normal_search = kwargs.get("max_normal_search") or 10
        multiplicity = kwargs.get("multiplicity") or (1, 1, 1)
        slabgen = SlabGenerator(structure, miller_index, min_slab_size=min_slab_size, min_vacuum_size=min_vacuum_size, center_slab=center_slab, primitive=primitive, max_normal_search=max_normal_search,)
        slabs = slabgen.get_slabs()
        struc = slabs[0]
        struc_supercell = struc.make_supercell(multiplicity)
        ads_site_finder = AdsorbateSiteFinder(struc_supercell)
        sites = ads_site_finder.find_adsorption_sites()
        adsorption_sites = ads_site_finder.find_adsorption_sites()
        if not adsorption_sites:
            raise ValueError(f"{adsorption_sites} site not available. Available: {list(sites.keys())}")
        adsorption_site = adsorption_sites[0]
        adsorbate_mol = process_mol(adsorbate_molecule)
        transformation = AddAdsorbateTransformation(adsorbate_mol, adsorption_site, height=height)
        adsorbed_structure = transformation.apply_transformation(struc_supercell)
        poscar = Poscar(adsorbed_structure, comment=f"adsorption of {adsorbate_molecule} on {miller_index} surface")
        poscar.write_file(filename)
        return filename