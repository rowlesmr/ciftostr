# -*- coding: utf-8 -*-
"""
Copyright © 2021 Matthew Rowles

Created on Mon Apr 19 17:10:25 2021

@author: Matthew Rowles
"""


# pylint: disable=line-too-long, invalid-name.

import math
import re
import os
import copy
import CifFile
from typing import List, Tuple, Union, Dict
import contextlib


class Vector3:
    def __init__(self, u: float, v: float, w: float):
        self.u = u
        self.v = v
        self.w = w

    def cross_product(self, other: 'Vector3') -> 'Vector3':
        return Vector3(self.v * other.w - self.w * other.v,
                       self.w * other.u - self.u * other.w,
                       self.u * other.v - self.v * other.u)

    def cross(self, other: 'Vector3') -> 'Vector3':
        return self.cross_product(other)

    def dot_product(self, other: 'Vector3') -> float:
        return self.u * other.u + self.v * other.v + self.w * other.w

    def dot(self, other: 'Vector3') -> float:
        return self.dot_product(other)

    def square_magnitude(self) -> float:
        return self.u ** 2 + self.v ** 2 + self.w ** 2

    def magnitude(self) -> float:
        return math.sqrt(self.square_magnitude())

    def __abs__(self) -> float:
        return self.magnitude()

    def __add__(self, other: 'Vector3') -> 'Vector3':
        if isinstance(other, Vector3):
            return Vector3(self.u + other.u, self.v + other.v, self.w + other.w)
        else:
            return NotImplemented

    def __sub__(self, other: 'Vector3') -> 'Vector3':
        if isinstance(other, Vector3):
            return Vector3(self.u - other.u, self.v - other.v, self.w - other.w)
        else:
            return NotImplemented

    def __mul__(self, other: (int, float)) -> 'Vector3':
        if isinstance(other, (int, float)):
            return Vector3(self.u * other, self.v * other, self.w * other)
        else:
            return NotImplemented

    __rmul__ = __mul__

    def __truediv__(self, other: (int, float)) -> 'Vector3':
        if isinstance(other, (int, float)):
            return Vector3(self.u / other, self.v / other, self.w / other)
        else:
            return NotImplemented

    def __floordiv__(self, other: (int, float)) -> 'Vector3':
        if isinstance(other, (int, float)):
            return Vector3(self.u // other, self.v // other, self.w // other)
        else:
            return NotImplemented

    def __repr__(self):
        return f"Vector3({self.u}, {self.v}, {self.w})"


class UnitCell:
    def __init__(self, a: str, b: str, c: str, al: str, be: str, ga: str):
        self.a = a
        self.b = b
        self.c = c
        self.al = al
        self.be = be
        self.ga = ga

        self.a_f = float(a)
        self.b_f = float(b)
        self.c_f = float(c)
        self.al_f = float(al)
        self.be_f = float(be)
        self.ga_f = float(ga)

        self.symmetry = self.find_symmetry()

    def __str__(self):
        return f"a {self.a} b {self.b} c {self.c}\nal {self.al} be {self.be} ga {self.ga}"

    def __getitem__(self, item: str) -> str:
        if item == "a":
            return self.a
        elif item == "b":
            return self.b
        elif item == "c":
            return self.c
        elif item == "al":
            return self.al
        elif item == "be":
            return self.be
        elif item == "ga":
            return self.ga
        else:
            raise ValueError(f"Got {item=}. Was expecting 'a', 'b', 'c', 'al', 'be', or 'ga'.")

    def find_symmetry(self) -> str:

        if all_close(self.a_f, self.b_f, self.c_f) and all_close(self.al_f, self.be_f, self.ga_f):
            return "Cubic" if all_close(self.al_f, float("90.0")) else "Rhombohedral"

        elif all_close(self.a_f, self.b_f) and all_close(self.al_f, self.be_f, float("90.0")):
            if all_close(self.ga_f, float("90.0")):
                return "Tetragonal"
            if all_close(self.ga_f, float("120.0")):
                return "Hexagonal"  # or trigonal, but I can't tell just from cell prms, and this matches the TOPAS macro

        elif not all_close(self.a_f, self.b_f, self.c_f) and all_close(self.al_f, self.be_f, self.ga_f, float("90.0")):
            return "Orthorhombic"

        if all_close(self.al_f, self.ga_f, float("90.0")) and not all_close(self.be_f, float("90.0")):
            return "Monoclinic_be"
        elif all_close(self.al_f, self.be_f, float("90.0")) and not all_close(self.ga_f, float("90.0")):
            return "Monoclinic_ga"
        elif all_close(self.ga_f, self.be_f, float("90.0")) and not all_close(self.al_f, float("90.0")):
            return "Monoclinic_al"

        return "Triclinic"

    def topas_str(self, sym: str = None) -> str:
        sym = self.symmetry if sym is None else sym

        if sym == "Cubic":
            return f"\t\tCubic({self.a})\t'{self.a}"
        elif sym == "Rhombohedral":
            return f"\t\tRhombohedral({self.a}, {self.al})\t'{self.a}, {self.al}"
        elif sym == "Tetragonal":
            return f"\t\tTetragonal({self.a}, {self.c})\t'{self.a}, {self.c}"
        elif sym == "Hexagonal":
            return f"\t\tHexagonal({self.a}, {self.c})\t'{self.a}, {self.c}"
        elif sym == "Orthorhombic":
            return f"\t\ta {self.a}\t'{self.a}\n\t\tb {self.b}\t'{self.b}\n\t\tc {self.c}\t'{self.c}"
        elif sym[-3:] in ["_al", "_be", "_ga"]:
            ang = sym[-2:]
            return f"\t\ta  {self.a}\t'{self.a}\n\t\tb  {self.b}\t'{self.b}\n\t\tc  {self.c}\t'{self.c}\n\t\t{ang} {self[ang]}\t'{self[ang]}"
        else:
            return f"\t\ta {self.a}\t'{self.a}\n\t\tb {self.b}\t'{self.b}\n\t\tc {self.c}\t'{self.c}\n\t\tal {self.al}\t'{self.al}\n\t\tbe {self.be}\t'{self.be}\n\t\tga {self.ga}\t'{self.ga}"


def all_close(*vals: float, rel_tol=1e-09, abs_tol=0.0):
    return all(math.isclose(vals[0], vals[i], rel_tol=rel_tol, abs_tol=abs_tol) for i in range(1, len(vals)))


def write_str(cif_file: str, str_file=None, data: str = "all", work: bool=False):
    """
    Writes all structures in a given CIF to one file each from a given cif file..

    Args:
        cif_file: path to a CIF file as a string
        str_file: give an explicit name for the resultant STR.
                  If the str_file name includes a full path, then that is used.
                  If the str_file is a relative path, then this is joined to the
                    path of the CIF file
        data : do you want 'all', the 'first', or the 'last' data blocks in the cif?
               if "append", then all the structures from the cif are written to
                 one file
        work : add the special requirements we figured out for work-related strs
    Returns:
        None. Writes file to disk.
     """
    print(f"Now reading {cif_file}.")
    cif = CifFile.ReadCif(cif_file, permissive=False)

    if data == "first":
        data_keys = [cif.keys()[0]]
    elif data == "last":
        data_keys = [cif.keys()[-1]]
    else:
        data_keys = cif.keys()

    path = os.path.dirname(cif_file)

    get_new_output_filename = str_file is None  # do I need to update the filename on each loop through the datakeys?
    append_data = (data == "append")  # am I appending all strs to the same output file?
    update_output_filename = True  # by default, I want to update the output filename

    for d in data_keys:
        # check if there are unit cell prms in the data. If not, it probably isn't a structure
        #  and we can skip it, and try the next data, if any.

        if "_cell_length_a" not in cif[d]:
            print(f"No structure detected in datablock {d}. Trying the next block...")
            continue

        if get_new_output_filename and update_output_filename:
            str_file = f"{clean_filename(get_phasename(cif, d))}.str"

        update_output_filename = not append_data  # if I'm appending, I don't want to update the filename each time

        filename = str_file if os.path.isabs(str_file) else os.path.join(path, str_file)
        write_type = "a" if append_data else "w"
        s = create_str(cif, d, work)

        with open(filename, write_type) as f:
            print(f"Now writing {filename}.")
            f.write(s)

    print("Done.")


def create_str(cif: Dict, data: str, work=False) -> str:
    """
    Creates a single string formatted as an STR representing the given data block in the cif.

    Args:
        cif: PyCifRW representation of a cif file
        data: the key of the data block to convert to a STR
        work: add the work-related bits and bobs
    Returns:
        A single string formatted as a TOPAS STR.
     """
    s = "\tstr\n"
    s += f'\t\tphase_name "{get_phasename(cif, data)}"\n'
    s += "\t\tPhase_Density_g_on_cm3( 0)\n"

    if work:
        s += "\t\tCS_L_local(!csL ,200)\n"
        s += "\t\tCS_G_local(!csG ,200)\n"
        s += "\t\tStrain_L_local(!strL ,0.001)\n"
        s += "\t\tStrain_G_local(!strG ,0.001)\n"
    else:
        s += "\t\tCS_L( ,200)\n"
    s += "\t\tMVW(0,0,0)\n"
    s += "\t\tscale @ 0.0001\n"
    s += get_unitcell(cif, data) + "\n"
    s += f'\t\tspace_group "{get_spacegroup(cif, data)}"\n'
    s += get_atom_sites_string(cif, data) + '\n'

    return s


def strip_brackets(s: str) -> Union[str, None]:
    """
    Strips the brackets at the end of a string.
    The main use case is to remove the error values, eg '0.123(45)' --> '0.123'.
    None is returned as None

    Args:
        s: A single string.
    Returns:
        A string, with no brackets at the end of each string
     """
    if s is None:
        return None

    return s[:m] if (m := s.find("(")) != -1 else s


def change_NA_value(s: str) -> Union[str, None]:
    """
    If a string consists of a single question mark ('?') or full stop ('.'),
    it is replaced with a None
    These are normally markers of "no value recorded"

    None is returned as None

    Args:
        s: A string.
    Returns:
        A string, with None in place of a single question mark or full stop
     """
    return None if s in {"?", "."} else s


def val_to_frac(s: str) -> Union[str, None]:
    """
    Checks a single string consisting of numeric values for
    values which are consistent with fractional values.

    If those values are consistent with 1/6, 1/3, 2/3, or 5/6, then the decimal
    representation is altered to be a fractional representation.

    TOPAS requires atoms on special positions with values that have non-terminating
    decimal representations to be represented as fractions.

    None is returned as None

    Args:
        s: A  single string represent numerical values.
    Returns:
        A single string, with fractions representing 1/6, 1/3, 2/3, or 5/6
    """
    if s is None:
        return None

    fraction_patterns_strings = ((re.compile(r"^([-+]?)0?\.16{2,}[67]$"), "1/6"),
                                 (re.compile(r"^([-+]?)0?\.3{4,}$"),      "1/3"),
                                 (re.compile(r"^([-+]?)0?\.6{3,}[67]$"),  "2/3"),
                                 (re.compile(r"^([-+]?)0?\.83{3,}$"),     "5/6"))

    for pattern, fraction in fraction_patterns_strings:
        if m := pattern.match(s):
            sign = m.group(1)
            return f"={sign}{fraction};"
    return s  # if we get here, it didn't match anything, so return what we were given


def get_dict_entry_copy(dct: Dict, *keys: str, default=None) -> Union[Dict, List[str], str]:
    """
    Returns a deepcopy of the first dictionary value to match from an arbitrary number of given keys.
    If there is no match, the default value is returned.
    #https://stackoverflow.com/a/67187811/36061

    Args:
        dct: A dictionary
        *keys: an arbitrary number of keys in order of preference for return.
        default: The value to return if no matching entry is found in the dictionary.
    Returns:
        A deepcopy of the dictionary value associated with the first key to match, or the default value.
     """
    for key in keys:
        with contextlib.suppress(KeyError):
            return get_dict_entry_copy_throw_error(dct, key)
    return default


def get_dict_entry_copy_throw_error(dct: Dict, key: str) -> Union[Dict, List[str], str]:
    """
    Returns a deepcopy of the dictionary value from a matching key.

    Args:
        dct: A dictionary
        key: a key.
    Returns:
        A deepcopy of the dictionary value associated with the key.

    Raises:
        KeyError if the key is not present
     """
    return copy.deepcopy(dct[key])


def pad_string_list(lst: List[str], pad="post") -> Union[List[str], None]:
    """
    Given a list of string, pad the length of the strings
    such that they are all the same

    Parameters
    ----------
    lst : list of strings
    pad: is it a "post" (after the string) pad, or a "pre" (before the string) pad?

    Returns
    -------
    List of strings, all of the same length

    """
    if lst is None:
        return None

    if any(s.startswith("-") for s in lst):  # Prepend the string if it looks like a negative number
        lst = [s if s.startswith("-") else f" {s}" for s in lst]

    # what is the max str len?
    max_len = 0
    for s in lst:
        if len(s) > max_len:
            max_len = len(s)

    lst = [pad_string(s, max_len, pad) for s in lst]

    return lst


def pad_string(s: str, d: int, pad: str) -> Union[str, None]:
    """
    Add spaces to the end or start  of a string until it is the length given by d
    If len is already more than d, it just returns s

    Parameters
    ----------
    s : string you want to pad at the end or start with spaces
    d : Integer that you want the final string length to be
    pad: is it a "post" (after the string) pad, or a "pre" (before the string) pad?

    Returns
    -------
    s : String of length d. If len(s) < d originaly, then s is returned unchanged.

    """
    if s is None:
        return None

    if pad == "post":
        return f"{s:{d}}"
    elif pad == "pre":
        return f"{s:>{d}}"
    else:
        return f"{s}"


def clean_phasename(s: str) -> str:
    """
    Removes new lines, carriage returns, and leading/trailing whitespace from a
    string representing a phase name.
    Also removes single quote characters "'", as they are comment characters in TOPAS

    Args:
        s: a string

    Returns:
        A string with no newlines, carriage returns, or leading/trailing whitespace
     """
    return s.strip().replace('\n', '').replace('\r', '').replace("'", "_")


def clean_filename(s: str) -> str:
    """
    Replaces illegal characters from a string which is to be used as a filename.
    This doesn't deal with the path, just the name. Also, not the extension.

    eg 'string' could be used to construct 'C:\\important\\string.cif' at a later
    point in the program.

    / \\ | : * ? " < >  are replaced by "_"

    Whitespace is also replaced by "_"

    Also checks for illegal names like 'CON', and 'PRN'.

    Yes, this is Windows-centric, but so is TOPAS.

    Trying to follow guidelines in:
    https://docs.microsoft.com/en-us/windows/win32/fileio/naming-a-file#naming-conventions

    Args:
        s: a string

    Returns:
        A string with no illegal characters
     """
    ILLEGAL_NAMES = ("CON", "PRN", "AUX", "NUL", "COM1", "COM2", "COM3", "COM4",
                     "COM5", "COM6", "COM7", "COM8", "COM9", "LPT1", "LPT2", "LPT3",
                     "LPT4", "LPT5", "LPT6", "LPT7", "LPT8", "LPT9", ".", "..")

    if s.upper() in ILLEGAL_NAMES:
        s = f"{s}_"

    if s.endswith("."):
        s = f"{s[:-1]}_"

    s = re.sub(r"[/\\:*?\"'<>|\s]", "_", s)  # remove all the illegal chars
    s = re.sub(r"__+", "_", s)  # remove multiple underscores
    s = re.sub(r"^_", "", s)  # remove a leading underscore

    return s


def get_phasename(cif: Dict, data: str) -> str:
    """
    Returns a string to be used as the name of the phase.

    The name is taken from (in order of preference):
        "_pd_phase_name", "_chemical_name_mineral", "_chemical_name_common",
        "_chemical_name_systematic", or "_chemical_name_structure_type"
    The value of the data block is then appended to this value.

    Args:
        cif: a PyCifRW dictionary
        data: the key value for a data block in the CIF
    Returns:
        A string denoting the name of the phase
     """
    phasename = get_dict_entry_copy(cif[data], "_pd_phase_name",
                                    "_chemical_name_mineral",
                                    "_chemical_name_common",
                                    "_chemical_name_systematic",
                                    "_chemical_name_structure_type",
                                    default="")
    phasename = clean_phasename(phasename)

    return data if phasename in ("", ".", "?", None) else f"{phasename}_{data}"


def get_spacegroup(cif: Dict, data: str) -> str:
    """
    Returns a string to be used as the spacegroup.

    The space group  is taken from (in order of preference):
        "_symmetry_space_group_name_H-M", "_symmetry_Int_Tables_number",
        "_space_group_name_H-M_alt", or "_space_group_IT_number",

    No validation is done on the space group.

    Args:
        cif: a PyCifRW dictionary
        data: the key value for a data block in the CIF
    Returns:
        A string denoting the space group of the phase, as given in the CIF.
     """
    spacegroup = get_dict_entry_copy(cif[data], "_symmetry_space_group_name_H-M",
                                     "_space_group_name_H-M_alt",
                                     "_symmetry_Int_Tables_number",
                                     "_space_group_IT_number",
                                     default="")
    if spacegroup.isdigit():
        print("Spacegroup given by number. Check that the SG setting matches that of the atom coordinates.")

    return spacegroup


def get_unitcell(cif: Dict, data: str, sym: str = None) -> str:
    """
    Returns a string giving the unit cell parameters of the structure.
    Some pattern matching is carried out to return TOPAS macros for
    unit cell prms, but no validation or space group consistency is done.

    Args:
        cif: a PyCifRW dictionary
        data: the key value for a data block in the CIF
        sym: Uses this symmetry to print out the results. Use "Triclinic" to print the general case.
              Use None to let the symmetry be determined by the unit cell.
    Returns:
        A string containing the unit cell parameters in STR format.
    Raises:
        KeyError: if any of the unit cell parameters are not present in the datablock.
    """
    a_s = strip_brackets(get_dict_entry_copy_throw_error(cif[data], "_cell_length_a"))
    b_s = strip_brackets(get_dict_entry_copy_throw_error(cif[data], "_cell_length_b"))
    c_s = strip_brackets(get_dict_entry_copy_throw_error(cif[data], "_cell_length_c"))
    al_s = strip_brackets(get_dict_entry_copy_throw_error(cif[data], "_cell_angle_alpha"))
    be_s = strip_brackets(get_dict_entry_copy_throw_error(cif[data], "_cell_angle_beta"))
    ga_s = strip_brackets(get_dict_entry_copy_throw_error(cif[data], "_cell_angle_gamma"))

    unit_cell = UnitCell(a_s, b_s, c_s, al_s, be_s, ga_s)

    return unit_cell.topas_str(sym)


def get_atom_sites_string(cif: Dict, data: str) -> str:
    """
    Returns a string giving the site parameters for all of the sites in the CIF.

    They are formatted as necessary for TOPAS. A single site consists of:
        "site _1_ num_posns 0 x _2_ y _3_ z _4_ occ _5_ _6_ beq _7_"

    1 - the site label is taken directly from "_atom_site_label"
    2 - the fractional x ordinate is taken from "_atom_site_fract_x"
    3 - the fractional y ordinate is taken from "_atom_site_fract_y"
    4 - the fractional z ordinate is taken from "_atom_site_fract_z"
    5 - the atom type is taken from "_atom_site_type_symbol" or inferred from the site label
    6 - the occupancy is taken from "_atom_site_occupancy" or given the value "1"
    7 - the isotropic displacement parameter,B, is taken from "_atom_site_B_iso_or_equiv"
        or 8*Pi^2*"_atom_site_U_iso_or_equiv". If anisotropic values are present as
        "_atom_site_aniso_B_11", "_22", "_33", then their average is used. If
        "_atom_site_aniso_U_11", "_22", "_33", then 8*Pi^2 times their average is used. If
        no values are present the the value "1" is assigned.

    If the decimal values of the fractional coordinates are consistent with the fractions
    1/6, 1/3, 2/3, or 5/6, then the decimal value is replaced by the fractional representation.

    Args:
        cif: a PyCifRW dictionary
        data: the key value for a data block in the CIF
    Returns:
        A string containing the atomic site parameters in STR format.
    Raises:
        KeyError: if any of the site label or fractional coordinate keys are not present.
    """
    labels = pad_string_list(get_dict_entry_copy_throw_error(cif[data], "_atom_site_label"))

    # to val_to_frac first, as a real match wouldn't have any errors to strip
    x = pad_string_list([strip_brackets(val_to_frac(i)) for i in cif[data]["_atom_site_fract_x"]])
    y = pad_string_list([strip_brackets(val_to_frac(i)) for i in cif[data]["_atom_site_fract_y"]])
    z = pad_string_list([strip_brackets(val_to_frac(i)) for i in cif[data]["_atom_site_fract_z"]])

    # type of atom
    try:
        atom_symbols = get_dict_entry_copy_throw_error(cif[data], "_atom_site_type_symbol")
        atom_labels = get_dict_entry_copy_throw_error(cif[data], "_atom_site_label")
        atoms = [fix_atom_type(symbol, label) for symbol, label in zip(atom_symbols, atom_labels)]
    except KeyError:
        print("Warning! Atom types inferred from site labels. Please check for correctness.")
        atoms = [convert_site_label_to_atom(label) for label in labels]
    atoms = pad_string_list(atoms)

    # occupancy
    try:
        occ = [strip_brackets(i) for i in cif[data]["_atom_site_occupancy"]]
    except KeyError:
        occ = ["1"] * len(labels)
    occ = pad_string_list(occ)

    # ADPs
    b_iso = pad_string_list(get_beq(cif, data))

    r = [make_atom_site_string(label, xf, yf, zf, atom, oc, b) for label, xf, yf, zf, atom, oc, b in zip(labels, x, y, z, atoms, occ, b_iso)]

    return "\n".join(r)


def make_atom_site_string(label: str, x: str, y: str, z: str, atom: str, occ: str, beq: str) -> str:
    """
    Returns a string giving the site parameters necessary for TOPAS:
        "site _1_ num_posns 0 x _2_ y _3_ z _4_ occ _5_ _6_ beq _7_"

    This is used by get_atom_site_string().

    Args:
        label: a string to describe the site (1)
        x: the fractional x ordinate (2)
        y: the fraction y ordinate (3)
        z: the fractional z ordinate (4)
        atom: the scattering factors to use for this site (5)
        occ: the statistical occupancy of this site (6)
        beq: the isotropic atomic displacement factor for this site (B) (7)

    Returns:
        A string containing the atomic site parameters in STR format.
    """
    return f"\t\tsite {label} num_posns 0\tx {x} y {y} z {z} occ {atom} {occ} beq {beq}"


def convert_site_label_to_atom(s: str) -> str:
    """
    Given a site label string, the first 2 or 1 characters are checked
    (case-sensitive) to see if they match an element symbols. If they do,
    this symbol is returned. If not, the entire site label is returned.

    If "W" is detected, then a warning is printed re tungsten or water.

    Args:
        s: a string denoting the site label

    Returns:
        A string containing the element symbol, or the entire label.
    """
    ELEMENTS = ("Ac", "Ag", "Al", "Am", "Ar", "As", "At", "Au", "Ba", "B",
                "Be", "Bi", "Br", "Ca", "Cd", "C", "Ce", "Cl", "Co", "Cr",
                "Cs", "Cu", "Dy", "Er", "Eu", "F", "Fe", "Fr", "Ga", "Gd",
                "Ge", "H", "He", "Hg", "Hf", "Ho", "I", "In", "Ir", "K",
                "Kr", "La", "Li", "Lu", "Mg", "Mn", "Mo", "N", "Na", "Nb",
                "Nd", "Ne", "Ni", "O", "Os", "Pb", "P", "Pa", "Pd", "Po",
                "Pr", "Pm", "Pt", "Ra", "Rb", "Re", "Rh", "Rn", "Ru", "S",
                "Sb", "Sc", "Sm", "Se", "Si", "Sn", "Sr", "Ta", "Tb", "Tc",
                "Te", "Th", "Tl", "Ti", "Tm", "W", "U", "V", "Xe", "Y",
                "Yb", "Zn", "Zr", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "D")

    WATER = ("Wat", "WAT", "wat")

    r = s[:3]
    if r in WATER:
        print(f"Site label '{s}' probably means 'water'. Please check that this atom really is oxygen.")
        return "O"

    r = s[:2]
    if r in ELEMENTS:
        return r

    r = s[:1]
    if r in ELEMENTS:
        if r == "W":
            print(f"W detected for site '{s}'. Do you mean tungsten or oxygen from a water molecule? Please check.")
        return r

    print(f"Can't decide what atom the site label '{s}' should be. Please check it.")
    return s  # if everything fails, just return what the label was


def fix_atom_type(a: str, orig_site_label=None) -> str:
    """
    Given an atom type from "_atom_site_type_symbol", any charge given is sometimes the wrong
    way around for TOPAS to handle. This function uses a regex to switch and given charge
    to the correct format, whilst leaving neutral atoms as they were.

    For example
    Cu   --> Cu
     I   -->  I
     K+  -->  K+1
     V3+ -->  V+3
    Fe2+ --> Fe+2
     O1- -->  O-1

    Furthermore, the charge is checked against the list of scattering factors supported
    by TOPAS, and if a mismatch is detected, a charged atom is replaced by the neutral one.
    eg: H+1 --> H

    Args:
        a: a string to denoting the atom type, with or without a charge
        orig_site_label: the site label for this atom, so I can use it to print error msgs, if needed.

    Returns:
        A string containing the atom type with the charge in the correct order, if one is given.
    """
    regex = re.search(r"^([A-Z][a-z]?)(\d{0,2})([+-]?)(\d{0,2})$", a)

    symbol = regex[1]
    charge = regex[2]
    sign = regex[3]
    digit = regex[4]  # if there is a trailing digit, then that is probably the charge

    if charge == "0":  # ie Si0+ is a valid symbol, but it is a neutral atom
        charge = ""
        sign = ""

    # check for a single sign with no charge. eg F-. Needs to return F-1
    if len(sign) == 1 and len(charge) == 0:
        if len(digit) == 0:
            charge = "1"
        else:  # the atom was probably the right way around to begin with
            return a

    return convert_atom_type_to_topas(symbol + sign + charge, orig_site_label)


def convert_atom_type_to_topas(a: str, orig_site_label=None) -> str:
    """
    Given an atom type in the correct TOPAS site format, the CIF dictionary can allow charges
    for which there are no scattering factors..

    The charge is checked against the list of scattering factors supported by TOPAS, and if
    a mismatch is detected, a charged atom is replaced by the neutral one.
    eg: H+1 --> H

    Args:
        a: a string to denoting the atom type, with or without a charge
        orig_site_label: the site label for this atom, so I can use it to print error msgs, if needed.

    Returns:
        A string containing the atom type with the charge in the correct order, if one is given.
    """
    ALLOWED = ("D", "H", "H-1", "D-1", "He", "Li", "Li+1", "Be", "Be+2", "B", "C", "N",
               "O", "O-1", "O-2", "F", "F-1", "Ne", "Na", "Na+1", "Mg", "Mg+2", "Al",
               "Al+3", "Si", "Si+4", "P", "S", "Cl", "Cl-1", "Ar", "K", "K+1", "Ca",
               "Ca+2", "Sc", "Sc+3", "Ti", "Ti+2", "Ti+3", "Ti+4", "V", "V+2", "V+3",
               "V+5", "Cr", "Cr+2", "Cr+3", "Mn", "Mn+2", "Mn+3", "Mn+4", "Fe", "Fe+2",
               "Fe+3", "Co", "Co+2", "Co+3", "Ni", "Ni+2", "Ni+3", "Cu", "Cu+1", "Cu+2",
               "Zn", "Zn+2", "Ga", "Ga+3", "Ge", "Ge+4", "As", "Se", "Br", "Br-1", "Kr",
               "Rb", "Rb+1", "Sr", "Sr+2", "Y", "Y+3", "Zr", "Zr+4", "Nb", "Nb+3", "Nb+5",
               "Mo", "Mo+3", "Mo+5", "Mo+6", "Tc", "Ru", "Ru+3", "Ru+4", "Rh", "Rh+3",
               "Rh+4", "Pd", "Pd+2", "Pd+4", "Ag", "Ag+1", "Ag+2", "Cd", "Cd+2", "In",
               "In+3", "Sn", "Sn+2", "Sn+4", "Sb", "Sb+3", "Sb+5", "Te", "I", "I-1", "Xe",
               "Cs", "Cs+1", "Ba", "Ba+2", "La", "La+3", "Ce", "Ce+3", "Ce+4", "Pr", "Pr+3",
               "Pr+4", "Nd", "Nd+3", "Pm", "Pm+3", "Sm", "Sm+3", "Eu", "Eu+2", "Eu+3", "Gd",
               "Gd+3", "Tb", "Tb+3", "Dy", "Dy+3", "Ho", "Ho+3", "Er", "Er+3", "Tm", "Tm+3",
               "Yb", "Yb+2", "Yb+3", "Lu", "Lu+3", "Hf", "Hf+4", "Ta", "Ta+5", "W", "W+6",
               "Re", "Os", "Os+4", "Ir", "Ir+3", "Ir+4", "Pt", "Pt+2", "Pt+4", "Au", "Au+1",
               "Au+3", "Hg", "Hg+1", "Hg+2", "Tl", "Tl+1", "Tl+3", "Pb", "Pb+2", "Pb+4",
               "Bi", "Bi+3", "Bi+5", "Po", "At", "Rn", "Fr", "Ra", "Ra+2", "Ac", "Ac+3",
               "Th", "Th+4", "Pa", "U", "U+3", "U+4", "U+6", "Np", "Np+3", "Np+4", "Np+6",
               "Pu", "Pu+3", "Pu+4", "Pu+6", "Am", "Cm", "Bk", "Cf")

    if a in ALLOWED:
        return a

    # if we get here, the atom doesn't exist in the tuple.
    regex = re.search(r"^([A-Z][a-z]?)([+-]?)(\d{0,2})$", a)  # Cu+2
    symbol = regex[1]

    insert_text = "" if orig_site_label is None else f" in site {orig_site_label}"

    print(f"{a} is not a legal TOPAS scattering factor. Atom{insert_text} replaced with {symbol}.")

    return symbol


def get_b_iso(cif: Dict, data: str) -> List[str]:
    """
    Returns a list of strings giving the isotropic atomic displacement parameters, B, of
    all atomic sites taken from "_atom_site_B_iso_or_equiv".

    Used by get_beq().

    Args:
        cif: a PyCifRW dictionary
        data: the key value for a data block in the CIF
    Returns:
        A list containing the isotropic atomic displacement parameters, B, of all atomic sites.
    Raises:
        KeyError: if the key "_atom_site_B_iso_or_equiv" is not present.
    """
    return [change_NA_value(strip_brackets(i)) for i in cif[data]["_atom_site_B_iso_or_equiv"]]


def convert_u_to_b(s: str) -> Union[str, None]:
    """
    Take a string representing a U value and return a string representing a B value
    B = 8*Pi^2 * U

    Parameters
    ----------
    s : a string representing a U value. Could be None

    Returns
    -------
    A string representing a B value. Could be None.
    """
    if s is None:
        return None

    s = float(s)
    s *= 8 * math.pi ** 2
    s = str(round(s, 3))  # round B value to 3 d.p.
    return s


def get_u_iso(cif: Dict, data: str) -> List[str]:
    """
    Returns a list of strings giving the isotropic atomic displacement parameters, B, of
    all atomic sites taken from "_atom_site_U_iso_or_equiv" found by multiplying their
    value by 8*Pi^2

    Used by get_beq().

    Args:
        cif: a PyCifRW dictionary
        data: the key value for a data block in the CIF
    Returns:
        A list containing the isotropic atomic displacement parameters, B, of all atomic sites.
    Raises:
        KeyError: if the key "_atom_site_U_iso_or_equiv" is not present.
    """

    return [convert_u_to_b(change_NA_value(strip_brackets(i))) for i in cif[data]["_atom_site_U_iso_or_equiv"]]


def get_beq(cif: Dict, data: str) -> List[str]:
    """
    Looks the Biso, Uiso, Baniso, Uaniso (in that order) to get beq values. If they are missing,
    then a default value of "1." is used.

    Args:
        cif: a PyCifRW dictionary
        data: the key value for a data block in the CIF
    Returns:
        A list containing the isotropic atomic displacement parameters, B, of all atomic sites.
    """
    labels = get_dict_entry_copy(cif[data], "_atom_site_label")  # I want everything to be in the same order as the site_label

    # get all potential ADPs from the str as dictionaries. The key is the atom label, the value is the ADP value
    beq_types = ["b_iso", "u_iso", "b_aniso", "u_aniso", "be_aniso"]  # the order here matters. If it matches one (in the loop below), then the others aren't tested.
    b_dicts = [make_b_dict(cif, data, beq_type) for beq_type in beq_types]

    r = []
    for label in labels:
        b = None
        for d in b_dicts:
            if (b := d.get(label)) is not None:  # I've found the value
                break

        if b is None or float(b) == 0.0:  # missing or zero value
            print(f"Warning! beq value missing or zero for site {label}! Default value of 1 entered")
            b = "1."
        elif b.startswith("-"):
            print(f"Warning! Negative atomic displacement parameter detected for site {label}! Please check.")
        r.append(b)  # the ADP values are now in the same order as the site label

    return r


def make_b_dict(cif: Dict, data: str, b_type: str) -> Dict:
    """
    A helper function for get_beq() to get the site labels and associated ADP values, and
    return a dictionary with labels as keys and values as ADP values. If the value is None,
    or 0, then that key is removed.

    Parameters
    ----------
    cif : PyCifRW dictionary

    data : the data block you're looking at

    b_type: "b_iso", "u_iso", "b_aniso", "u_aniso", or "be_aniso"

    Returns
    -------
    Dictionary containing site labels as keys and ADPs as values.
    If a value is None, it is removed from the dictionary.
    It is possible to return an empty dictionary.
    If the b_type selected is not in the cif, an empty dictionary is returned

    Raises
    ------
    ValueError if b_type is not one of "b_iso", "u_iso", "b_aniso", "u_aniso", or "be_aniso"
    """

    b_info = {"b_iso": (get_b_iso, True), "u_iso": (get_u_iso, True),
              "b_aniso": (get_b_aniso, False), "u_aniso": (get_u_aniso, False),
              "be_aniso": (get_beta_aniso_as_b, False)}

    try:
        f, iso = b_info[b_type]
    except KeyError as e:
        raise ValueError(f'Invalid choice. Got "{b_type}", expecting "b_iso", "u_iso", "b_aniso", "u_aniso", or "be_aniso".') from e

    try:
        b_list = f(cif, data)
        tag = "_atom_site_label" if iso else "_atom_site_aniso_label"
        b_labels = get_dict_entry_copy(cif[data], tag)
    except KeyError:
        return {}

    d = dict(zip(b_labels, b_list))

    # filter out any None/zero values
    d = {k: v for k, v in d.items() if v is not None or v != "0.0"}

    return d


def get_b_aniso(cif: Dict, data: str) -> List[str]:
    """
    Returns a list of strings giving the equivalent isotropic atomic displacement parameters, B,
    of all atomic sites by taking the average value of "_atom_site_aniso_B_11", "_22", and "_33".

    Used by getAniso().

    Args:
        cif: a PyCifRW dictionary
        data: the key value for a data block in the CIF
    Returns:
        A list containing the isotropic atomic displacement parameters, B, of all atomic sites.
    Raises:
        KeyError: if any of the keys "_atom_site_aniso_B_11", "_22", or "_33" is not present.
    """
    # convert the str lists to float lists
    B11 = [float(strip_brackets(i)) for i in cif[data]["_atom_site_aniso_B_11"]]
    B22 = [float(strip_brackets(i)) for i in cif[data]["_atom_site_aniso_B_22"]]
    B33 = [float(strip_brackets(i)) for i in cif[data]["_atom_site_aniso_B_33"]]

    b_equiv = [(b11 + b22 + b33) / 3.0 for b11, b22, b33 in zip(B11, B22, B33)]  # get the average of the three values
    b_equiv = [str(round(b, 3)) for b in b_equiv]  # round to 3 d.p.

    return b_equiv


def get_u_aniso(cif: Dict, data: str) -> List[str]:
    """
    Returns a list of strings giving the equivalent isotropic atomic displacement parameters, B,
    of all atomic sites by multiplying the average value of "_atom_site_aniso_U_11", "_22",
    and "_33" by 8*Pi^2.

    Used by getAniso().

    Args:
        cif: a PyCifRW dictionary
        data: the key value for a data block in the CIF
    Returns:
        A list containing the isotropic atomic displacement parameters, B, of all atomic sites.
    Raises:
        KeyError: if any the keys "_atom_site_aniso_U_11", "_22", or "_33" is not present.
    """
    # convert the str lists to float lists
    U11 = [float(strip_brackets(i)) for i in cif[data]["_atom_site_aniso_U_11"]]
    U22 = [float(strip_brackets(i)) for i in cif[data]["_atom_site_aniso_U_22"]]
    U33 = [float(strip_brackets(i)) for i in cif[data]["_atom_site_aniso_U_33"]]

    return [convert_u_to_b(str((u11 + u22 + u33) / 3.0)) for u11, u22, u33 in zip(U11, U22, U33)]


def get_beta_aniso_as_b(cif: Dict, data: str) -> List[str]:
    """
    Returns a list of strings giving the equivalent isotropic atomic displacement parameters, B.
    Converts beta values to B values, and then takes the average value of "_atom_site_aniso_beta_11", "_22", and "_33".

    beta_11 = B_11 * (a*^2/4)
    beta_22 = B_22 * (b*^2/4)
    beta_33 = B_33 * (c*^2/4)
    beta_23 = B_23 * ((b* c*)/4)
    beta_31 = B_31 * ((c* a*)/4)
    beta_23 = B_23 * ((a* b*)/4)

    a* == reciprocal cell length  - ie |ã*|

    Used by getAniso().

    Args:
        cif: a PyCifRW dictionary
        data: the key value for a data block in the CIF
    Returns:
        A list containing the isotropic atomic displacement parameters, B, of all atomic sites.
    Raises:
        KeyError: if any of the keys "_atom_site_aniso_B_11", "_22", or "_33" is not present, or any of the "_cell_length_" or "_cell_angle" keys.
    """
    # convert the str lists to float lists
    be11 = [float(strip_brackets(i)) for i in cif[data]["_atom_site_aniso_beta_11"]]
    be22 = [float(strip_brackets(i)) for i in cif[data]["_atom_site_aniso_beta_22"]]
    be33 = [float(strip_brackets(i)) for i in cif[data]["_atom_site_aniso_beta_33"]]

    a = strip_brackets(get_dict_entry_copy_throw_error(cif[data], "_cell_length_a"))
    b = strip_brackets(get_dict_entry_copy_throw_error(cif[data], "_cell_length_b"))
    c = strip_brackets(get_dict_entry_copy_throw_error(cif[data], "_cell_length_c"))
    al = strip_brackets(get_dict_entry_copy_throw_error(cif[data], "_cell_angle_alpha"))
    be = strip_brackets(get_dict_entry_copy_throw_error(cif[data], "_cell_angle_beta"))
    ga = strip_brackets(get_dict_entry_copy_throw_error(cif[data], "_cell_angle_gamma"))

    as_v, bs_v, cs_v = make_reciprocal_unit_cell_vectors(a, b, c, al, be, ga)

    # convert be values to b values
    be11 = [4 * v / as_v.square_magnitude() for v in be11]
    be22 = [4 * v / bs_v.square_magnitude() for v in be22]
    be33 = [4 * v / cs_v.square_magnitude() for v in be33]

    b_equiv = [(b11 + b22 + b33) / 3.0 for b11, b22, b33 in zip(be11, be22, be33)]  # get the average of the three values
    b_equiv = [str(round(b, 3)) for b in b_equiv]  # round to 3 d.p.

    return b_equiv


def make_reciprocal_unit_cell_vectors(a: str, b: str, c: str, al: str, be: str, ga: str) -> Tuple[Vector3, Vector3, Vector3]:
    a_v, b_v, c_v = make_unit_cell_vectors(a, b, c, al, be, ga)
    V = a_v.dot(b_v.cross(c_v))

    cas_v = b_v.cross(c_v)
    cbs_v = c_v.cross(a_v)
    ccs_v = a_v.cross(b_v)

    as_v = cas_v / V
    bs_v = cbs_v / V
    cs_v = ccs_v / V

    return as_v, bs_v, cs_v


def make_unit_cell_vectors(a: str, b: str, c: str, al: str, be: str, ga: str) -> Tuple[Vector3, Vector3, Vector3]:
    # convert strings to float
    a_f = float(strip_brackets(a))
    b_f = float(strip_brackets(b))
    c_f = float(strip_brackets(c))
    al_f = math.radians(float(strip_brackets(al)))
    be_f = math.radians(float(strip_brackets(be)))
    ga_f = math.radians(float(strip_brackets(ga)))

    #  this comes from eqn 50 in Trueblood et al Acta Cryst 1996 A52 770-781
    a_v = Vector3(a_f, 0.0, 0.0)
    b_v = Vector3(b_f * math.cos(ga_f), b_f * math.sin(ga_f), 0.0)
    c_v = Vector3(c_f * math.cos(be_f),
                  c_f * (math.cos(al_f) - (math.cos(ga_f) * math.cos(be_f))) / math.sin(ga_f),
                  c_f * math.sqrt(1 + 2 * math.cos(al_f) * math.cos(be_f) * math.cos(ga_f) - math.cos(al_f)**2 - math.cos(be_f)**2 - math.cos(ga_f)**2) / math.sin(ga_f))
    return a_v, b_v, c_v


if __name__ == "__main__":
    print("Hello, world!")

