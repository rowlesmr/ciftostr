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


# TODO: add test
def write_str(cif_file, str_file=None, data="all"):
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
    Returns:
        None. Writes file to disk.
     """
    print(f"Now reading {cif_file}.")
    cif = CifFile.ReadCif(cif_file)

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
            str_file = clean_filename(get_phasename(cif, d)) + ".str"

        update_output_filename = not append_data  # if I'm appending, I don't want to update the filename each time

        if os.path.isabs(str_file):  # if the given strfile is a full path, use it, don't change it.
            f = str_file
        else:  # otherwise, take the path from the cif file
            f = os.path.join(path, str_file)

        if append_data:
            file_to_write = open(f, "a")
        else:
            file_to_write = open(f, "w")

        s = create_str(cif, d)

        print(f"Now writing {file_to_write.name}.")
        file_to_write.write(s)
        file_to_write.close()

    print("Done.")


# TODO: add test
def create_str(cif, data):
    """
    Creates a single string formatted as an STR representing the given data block in the cif.

    Args:
        cif: PyCifRW representation of a cif file
        data: the key of the data block to convert to a STR
    Returns:
        A single string formatted as a TOPAS STR.
     """
    s = "\tstr\n"
    s += f'\t\tphase_name "{get_phasename(cif, data)}"\n'
    s += "\t\tPhase_Density_g_on_cm3( 0)\n"
    s += "\t\tCS_L( ,200)\n"
    s += "\t\tMVW(0,0,0)\n"
    s += "\t\tscale @ 0.0001\n"
    s += get_unitcell(cif, data) + "\n"
    s += f'\t\tspace_group "{get_spacegroup(cif, data)}"\n'
    s += get_atom_sites_string(cif, data) + '\n'

    return s


def strip_brackets(s):
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

    try:
        bracket = s.index("(")
    except ValueError:  # if there isn't a bracket, just keep going.
        bracket = len(s)

    return s[0:bracket]


def change_NA_value(s):
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
    if s in ("?", "."):
        s = None
    return s


def val_to_frac(s):
    """
    Checks a single string consisting of numeric values for
    values which are consistent with fractional values.

    If those values are consistent with 1/6, 1/3, 2/3, or 5/6, then the decimal
    representation is altered to be a fractional representation.

    TOPAS requires atoms on special positions with values that have non-terminating
    decimal representations to have be represented as fractions.

    None is returned as None

    Args:
        s: A  single string represent numerical values.
    Returns:
        A single string, with fractions representing 1/6, 1/3, 2/3, or 5/6
    """
    ONE_SIXTH = ["0.1666", "0.16666", "0.166666", "0.1666666",
                 "0.1667", "0.16667", "0.166667", "0.1666667"]
    ONE_THIRD = ["0.3333", "0.33333", "0.333333", "0.3333333"]
    TWO_THIRD = ["0.6666", "0.66666", "0.666666", "0.6666666",
                 "0.6667", "0.66667", "0.666667", "0.6666667"]
    FIVE_SIXTH = ["0.8333", "0.83333", "0.833333", "0.8333333"]
    NEG_ONE_SIXTH = [f"-{s}" for s in ONE_SIXTH]
    NEG_ONE_THIRD = [f"-{s}" for s in ONE_THIRD]
    NEG_TWO_THIRD = [f"-{s}" for s in TWO_THIRD]
    NEG_FIVE_SIXTH = [f"-{s}" for s in FIVE_SIXTH]

    ONE_SIXTH_FRAC = ["=1/6;"] * len(ONE_SIXTH)
    ONE_THIRD_FRAC = ["=1/3;"] * len(ONE_THIRD)
    TWO_THIRD_FRAC = ["=2/3;"] * len(TWO_THIRD)
    FIVE_SIXTH_FRAC = ["=5/6;"] * len(FIVE_SIXTH)
    NEG_ONE_SIXTH_FRAC = ["=-1/6;"] * len(ONE_SIXTH)
    NEG_ONE_THIRD_FRAC = ["=-1/3;"] * len(ONE_THIRD)
    NEG_TWO_THIRD_FRAC = ["=-2/3;"] * len(TWO_THIRD)
    NEG_FIVE_SIXTH_FRAC = ["=-5/6;"] * len(FIVE_SIXTH)

    ONE_SIXTH_d = dict(zip(ONE_SIXTH + NEG_ONE_SIXTH, ONE_SIXTH_FRAC + NEG_ONE_SIXTH_FRAC))
    ONE_THIRD_d = dict(zip(ONE_THIRD + NEG_ONE_THIRD, ONE_THIRD_FRAC + NEG_ONE_THIRD_FRAC))
    TWO_THIRD_d = dict(zip(TWO_THIRD + NEG_TWO_THIRD, TWO_THIRD_FRAC + NEG_TWO_THIRD_FRAC))
    FIVE_SIXTH_d = dict(zip(FIVE_SIXTH + NEG_FIVE_SIXTH, FIVE_SIXTH_FRAC + NEG_FIVE_SIXTH_FRAC))

    fracs = {}
    fracs.update(ONE_SIXTH_d)
    fracs.update(ONE_THIRD_d)
    fracs.update(TWO_THIRD_d)
    fracs.update(FIVE_SIXTH_d)

    r = None

    if s is None:
        return r

    try:
        r = fracs[s]
        print(f"Fractional atomic coordinate '{s}' detected and replaced with '{r}'.")
    except KeyError:
        r = s  # s isn't a fraction. probably.

    return r


def get_dict_entry_copy(dct, *keys, default=None):
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
        try:
            return get_dict_entry_copy_throw_error(dct, key)
        except KeyError:
            pass
    return default


def get_dict_entry_copy_throw_error(dct, key):
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


def pad_string_list(l, pad="post"):
    """
    Given a list of atomic ordinates, labels, or the like, pad the length of the strings
    such that they are all the same

    Parameters
    ----------
    l : list of strings
    pad: is it a "post" (after the string) pad, or a "pre" (before the string) pad?

    Returns
    -------
    List of strings, all of the same length

    """
    if l is None:
        return None

    # if any strings start with '-', it's probably a negative number, and I need to prepend
    #  a space to those that don't start with a '-'.
    negative_string = False

    for s in l:
        if s.startswith("-"):
            negative_string = True
            break  # only need to see if one is negative

    if negative_string:  # Prepend the string if it looks like a negative number
        l = [" " + s if not s.startswith("-") else s for s in l]

    # what is the max str len?
    max_len = 0
    for s in l:
        if len(s) > max_len:
            max_len = len(s)

    l = [pad_string(s, max_len, pad) for s in l]

    return l


def pad_string(s, d, pad):
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

    while len(s) < d:
        if pad == "post":
            s = s + " "
        elif pad == "pre":
            s = " " + s
        else:
            return s
    return s


def count_nones(l):
    """
    Counts the number of Nones in an iterable (eg list)

    Parameters
    ----------
    l : an iterable

    Returns
    -------
    num : integer with the number of Nones is the iterable

    """
    num = 0
    for s in l:
        if s is None:
            num += 1
    return num


def clean_phasename(s):
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


def clean_filename(s):
    """
    Replaces illegal characters from a string which is to be used as a filename.
    This doesn't deal with the path, just the name. Also, not the extension.

    eg 'string' could be used to construct 'C:\important\string.cif' at a later
    point in the program.

    / \ | : * ? " < >  are replaced by "_"

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
        s = s + "_"

    if s.endswith("."):
        s = s[:-1] + "_"

    return re.sub("[\\\\/:*?\"<>|\s]", "_", s)


def get_phasename(cif, data):
    """
    Returns a string to be used as the name of the phase.

    The name is taken from (in order of preference):
        "_chemical_name_mineral", "_chemical_name_common",
        "_chemical_name_systematic", or "_chemical_name_structure_type"
    The value of the data block is then appended to this value.

    Args:
        cif: a PyCifRW dictionary
        data: the key value for a data block in the CIF
    Returns:
        A string denoting the name of the phase
     """
    phasename = get_dict_entry_copy(cif[data], "_chemical_name_mineral",
                                    "_chemical_name_common",
                                    "_chemical_name_systematic",
                                    "_chemical_name_structure_type",
                                    default="")
    phasename = clean_phasename(phasename)

    r = ""
    # check that all the types of non-name names are accounted for.
    if phasename in ("", ".", "?", None):
        r = data
    else:
        r = f"{phasename}_{data}"

    return r


def get_spacegroup(cif, data):
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
    try:
        _ = int(spacegroup)
        print("Spacegroup given by number. Check that the SG setting matches that of the atom coordinates.")
    except ValueError:
        pass  # if it gets here, it wasn't a SG number...

    return spacegroup


def get_unitcell(cif, data):
    """
    Returns a string giving the unit cell parameters of the structure.
    Some pattern matching is carried out to return TOPAS macros for
    unit cell prms, but no validation or space group consistency is done.

    see also get_unitcell2() - it makes no attempt to pattern-match

    Args:
        cif: a PyCifRW dictionary
        data: the key value for a data block in the CIF
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

    a = float(a_s)
    b = float(b_s)
    c = float(c_s)
    al = float(al_s)
    be = float(be_s)
    ga = float(ga_s)

    s = ""
    if a == b and b == c:  # cubic or rhombohedral
        if al == be and be == ga and al == float("90"):  # cubic
            s = f"\t\tCubic({a_s})\t'{a_s}"
        if al == be and be == ga and al != float("90"):  # rhombohedral
            s = f"\t\tRhombohedral({a_s}, {al_s})\t'{a_s}, {al_s}"

    elif a == b and b != c:  # tetragonal or hexagonal/trigonal
        if al == be and be == ga and al == float("90"):  # tetragonal
            s = f"\t\tTetragonal({a_s}, {c_s})\t'{a_s}, {c_s}"
        if al == be and al == float("90") and ga == float("120"):  # hexagonal or trigonal
            s = f"\t\tHexagonal({a_s}, {c_s})\t'{a_s}, {c_s}"

    elif a != b and a != c and b != c:  # ortho, mono, tric
        if al == be and be == ga and al == float("90"):  # ortho
            s = f"\t\ta {a_s}\t'{a_s}\n\t\tb {b_s}\t'{b_s}\n\t\tc {c_s}\t'{c_s}"
        if al != be and al != ga and be != ga:  # tric
            s = f"\t\ta  {a_s}\t'{a_s}\n\t\tb  {b_s}\t'{b_s}\n\t\tc  {c_s}\t'{c_s}\n\t\tal {al_s}\t'{al_s}\n\t\tbe {be_s}\t'{be_s}\n\t\tga {ga_s}\t'{ga_s}"
        if al == be and al != ga and al == float("90"):  # mono_1
            s = f"\t\ta  {a_s}\t'{a_s}\n\t\tb  {b_s}\t'{b_s}\n\t\tc  {c_s}\t'{c_s}\n\t\tga {ga_s}\t'{ga_s}"
        if al == ga and al != be and al == float("90"):  # mono_2
            s = f"\t\ta  {a_s}\t'{a_s}\n\t\tb  {b_s}\t'{b_s}\n\t\tc  {c_s}\t'{c_s}\n\t\tbe {be_s}\t'{be_s}"
        if be == ga and be != al and be == float("90"):  # mono_3
            s = f"\t\ta  {a_s}\t'{a_s}\n\t\tb  {b_s}\t'{b_s}\n\t\tc  {c_s}\t'{c_s}\n\t\tal {al_s}\t'{al_s}"

    # to catch everything else
    else:
        s = f"\t\ta {a_s}\t'{a_s}\n\t\tb {b_s}\t'{b_s}\n\t\tc {c_s}\t'{c_s}\n\t\tal {al_s}\t'{al_s}\n\t\tbe {be_s}\t'{be_s}\n\t\tga {ga_s}\t'{ga_s}"

    return s


def get_unitcell2(cif, data):
    """
    Returns a string giving the unit cell parameters of the structure.
    see also get_unitcell()

    Args:
        cif: a PyCifRW dictionary
        data: the key value for a data block in the CIF
    Returns:
        A string containing the unit cell parameters in STR format.
    Raises:
        KeyError: if any of the unit cell parameters are not present in the datablock.
    """
    a = strip_brackets(get_dict_entry_copy_throw_error(cif[data], "_cell_length_a"))
    b = strip_brackets(get_dict_entry_copy_throw_error(cif[data], "_cell_length_b"))
    c = strip_brackets(get_dict_entry_copy_throw_error(cif[data], "_cell_length_c"))
    al = strip_brackets(get_dict_entry_copy_throw_error(cif[data], "_cell_angle_alpha"))
    be = strip_brackets(get_dict_entry_copy_throw_error(cif[data], "_cell_angle_beta"))
    ga = strip_brackets(get_dict_entry_copy_throw_error(cif[data], "_cell_angle_gamma"))

    return f"\t\ta {a} b {b} c {c}\n\t\tal {al} be {be} ga {ga}"


def get_atom_sites_string(cif, data):
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

    r = [make_atom_site_string(label, xf, yf, zf, atom, oc, b) + "\n" for label, xf, yf, zf, atom, oc, b in zip(labels, x, y, z, atoms, occ, b_iso)]

    return "".join(r)


def make_atom_site_string(label, x, y, z, atom, occ, beq):
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


def convert_site_label_to_atom(s):
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

    r = s[0:3]  # take the first three characters, maybe you mean water?
    if r in WATER:
        print(f"Site label '{s}' probably means 'water'. Please check that this atom really is oxygen.")
        return "O"

    r = s[0:2]  # take the first two characters
    if r in ELEMENTS:
        return r

    r = s[0:1]  # take the first character
    if r in ELEMENTS:
        if r == "W":
            print(f"W detected for site '{s}'. Do you mean tungsten or oxygen from a water molecule? Please check.")
        return r

    print(f"Can't decide what atom the site label '{s}' should be. Please check it.")
    return s  # if everything fails, just return what the label was


def fix_atom_type(a, orig_site_label=None):
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

    Returns:
        A string containing the atom type with the charge in the correct order, if one is given.
    """
    regex = re.search("([A-Za-z]{1,2})(\d{0,2})([+-]{0,1})(\d{0,2})", a)

    symbol = regex.group(1)
    charge = regex.group(2)
    sign = regex.group(3)
    digit = regex.group(4)  # if there is a trailing digit, then that is probably the charge

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


def convert_atom_type_to_topas(a, orig_site_label=None):
    """
    Given an atom type in the correct TOPAS site format, the CIF dictionary can allow charges
    for which there are no scattering factors..

    The charge is checked against the list of scattering factors supported by TOPAS, and if
    a mismatch is detected, a charged atom is replaced by the neutral one.
    eg: H+1 --> H

    Args:
        a: a string to denoting the atom type, with or without a charge

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
    regex = re.search("([A-Za-z]{1,2})([+-]{0,1})(\d{0,2})", a)  # Cu+2
    symbol = regex.group(1)

    if orig_site_label is not None:
        insert_text = f" in site {orig_site_label}"
    else:
        insert_text = ""

    print(f"{a} is not a legal TOPAS scattering factor. Atom{insert_text} replaced with {symbol}.")

    return symbol


def get_b_iso(cif, data):
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
    b_iso = [change_NA_value(strip_brackets(i)) for i in cif[data]["_atom_site_B_iso_or_equiv"]]

    num = count_nones(b_iso)
    if num > 0:
        print(f"{num} missing Biso values.")

    return b_iso


def convert_u_to_b(s):
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
    s = 8 * math.pi ** 2 * s
    s = str(round(float(s), 3))  # round B value to 3 d.p.
    return s


def get_u_iso(cif, data):
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
    b_iso = [convert_u_to_b(change_NA_value(strip_brackets(i))) for i in cif[data]["_atom_site_U_iso_or_equiv"]]

    num = count_nones(b_iso)
    if num > 0:
        print(f"{num} missing Uiso values.")

    return b_iso


def get_beq(cif, data):
    """
    Looks the Biso, Uiso, Baniso, Uaniso (in that order) to get beq values. If they are missing,
    then a default value of "1." is used.

    Args:
        cif: a PyCifRW dictionary
        data: the key value for a data block in the CIF
    Returns:
        A list containing the isotropic atomic displacement parameters, B, of all atomic sites.

    """
    # these will act as keys later on
    labels = get_dict_entry_copy(cif[data], "_atom_site_label")  # I want everything to be in the same order as the site_label

    # get all potential ADPs from the str as dictionaries. The key is the atom label, the value is the ADP value
    b_iso = make_b_dict(cif, data, "b_iso")
    u_iso = make_b_dict(cif, data, "u_iso")
    b_aniso = make_b_dict(cif, data, "b_aniso")
    u_aniso = make_b_dict(cif, data, "u_aniso")

    # check all the dictionaries in my order of preference to get the ADP value.
    r = []
    for key in labels:
        try:
            b = b_iso[key]
        except KeyError:
            try:
                b = u_iso[key]
            except KeyError:
                try:
                    b = b_aniso[key]
                    print(f"beq value for site {key} calculated from anisotropic B values.")
                except KeyError:
                    try:
                        b = u_aniso[key]
                        print(f"beq value for site {key} calculated from anisotropic U values.")
                    except KeyError:
                        b = None

        if b is None or float(b) == 0.0:  # missing or zero value
            print(f"Warning! beq value missing or zero for site {key}! Default value of 1 entered")
            b = "1."
        elif b.startswith("-"):
            print(f"Warning! Negative atomic displacement parameter detected for site {key}! Please check.")
        r.append(b)  # the ADP values are now in the same order as the site label

    return r


def make_b_dict(cif, data, b_type):
    """
    A helper function for get_beq() to get the site labels and associated ADP values, and
    return a dictionary with labels as keys and values as ADP values. If the value is None,
    then that key is removed

    Parameters
    ----------
    cif : PyCifRW dictionary

    data : the data block you're looking at

    b_type: "b_iso", "u_iso", "b_aniso", or "u_aniso"

    Returns
    -------
    Dictionary containing site labels as keys and ADPs as values.
    If a value is None, it is removed from the dictionary.
    It is possible to return an empty dictionary.
    If the b_type selected is not in the cif, an empty dictionary is returned

    Raises
    ------
    ValueError if b_type is not one of "b_iso", "u_iso", "b_aniso", or "u_aniso".

    """
    iso = True
    if b_type == "b_iso":
        f = get_b_iso
    elif b_type == "u_iso":
        f = get_u_iso
    elif b_type == "b_aniso":
        f = get_b_aniso
        iso = False
    elif b_type == "u_aniso":
        f = get_u_aniso
        iso = False
    else:
        raise ValueError(f'Invalid choice. Got "{b_type}", expecting "b_iso", "u_iso", "b_aniso", or "u_aniso".')

    try:
        b_list = f(cif, data)
        if iso:
            b_labels = get_dict_entry_copy(cif[data], "_atom_site_label")
        else:
            b_labels = get_dict_entry_copy(cif[data], "_atom_site_aniso_label")
    except KeyError:
        return {}
    d = dict(zip(b_labels, b_list))

    # filter out any None values
    d = {k: v for k, v in d.items() if v is not None}
    return d


def get_b_aniso(cif, data):
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
    b_equiv = [str(round(float(i), 3)) for i in b_equiv]  # round to 3 d.p.

    return b_equiv


def get_u_aniso(cif, data):
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

    # get the average of the three values
    b_equiv = [convert_u_to_b(str((u11 + u22 + u33) / 3.0)) for u11, u22, u33 in zip(U11, U22, U33)]

    return b_equiv
