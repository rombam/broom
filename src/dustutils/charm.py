# Tools to work with CHARM files
import numpy as np

from warnings import warn
from .polartable import PolarTable


def parse_rw(filename):
    """Parse a rotor wake (rw) CHARM file and return useful DUST parameters of the case.

    Parameters
    ----------
    filename : str, Path
        Path to the rotor wake (rw) file.

    Returns
    -------
    rw_dict : dict
        Dictionary containing useful rotor parameters.

    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Assign number of expected values and expected type
    n_val = {'NBLADE': [1, int],
             'OMEGA': [1, float],
             'IROTAT': [1, int],
             'XROTOR': [3, float],
             'X,Y,Z tilt': [3, float],
             'ITILT': [1, int],
             'ICOLL': [1, int],
             'COLL': [1, float],
             'CT': [1, float]}

    rw_dict = {}

    for idx, line in enumerate(lines):
        found = []
        for keyword in n_val.keys():
            if keyword in lines[idx-1]:
                found.append(keyword)
        if found:
            values = line.strip().replace('\n', '').split(' ')
            values = [float(val) for val in values if val]

            count = 0
            for key in found:
                n_values = n_val[key][0]
                type_values = n_val[key][1]
                if n_values > 1:
                    rw_dict[key] = [type_values(val)
                                    for val in values[count:count+n_values]]
                else:
                    rw_dict[key] = type_values(values[count])
                count += n_values
                n_val.pop(key)

    # Postprocess tilt angle name
    rw_dict['TILT'] = rw_dict.pop('X,Y,Z tilt')

    return rw_dict


def parse_main(filename):
    """Parse the main CHARM file and return useful DUST parameters.

    Parameters
    ----------
    filename : str, Path
        Path to the main CHARM file.

    Returns
    -------
    main_dict : dict
        Dictionary containing useful CHARM main input file parameters.

    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Assign number of expected values and expected type
    n_val = {'NROTOR': [1, int],
             'IMKS': [1, int],
             'SSPD': [1, float],
             'RHO': [1, float],
             'SFRAME': [1, int],
             'UCG': [1, float],
             'VCG': [1, float],
             'WCG': [1, float],
             'PCG': [1, float],
             'QCG': [1, float],
             'RCG': [1, float],
             'NPSI': [1, int],
             'NREV': [1, int]}

    main_dict = {}

    # Parse main parameters
    for idx, line in enumerate(lines):
        found = []
        for keyword in n_val.keys():
            if keyword in lines[idx-1]:
                found.append(keyword)
        if found:
            values = line.strip().replace('\n', '').split(' ')
            values = [float(val) for val in values if val]

            count = 0
            for key in found:
                n_values = n_val[key][0]
                type_values = n_val[key][1]
                if n_values > 1:
                    main_dict[key] = [type_values(val)
                                      for val in values[count:count+n_values]]
                else:
                    main_dict[key] = type_values(values[count])
                count += n_values
                n_val.pop(key)

    # Parse filename blocks
    main_dict['ROTOR_FILES'] = {}
    count = 0
    for idx, line in enumerate(lines):
        if "PATHNAME" in line:
            main_dict['PATHNAME'] = lines[idx+1].strip()
        elif "INPUT FILENAMES" in line:
            name = line.strip().split('for')[-1].split('#')[0].strip()
            main_dict['ROTOR_FILES'][name] = {'rw': lines[idx+1].strip(),
                                               'bg': lines[idx+2].strip(),
                                               'bd': lines[idx+3].strip(),
                                               'af': lines[idx+4].strip(),
                                               'cs': lines[idx+5].strip()}
            count += 1

    # Warn if numebr of filename blocks is not equal to NROTOR
    if count != main_dict["NROTOR"]:
        warn(f'Found {count} rotor filename entries and ' /
             +f'NROTOR={main_dict["NROTOR"]} in {filename}')

    return main_dict


def parse_bg(filename):
    """Parse a blade geometry CHARM file and return useful DUST parameters.

    Parameters
    ----------
    filename : str, Path
        Path to the blade geometry (bg) file.

    Returns
    -------
    bg_dict : dict
        Dictionary containing useful CHARM blade geometry file parameters.

    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    bg_dict = {}

    # Pre-parse number of expected segments
    for idx, line in enumerate(lines):
        if 'NSEG' in line:
            bg_dict['NSEG'] = nseg = int(lines[idx+1].strip())
            break

    # Assign number of expected values and expected type for single-value keys
    n_val = {'CUTOUT': [1, float],
             'TWRD': [1, float],
             'NCHORD': [1, int],
             'NSPAN': [1, int],
             'ICOS': [1, int]}

    # Parse single parameters
    for idx, line in enumerate(lines):
        found = []
        for keyword in n_val.keys():
            if keyword in lines[idx-1]:
                found.append(keyword)
        if found:
            values = line.strip().replace('\n', '').split(' ')
            values = [float(val) for val in values if val]

            count = 0
            for key in found:
                n_values = n_val[key][0]
                type_values = n_val[key][1]
                if n_values > 1:
                    bg_dict[key] = [type_values(val)
                                    for val in values[count:count+n_values]]
                else:
                    bg_dict[key] = type_values(values[count])
                count += n_values
                n_val.pop(key)

    distr_val = {'SL': [nseg, float],
                 'CHORD': [nseg+1, float],
                 'SWEEPD': [nseg, float],
                 'TWSTGD': [nseg, float],
                 'ANHD': [nseg, float],
                 'THIKND': [nseg+1, float]}

    # Parse multiple parameters
    for idx, line in enumerate(lines):
        found = []
        for keyword in distr_val.keys():
            if keyword in lines[idx-1]:
                found.append(keyword)
        if found:
            values = line.strip().replace('\n', '').split(' ')
            values = [val for val in values if val]
            if len(values) == 1:
                values = values[0].split('*')
                values = int(values[0])*[float(values[1])]
            else:
                values = [float(val) for val in values]

            count = 0
            for key in found:
                n_values = distr_val[key][0]
                type_values = distr_val[key][1]
                if n_values > 1:
                    bg_dict[key] = [type_values(val)
                                    for val in values[count:count+n_values]]
                else:
                    bg_dict[key] = type_values(values[count])
                count += n_values
                distr_val.pop(key)

    return bg_dict


def parse_af(filename):
    """Parse a CHARM airfoil table file and return useful DUST parameters.

    Note
    ----
    Reynolds number is not explicitly included in the CHARM C81 tables. Thus, it is set to
    0.0 for all airfoils. This should be adjusted afterwards.

    Parameters
    ----------
    filename : str, Path
        Path to the airfoil table (af) file.

    Returns
    -------
    af_dict : dict
        Dictionary containing useful CHARM airfoil table parameters.

    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    af_dict = {}

    # Postprocess every line
    lines = [line.strip().replace('\n', '').split(' ') for line in lines]
    lines = [[item for item in line if item] for line in lines]

    # Parse header and remove the initial lines
    af_dict['N_AF'] = int(lines[0][0])
    af_dict['SEC'] = {}
    r_R = [float(item) for item in lines[1]]
    thickd = [float(item) for item in lines[2]]
    lines = lines[3:]

    # Parse body
    for af in range(af_dict['N_AF']):

        # Comment header
        com_1 = ' '.join(lines[0])
        com_2 = ' '.join(lines[1])
        desc = f'{com_1} -- {com_2}'
        name = lines[2][0]
        lpoints, dpoints, mpoints = lines[2][1:4]
        _, al = (int(lpoints[:-2]), int(lpoints[-2:]))
        _, ad = (int(dpoints[:-2]), int(dpoints[-2:]))
        _, am = (int(mpoints[:-2]), int(mpoints[-2:]))
        lines = lines[3:]

        # NOTE: as a stopgap solution, Mach numbers will be assumed to be the same for
        # Cl, Cd and Cm. Later on, PolarTable should be generalized to include all values.
        # NOTE: similarly, alpha is assumed to be the same for Cl, Cd and Cm. This is most
        # likely going to always be the case
        # NOTE: since CHARM files do not include structured Reynolds number information,
        # it will be set to 0.0

        # CL
        mach_cl = [float(item) for item in lines.pop(0)]
        cl_table = np.array(lines[0:al]).astype(float)
        alpha_cl = cl_table[:, 0]
        val_cl = cl_table[:, 1:]
        lines = lines[al:]

        # CD
        mach_cd = [float(item) for item in lines.pop(0)]
        cd_table = np.array(lines[0:ad]).astype(float)
        alpha_cd = cd_table[:, 0]
        val_cd = cd_table[:, 1:]
        lines = lines[ad:]

        # CM
        mach_cm = [float(item) for item in lines.pop(0)]
        cm_table = np.array(lines[0:am]).astype(float)
        alpha_cm = cm_table[:, 0]
        val_cm = cm_table[:, 1:]
        lines = lines[am:]

        # Build the polar object
        polar = PolarTable(name=name, desc=desc)
        for idx, m in enumerate(mach_cl):
            polar.append_point(re=0.0, m=m, alpha=alpha_cl, cl=val_cl[:, idx])
        for idx, m in enumerate(mach_cd):
            polar.append_point(re=0.0, m=m, alpha=alpha_cd, cd=val_cd[:, idx])
        for idx, m in enumerate(mach_cm):
            polar.append_point(re=0.0, m=m, alpha=alpha_cm, cm=val_cm[:, idx])
        af_dict['SEC'][af] = {'name': name,
                              'rR': r_R[af],
                              'thick': thickd[af],
                              'polar': polar}

    return af_dict
