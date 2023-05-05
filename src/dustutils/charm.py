# Tools to work with CHARM files
import numpy as np

from pathlib import Path
from warnings import warn
from copy import deepcopy
from .polartable import PolarTable
from .solver import Settings, TimeOpts, WakeOpts, FMMOpts, FlowOpts, ModelOpts
from .mesh import Point, Line, Pointwise, MeshMirror, MeshSymmetry
from .reference import Reference, RotorMulti, RotorDOF

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
    # TODO: dependency on SFRAME. If SFRAME=1, then UCG, VCG, WCG, if SFRAME=2, then
    # WIND
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
                 'THCKND': [nseg+1, float]}

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


def mesh_charm(bg, af, rw):
    """Create a Pointwise mesh object from CHARM case files.

    Parameters
    ----------
    bg : dict
        Blade geometry dictionary from CHARM.
    af : dict
        Airfoil dictionary from CHARM.
    rw : dict
        Rotor wake dictionary from CHARM.

    """
    R = sum(bg['SL'])+bg['CUTOUT']

    if rw['OMEGA'] > 0.0:
        if rw['IROTAT'] == -1:
            mesh_mirror = MeshMirror(mesh_mirror=True)
            mult = -1
        else:
            mesh_mirror = MeshMirror()
            mult = 1
        mesh_symmetry = MeshSymmetry()
    else:
        mesh_mirror = MeshMirror()
        mesh_symmetry = MeshSymmetry(mesh_symmetry=True)
        mult = 1

    prop_data = {}
    prop_data['af'] = {}
    prop_data['af']['rR'] = np.array([af['SEC'][i]['rR'] for i in af['SEC'].keys()])
    prop_data['af']['table'] = [af['SEC'][i]['polar'] for i in af['SEC'].keys()]
    prop_data['chord'] = {}
    prop_data['chord']['span'] = np.array(bg['SL'])
    prop_data['chord']['sweep'] = np.array(bg['SWEEPD'])
    prop_data['chord']['anhed'] = np.array(bg['ANHD'])
    tmptwst = [bg['TWRD']]+bg['TWSTGD']
    prop_data['chord']['twist'] = np.array([sum(tmptwst[0:i+1]) for i in range(len(tmptwst))])
    prop_data['chord']['chord'] = np.array(bg['CHORD'])
    prop_data['chord']['rR'] = [round(bg['CUTOUT']/R+sum(prop_data['chord']['span'][0:i]/R), 3) for i in range(len(prop_data['chord']['span'])+1)]

    # Sectional properties
    prop_data['interp'] = {'rR': sorted(set(prop_data['af']['rR']).union(set(prop_data['chord']['rR'])))}
    prop_data['interp']['r'] = [r*R for r in prop_data['interp']['rR']]
    prop_data['interp']['chord'] = np.interp(prop_data['interp']['rR'], prop_data['chord']['rR'], prop_data['chord']['chord'])
    prop_data['interp']['twist'] = np.interp(prop_data['interp']['rR'], prop_data['chord']['rR'], prop_data['chord']['twist'])
    prop_data['interp']['af'] = [prop_data['af']['table'][list(prop_data['af']['rR']).index(rr)] if rr in prop_data['af']['rR'] else 'interp' for rr in prop_data['interp']['rR']]
    prop_data['interp']['af_name'] = [f"airfoil/{prop_data['interp']['af'][i].name}.c81" if prop_data['interp']['af'][i] != 'interp' else 'interp' for i in range(len(prop_data['interp']['af']))]

    # Strip properties
    prop_data['interp']['sweep'] = np.interp(prop_data['interp']['rR'][1:], prop_data['chord']['rR'][1:], prop_data['chord']['sweep'])
    prop_data['interp']['sweep'] = np.insert(prop_data['interp']['sweep'], 0, 0.0)
    prop_data['interp']['anhed'] = np.interp(prop_data['interp']['rR'][1:], prop_data['chord']['rR'][1:], prop_data['chord']['anhed'])
    prop_data['interp']['anhed'] = np.insert(prop_data['interp']['anhed'], 0, 0.0)
    prop_data['interp']['span'] = np.array([0.0] + [prop_data['interp']['r'][i]-prop_data['interp']['r'][i-1]
                                                    for i in range(1, len(prop_data['interp']['r']))])

    nsec = len(prop_data['interp']['af_name'])
    airfoils = prop_data['interp']['af_name']

    # Parse section positions
    prop_data['interp']['dx'] = np.array([np.tan(np.deg2rad(prop_data['interp']['sweep'][i]))\
                                          * mult * prop_data['interp']['span'][i]
                                          for i in range(nsec)])
    prop_data['interp']['dz'] = np.array([-np.tan(np.deg2rad(prop_data['interp']['anhed'][i]))\
                                          * mult * prop_data['interp']['span'][i]
                                          for i in range(nsec)])

    nsec = len(prop_data['interp']['af_name'])
    airfoils = prop_data['interp']['af_name']
    positions = [[np.sum(prop_data['interp']['dx'][0:i+1]),
                  mult*prop_data['interp']['r'][i],
                  np.sum(prop_data['interp']['dz'][0:i+1])] for i in range(nsec)]

    points = [Point(i, positions[i-1], prop_data['interp']['chord'][i-1], prop_data['interp']['twist'][i-1], airfoil_table=airfoils[i-1]) for i in range(1, nsec+1)]
    lines = [Line('Straight', points[i], points[i+1], 1) for i in range(0, nsec-1)]
    pointgeom = Pointwise(el_type='l',
                          nelem_chord=1,
                          points=points,
                          lines=lines,
                          type_chord='cosineLE',
                          mesh_symmetry=mesh_symmetry,
                          mesh_mirror=mesh_mirror)

    return pointgeom


def ref_charm(rw, tag='propsys', parent='ac'):
    """Create a Reference object from CHARM case files.

    Parameters
    ----------
    rw : dict
        Rotor wake parameter dictionary.
    tag : str, optional
        Reference tag. Default: 'propsys'.
    parent : str, optional
        Parent tag. Default: 'ac'.

    """
    rpm = rw['OMEGA']*60/(2*np.pi)
    origin = rw['XROTOR']
    yaw = rw['TILT'][2]
    pitch = rw['TILT'][1]
    roll = rw['TILT'][0]
    nblades = rw['NBLADE']
    irotat = rw['IROTAT']
    coll = rw['COLL']
    if irotat == 1:
        rot_axis = np.array([0.0, 0.0, -1.0])
    elif irotat == -1:
        rot_axis = np.array([0.0, 0.0, 1.0])
    if nblades > 1:
        propref = Reference(reference_tag=tag, parent_tag=parent,
                            origin=np.array([0.0, 0.0, 0.0]),
                            orientation=np.array([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]),
                            multiple=True,
                            multiplicity=RotorMulti(n_blades=nblades, rot_axis=rot_axis,
                                                    rot_rate=rpm/60*2*np.pi,
                                                    dofs=[RotorDOF('Flap'),
                                                          RotorDOF('Lag'),
                                                          RotorDOF('Pitch',
                                                                   collective=coll)]))
    else:
        # Modify the yaw (azimuth) and pitch angles in the CHARM reference frame
        yaw = 0.0
        pitch += 180.0
        propref = Reference(reference_tag=tag, parent_tag=parent,
                            origin=np.array([0.0, 0.0, 0.0]),
                            orientation=np.array([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]))

    propref.transform(origin, yaw, pitch, roll)

    return propref


def geom_charm(bg, af, rw, name):
    """Create a Geom object from CHARM case files.

    Parameters
    ----------
    bg : dict
        Dictionary of blade geometry case file parameters.
    af : dict
        Dictionary of airfoil case file parameters.
    rw : dict
        Dictionary of rotor wake case file parameters.
    name : str
        Name of the Geom object.

    """
    geom = mesh_charm(bg, af, rw)
    ref = ref_charm(rw, tag=name)

    return deepcopy(geom), deepcopy(ref)

def opts_charm(main, rw, parts=200000, part_box_min=np.array([-15.0, -15.0, -15.0]),
               part_box_max=np.array([15.0, 15.0, 15.0])):
    """Create a Settings object from CHARM case files.

    Parameters
    ----------
    main : dict
        Dictionary of main case file parameters.
    rw : dict
        Dictionary of rotor wake case file parameters.
    parts : int, optional
        Number of wake particles to use. Default: 200000.
    part_box_min : np.ndarray, optional
        Minimum corner of the wake particle box. Default: np.array([-15.0, -15.0, -15.0]).
    part_box_max : np.ndarray, optional
        Maximum corner of the wake particle box. Default: np.array([15.0, 15.0, 15.0]).

    """
    rpm = rw['OMEGA']*60/(2*np.pi)
    n_rev = main['NREV']
    step = 360/main['NPSI']

    if rpm == 0:
        topts = TimeOpts(tstart=0.0, tend=1.0, timesteps=10)
    else:
        topts = TimeOpts.from_rpm(rpm=rpm, step=step, nrev=n_rev)

    flopts = FlowOpts(u_inf=np.array([main['UCG'], -main['VCG'], main['WCG']]))
    fmmopts = FMMOpts(fmm=True, box_length=10,
                      n_box=[3, 3, 3], octree_origin=part_box_min,
                      n_octree_levels=6,
                      min_octree_part=5,
                      multipole_degree=2)
    wakeopts = WakeOpts(n_wake_particles=parts,
                        particles_box_min=part_box_min,
                        particles_box_max=part_box_max)
    modelopts = ModelOpts(penetration_avoidance=True)
    settings = Settings(time=topts, flow=flopts, fmm=fmmopts, wake=wakeopts,
                        model=modelopts)

    return settings


def write_airfoils(main, target):
    """Write the case airfoils to a target folder.

    Parameters
    ----------
    main : str, Path
        Path to the main file.
    target : str, Path
        Path to the target folder.

    """
    mainpath = Path(main)
    root = mainpath.parent.absolute()
    main = parse_main(mainpath)
    afolder = Path(target)
    afolder.mkdir(exist_ok=True, parents=True)
    for rdata in main['ROTOR_FILES'].values():
        af_path = root / Path(rdata['af'])
        af = parse_af(af_path)
        for foil in af['SEC'].values():
            savename = afolder / Path(foil['name'] + '.c81')
            if not Path(savename).exists():
                foil['polar'].to_dust(savename, re=[0.0])
            else:
                pass
