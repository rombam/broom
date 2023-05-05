import os
import numpy as np

from dataclasses import dataclass
from typing import List, Union
from pathlib import Path
from copy import deepcopy

from .utils import Printable
from .mesh import CGNS, Parametric, Pointwise
from .reference import Reference
from .solver import Settings
from .post import Post, basic_post
from .charm import parse_main, parse_af, parse_bg, parse_rw, geom_charm, opts_charm

@dataclass
class Geom(Printable):
    """DUST geometry class.

    Attributes
    ----------
    comp_name : str
        Component name.
    geom : pdust.mesh.CGNS, pdust.mesh.Parametric, pdust.mesh.Pointwise
        Geometry mesh object.
    ref : pdust.reference.Reference, optional
        Reference frame object. If None, will use the default global reference frame.
    save_name : str, optional
        Name of the geometry file to be saved. If None, will use the component name.

    """
    comp_name: str
    geom: Union[CGNS, Parametric, Pointwise]
    ref: Reference = None
    save_name: str = None

    def to_fort(self):
        """Modified to_fort superclass method.

        Returns
        -------
        dict : dict
            Dictionary containing the preprocessing, geometry and reference frame Fortran
            string representations.

        """
        if not self.save_name:
            save_name = self.comp_name
        else:
            save_name = self.save_name

        pre_string = '\n'.join([f'comp_name = {self.comp_name}',
                                f'geo_file = {save_name}.in',
                                f'ref_tag = {self.ref.reference_tag}'])

        return {'pre': pre_string, 'geom': deepcopy(self.geom).to_fort(), 'ref': deepcopy(self.ref).to_fort()}


@dataclass
class Case(Printable):
    """DUST simulation object.

    Attributes
    ----------
    name : str
        Name of the simulation case.
    geoms : Geom, List[Geom]
        Geometry or list of geometries to be included in the simulation.
        DustGeom objects include information about the mesh and reference system of each
        geometry.
    settings : Settings
        Settings object containing the simulation settings and solver options.
        TODO: implement settings templates to make it easier for users.
    references : List[Reference], optional
        List of extra reference frames to be used in the simulation/postprocessing, apart
        from the geometry-bound ones. This is useful for adding extra reference frames
        to postprocess loads in.
    post : Post, optional
        [WIP] Postprocessing settings object.

    """
    name: str
    geoms: Union[Geom, List[Geom]]
    settings: Settings
    references: Union[Reference, list[Reference]] = None
    post: Post = None

    def __post_init__(self):
        """Post-initialization method."""
        if not isinstance(self.geoms, list):
            self.geoms = [self.geoms]
        if not isinstance(self.references, list):
            self.references = [self.references]

    @classmethod
    def from_charm(cls, main, name):
        """Create a DUST case from a CHARM case.

        Notes
        -----
        By default, returns a case with lifting-line elements, n_part=200000, a bounding box between
        (-15, -15, -15) and (15, 15, 15), and a timestep set by using the RPM and NPSI from the CHARM
        case. The user can modify these settings by changing the case properties after creation.

        Parameters
        ----------
        main : str, Path
            Path to the CHARM main file.
        name : str
            Name of the case.

        """
        mainpath = Path(main)
        root = mainpath.parent.absolute()
        main = parse_main(mainpath)
        geoms = []
        omegas = {}
        for rname, rdata in main['ROTOR_FILES'].items():
            bg_path = root / Path(rdata['bg'])
            af_path = root / Path(rdata['af'])
            rw_path = root / Path(rdata['rw'])
            bg = parse_bg(bg_path)
            af = parse_af(af_path)
            rw = parse_rw(rw_path)
            omegas[rname] = rw['OMEGA']
            geom, ref = geom_charm(bg, af, rw, rname)
            geoms.append(Geom(rname, geom, ref))

        rw = parse_rw(root / Path(main['ROTOR_FILES'][max(omegas, key=omegas.get)]['rw']))
        flowsets = opts_charm(main, rw)
        last_step = int((flowsets.time.tend-flowsets.time.tstart)/flowsets.time.dt)
        post = basic_post(res=(1, last_step), name=name)
        ac_ref = Reference(reference_tag='ac', parent_tag='0', origin=np.array([0.0, 0.0, 0.0]),
                           orientation=np.array([-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0]))

        return cls(name, geoms, flowsets, references=[ac_ref], post=post)

    def to_fort(self, geom_name=None, res_name=None, post_name=None):
        """Modified to_fort superclass method.

        Parameters
        ----------
        geom_name : str, optional
            Name of the geometry file to be saved. If None, will use 'geo_input.h5'.
            Should be relative to the case folder.
        res_name : str, optional
            Name of the results file to be saved. Should include the path and a preffix
            for the results. If None, will use 'Output/<self.name>'.
            Should be relative to the case folder.
        post_name : str, optional
            Name of the post-processed files to be saved. Should include the path and a
            preffix for the results. If None, will use 'Postpro/<self.name>'.
            Should be relative to the case folder.

        Returns
        -------
        dict : dict
            Dictionary containing the preprocessing, references, geometries, solver sett-
            ings and postprocessing settings Fortran string representations.

        """
        geom_fort = [geomobj.to_fort() for geomobj in deepcopy(self.geoms)]
        geom_names = [geomobj.comp_name for geomobj in deepcopy(self.geoms)]

        pre_str = '\n\n'.join([geomd['pre'] for geomd in geom_fort]
                              + [f'file_name = {geom_name}'])
        geom_dict = {geom_names[i]: geom_fort[i]['geom'] for i in range(len(geom_names))}

        if self.references is not None:
            ref_str = '\n\n'.join([refobj.to_fort() for refobj in deepcopy(self.references)])
        else:
            ref_str = ''

        ref_str += '\n\n' + '\n\n'.join([geomd['ref'] for geomd in geom_fort])

        set_str = f'basename = {res_name}\n' + deepcopy(self.settings).to_fort()\
            + f'\n\ngeometry_file = {geom_name}'\
            + '\nreference_file = references.in'

        if self.post:
            post_str = deepcopy(self.post).to_fort()
        else:
            post_str = ''

        return {'pre': pre_str, 'geom': geom_dict,
                'ref': ref_str, 'settings': set_str,
                'post': post_str}

    def write_case(self, folder, geom_name=None, res_name=None, post_name=None):
        """Write the DUST case files to the specified folder.

        Parameters
        ----------
        folder : str, Path
            Path to the folder where the case files will be written. The folder will be
            created if it does not exist. Inside the folder, another folder named
            self.name will be created and the case files will be written there.
        geom_name : str, optional
            Name of the geometry file to be saved. If None, will use 'geo_input.h5'.
            Should be relative to the case folder.
        res_name : str, optional
            Name of the results file to be saved. Should include the path and a preffix
            for the results. If None, will use 'Output/<self.name>'.
            Should be relative to the case folder.
        post_name : str, optional
            Name of the post-processed files to be saved. Should include the path and a
            preffix for the results. If None, will use 'Postpro/<self.name>'.
            Should be relative to the case folder.

        """
        # Create base directory
        folder = Path(folder)
        os.makedirs(folder, exist_ok=True)

        if not geom_name:
            geom_name = 'geo_input.h5'
        if not res_name:
            res_name = f'Output/{self.name}'
        if not post_name:
            post_name = f'Postpro/{self.name}'

        # Create subdirectories
        res_folder = Path(res_name).parent if Path(res_name).parent != '.' else ''
        post_folder = Path(post_name).parent if Path(post_name).parent != '.' else ''
        res_path = Path.joinpath(folder, res_folder)
        post_path = Path.joinpath(folder, post_folder)
        os.makedirs(res_path, exist_ok=True)
        os.makedirs(post_path, exist_ok=True)

        # Write case files
        str_dict = self.to_fort(geom_name=geom_name, res_name=res_name,
                                post_name=post_name)
        with open(Path.joinpath(folder, 'dust_pre.in'), 'w+') as f:
            f.writelines(str_dict['pre'])
        with open(Path.joinpath(folder, 'dust.in'), 'w+') as f:
            f.writelines(str_dict['settings'])
        with open(Path.joinpath(folder, 'dust_post.in'), 'w+') as f:
            f.writelines(str_dict['post'])
        with open(Path.joinpath(folder, 'references.in'), 'w+') as f:
            f.writelines(str_dict['ref'])
        for geom in str_dict['geom'].keys():
            with open(Path.joinpath(folder, f'{geom}.in'), 'w+') as f:
                f.writelines(str_dict['geom'][geom])
