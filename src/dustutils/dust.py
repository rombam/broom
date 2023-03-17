import numpy as np

from dataclasses import dataclass, field
from typing import List, Union, Literal
from pathlib import Path
from copy import deepcopy

from dustutils.utils import Printable
from dustutils.mesh import CGNS, Parametric, Pointwise
from dustutils.reference import Reference
from dustutils.solver import Settings

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
                                f'ref_tag={self.ref.reference_tag}'])

        return {'pre': pre_string, 'geom': self.geom.to_fort(), 'ref': self.ref.to_fort()}

@dataclass
class Analysis(Printable):
    """DUST postprocessing analysis."""
    type: Literal['viz', 'integral_loads', 'hinge_loads', ]


@dataclass
class Post(Printable):
    """DUST postprocessing class."""
    data_basename: Union[str, Path]
    basename: Union[str, Path]
    analyses: List[Analysis]

@dataclass
class Case(Printable):
    """DUST simulation object.

    Attributes
    ----------
    geoms : Geom, List[Geom]
        Geometry or list of geometries to be included in the simulation.
        DustGeom objects include information about the mesh and reference system of each
        geometry.
    settings : Settings
        Settings object containing the simulation settings and solver options.
        TODO: implement settings templates to make it easier for users.
    post : Post, optional
        [WIP] Postprocessing settings object.

    """
    geoms: Union[Geom, List[Geom]]
    settings: Settings
    post: Post = None

    def to_fort(self):
        """Modified to_fort superclass method.

        Returns
        -------
        dict : dict
            Dictionary containing the preprocessing, references, geometries, solver sett-
            ings and postprocessing settings Fortran string representations.

        """
        if not isinstance(self.geoms, list):
            self.geoms = [self.geoms]
        pre = '\n\n'.join([geomobj.to_fort()['pre'] for geomobj in self.geoms])
        geom_dict = {geomobj.comp_name: geomobj.geom.to_fort() for geomobj in self.geoms}
        ref_str = '\n\n'.join([refobj.ref.to_fort() for refobj in self.geoms])
        return {'pre': pre, 'geom': geom_dict,
                'ref': ref_str, 'settings': self.settings.to_fort(),
                'post': self.post.to_fort()}
