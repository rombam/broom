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

    """
    comp_name: str
    geom: Union[CGNS, Parametric, Pointwise]
    ref: Reference = None

    def to_fort(self):
        """Modified to_fort superclass method.

        Returns
        -------
        dict : dict
            Dictionary containing the geometry and reference frame Fortran string repre-
            sentations.

        """
        return {'geom': self.geom.to_fort(), 'ref': self.ref.to_fort()}

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
    post : Postpro, optional
        [WIP] Postprocessing settings object. Not yet implemented.

    """
    geoms: Union[Geom, List[Geom]]
    settings: Settings
    post: None

    def to_fort(self):
        """Modified to_fort superclass method.

        Returns
        -------
        dict : dict
            Dictionary containing the preprocessing, references, geometries, solver sett-
            ings and postprocessing settings Fortran string representations.

        """
        geoms = {}
        return {'geom': self.geom.to_fort(), 'ref': self.ref.to_fort()}
