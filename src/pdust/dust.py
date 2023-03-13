import numpy as np

from dataclasses import dataclass, field
from typing import List, Union, Literal
from pathlib import Path
from copy import deepcopy

from pdust.utils import Printable
from pdust.mesh import CGNS, Parametric, Pointwise
from pdust.reference import Reference

@dataclass
class DustGeom(Printable):
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

class DustCase(Printable):
    """DUST simulation object.

    Attributes
    ----------
    geoms : List[DustGeom]
        List of geometries to be included in the simulation. DustGeom objects include
        information about the mesh and reference system of each geometry.
    

    """
    pass

