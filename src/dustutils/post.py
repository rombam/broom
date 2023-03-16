import numpy as np

from dataclasses import dataclass, field
from typing import List, Union, Literal
from pathlib import Path
from copy import deepcopy

from dustutils.utils import Printable

@dataclass
class Analysis(Printable):
    """DUST postprocessing analysis."""
    type: Literal['viz', 'integral_loads', 'hinge_loads']

@dataclass
class Viz(Printable):
    """A visualization analysis."""
    name: str
    start_res: int
    end_res: int
    step_res: int = 1
    format: Literal['vtk', 'tecplot'] = 'vtk'
    wake: bool = True
    average: bool = False
    avg_res: int = None
    variable: Union[Literal['vorticity', 'vorticity_vector', 'velocity',
                            'surface_velocity','pressure', 'cp', 'turbulent_viscosity',
                            'vortex_rad'],
                    List[str]]
    component: Union[str, List[str]] = 'all'

    def __post_init__(self):
        """Post-constructor method to assign fixed values."""
        # Prepare the target attribute list in correct order
        attlist = deepcopy(list(self.__dict__.keys()))
        attlist.insert(0, 'type')

        # Reconstruct and reorder attributes
        self.type = 'viz'
        self.__dict__ = {attlist[i]: self.__dict__[attlist[i]]
                         for i in range(len(attlist))}

@dataclass
class IntLoads(Printable):
    """A visualization analysis."""
    name: str
    start_res: int
    end_res: int
    step_res: int = 1
    format: Literal['dat', 'tecplot'] = 'dat'
    average: bool = False
    component: Union[str, List[str]] = 'all'
    reference_tag: str = '0'

    def __post_init__(self):
        """Post-constructor method to assign fixed values."""
        # Prepare the target attribute list in correct order
        attlist = deepcopy(list(self.__dict__.keys()))
        attlist.insert(0, 'type')

        # Reconstruct and reorder attributes
        self.type = 'integral_loads'
        self.__dict__ = {attlist[i]: self.__dict__[attlist[i]]
                         for i in range(len(attlist))}
    pass

@dataclass
class HingeLoads(Printable):
    """A visualization analysis."""
    pass

@dataclass
class Probe(Printable):
    """A visualization analysis."""
    pass

@dataclass
class FlowField(Printable):
    """A visualization analysis."""
    pass

@dataclass
class Sectional(Printable):
    """A visualization analysis."""
    pass

@dataclass
class Chordwise(Printable):
    """A visualization analysis."""
    pass

@dataclass
class Post(Printable):
    """DUST postprocessing class.

    Attributes
    ----------
    data_basename : str, Path
        Base name (with path) of the data which must be analysed.
    basename : str, Path
        Base name (with path) of the postprocessing results.
    analyses : list
        List of analyses to be performed.

    """
    data_basename: Union[str, Path]
    basename: Union[str, Path]
    analyses: List[Union[Viz, IntLoads, HingeLoads, Probe, FlowField, Sectional,
                         Chordwise]]
