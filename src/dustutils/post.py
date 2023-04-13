import numpy as np

from dataclasses import dataclass
from typing import List, Union, Literal
from pathlib import Path
from copy import deepcopy

from dustutils.utils import Printable


@dataclass
class PrintAnalysis(Printable):
    """Printable analysis helper class.

    Note
    ----
    The sole purpose of this class is to modify the to_fort method of the Printable class
    for the analysis classes.

    """
    def _fort_strs(self, ignore=[], indent=4):
        return ['\nanalysis = {'] + super()._fort_strs(ignore=ignore, indent=indent)\
               + ['}']


@dataclass
class Viz(PrintAnalysis):
    """Visualization analysis.

    Attributes
    ----------
    name : str
        Name of the analyis, will be appended as a suffix to the basename.
    start_res : int
        First result to be visualized.
    end_res : int
        Last result to be visualized.
    variable : str, list
        Variable(s) to be visualized.
    step_res : int, optional
        Step between results to be visualized. Default is 1.
    format : str, optional
        Format of the output file. Default is 'vtk'.
    wake : bool, optional
        If True, will include the wake in the postprocessing. Default is True.
    average : bool, optional
        If True, will average the results. Default is False.
    avg_res : int, optional
        Number of results to be averaged. Default is None.
    component : str, list, optional
        Name(s) of the component(s) to be postprocessed. Default is 'all'.

    """
    name: str
    start_res: int
    end_res: int
    variable: Union[Literal['vorticity', 'vorticity_vector', 'velocity',
                            'surface_velocity', 'pressure', 'cp', 'turbulent_viscosity',
                            'vortex_rad'],
                    List[str]]
    step_res: int = 1
    format: Literal['vtk', 'tecplot'] = 'vtk'
    wake: bool = True
    average: bool = False
    avg_res: int = None
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
class Integral(PrintAnalysis):
    """Integral loads analysis.

    Attributes
    ----------
    name : str
        Name of the analyis, will be appended as a suffix to the basename.
    start_res : int
        First result to be visualized.
    end_res : int
        Last result to be visualized.
    step_res : int, optional
        Step between results to be visualized. Default is 1.
    format : str, optional
        Format of the output file. Default is 'dat'.
    average : bool, optional
        If True, will average the results. Default is False.
    component : str, list, optional
        Name(s) of the component(s) to be postprocessed. Default is 'all'.
    reference_tag : str, optional
        Tag of the reference frame to compute the loads in. Default is '0', which corres-
        ponds to the global reference frame.

    """
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


@dataclass
class Hinge(PrintAnalysis):
    """Hinge loads analysis."""
    pass


@dataclass
class Probe(PrintAnalysis):
    """Probe analysis."""
    pass


@dataclass
class FlowField(PrintAnalysis):
    """Flowfield analysis.

    Attributes
    ----------
    name : str
        Name of the analyis, will be appended as a suffix to the basename.
    start_res : int
        First result to be visualized.
    end_res : int
        Last result to be visualized.
    variable : str, list
        Variable(s) to be visualized.
    n_xyz : list, np.ndarray
        Number of points in each direction. For a plane or a line, insert 1 point in the
        applicable direction.
    min_xyz : list, np.ndarray
        Minimum coordinates of the flowfield.
    max_xyz : list, np.ndarray
        Maximum coordinates of the flowfield.
    step_res : int, optional
        Step between results to be visualized. Default is 1.
    format : str, optional
        Format of the output file. Default is 'vtk', can be 'vtk' or 'tecplot'.
    average : bool, optional
        If True, will average the results. Default is False.
    component : str, list, optional
        Name(s) of the component(s) to be postprocessed. Default is 'all'.
    reference_tag : str, optional
        Tag of the reference frame to compute the loads in. Default is '0', which corres-
        ponds to the global reference frame.

    """
    name: str
    start_res: int
    end_res: int
    variable: Union[Literal['vorticity', 'vorticity_vector', 'velocity',
                            'surface_velocity', 'pressure', 'cp', 'turbulent_viscosity',
                            'vortex_rad'],
                    List[str]]
    n_xyz: Union[List[int], np.ndarray]
    min_xyz: Union[List[Union[int, float, np.number]], np.ndarray]
    max_xyz: Union[List[Union[int, float, np.number]], np.ndarray]
    step_res: int = 1
    format: Literal['vtk', 'tecplot'] = 'vtk'
    average: bool = False
    component: Union[str, List[str]] = 'all'
    reference_tag: str = '0'

    def __post_init__(self):
        """Post-constructor method to assign fixed values."""
        # Prepare the target attribute list in correct order
        attlist = deepcopy(list(self.__dict__.keys()))
        attlist.insert(0, 'type')

        # Reconstruct and reorder attributes
        self.type = 'flow_field'
        self.__dict__ = {attlist[i]: self.__dict__[attlist[i]]
                         for i in range(len(attlist))}


@dataclass
class Sectional(PrintAnalysis):
    """Sectional loads analysis."""
    pass


@dataclass
class Chordwise(PrintAnalysis):
    """Chordwise loads analysis."""
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
    analyses: List[Union[Viz, Integral, Hinge, Probe, FlowField, Sectional,
                         Chordwise]]
