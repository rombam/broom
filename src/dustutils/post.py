from dataclasses import dataclass
from typing import List, Union, Literal
from pathlib import Path
from copy import deepcopy

from .utils import Printable


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
    separate_wake : bool, optional
        If True, will output the wake and the body in separate files. Default is False.
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
    separate_wake: bool = False
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
    """Flow analysis."""
    pass


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


def basic_post(res, name):
    """Create a basic postprocessing object for a DUST case.

    Note
    ----
    Includes a visualization and integral analysis.

    Parameters
    ----------
    res : tuple, list
        (start_res, end_res) for the postprocessing object.
    name : str
        Name of the case.

    """
    viz_post = Viz(name='viz',
                   start_res=res[0],
                   end_res=res[1],
                   variable=['vorticity', 'surface_velocity', 'velocity', 'pressure', 'cp', 'vorticity_vector'],
                   step_res=1)

    int_post = Integral(name='int_wind',
                        start_res=res[0],
                        end_res=res[1],
                        step_res=1)

    int_post = Integral(name='int_ac',
                        start_res=res[0],
                        end_res=res[1],
                        step_res=1,
                        reference_tag='ac')

    post_obj = Post(data_basename=f'Output/{name}',
                    basename=f'Postpro/{name}',
                    analyses=[viz_post, int_post])

    return post_obj
