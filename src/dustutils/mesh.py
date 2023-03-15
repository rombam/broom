import numpy as np

from dataclasses import dataclass
from typing import List, Union, Literal
from pathlib import Path
from copy import deepcopy

from dustutils.utils import Printable


@dataclass
class MeshSymmetry(Printable):
    """Mesh symmetry information.

    Attributes
    ----------
    mesh_symmetry : bool, optional
        Whether the mesh is symmetric with respect to a plane defined by symmetry_normal
        and symmetry_point.
    symmetry_point : List[number], np.ndarray
        (3, ) Array containing the coordinates of the symmetry point.
    symmetry_normal : List[number], np.ndarray
        (3, ) Array containing the coordinates of the symmetry plane normal vector.

    """
    mesh_symmetry: bool = False
    symmetry_point: Union[List[Union[int, float, np.number]], np.ndarray]\
        = np.zeros(3)
    symmetry_normal: Union[List[Union[int, float, np.number]], np.ndarray]\
        = np.array([0.0, 1.0, 0.0])


@dataclass
class MeshMirror(Printable):
    """Mesh mirror information.

    Note
    ----
    As per the DUST documentation, the mirroring operation is the same as the symmetry
    operation, except that it does not keep the original half of the mesh.

    Attributes
    ----------
    mesh_mirror : bool, optional
        Whether the mesh is symmetric with respect to a plane defined by symmetry_normal
        and symmetry_point.
    mirror_point : List[number], np.ndarray
        (3, ) Array containing the coordinates of the mirror point.
    mirror_normal : List[number], np.ndarray
        (3, ) Array containing the coordinates of the mirror plane normal vector.

    """
    mesh_mirror: bool = False
    mirror_point: Union[List[Union[int, float, np.number]], np.ndarray]\
        = np.zeros(3)
    mirror_normal: Union[List[Union[int, float, np.number]], np.ndarray]\
        = np.array([0.0, 1.0, 0.0])

@dataclass
class Point(Printable):
    """Pointwise mesh point.

    Attributes
    ----------
    id : int
        Identifier of the point. Should be a unique integer in increasing value.
        Overrides the id method, although it only applies to the current namespace.
    coordinates : List[number], np.ndarray
        (3, ) Array of coordinates of the point.
    chord : number
        Chord of the section airfoil. Units: length.
    twist : number
        Twist angle of the point section. Units: degrees.
    airfoil : str, Path, optional
        Path to the airfoil coordinate file or airfoil identifier. 
        If a filename, it will assume leading edge to trailing edge along the upper
        surface then back to trailing edge format (Selig).
        If NACA code, will take aNACA 4-digit code (e.g. 'NACA0012') and some NACA 5-digit
        codes.
        This entry will not be used for lifting line geometries ('l').
    airfoil_table : str, Path, optional
        Path to the airfoil polar table in C81 format. Will not be used if using a surface
        panel mesh ('p') or a vortex lattice mesh without airfoil table correction on.
        This entry is required for lifting line geometries ('l').
    section_normal : str, optional
        Type of normal to each point section. Can be 'referenceLine', 'y_axis',
        'y_axis_neg' or 'vector'.
    section_normal_vector : List[number], np.ndarray, optional
        Normal vector to the section. Only applies if section_normal = 'vector'.
    flip_section : bool, optional
        Whether to flip the section upside down. Useful for looping geometries. 

    """
    id: int
    coordinates: Union[List[Union[int, float, np.number]], np.ndarray]
    chord: Union[int, float, np.number]
    twist: Union[int, float, np.number]
    airfoil: Union[str, Path] = None
    airfoil_table: Union[str, Path] = None
    section_normal: Literal['referenceLine', 'y_axis', 'y_axis_neg', 'vector']\
        = 'referenceLine'
    section_normal_vector: Union[List[Union[int, float, np.number]], np.ndarray]\
        = np.array([0.0, 1.0, 0.0])
    flip_section: bool = False

    def _fort_strs(self, ignore=[], indent=4):
        if self.section_normal != 'vector':
            ignore.extend(['section_normal_vector'])
        return ['\npoint = {'] + super()._fort_strs(ignore=ignore, indent=indent) + ['}']

@dataclass
class Line(Printable):
    """Pointwise mesh line. Should be defined by two end Points.

    Attributes
    ----------
    type : str
        Type of the line. Can be 'Straight' or 'Spline'.
    start_point : Point
        Starting point of the Line.
    end_point : Point
        End point of the Line.
    nelems : int
        Number of spanwise elements of the Line.
    type_span : str, optional
        Type of spanwise mesh element distribution. Can be 'uniform', 'cosine',
        'cosineLE' or 'cosineTE'.
    tangent_vec1 : List[number], np.ndarray, optional
        Inboard tangent vector at the start point of a spline. This only applies if the
        line is a spline. Can be ommitted if the surrounding lines are straight, as the
        spline will inherit the tangent vector from the neighboring lines.
    tangent_vec2 : List[number], np.ndarray, optional
        Inboard tangent vector at the end point of a spline. This only applies if the line
        is a spline. Can be ommitted if the surrounding lines are straight, as the spline
        will inherit the tangent vector from the neighboring lines.
    tension : number
        Tension parameter of the Hermitian spline.
    bias : number
        Bias parameter of the Hermitian spline.
    end_points : List[int]
        (2, ) Start and end point IDs. For printing purposes.

    """
    type: Literal['Straight', 'Spline']
    start_point: Point
    end_point: Point
    nelems: int
    type_span: Literal['uniform', 'cosine', 'cosineLE', 'cosineTE'] = 'uniform'
    tangent_vec1: Union[List[Union[int, float, np.number]], np.ndarray]\
        = None
    tangent_vec2: Union[List[Union[int, float, np.number]], np.ndarray]\
        = None
    tension: Union[int, float, np.number] = 0.0
    bias: Union[int, float, np.number] = 0.0

    def __post_init__(self):
        """Post-constructor method to assign fixed values."""
        # Prepare the target attribute list in correct order
        attlist = deepcopy(list(self.__dict__.keys()))
        attlist.insert(1, 'end_points')

        # Reconstruct and reorder attributes
        self.end_points = [self.start_point.id, self.end_point.id]
        self.__dict__ = {attlist[i]: self.__dict__[attlist[i]]
                         for i in range(len(attlist))}

    def _fort_strs(self, ignore=[], indent=4):
        ignore.extend(['start_point', 'end_point'])
        if self.type.upper() == 'STRAIGHT':
            ignore.extend(['tangent_vec1', 'tangent_vec2', 'tension', 'bias'])
        return ['\nline = {'] + super()._fort_strs(ignore=ignore, indent=indent) + ['}']


@dataclass
class Section(Printable):
    """Parametric mesh cross-section.

    Attributes
    ----------
    chord: number
        Chord of the section. Units: length.
    twist: number
        Twist of the section. Units: degrees.
    airfoil: str, Path, optional
        Path to the airfoil coordinate file or airfoil identifier. 
        If a filename, it will assume leading edge to trailing edge along the upper
        surface then back to trailing edge format (Selig).
        If NACA code, will take aNACA 4-digit code (e.g. 'NACA0012') and some NACA 5-digit
        codes.
        This entry will not be used for lifting line geometries.
    airfoil_table: str, Path, optional
        Path to the airfoil polar table in C81 format. Will not be used if using a surface
        panel mesh ('p') or a vortex lattice mesh without airfoil table correction on.
        This entry is required for lifting line geometries ('l').

    """
    chord: Union[int, float, np.number]
    twist: Union[int, float, np.number]
    airfoil: Union[str, Path] = None
    airfoil_table: Union[str, Path] = None

    def _fort_strs(self, ignore=[], indent=0):
        return [''] + super()._fort_strs(ignore=ignore, indent=indent)

@dataclass
class Region(Printable):
    """Parametric mesh region. Should be between two Sections.

    Attributes
    ----------
    span : number
        Span of the region. Units: length.
    sweep : number
        Sweep angle of the region. Units: degrees.
    dihed : number
        Dihedral angle of the region. Units: degrees.
    nelem_span : int
        Number of spanwise elements in the region.
    type_span : str, optional
        Type of spanwise mesh element distribution. Can be 'uniform', 'cosine',
        'cosineLE' or 'cosineTE'.

    """
    span: Union[int, float, np.number]
    sweep: Union[int, float, np.number]
    dihed: Union[int, float, np.number]
    nelem_span: int
    type_span: Literal['uniform', 'cosine', 'cosineLE', 'cosineTE'] = 'uniform'

    def _fort_strs(self, ignore=[], indent=0):
        return [''] + super()._fort_strs(ignore=ignore, indent=indent)

@dataclass
class BaseOpts(Printable):
    """General geometry mesh options.

    Attributes
    ----------
    tol_se_wing : number, optional
        Tolerance in trailing edge merging of nodes.
    inner_product_te : number, optional
        Tolerance for the identification of trailing edges using the inner product of the
        normals.
    proj_te : bool, optional
        Whether to force the projection of the trailing edge in a specific
        direction.
    proj_te_dir : str, optional
        Projection direction of the trailing edge. Can be 'parallel' or 'normal', and will
        apply only if proj_te = True.
    proj_te_vector : List[number], np.ndarray, optional
        (3, ) Vector to specify the direction declared in proj_te_dir.
    suppress_te : bool, optional
        Whether to suppress the trailing edge from the component. Even if it's found, it
        will not release wakes during the simulation.
    offset : List[number], np.ndarray, optional
        (3, ) Offset to apply to the loaded points.
    scaling_factor : number, optional
        Scaling factor applied to the loaded points. Scaling is applied after the offset,
        refer to DUST documentation.
    mesh_flat : bool, optional
        Whether to model lifting lines as a flat surface twisted according to the input
        twist distribution. Applicable only to lifting line geometries.

    """
    tol_se_wing: Union[int, float, np.number] = 0.001
    inner_product_te: Union[int, float, np.number] = -0.5
    proj_te: bool = False
    proj_te_dir: Literal['parallel', 'normal'] = 'parallel'
    proj_te_vector: Union[List[Union[int, float, np.number]], np.ndarray]\
        = np.array([1.0, 0.0, 0.0])
    suppress_te: bool = False
    offset: Union[List[Union[int, float, np.number]], np.ndarray]\
        = np.array([0.0, 0.0, 0.0])
    scaling_factor: Union[int, float, np.number] = 1.0
    mesh_flat: bool = False


@dataclass
class CGNS(Printable):
    """CGNS mesh class.

    Attributes
    ----------
    mesh_file : str, Path
        Path to the mesh file.
    section_name : List[str], optional
        List of CGNS component section names to include in the simulation. If None, will
        include all of the components found.

    """
    mesh_file: Union[str, Path]
    section_name: List[str] = None
    options: BaseOpts = BaseOpts()

    def __post_init__(self):
        """Post-constructor method to assign fixed values."""
        # Prepare the target attribute list in correct order
        attlist = deepcopy(list(self.__dict__.keys()))
        attlist.insert(1, 'mesh_file_type')
        attlist.insert(2, 'el_type')
        attlist.insert(3, 'mesh_flat')

        # Reconstruct and reorder attributes
        self.mesh_file_type = 'cgns'
        self.el_type = 'p'
        self.mesh_flat = False
        self.__dict__ = {attlist[i]: self.__dict__[attlist[i]]
                         for i in range(len(attlist))}


@dataclass
class Parametric(Printable):
    """Parametric type mesh.

    Attributes
    ----------
    el_type : str
        Element type. Can be 'p' for surface panels, 'v' for vortex
        lattice elements and 'l' for lifting line elements.
    nelem_chord : int
        Number of chordwise mesh elements. If the mesh type is surface panels, this
        number of elements will be applied to the upper and lower surfaces, thus producing
        2*nelem_chord total chordwise elements. If the mesh type is lifting line, this
        parameter is ignored.
    sections : List[Section]
        List of sections to be included in the parametric mesh.
    regions : List[Region]
        List of regions to be included in the parametric mesh. Each region is
        contained between two sections, in sequential order.
    airfoil_table_correction : bool, optional
        Whether to correct a VLM mesh using C81 formatted polar tables. Only applies if
        el_type = 'v'.
    twist_linear_interpolation : bool, optional
        Whether to linearly interpolate twist between sections, instead of interpolating
        the point coordinates between two neighboring Sections.
    starting_point : List[number], np.ndarray, optional
        (3, ) Coordinates of the starting point of the mesh, relative to the component
        frame of reference.
    reference_chord_fraction : number, optional
        Fraction of the chord at which to place the axis which will be rotated according
        to sweep and dihedral angles, and around which airfoils are twisted. Ignored
        if element type is lifting line, since this axis will automatically be the
        1/4th chord point.
    type_chord : str, optional
        Type of chordwise mesh element distribution.

    """
    el_type: Literal['p', 'v', 'l']
    nelem_chord: int
    sections: List[Section]
    regions: List[Region]
    airfoil_table_correction: bool = False
    twist_linear_interpolation: bool = False
    starting_point: Union[List[Union[int, float, np.number]], np.ndarray]\
        = np.array([0.0, 0.0, 0.0])
    reference_chord_fraction: Union[int, float, np.number] = 0.0
    type_chord: Literal['uniform', 'cosine', 'cosineLE', 'cosineTE'] = 'uniform'
    options: BaseOpts = BaseOpts()

    def __post_init__(self):
        """Post-constructor method to assign fixed values."""
        # Prepare the target attribute list in correct order
        attlist = deepcopy(list(self.__dict__.keys()))
        attlist.insert(0, 'mesh_file_type')
        attlist.remove('sections')
        attlist.remove('regions')
        attlist.remove('nelem_chord')
        attlist.extend(['nelem_chord', 'sections', 'regions', '_sectlines'])

        # Reconstruct and reorder attributes
        self.mesh_file_type = 'parametric'
        self._sectlines = [self.sections[0]]
        for i in range(len(self.regions)):
            self._sectlines.extend([self.regions[i], self.sections[i+1]])
        self.__dict__ = {attlist[i]: self.__dict__[attlist[i]]
                         for i in range(len(attlist))}

    def _fort_strs(self, ignore=[], indent=0):
        ignore.extend(['sections', 'regions'])
        return super()._fort_strs(ignore=ignore, indent=indent)

    def to_dict(self, ignore=[]):
        ignore.extend(['_sectlines'])
        return super().to_dict(ignore=ignore)

@dataclass
class Pointwise(Printable):
    """Pointwise type mesh.

    Attributes
    ----------
    el_type : str
        Element type. Can be 'p' for surface panels, 'v' for vortex
        lattice elements and 'l' for lifting line elements.
    airfoil_table_correction : bool, optional
        Whether to correct a VLM mesh using C81 formatted polar tables. Only applies if
    twist_linear_interpolation : bool, optional
        Whether to linearly interpolate twist between sections, instead of interpolating
        the point coordinates between two neighboring Sections.
    starting_point : List[number], np.ndarray, optional
        (3, ) Coordinates of the starting point of the mesh, relative to the component
        frame of reference.
    reference_chord_fraction : number, optional
        Fraction of the chord at which to place the axis which will be rotated according
        to sweep and dihedral angles, and around which airfoils are twisted. Ignored
        if element type is lifting line, since this axis will automatically be the
        1/4th chord point.
    nelem_chord : int
        Number of chordwise mesh elements. If the mesh type is surface panels, this
        number of elements will be applied to the upper and lower surfaces, thus producing
        2*nelem_chord total chordwise elements. If the mesh type is lifting line, this
        parameter is ignored.
    type_chord : str, optional
        Type of chordwise mesh element distribution.

    """
    el_type: Literal['p', 'v', 'l']
    nelem_chord: int
    points: List[Point]
    lines: List[Line]
    mesh_symmetry: MeshSymmetry = MeshSymmetry()
    mesh_mirror: MeshMirror = MeshMirror()
    airfoil_table_correction: bool = False
    reference_chord_fraction: Union[int, float, np.number] = 0.0
    type_chord: Literal['uniform', 'cosine', 'cosineLE', 'cosineTE'] = 'uniform'
    options: BaseOpts = BaseOpts()

    def __post_init__(self):
        """Post-constructor method to assign fixed values."""
        # Prepare the target attribute list in correct order
        attlist = deepcopy(list(self.__dict__.keys()))
        attlist.insert(0, 'mesh_file_type')
        attlist.remove('points')
        attlist.remove('lines')
        attlist.remove('nelem_chord')
        attlist.extend(['nelem_chord', 'points', 'lines'])

        # Reconstruct and reorder attributes
        self.mesh_file_type = 'pointwise'
        self.__dict__ = {attlist[i]: self.__dict__[attlist[i]]
                         for i in range(len(attlist))}

@dataclass
class ActuatorDisk(Printable):
    """Actuator disk element.

    Attributes
    ----------
    radius : number
        Radius of the actuator disk. Units: length.
    nstep : int
        Number of radial discretization points.
    axis : int
        Index of the axis the disk is aligned to. Can be 1 for X, 2 for Y or 3 for Z.
        These are the associated reference frame axes.
    traction : number
        Thrust of the actuator disk. Units: N.

    """
    radius: Union[int, float, np.number]
    nstep: int
    axis: int
    traction: Union[int, float, np.number]

    def __post_init__(self):
        """Post-constructor method to assign fixed values."""
        # Prepare the target attribute list in correct order
        attlist = deepcopy(list(self.__dict__.keys()))
        attlist.insert(0, 'mesh_file_type')
        attlist.insert(1, 'el_type')

        # Reconstruct and reorder attributes
        self.mesh_file_type = 'parametric'
        self.el_type = 'a'
        self.__dict__ = {attlist[i]: self.__dict__[attlist[i]]
                         for i in range(len(attlist))}
