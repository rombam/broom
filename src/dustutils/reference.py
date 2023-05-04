import numpy as np

from dataclasses import dataclass, field
from typing import List, Union, Literal
from pathlib import Path
from copy import deepcopy

from dustutils.utils import Printable, inplacify


@dataclass
class RotFunc(Printable):
    """Rotation motion input as a simple function.

    Attributes
    ----------
    func : int
        Identifier of the function type. Possible values are 0 for
        a constant function and 1 for a sinusoidal function.
    axis : List[number], np.ndarray
        (3, ) Axis of rotation.
    amplitude : number
        Collective amplitude of the rotation. Units: deg / deg/s.
    omega : number
        Pulsation of the sinusoidal motion. Considered only if func = 1.
    phase : number
        Phase of the sinusoidal motion. Considered only if func = 1.
    offset : number
        Constant offset or rotation rate. Units: deg / deg/s.
    psi_0 : number
        Starting angle of the rotation, considered only if the associated RotType is
        'velocity'.
    """
    func: int
    axis: Union[List[Union[int, float, np.number]], np.ndarray]
    amplitude: Union[int, float, np.number] = 1.0
    omega: Union[int, float, np.number] = 1.0
    phase: Union[int, float, np.number] = 0.0
    offset: Union[int, float, np.number] = 0.0
    psi_0: Union[int, float, np.number] = 0.0

    def __post_init__(self):
        """Post-constructor method to assign fixed values."""
        self.input_type = 'simple_function'


@dataclass
class RotFile(Printable):
    """Rotation motion input from a file.

    Attributes
    ----------
    file : str, Path
        Path to the file containing the rotation motion coordinates.
        It should have two columns, the first containing a series of time values
        adequate for the simulation time, and the second one containing the angles or
        rotation angles.

    """
    file: Union[str, Path]

    def __post_init__(self):
        """Post-constructor method to assign fixed values."""
        self.input_type = 'from_file'

    def __post_init__(self):
        """Post-constructor method to assign fixed values."""
        # Prepare the target attribute list in correct order
        attlist = deepcopy(list(self.__dict__.keys()))
        attlist.insert(0, 'input_type')

        # Reconstruct and reorder attributes
        self.input_type = 'from_file'
        self.__dict__ = {attlist[i]: self.__dict__[attlist[i]]
                         for i in range(len(attlist))}

@dataclass
class PoleFunc:
    """Pole motion input as a simple function.

    Attributes
    ----------
    func : List[int], np.ndarray
        (3, ) Identifier of the function type. Possible values are 0 for
        a constant function and 1 for a sinusoidal function.
        Must be a three element iterable of integers.
    amplitude : number
        Amplitude of the pole motion. It will be applied to all the components.
    vector : List[number], np.ndarray
        (3, ) Relative amplitude of the motion for each component.
    omega : List[number], np.ndarray
        (3, ) Pulsation of each sinusoidal motion. Considered only for the components with
        func[i] = 1.
    phase : List[number], np.ndarray
        (3, ) Phase of each sinusoidal motion. Considered only for the components with
        func[i] = 1.
    offset : List[number], np.ndarray
        (3, ) Constant offset for the components of the pole motion.
    position_0 : List[number], np.ndarray
        (3, ) Starting position of the pole, considered only if the associated PoleType is
        'velocity'.
    """
    func: Union[List[int], np.ndarray]
    amplitude: Union[int, float, np.number] = 1.0
    vector: Union[List[Union[int, float, np.number]], np.ndarray]\
        = np.ones(3)
    omega: Union[List[Union[int, float, np.number]], np.ndarray]\
        = np.ones(3)
    phase: Union[List[Union[int, float, np.number]], np.ndarray]\
        = np.zeros(3)
    offset: Union[List[Union[int, float, np.number]], np.ndarray]\
        = np.zeros(3)
    position_0: Union[List[Union[int, float, np.number]], np.ndarray]\
        = np.zeros(3)

    def __post_init__(self):
        """Post-constructor method to assign fixed values."""
        self.input_type = 'simple_function'


@dataclass
class PoleFile(Printable):
    """Pole motion input from a file.

    Attributes
    ----------
    file : str, Path
        Path to the file containing the pole motion coordinates.
        It should have four columns, the first containing a series of time values
        adequate for the simulation time, and the following three the three components
        of the pole motion.

    """
    file: Union[str, Path]

    def __post_init__(self):
        """Post-constructor method to assign fixed values."""
        # Prepare the target attribute list in correct order
        attlist = deepcopy(list(self.__dict__.keys()))
        attlist.insert(0, 'input_type')

        # Reconstruct and reorder attributes
        self.input_type = 'from_file'
        self.__dict__ = {attlist[i]: self.__dict__[attlist[i]]
                         for i in range(len(attlist))}

    def _fort_strs(self, ignore=[], indent=2):
        return super()._fort_strs(ignore=ignore, indent=indent)

@dataclass
class Pole(Printable):
    """Pole motion of the Reference frame.

    Attributes
    ----------
    variable : PoleType
        Indicator of the motion variable (position/velocity).
    values : PoleFunc, PoleFile
        Motion values. Can be a simple built-in function (PoleFunc) or file-based
        (PoleFile).

    """
    variable: Literal['position', 'velocity']
    values: Union[PoleFunc, PoleFile]

@dataclass
class RotorDOF(Printable):
    """Rotor degree of freedom.

    Attributes
    ----------
    hinge_type : str
        Hinge type. Can be initialized through a string, which can take the values
        'Flap', 'Lag' and 'Pitch'.
    hinge_offset : np.ndarray, optional
        (3, ) Hinge point offset from the reference frame.
    collective : int, float, np.number, optional
        Collective angle of the rotor blade. Unit: degrees.
    cyclic_ampl : int, float, np.number, optional
        Cyclic amplitude of the periodic motion. Unit: degrees.
    cyclic_phas : int, float, np.number, optional
        Cyclic phase of the periodic motion. Unit: degrees.

    """
    hinge_type: Literal['Flap', 'Lag', 'Pitch']
    hinge_offset: Union[List[Union[int, float, np.number]], np.ndarray]\
        = np.zeros(3)
    collective: Union[int, float, np.number] = 0.0
    cyclic_ampl: Union[int, float, np.number] = 0.0
    cyclic_phas: Union[int, float, np.number] = 0.0

    def _fort_strs(self, ignore=[], indent=4):
        return ['\n  dof = {'] + super()._fort_strs(ignore=ignore, indent=indent) + ['  }']

@dataclass
class RotorMulti(Printable):
    """Rotor multiplicity class.

    Attributes
    ----------
    n_blades : int
        Number of rotor blades.
    rot_axis : List[number], np.ndarray
        (3, ) Axis of rotation.
    rot_rate : number
        Rotation rate of the blades around the rotation axis. Units: rad/s.
    psi_0 : number, optional
        Starting angle of the rotor at the beginning of the simulation.
    hub_offset : number, optional
        Offset from the rotation pole of the beginning of the chain of reference frames
        for each blade.
    n_dofs : int, optional
        Number of degrees of freedom for each rotor blade.
    dofs : List[RotorDOF], optional
        List of degrees of freedom for each rotor blade.
        # TODO: Should evaluate whether having a 0.0-value DOF affects performance. By
        default, all three Flap, Lag and Pitch DOFs are created.

    """
    n_blades: int
    rot_axis: Union[List[Union[int, float, np.number]], np.ndarray]
    rot_rate: Union[int, float, np.number]
    psi_0: Union[int, float, np.number] = 0.0
    hub_offset: Union[int, float, np.number] = 0.0
    n_dofs: int = 3
    dofs: List[RotorDOF] = field(default_factory=lambda: [RotorDOF('Flap'),
                                                          RotorDOF('Lag'),
                                                          RotorDOF('Pitch')])

    def __post_init__(self):
        """Post-constructor method to assign fixed values."""
        # Prepare the target attribute list in correct order
        attlist = deepcopy(list(self.__dict__.keys()))
        attlist.insert(0, 'mult_type')

        # Reconstruct and reorder attributes
        self.mult_type = 'rotor'
        self.__dict__ = {attlist[i]: self.__dict__[attlist[i]]
                         for i in range(len(attlist))}

    def _fort_strs(self, ignore=[], indent=2):
        return ['\nmultiplicity = {'] + super()._fort_strs(ignore=ignore, indent=indent) + ['  }']

@dataclass
class Rotation(Printable):
    """Rotation motion of the Reference frame."""
    variable: Literal['angle', 'velocity']
    values: Union[RotFunc, RotFile]

@dataclass
class Motion(Printable):
    """Motion of a Reference frame."""
    pole: Pole = None
    rotation: Rotation = None


@dataclass
class Reference(Printable):
    """Reference frame class.

    Attributes
    ----------
    reference_tag : str
        Reference frame identifier.
    parent_tag : str
        Parent reference frame identifier. If '0', global frame.
    origin : List[number], np.ndarray
        (3, ) X-Y-Z coordinates of the origin, in parent reference frame coordinates.
    orientation : List[number], np.ndarray
        (9, ) i-j-k vectors of the frame.
        # TODO: make this accept 2D 3x3 arrays to make it more user friendly.
    multiplicity : RotorMulti
        Multiplicity object. For now, DUST only supports RotorMulti. Only used if multiple
        is not False.
    motion : Motion
        Motion object. Only used if moving is not False.

    """
    reference_tag: str
    parent_tag: str
    origin: Union[List[Union[int, float, np.number]], np.ndarray]
    orientation: Union[List[Union[int, float, np.number]], np.ndarray]
    multiple: bool = False
    moving: bool = False
    multiplicity: RotorMulti = None
    motion: Motion = None

    @property
    def attitude(self):
        """Return the attitude matrix of the reference frame."""
        return np.array(self.orientation).reshape((3, 3))

    def _translate(self, r):
        """Translate the reference frame.

        Parameters
        ----------
        r : List[number], np.ndarray
            (3, ) Translation vector with respect to the parent reference frame.

        """
        self.origin = np.array(self.origin) + r

    def _rotate(self, R):
        """Rotate the reference frame.

        Parameters
        ----------
        R : np.ndarray
            (3, 3) Rotation matrix.

        """
        self.orientation = np.dot(R, self.attitude).flatten()

    def _rotzyx(self, yaw, pitch, roll):
        """Return the intrinsic rotation matrix for the given Tait-Bryan angles about axes
        Z-Y-X (yaw-pitch-roll).

        Parameters
        ----------
        yaw : number
            Yaw angle. Unit: degrees.
        pitch : number
            Pitch angle. Unit: degrees.
        roll : number
            Roll angle. Unit: degrees.

        Returns
        -------
        np.ndarray
            (3, 3) Rotation matrix.

        """
        yaw = np.deg2rad(yaw)
        pitch = np.deg2rad(pitch)
        roll = np.deg2rad(roll)

        Rz = np.array([[np.cos(yaw), -np.sin(yaw), 0],
                       [np.sin(yaw), np.cos(yaw), 0],
                       [0, 0, 1]])

        Ry = np.array([[np.cos(pitch), 0, np.sin(pitch)],
                       [0, 1, 0],
                       [-np.sin(pitch), 0, np.cos(pitch)]])

        Rx = np.array([[1, 0, 0],
                       [0, np.cos(roll), -np.sin(roll)],
                       [0, np.sin(roll), np.cos(roll)]])

        return np.dot(Rz, np.dot(Ry, Rx))

    @inplacify
    def transform(self, r, yaw, pitch, roll):
        """Transform the reference frame. Translation is performed first, then rotation.
        Rotation is performed about the Z-Y-X axes (yaw-pitch-roll).

        Parameters
        ----------
        r : List[number], np.ndarray
            (3, ) Translation vector with respect to the parent reference frame.
        yaw : number
            Yaw angle. Unit: degrees.
        pitch : number
            Pitch angle. Unit: degrees.
        roll : number
            Roll angle. Unit: degrees.
        inplace : bool, optional
            If True, perform the transformation in place. If False, return a new object.
            Default: True.

        """
        self._translate(r)
        self._rotate(self._rotzyx(yaw, pitch, roll))
        return self
