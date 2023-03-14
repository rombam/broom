import numpy as np

from dataclasses import dataclass, field
from typing import List, Union, Literal
from pathlib import Path
from copy import deepcopy

from pdust.utils import Printable

@dataclass
class Gust(Printable):
    """Gust settings.

    Attributes
    ----------
    gust_type : Literal['ACM', 'linear']
        Type of gust.
    gust_origin : List[int, float, np.number], np.ndarray
        (3, ) Position of the point whose airstream velocity is being computed.
    gust_front_direction : List[int, float, np.number], np.ndarray
        (3, ) Unit vector that defines the direction of propagation of the front.
        Units: m/s.
    gust_front_speed : int, float, np.number
        Velocity of propagation of the front in direction gust_front_direction.
        Units: m/s.
    gust_u_des : int, float, np.number
        Gust amplitude. Units: m/s.
    gust_perturbation_direction : List[int, float, np.number], np.ndarray, optional
        (3, ) Direction of perturbation. Units: m/s.
    gust_gradient : int, float, np.number, optional
        Gradient of the gust. Units: m/s.
    gust_start_time : int, float, np.number, optional
        Time at which the gust starts. Units: seconds.

    """
    gust_type: Literal['ACM', 'linear']
    gust_origin: Union[List[int, float, np.number], np.ndarray]
    gust_front_direction: Union[List[int, float, np.number], np.ndarray]
    gust_front_speed: Union[int, float, np.number]
    gust_u_des: Union[int, float, np.number]
    gust_perturbation_direction: Union[List[int, float, np.number], np.ndarray]\
        = np.array([0.0, 0.0, 1.0])
    gust_gradient: Union[int, float, np.number] = 1.0
    gust_start_time: Union[int, float, np.number] = 0.0

    def __post_init__(self):
        """Post-constructor method to assign fixed values."""
        # Prepare the target attribute list in correct order
        attlist = deepcopy(list(self.__dict__.keys()))
        attlist.insert(0, 'gust')

        # Reconstruct and reorder attributes
        self.gust = True
        self.__dict__ = {attlist[i]: self.__dict__[attlist[i]]
                         for i in range(len(attlist))}


@dataclass
class TimeOpts(Printable):
    """Solver time-related settings.

    Attributes
    ----------
    tstart : int, float, np.number
        Start time. Units: seconds.
    tend : int, float, np.number
        End time. Units: seconds.
    timesteps : int, optional
        Number of time steps. Only used if dt is not specified.
    dt : int, float, np.number, optional
        Time step size. Units: seconds.
    dt_out : int, float, np.number, optional
        Output time step size. Units: seconds.
    dt_debug_out : int, float, np.number, optional
        Debug output time step size. Units: seconds.
    restart_from_file : bool, optional
        Whether to restart from a file.
    restart_file : str, Path, optional
        Path to restart file.
    reset_time : bool, optional
        Whether to reset the time to tstart after a restart.   

    """
    tstart: Union[int, float, np.number]
    tend: Union[int, float, np.number]
    timesteps: int = 1
    dt: Union[int, float, np.number] = None
    dt_out: Union[int, float, np.number] = None
    dt_debug_out: Union[int, float, np.number] = None
    restart_from_file: bool = False
    restart_file: Union[str, Path] = None
    reset_time: bool = False
    
    def __post_init__(self):
        if self.dt:
            print(f'{self.dt = } is specified, ignoring {self.timesteps = }')
            self.timesteps = None
    

@dataclass
class FlowOpts(Printable):
    """Flow settings.

    Attributes
    ----------
    u_inf : List[int, float, np.number], np.ndarray
        (3, ) Freestream velocity vector. Units: m/s.
    u_ref : int, float, np.number, optional
        Reference velocity. Units: m/s.
    gust : pdust.gust.Gust, optional
        Gust object.
    p_inf : int, float, np.number, optional
        Freestream pressure. Units: Pa.
    a_inf : int, float, np.number, optional
        Freestream speed of sound. Units: m/s.
    mu_inf : int, float, np.number, optional
        Freestream dynamic viscosity. Units: Pa*s.
    rho_inf : int, float, np.number, optional
        Freestream density. Units: kg/m^3.

    """
    u_inf: Union[List[int, float, np.number], np.ndarray]
    u_ref: Union[int, float, np.number] = np.linalg.norm(u_inf)
    gust: Gust = None
    p_inf: Union[int, float, np.number] = 101325.0
    a_inf: Union[int, float, np.number] = 340.0
    mu_inf: Union[int, float, np.number] = 1.8e-5
    rho_inf: Union[int, float, np.number] = 1.225


@dataclass
class WakeOpts(Printable):
    """Wake settings.

    Attributes
    ----------
    n_wake_panels : int, optional
        Number of wake panels.
    n_wake_particles : int, optional
        Number of wake particles.
    particles_box_min : List[int, float, np.number], np.ndarray, optional
        (3, ) Minimum coordinates of the wake particles box. Units: m.
    particles_box_max : List[int, float, np.number], np.ndarray, optional
        (3, ) Maximum coordinates of the wake particles box. Units: m.
    implicit_panel_scale : int, float, np.number, optional
        Scaling factor for the first panel of the wake, the one which enforces the Kutta
        condition. Scales the second row of the first panel, whose size is determined
        by the local velocity and timestep.
    implicit_panel_min_vel : int, float, np.number, optional
        Minimum velocity for the implicit panel size, in order to avoid an implicit panel
        of zero length. Units: m/s.
    rigid_wake : bool, optional
        Whether to use a rigid wake that evolves according to a given global velocity
        rather than the local velocity field.
    rigid_wake_vel : List[int, float, np.number], np.ndarray, optional
        (3, ) Velocity of the rigid wake. Units: m/s.
    join_te : bool, optional
        Whether to use trailing edge joining for close trailing edges.
    join_te_factor : int, float, np.number, optional
        Factor to join the trailing edge of the wake panels.

    """
    n_wake_panels: int = 1
    n_wake_particles: int = 10000
    particles_box_min: Union[List[int, float, np.number], np.ndarray]\
        = np.array([-10.0, -10.0, -10.0])
    particles_box_max: Union[List[int, float, np.number], np.ndarray]\
        = np.array([10.0, 10.0, 10.0])
    implicit_panel_scale: Union[int, float, np.number] = 0.3
    implicit_panel_min_vel: Union[int, float, np.number] = 1e-8
    rigid_wake: bool = False
    rigid_wake_vel: Union[List[int, float, np.number], np.ndarray]\
        = None
    join_te: bool = False
    join_te_factor: Union[int, float, np.number] = 1.0

@dataclass
class ModelOpts(Printable):
    pass

@dataclass
class FMMOpts(Printable):
    pass

@dataclass
class LLOpts(Printable):
    pass

@dataclass
class VLMOpts(Printable):
    pass

@dataclass
class Settings(Printable):
    """DUST solver settings.

    Attributes
    ----------
    time : pdust.dust.TimeOpts
        Time settings.
    flow : pdust.dust.FlowOpts, optional
        Flow settings.
    wake : pdust.dust.WakeOpts, optional
        Wake settings.
    model : pdust.dust.ModelOpts, optional
        Model settings.
    fmm : pdust.dust.FMMOpts, optional
        Fast multipole method settings.
    liftl : pdust.dust.LLOpts, optional
        Lifting line settings.
    vlm : pdust.dust.VLMOpts, optional
        Vortex lattice method settings.

    """
    time: TimeOpts
    flow: FlowOpts
    wake: WakeOpts
    model: ModelOpts
    fmm: FMMOpts
    liftl: LLOpts
    vlm: VLMOpts
