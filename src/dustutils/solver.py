import numpy as np

from dataclasses import dataclass, field
from typing import List, Union, Literal
from pathlib import Path
from copy import deepcopy

from dustutils.utils import Printable

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
    u_ref: Union[int, float, np.number] = u_inf
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
    """Model numerical settings.

    Attributes
    ----------
    far_field_ratio_doublet : int, float, np.number, optional
        Ratio with respect to element length to set the thresholds for far field
        approximation. When evaluating the influence of the doublets of an element, if the
        evaluation point is distant more than far_field_ratio_doublet times the character-
        istic length of the element, simplified cheaper far field approximated formulae
        are emplozed instead of the standard ones. The characteristic length of the elem-
        ent is taken as the maximum length of all the element edges.
    far_field_ratio_source : int, float, np.number, optional
        Threshold to use far field approximations, similarly to far_field_ratio_doublet.
        Only applies to three dimensional surface panels.
    doublet_threshold : int, float, np.number, optional
        Distance threshold under which the evaluation point, with respect to a panel, is
        considered inside the plane of the panel.
    rankine_rad : int, float, np.number, optional
        Parameter which sets the radius under which the Rankine approximation of vortexes 
        cores is employed. Used for aerodynamic elements and panels (i.e. everything
        except vortex particles).
    vortex_rad : int, float, np.number, optional
        Parameter which uniformly sets the radius of the vortex particles, only used if
        k_vortex_rad is disabled.
    k_vortex_rad : int, float, np.number, optional
        Coefficient for the automatic computation of the radius of each vortex particle;
        set a negative number to disable the feature and revert to uniform vortex radius.
    cutoff_rad : int, float, np.number, optional
        Parameter which sets the radius under which the vortexes interaction is completely
        set to zero.
    vortstretch : bool, optional
        Calculate the evolution of vorticity of the particles considering the vortex
        stretching.
    vortstretch_from_elems : bool, optional
        Compute also the contribution to the vortex stretching of the particles due to the
        model elements.
    diffusion : bool, optional
        Calculate the evolution of vorticity of the particles considering the vorticity
        diffusion.
    turbulent_viscosity : bool, optional
        Introduce additional turbulent viscosity (Smagorinsky style) to the vorticity di-
        ffusion term. Working only when fmm is True. Experimental.
    penetration_avoidance : bool, optional
        Apply the penetration avoidance algorithm to avoid the penetration of particles
        inside the solid bodies.
    penetration_avoidance_check_radius : int, float, np.number, optional
        Radius multiplication factor of the check radius. All the particles within a dis-
        tance d ≤ PrUrefΔt from each element are checked for potential penetration, where 
        Pr is the multiplication factor. A bigger factor minimizes the risk of penetration
        of extremely fast particles, but affects the performance.
    penetration_avoidance_element_radius : int, float, np.number, optional
        Surface correction element radius multiplication factor. See DUST official docu-
        mentation for more information.
    divergence_filtering : bool, optional
        Whether to employ the divergence filtering to keep the vorticity field divergence-
        free.
    filter_time_scale : int, float, np.number, optional
        Timescale of the time filter to filter the divergence. The input is not the actual
        timescale but the number of simulation timesteps of which the timescale is
        consisting.

    """
    far_field_ratio_doublet: Union[int, float, np.number] = 10.0
    far_field_ratio_source: Union[int, float, np.number] = 10.0
    doublet_threshold: Union[int, float, np.number] = 1e-6
    rankine_rad: Union[int, float, np.number] = 0.1
    vortex_rad: Union[int, float, np.number] = 0.1
    k_vortex_rad: Union[int, float, np.number] = 1.0
    cutoff_rad: Union[int, float, np.number] = 0.001
    vortstretch: bool = True
    vortstretch_from_elems: bool = False
    diffusion: bool = True
    turbulent_viscosity: bool = False
    penetration_avoidance: bool = False
    penetration_avoidance_check_radius: Union[int, float, np.number] = 5.0
    penetration_avoidance_element_radius: Union[int, float, np.number] = 1.5
    divergence_filtering: bool = True
    filter_time_scale: Union[int, float, np.number] = 40.0


@dataclass
class FMMOpts(Printable):
    """Fast multipole method (FMM) settings."""
    fmm: bool
    box_length: Union[int, float, np.number] = None
    n_box: Union[List[int], np.ndarray] = None
    octree_origin: Union[List[int, float, np.number], np.ndarray] = None
    n_octree_levels: int = None
    min_octree_part: int = None
    multipole_degree: int = None
    dynamic_layers: bool = False
    nmax_octree_levels: int = None
    leaves_time_ratio: Union[int, float, np.number] = None
    fmm_panels: bool = False
    viscosity_effects: bool = False
    particles_redistribution: bool = False
    particles_redistribution_ratio: Union[int, float, np.number] = 3.0
    octree_level_solid: int = None
    turbulent_viscosity: bool = False
    hcas: bool = False
    hcas_time: Union[int, float, np.number] = None
    hcas_velocity: Union[List[int, float, np.number], np.ndarray]\
        = None
    refine_wake: bool = True
    k_refine: int = 1


@dataclass
class LLOpts(Printable):
    """Lifting line model options."""
    ll_solver: Literal['GammaMethod', 'AlphaMethod'] = 'GammaMethod'
    ll_reynolds_corrections: bool = False
    ll_max_iter: int = 100
    ll_tol: Union[int, float, np.number] = 1e-6
    ll_damp: Union[int, float, np.number] = 25.0
    ll_stall_regularisation: bool = True
    ll_stall_regularisation_nelems: int = 1
    ll_stall_regularisation_niters: int = 1
    ll_stall_regularisation_alpha_stall: Union[int, float, np.number] = 15.0
    ll_loads_avl: bool = False
    ll_artificial_viscosity: Union[int, float, np.number] = 0.0
    ll_artificial_viscosity_adaptive: bool = False
    ll_artificial_viscosity_adaptive_alpha: Union[int, float, np.number] = None
    ll_artificial_viscosity_adaptive_dalpha: Union[int, float, np.number] = None


@dataclass
class VLMOpts(Printable):
    """Vortex lattice model options."""
    vl_tol: Union[int, float, np.number] = 1e-4
    vl_relax: Union[int, float, np.number] = 0.3
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
