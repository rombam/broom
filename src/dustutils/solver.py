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
    gust_type : str
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
    gust_origin: Union[List[Union[int, float, np.number]], np.ndarray]
    gust_front_direction: Union[List[Union[int, float, np.number]], np.ndarray]
    gust_front_speed: Union[int, float, np.number]
    gust_u_des: Union[int, float, np.number]
    gust_perturbation_direction: Union[List[Union[int, float, np.number]], np.ndarray]\
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

    @classmethod
    def from_rpm(cls, rpm, step, nrev, tstart=0.0, dt_out=None):
        """Construct a TimeOpts object from RPM value, step size, and number of revolu-
        tions. Mainly aimed at systems with propellers.

        Parameters
        ----------
        rpm : int, float, np.number
            Rotational speed. Units: rpm.
        step : int, float, np.number
            Step size. Units: degrees.
        nrev : int
            Number of revolutions.
        tstart : int, float, np.number, optional
            Start time. Units: seconds.
        dt_out : int, float, np.number, optional
            Output time step size. Units: seconds.

        Returns
        -------
        cls : TimeOpts
            TimeOpts object.

        """
        degs = rpm/60*360
        dt = step/degs
        tend = nrev/(rpm/60.0)

        if dt_out is None:
            dt_out = dt

        return cls(tstart=tstart, tend=tend, dt=dt, dt_out=dt_out)


@dataclass
class FlowOpts(Printable):
    """Flow settings.

    Attributes
    ----------
    u_inf : List[int, float, np.number], np.ndarray
        (3, ) Freestream velocity vector. Units: m/s.
    u_ref : int, float, np.number, optional
        Reference velocity. By default is None and the solver will assume u_ref = u_inf.
        Units: m/s.
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
    u_inf: Union[List[Union[int, float, np.number]], np.ndarray]
    u_ref: Union[int, float, np.number] = None
    gust: Gust = None
    p_inf: Union[int, float, np.number] = 101325.0
    a_inf: Union[int, float, np.number] = 340.0
    mu_inf: Union[int, float, np.number] = 1.8e-5
    rho_inf: Union[int, float, np.number] = 1.225

    def __post_init__(self):
        if np.linalg.norm(self.u_inf) < 1e-6:
            print('Warning: u_inf is close to 0. Setting u_ref to 1.0.')
            self.u_ref = 1.0
        elif self.u_ref is None:
            self.u_ref = np.linalg.norm(self.u_inf)


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
        by the local velocity and timestep. Default: 0.3.
    implicit_panel_min_vel : int, float, np.number, optional
        Minimum velocity for the implicit panel size, in order to avoid an implicit panel
        of zero length. Units: m/s. Default: 1e-8.
    rigid_wake : bool, optional
        Whether to use a rigid wake that evolves according to a given global velocity
        rather than the local velocity field. Default: False.
    rigid_wake_vel : List[int, float, np.number], np.ndarray, optional
        (3, ) Velocity of the rigid wake, required if rigid_wake is True. Units: m/s.
    join_te : bool, optional
        Whether to use trailing edge joining for close trailing edges. Default: False.
    join_te_factor : int, float, np.number, optional
        Factor to join the trailing edge of the wake panels. Default: 1.0.

    """
    n_wake_panels: int = 1
    n_wake_particles: int = 10000
    particles_box_min: Union[List[Union[int, float, np.number]], np.ndarray]\
        = np.array([-10.0, -10.0, -10.0])
    particles_box_max: Union[List[Union[int, float, np.number]], np.ndarray]\
        = np.array([10.0, 10.0, 10.0])
    implicit_panel_scale: Union[int, float, np.number] = None
    implicit_panel_min_vel: Union[int, float, np.number] = None
    rigid_wake: bool = None
    rigid_wake_vel: Union[List[Union[int, float, np.number]], np.ndarray]\
        = None
    join_te: bool = None
    join_te_factor: Union[int, float, np.number] = None

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
        Default: 1.0.
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
    k_vortex_rad: Union[int, float, np.number] = None
    cutoff_rad: Union[int, float, np.number] = 0.001
    vortstretch: bool = True
    vortstretch_from_elems: bool = False
    diffusion: bool = True
    penetration_avoidance: bool = False
    penetration_avoidance_check_radius: Union[int, float, np.number] = 5.0
    penetration_avoidance_element_radius: Union[int, float, np.number] = 1.5
    divergence_filtering: bool = True
    filter_time_scale: Union[int, float, np.number] = 40.0


@dataclass
class FMMOpts(Printable):
    """Fast multipole method (FMM) settings.

    Attributes
    ----------
    fmm : bool
        Whether to use the fast multipole method for particles evolution.
    fmm_panels : bool, optional
        Whether to use the fast multipole method also on panels.
    box_length : int, float, np.number, optional
        Length of the level 1 cubic boxes composing the octree.
    n_box : List[int], np.ndarray, optional
        (3, ) Number of base level 1 cubic boxes in each spatial direction.
    octree_origin : List[int, float, np.number], np.ndarray, optional
        (3, ) Origin of the octree. Starting from the origin, the octree extends in each
        direction of n_box times box_length.
    n_octree_levels : int, optional
        Number of levels in which the base boxes are divided. At each level the upper
        level boxes are divided into eight half sized boxes.
    min_octree_part : int, optional
        Minimum number of particles contained in an octree box in order to consider it a
        leaf (lowest level box). If not enough particles are contained in a box, the box
        is not considered and the particles are gathered at the higher level parent box.
    multipole_degree : int, optional
        Degree of the expansions in the multipole method.
    dynamic_layers : bool, optional
        Use dynamic octree layers, i.e. a further division layer in the octree is added
        everytime the time spent in the particle to particle calculations is greater than
        the one spent in the fast multipole part.
    nmax_octree_levels : int, optional
        Maximum number of divisions allowed during dynamic layers. The number of starting
        layers is still n_octree_levels, which however might increase during the simula-
        tion, but are always kept under a maximum number.
    leaves_time_ratio : int, float, np.number, optional
        Ratio between the time spent in the particle to particle computations in the lea-
        ves of the octree with respect to the rest of the fast multipole computations that
        triggers the increase of the octree levels.
    viscosity_effects : bool, optional
        Whether to take into account viscosity effects on the geometry surface, enabling
        the release of vortex particles from different points on the geometry surface.
        Experimental feature.
    particles_redistribution : bool, optional
        Redistribute the particles having a small intensity to the neighboring ones. Acti-
        ve only if fmm is True. Redistributed particles are then erased to reduce the tot-
        al number of particles.
    particles_redistribution_ratio : int, float, np.number, optional
        Ratio that controls the intensity threshold for particle redistribution. Refer to
        DUST documentation for more information.
    octree_level_solid : int, optional
        Controls the resolution of the octree levels containing solid boundaries when
        marking them as solid. Refer to DUST documentation for more information.
    turbulent_viscosity : bool, optional
        Employ an additional turbuleent viscosity, using a Smagorinsky-like model to take
        into account the dissipation of energy in turbulent conditions towards small, not
        resolved turbulent scales. Only applicable if fmm is True.
    hcas : bool, optional
        Hover Convergence Augmentation System. Refer to DUST documentation for more infor-
        mation.
    hcas_time : int, float, np.number, optional
        Duration of HCAS application. Units: seconds.
    hcas_velocity : List[int, float, np.number], np.ndarray
        Velocity to apply to the particles during the HCAS use.
    refine_wake : bool, optional
        Allow splitting the wake panel into subparticles with a division that is sub-mul-
        tiple of the shortest panel edge.
    k_refine : int, optional
        Number of subdivision of the shortest panel edge if refine_wake is True.

    """
    fmm: bool
    fmm_panels: bool = False
    box_length: Union[int, float, np.number] = None
    n_box: Union[List[int], np.ndarray] = None
    octree_origin: Union[List[Union[int, float, np.number]], np.ndarray] = None
    n_octree_levels: int = None
    min_octree_part: int = None
    multipole_degree: int = None
    dynamic_layers: bool = False
    nmax_octree_levels: int = None
    leaves_time_ratio: Union[int, float, np.number] = None
    viscosity_effects: bool = False
    particles_redistribution: bool = False
    particles_redistribution_ratio: Union[int, float, np.number] = 3.0
    octree_level_solid: int = None
    turbulent_viscosity: bool = False
    hcas: bool = False
    hcas_time: Union[int, float, np.number] = None
    hcas_velocity: Union[List[Union[int, float, np.number]], np.ndarray]\
        = None
    refine_wake: bool = False
    k_refine: int = 1


@dataclass
class LLOpts(Printable):
    """Lifting line model options.

    Attributes
    ----------
    ll_solver : str, optional
        Lifting line solver to be employed. Can be GammaMethod or AlphaMethod. Refer to
        DUST documentation for more information.
    ll_reynolds_corrections : bool, optional
        Employ a Reynolds number correction to obtain an extrapolation of the lifting li-
        ne tables at the simulation conditions Reynolds number if different than the ones
        provided in the lookup table.
    ll_reynolds_corrections_nfact : int, float, np.number, optional
        Power law exponent of the Reynolds correction extrapolation formula. Check the
        DUST documentation for more information.
    ll_max_iter : int, optional
        Number of iterations of the fixed point non-linear solver used to obtain the lift-
        ing lines solution.
    ll_tol : int, float, np.number, optional
        Relative tolerance at which the fixed point lifting linse solver stops.
    ll_damp : int, float, np.number, optional
        Value of the damping coefficient employed during fixed point iterations, to sup-
        press possible oscillations.
    ll_stall_regularisation : bool, optional
        Avoid certain elements converging to a stalled state during the first timesteps of
        a simulation while the rest are not stalled. Refer to the DUST documentation for
        more information.
    ll_stall_regularisation_nelems : int, optional
        Number of lifting line elements to correct in case of isolated stall among non
        stalled elements. Cannot be higher than 1 at the moment.
    ll_stall_regularisation_niters : int, optional
        Number of lifting line iterations between two regularisation processes.
    ll_stall_regularisation_alpha_stall : int, float, np.number, optional
        Stall angle used as a threshold for regularisation process. Units: degrees.
    ll_loads_avl : bool, optional
        Use AVL expression for the inviscid load contributions of lifting line elements,
        in the same way of load computation used for vortex lattice elements.
    ll_artificial_viscosity: int, float, np.number, optional
        Artificial viscosity used to spatially regularize the solution with a gaussian
        kernel, to be used if ll_solver = 'AlphaMethod' to regularize post-stall situa-
        tions. Refer to DUST documentation for more information.
    ll_artificial_viscosity_adaptive : bool, optional
        Whether to use an adaptive strategy to introduce artificial viscosity for regular-
        ization, in order to regularize post stall configuration while not influencing non
        stalled configurations.
    ll_artificial_viscosity_adaptive_alpha : int, float, np.number, optional
        Angle of attack after which the full artificial viscosity is introduced, thus
        after which the maximum regularization is operated. Should be set around or over
        the stall angle. Only applicable if ll_artificial_viscosity_adaptive is True.
        Units: degrees.
    ll_artificial_viscosity_adaptive_dalpha : int, float, np.number, optional
        Angle of attack range before ll_artificial_viscosity_adaptive_alpha where the
        artificial viscosity is gradually introduced from zero to the maximum value set in
        ll_artificial_viscosity.

    """
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
    """Vortex lattice model options.

    Attributes
    ----------
    vl_tol : int, float, np.number, optional
        Tolerance for the absolute error on lift coefficient in fixed point iteration
        for corrected vortex lattice.
    vl_relax : int, float, np.number, optional
        Constant relaxation factor for RHS update.
    aitken_relaxation : bool, optional
        Activate the Aitken acceleration and stabilization method involved during the
        vortex lattice fixed point iterations. The initial relaxation is the one set in
        vl_relax.
    vl_maxiter : int, optional
        Maximum number of iterations for the correction of the vortex lattice.
    vl_start_step : int, optional
        Time step at which the vortex lattice correction starts.
    vl_dynstall : bool, optional
        Activate the Boeing dynamic stall correction for the vortex lattice. Experimental.
    vl_average : bool, optional
        Activate the averaging of the vortex lattice correction over a number of time
        steps, set in vl_average_iter. This may further stabilize the solution in case of
        stalled conditions.
    vl_average_iter : int, optional
        Number of time steps over which the vortex lattice correction is averaged.

    """
    vl_tol: Union[int, float, np.number] = 1e-4
    vl_relax: Union[int, float, np.number] = 0.3
    aitken_relaxation: bool = True
    vl_maxiter: int = 100
    vl_start_step: int = 1
    vl_dynstall: bool = False
    vl_average: bool = False
    vl_average_iter: int = 10


@dataclass
class Settings(Printable):
    """DUST solver settings.

    Attributes
    ----------
    time : pdust.dust.TimeOpts
        Time settings.
    flow : pdust.dust.FlowOpts
        Flow settings.
    fmm : pdust.dust.FMMOpts, optional
        Fast multipole method settings.
        NOTE: if not specified, fmm will not be used, in contrast to the regular DUST
        logic. This is for ease of implementation and may change.
    wake : pdust.dust.WakeOpts, optional
        Wake settings.
    model : pdust.dust.ModelOpts, optional
        Model settings.
    liftl : pdust.dust.LLOpts, optional
        Lifting line settings.
    vlm : pdust.dust.VLMOpts, optional
        Vortex lattice method settings.

    """
    time: TimeOpts
    flow: FlowOpts
    fmm: FMMOpts = FMMOpts(fmm=False)
    wake: WakeOpts = WakeOpts()
    model: ModelOpts = ModelOpts()
    liftl: LLOpts = LLOpts()
    vlm: VLMOpts = VLMOpts()
