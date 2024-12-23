/******************************************************************************
 *
 * CitcomSVE-3.0 is a finite element package that solves for dynamic
 * deformation of planetary mantle and crust in a spherical shell under either
 * surface or tidal loads. The mantle and crust may assume a fully compressible
 * or incompressible medium with 3-D viscoelastic structures. The package works
 * with parallel computers with the number of cores or CPUs greater than 12
 * (test runs with >6,000 cores). This version was built on CitcomS-3.1.1,
 * which is a code for planetary mantle convection of a purely viscous flow.
 * This heading was added to the files in CitcomSVE-3.0 that has been modified
 * extensively from or new to CitcomS. Following CitcomS, CitcomSVE has been
 * publicly available as an open source package at github since 2022.
 *
 *               Shijie Zhong, University of Colorado
 ******************************************************************************/


#define NCS      14
	
struct VE_DATA_CONT {
	
	int KERNEL;
	int DIRECT;
	int Heaviside; 
	int variable_te; 
	int read_ice_format; 
	int apply_potential; 
	int plate_margins; 
	int change_of_load; 
	int iteration; // from inputfile 'iteration'
	int SELFG;	   // from inputfile 'self_gravitation'
	int SLE;	   // from inputfile 'apply_SLE'
	int SLE_lmax;  // from inputfile 'SLE_lmax'
	int polar_wander;  // from inputfile 'polar_wander'
	double CApert;	   // from inputfile 'polar_wander_CApert'
	double kf;	   // from inputfile 'polar_wander_kf'
	int pure_visc; // from inputfile 'pure_visc' (0 => visco-elastic)
	double damping;
	
	int  save_binary;	   /* added by Archie, March 2005 */
	char ice_file[100];
	char ocean_file[100];  // filename for ocean function, for SLE
	char timestage_file[100];  // filename for time stages
	char slabgeoid_file[100];
	char visc_file[100];
	char te_file[100];
	char platemargins_file[100];
	int compressible;
	int OneDmodel_read;  // was 1Dmodel_read, modify by tao. should not start with num,
	int OneDmodel_read_g; // whether to read g from model, or compute from density.
	char OneDmodel_file[100];
	
	int stages; // number of stages in ice loading, from inputfile 'stages'
	int stage;	// current stage
	int stages_step[400];  // timestep at end of each stage, from inputfile 'step'
	double stages_timestep[400]; // time per step (ie, stages_time/delta-stages_step)
	double stages_time[400]; // duration of each stage, inputfile 'timestep'
	double stages_ampt[400];

	int NonZero_Init_ICELoad; // added by tao  1: step 0 has non-zero ice loading
	
	int num_icecaps;
	double max_iceV[40];
	double ice_center_t[40];
	double ice_center_f[40];
	double ice_radius[40];
	
	double inc_iceV[5][40];
	
	double potential_vary_PW;  // mult this by (eg) E->init_potential to get
	
	double potential_scaling;  // mult this by (eg) E->init_potential to get
						   // MKS units of potential (J/kg)
	double rotation_rate;
	
	double	 youngs_mod;
	double	 shear_mod;   // from inputfile 'shearmodulus' (Pa)
	double	 tau;	// maxwell time in sec
	double	 tau_in_years;	// maxwell time in yrs (t*tau = dimensional time)
	double	 Te;
	
	double surf_scaling;  // scales nondim height -> nondim stress	(rho.g.R/mu)
	double botm_scaling;   // like surf_scaling but for CMB   (drho.g.R_bot/mu)
	double ice_stress_scale; // like surf_scaling but for ice (uses density_ice)
	double Rsg;
	   // the following are the SLE (Sea Level Equation):
	double incr_ice_volume;  // dimless incremental volume of glaciers
								 // from ice3g . for SLE
	double ocean_area; // ocean_area (dim-less) (for SLE)
	double incr_ocean_area; // incremental
	double eustatic_sea_level;
	double barystatic_sea_level;
	   // the following are for CM (center of mass) calculations:
	double CM[3];	   // Center of Mass (for degree-1)
	double CM_incr[3]; //
	double CM_pot[3]; //
	double CM_incr_ice_static_ocean[3]; // CM motion induced by ice and static ocean load; see apply_new_loads
	double PW[3];	   // polar motion 
	double PW_incr[3];		// polar motion incr
	double PW_pot[3];	   // for 2,1 grav. 
	double Omega_surface[3]; // net rotation of surface (relative to whole earth); see remove_rigid_rot

	int ELASTIC_FROM_FILE; // option for read 3D mu and lambda
	char elastic_file[100];
	  } ve_data_cont;
	
struct SLICE_VE {
	
	double *surf[4][NCS]; // nondim height at surface
						// 1 has instantaneous topography
						// 3 has cumulative topography
	double *botm[4][NCS]; // nondim height at surface
						// 1 has instantaneous topography
						// 3 has cumulative topography
	// note: slice.surf/botm[2] is used for potential calculation
	//		 and for center-of-mass calculations
	double *load[6][NCS];  // nondim stress applied at surfaces
					   // 0 has cumulative surface ice+topo stress
					   // 2 has cumulative CMB topo stress
					   // 0,1 for surface;	2,3 for CMB
	double *elast_thick[NCS];
	double *plate_margins[NCS];
	double *all_load[NCS];
	double *total_load[NCS];
	double *all_potential[NCS];
	double *all_load_potential[NCS];
	double *total_potential[NCS];
	double *init_total_potential[NCS];
	
	double *iceload[6][NCS];
	double *ice_height_prev[NCS];
	double *ice_height_curr[NCS];
   // the following are for the SLE:
	double *ocean_fcn[NCS];
	double *incr_ocean_fcn[NCS];
	double *ocean_fcn_prev[NCS];
	double *ocean_fcn_next[NCS];
	double *static_oceanload[NCS];
	double *dynamic_oceanload[NCS];
	double *total_static_oceanload[NCS];
	double *total_dynamic_oceanload[NCS];
	double *init_oceanload[NCS];

		// modify by tao
	double *all_stress[NCS];  // similiar to surface/botm[2], used in potential & CM
	
  } slice_ve;

