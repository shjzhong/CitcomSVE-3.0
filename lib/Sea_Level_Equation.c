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



#include "global_defs.h"
#include "element_definitions.h"
#include <math.h>
#include <stdio.h>

/*============================================================================
 * get_dynamic_oceanload
 * =====================
 * Sets E->slice.oceanload and E->slice.incr_oceanload based on
 *     E->potential[0] (geoid) and E->slice.surf[3] (topo).
 * Calculates the new ocean load due to the current geoid and topography
 *     according to the Sea Level Equation (SLE). The average is removed.
 *     Results have units of nondimensional stress (like iceload and 
 *     slice.load[0]). Also sets slice.incr_oceanload (the change from the
 *     previous oceanload calculation).
 * Updates these variables:
 *     slice.dynamic_oceanload (incremental, used for loading the earth)
 *     slice.total_dynamic_oceanload (cumulative, only used locally)
 *============================================================================*/
void  get_dynamic_oceanload(E,count)
    struct All_variables *E;
    int count;
{
    static double scale_geoid;
    static int been_here=0;
    double c, nu_volume, tmp_oc;
    double tmp,tmp1,tmp2,tmp3;
    int m,i,j ; 
    FILE *fp;
    double t;
    double topo, geoid;
    double total_surface_integral();
    void remove_average();
    void sphere_expansion_output();
    void truncate_Ylm_expansion();
    //
    int verbose=1;   // write out data from SLE calculations
    //


    if (!been_here) { // business upon first invocation of function
        been_here = 1;
        // conversion of potential -> (non-dim) geoid
        scale_geoid =  4.0*M_PI * E->data.grav_const * E->data.density
                     * E->sphere.dradius * E->sphere.dradius ; // dim'ful pot.
        scale_geoid /= E->data.grav_acc;                       // dim'ful geoid
        scale_geoid /= E->sphere.dradius ;                     // dimless geoid
        // allocate and initialize init_ and total_dynamic_oceanload
        for (m=1;m<=E->sphere.caps_per_proc;m++) {
            E->slice_ve.total_dynamic_oceanload[m] = (double *)
                              malloc((E->lmesh.nsf+2)*sizeof(double));
            E->slice_ve.init_oceanload[m] = (double *)
                              malloc((E->lmesh.nsf+2)*sizeof(double));
            for (j=1;j<=E->lmesh.nsf;j++) {
                E->slice_ve.total_dynamic_oceanload[m][j] = 0.0 ;
                E->slice_ve.init_oceanload[m][j] = 0.0 ;
            }
        }
    }

if (E->parallel.me_loc[3] == E->parallel.nprocz-1) {
    if (count==0) {   // business at the beginning of the timestep
        // dynamic_oceanload is calculated (below) as the change with 
        // respect to init_oceanload, so put the previous total into init_oc
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)    
          E->slice_ve.init_oceanload[m][j]= E->slice_ve.total_dynamic_oceanload[m][j];
    }
    
    // Calculate N-U, where N is geoid and U is topography
    //     note: both N and U here should have no l=0 component.
    if (count==0) { // beginning (or not performing) a self-grav iteration
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)  {
            topo = E->slice_ve.surf[3][m][j] ;
            //geoid = 0.0;
            geoid = scale_geoid*E->potential[0][m][j] ;
            E->Xsurf[3][m][j] = (geoid - topo) * E->slice_ve.ocean_fcn[m][j] ;
        }
    } else {       // within a self-grav iteration (includes U increment)
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)  {
            i = j*E->lmesh.noz;
            topo = E->slice_ve.surf[3][m][j] + E->U[m][E->id[m][i].doff[3]] ;
            geoid = scale_geoid*E->potential[0][m][j] ;
            //geoid = 0.0;
            E->Xsurf[3][m][j] = geoid - topo ;
        }

        remove_average(E,E->Xsurf[3],1); // should have no average in order
                                         // to agree with count=0 calculation
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)  
            E->Xsurf[3][m][j] *= E->slice_ve.ocean_fcn[m][j] ;
    }

    // find (dynamic) eustatic sea level, c = [ - integral((N-U).Oc.dA) ] / Ao
    nu_volume = total_surface_integral( E, E->Xsurf[3], 1) ; 
    c =  -nu_volume / E->ve_data_cont.ocean_area ;

    // calculate total_dynamic_oceanload = ( N-U+c )* Oc
    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (j=1;j<=E->lmesh.nsf;j++)
        E->slice_ve.total_dynamic_oceanload[m][j] = 1000.0/917.4*E->ve_data_cont.ice_stress_scale *
                       ( E->Xsurf[3][m][j] + ( c * E->slice_ve.ocean_fcn[m][j] )) ;
    // truncate oceanload, if desired (SLE_lmax==0 => no truncation):
    if (E->ve_data_cont.SLE_lmax)  
        truncate_Ylm_expansion(E, E->slice_ve.total_dynamic_oceanload, 
                                  E->ve_data_cont.SLE_lmax );
    // calculate incremental dynamic_oceanload:
    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (j=1;j<=E->lmesh.nsf;j++)
        E->slice_ve.dynamic_oceanload[m][j]= E->slice_ve.total_dynamic_oceanload[m][j]
                                        - E->slice_ve.init_oceanload[m][j] ;
    // remove average from incremental oceanload:
    remove_average(E,E->slice_ve.dynamic_oceanload,1);

	E->ve_data_cont.barystatic_sea_level = E->ve_data_cont.eustatic_sea_level + c * E->sphere.dradius;
    // print some messages
    if (verbose && E->parallel.me==E->parallel.nprocz-1)
    for (i=0;i<2;i++) {
        fp = (i==0)? stderr : E->fp_out ;
        fprintf(fp,"Dynamic oceanload calculation at %.1f yrs: ",
                    E->ve_data_cont.tau_in_years*E->monitor.elapsed_time );    // time
        fprintf(fp,"ESL(meters)= %g nu_volume= %g Ao= %g \n", 
                    c * E->sphere.dradius,                // eust. SL (meters)
                    nu_volume , E->ve_data_cont.ocean_area );
        fflush(fp);
    }

  }

    return;
}

/*============================================================================= 
 * get_static_oceanload(E)
 * =======================
 * Sets slice.static_oceanload, the incremental static (ie, non-dynamic, 
 * having no dependence on topography or geoid) part of the oceanload.
 * The static oceanload is given by
 *    Lo_stat = ( - M_ice / Ao.rho_ice )* Oc(t)
 *    where M_ice is mass of ice, Ao is ocean area, rho_ice is density of ice, 
 *          and Oc is the ocean function.
 * Updates these variables:
 *     data.ocean_area (adds data.incr_ocean_area)
 *     slice.ocean_fcn (adds slice.incr_ocean_fcn)
 *     slice.static_oceanload (incremental, used for loading the earth)
 *     slice.total_static_oceanload (cumulative, only used locally)
 *============================================================================*/

void get_static_oceanload(E)
    struct All_variables *E;
{
    static int been_here=0;
    static double ice_volume=0.0;
    double c, tmp_oc;
    int n,m,i,j ; 
    FILE *fp;
    double total_surface_integral();
    void truncate_Ylm_expansion();
    void remove_average();
    void get_ocean_fcn();
    // debug flags:
    int debug_const_oc=0; // use constant (not time-dep) ocean fcn
    int verbose=1;   // write out data from SLE calculations
    //

    // set data.incr_ocean_fcn, data.incr_ocean_area, data.incr_ice_volume:
    get_ocean_fcn(E);

    if (!been_here) { // business upon first invocation of function
        been_here = 1;
        // calculate Ao = (non-dim) area of ocean
        E->ve_data_cont.ocean_area = total_surface_integral( E, E->slice_ve.ocean_fcn, 1) ;
        // allocate and initialize total_static_oceanload
        for (m=1;m<=E->sphere.caps_per_proc;m++)   {
            E->slice_ve.total_static_oceanload[m] = (double *)
                              malloc((E->lmesh.nsf+2)*sizeof(double));
            for (n=1;n<=E->lmesh.nsf;n++)
                E->slice_ve.total_static_oceanload[m][n] = 0.0;
        }
    }

    // get new water volume (non-dim) taken from oceans
    ice_volume += E->ve_data_cont.incr_ice_volume ;
    // get new ocean area and ocean_fcn by adding incrementals
    if (!debug_const_oc) {
        E->ve_data_cont.ocean_area += E->ve_data_cont.incr_ocean_area ;
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)    
            E->slice_ve.ocean_fcn[m][j] += E->slice_ve.incr_ocean_fcn[m][j] ;
    }

    c =  -ice_volume / E->ve_data_cont.ocean_area * E->ve_data_cont.ice_stress_scale ;

    // calculate new total_static_oceanload = c * Oc
    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (j=1;j<=E->lmesh.nsf;j++)
        E->Xsurf[3][m][j] = c * E->slice_ve.ocean_fcn[m][j] ;  
    // truncate oceanload, if desired (SLE_lmax==0 => no trucation):
    if (E->ve_data_cont.SLE_lmax) 
        truncate_Ylm_expansion(E, E->Xsurf[3], E->ve_data_cont.SLE_lmax );
    // calculate incremental static_oceanload
    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (j=1;j<=E->lmesh.nsf;j++) {
        E->slice_ve.static_oceanload[m][j] = E->Xsurf[3][m][j]
                                        - E->slice_ve.total_static_oceanload[m][j];
        E->slice_ve.total_static_oceanload[m][j] = E->Xsurf[3][m][j] ;
    }
    // remove average from incremental oceanload:
    remove_average(E,E->slice_ve.static_oceanload,1);

	E->ve_data_cont.eustatic_sea_level = 917.4/1000 * c/E->ve_data_cont.ice_stress_scale * E->sphere.dradius;

    // print some messages
    if (verbose && E->parallel.me==E->parallel.nprocz-1)
    for (i=0;i<2;i++) {
       fp = (i==0)? stderr : E->fp_out ;
       fprintf(fp,"Static oceanload calculation at %.1f yrs: ESL(meters)= %g\n",
               E->ve_data_cont.tau_in_years*E->monitor.elapsed_time ,          // time
               c/E->ve_data_cont.ice_stress_scale * E->sphere.dradius // Eus SL (meters)
              );
       fflush(fp);
    }

    return;
}

/*============================================================================= 
 * get_ocean_fcn(E)
 * ================
 * Sets the variables that store the time-dependent ocean function.
 * The ocean funcion is 1 over ocean, 0 over land. The time dependence comes 
 *     from the ice loading---ie, those parts of the ocean where there was 
 *     glacier coverage are set to 0 (land) in the ocean function. Thus the
 *     ocean fcn depends on the glaciation model (iceNg).
 * Note: an 'epoch' here is either the timespan of glaciation (if stage==0)
 *                              or a 1000 year interval       (if stage==1)
 * This function sets the following variables:
 *     data.incr_ocean_fcn
 *     data.incr_ocean_area
 *     data.incr_ice_volume
 *============================================================================*/

void get_ocean_fcn_FE(E)
    struct All_variables *E;
{
	void parallel_process_termination();
    FILE *fp;
    int m,i,j,n;
    static int ifile=0;  // epoch-index for ocean function files
    static int been_here=0;
    void remove_average();
    double total_surface_integral();
    double total_time, steps_per_epoch;
    int temp0, start_epoch;
    char outfile[255],input_s[200];
    static double ice_volume[16] = { 0.000000000e+00,   // non-dim total volume
                                     1.712407455e-04,   // of ice model glaciers
                                     1.629934092e-04,   // for each epoch
                                     1.540671424e-04,
                                     1.442003899e-04,
                                     1.343447674e-04,
                                     1.205491219e-04,
                                     1.044328696e-04,
                                     8.725370150e-05,
                                     7.192767982e-05,
                                     5.665730818e-05,
                                     3.550639656e-05,
                                     2.450215717e-05,
                                     1.454970359e-05,
                                     5.325542109e-06,
                                     0.000000000e+00 };
    // debug flags:
    int debug_start_w_ice=0; // use the epoch 1 (maximum ice) oc-fcn for start
    int verbose=1;           // not used
    //
   
    if (!been_here) {   // first visit
        been_here = 1;
        // allocate and initialize some arrays for ocean_fcn:
        for (m=1;m<=E->sphere.caps_per_proc;m++)   {
            E->slice_ve.ocean_fcn[m] = (double *)
                          malloc((E->lmesh.nsf+2)*sizeof(double));
            E->slice_ve.incr_ocean_fcn[m] = (double *)
                          malloc((E->lmesh.nsf+2)*sizeof(double));
            E->slice_ve.ocean_fcn_prev[m] = (double *)
                          malloc((E->lmesh.nsf+2)*sizeof(double));
            E->slice_ve.ocean_fcn_next[m] = (double *)
                          malloc((E->lmesh.nsf+2)*sizeof(double));
            for (n=1;n<=E->lmesh.nsf;n++)  {  
                E->slice_ve.ocean_fcn[m][n]      = 0.0;
                E->slice_ve.incr_ocean_fcn[m][n] = 0.0;
                E->slice_ve.ocean_fcn_prev[m][n] = 0.0;
                E->slice_ve.ocean_fcn_next[m][n] = 0.0;
            }
        } // end allocate and initialize

        // read starting ocean fcn:
        if (debug_start_w_ice) start_epoch = 1;
        else                   start_epoch = 0;
        sprintf(outfile,"%s.%d.time%d", 
                        E->ve_data_cont.ocean_file, E->parallel.me, start_epoch);
        fp = fopen(outfile,"r");
        if (fp==NULL) {
            fprintf(stderr,"ERROR: cannot find file %s\n",outfile); 
            parallel_process_termination();
        }
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)  {
            fgets(input_s,200,fp);
            sscanf(input_s,"%d",&(temp0));
            E->slice_ve.ocean_fcn_prev[m][j] = (float)temp0 ;
            E->slice_ve.ocean_fcn[m][j] = E->slice_ve.ocean_fcn_prev[m][j] ;
        }
        fclose(fp);
        ifile = 1;  // index for next-epoch ocean function
    }
    else  { // not the first visit, so check if the elapsed_time is such
            // that we are entering the next epoch
        total_time = E->ve_data_cont.tau_in_years*E->monitor.elapsed_time;
        i = (total_time - E->ve_data_cont.stages_time[0]);
        // note: if E->control.stage==0 (during glaciation),   i<0 (returns)
        //       if E->control.stage==1 (during deglaciation), check ifile
        //       if E->control.stage==2 (after deglaciation),  see next block
        if ( i>((ifile-1)*1000+0.01) ) ifile += 1;    // entering next epoch
        else return ;                                 // remaining in same epoch
    }
    if (E->ve_data_cont.stage==2)  {                 // completed deglaciation 
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++) 
            E->slice_ve.incr_ocean_fcn[m][j] = 0.0;    // 0 => no more change
        E->ve_data_cont.incr_ocean_area = 0.0 ;
        E->ve_data_cont.incr_ice_volume = 0.0  ;
        return;
    }

    // the remaining code reads in data for new epoch, and resets variables
    //     incr_ocean_fcn, incr_ocean_area, and incr_ice_volume

    if (ifile>1)  // if this is not the first visit, curr->prev
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)  
            E->slice_ve.ocean_fcn_prev[m][j] = E->slice_ve.ocean_fcn_next[m][j];   

    // read new epoch ocean fcn into slice.ocean_fcn_next
    sprintf(outfile,"%s.%d.time%d", E->ve_data_cont.ocean_file,E->parallel.me,ifile);
    fp = fopen(outfile,"r");
    if (fp==NULL) {
        fprintf(stderr,"ERROR: cannot find file %s\n",outfile); 
        parallel_process_termination();
    }
    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (j=1;j<=E->lmesh.nsf;j++)  {
        fgets(input_s,200,fp);
        sscanf(input_s,"%d",&(temp0));
        E->slice_ve.ocean_fcn_next[m][j] = (float)temp0 ;
    }
    fclose(fp);

    // steps_per_epoch is the number of timesteps in this epoch
    if (ifile==1) steps_per_epoch = E->ve_data_cont.stages_step[0];
    else  steps_per_epoch =  1000.0/E->ve_data_cont.tau_in_years   // nondim 1kyr
                             / E->ve_data_cont.stages_timestep[E->ve_data_cont.stage] ;

    // set incr_ocean_fcn
    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (j=1;j<=E->lmesh.nsf;j++) 
        E->slice_ve.incr_ocean_fcn[m][j] = ( E->slice_ve.ocean_fcn_next[m][j]
                                         -E->slice_ve.ocean_fcn_prev[m][j]) 
                                       /steps_per_epoch;
    // set incr_ocean_area
    E->ve_data_cont.incr_ocean_area=total_surface_integral(E,E->slice_ve.incr_ocean_fcn,1);

    // set incr_ice_volume (for the incremental non-dim volume of meltwater)
    E->ve_data_cont.incr_ice_volume = (ice_volume[ifile] - ice_volume[ifile-1]) 
                             /steps_per_epoch;
                                 
    return;
}

void get_ocean_fcn(E)
    struct All_variables *E;
{
    FILE *fp;
    int m,i,j,n;
    static int ifile=0;  // epoch-index for ocean function files
    static int been_here=0;
    static int step_prev=0;
    void remove_average();
    void read_reg_grids();
    void parallel_process_termination();
    double total_surface_integral();
    double total_time, steps_per_epoch;
    double ice_volume_new, ice_volume_old,temp0;
    int start_epoch;
    char outfile[255],input_s[200];

    // debug flags:
    int debug_start_w_ice=0; // use the epoch 1 (maximum ice) oc-fcn for start
    int verbose=1;           // not used
    //
   
    if (!been_here) {   // first visit
        // allocate and initialize some arrays for ocean_fcn:
        for (m=1;m<=E->sphere.caps_per_proc;m++)   {
            E->slice_ve.ocean_fcn[m] = (double *)
                          malloc((E->lmesh.nsf+2)*sizeof(double));
            E->slice_ve.incr_ocean_fcn[m] = (double *)
                          malloc((E->lmesh.nsf+2)*sizeof(double));
            E->slice_ve.ocean_fcn_prev[m] = (double *)
                          malloc((E->lmesh.nsf+2)*sizeof(double));
            E->slice_ve.ocean_fcn_next[m] = (double *)
                          malloc((E->lmesh.nsf+2)*sizeof(double));
            for (n=1;n<=E->lmesh.nsf;n++)  {  
                E->slice_ve.ocean_fcn[m][n]      = 0.0;
                E->slice_ve.incr_ocean_fcn[m][n] = 0.0;
                E->slice_ve.ocean_fcn_prev[m][n] = 0.0;
                E->slice_ve.ocean_fcn_next[m][n] = 0.0;
            }
        } // end allocate and initialize

    }

   if (E->ve_data_cont.DIRECT == 0 && been_here!=0)
       return;

    been_here = 1;

        // read starting ocean fcn:

   ifile = E->ve_data_cont.stage + 1;

    // the remaining code reads in data for new epoch, and resets variables
    //     incr_ocean_fcn, incr_ocean_area, and incr_ice_volume

    if (ifile==1)   {   // if this is not the first visit, curr->prev
       sprintf(outfile,"%s.%d",E->ve_data_cont.ocean_file, ifile-1); 
       read_reg_grids(E,outfile,E->slice_ve.ocean_fcn_prev);

       for (m=1;m<=E->sphere.caps_per_proc;m++)
         for (j=1;j<=E->lmesh.nsf;j++)  
            E->slice_ve.ocean_fcn[m][j] = E->slice_ve.ocean_fcn_prev[m][j];   
  
       step_prev = 0;
       }
    else if (ifile>1)   {   // if this is not the first visit, curr->prev
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)  
            E->slice_ve.ocean_fcn_prev[m][j] = E->slice_ve.ocean_fcn_next[m][j];   
        step_prev = E->ve_data_cont.stages_step[E->ve_data_cont.stage-1];
        }

   sprintf(outfile,"%s.%d",E->ve_data_cont.ocean_file, ifile); 
   read_reg_grids(E,outfile,E->slice_ve.ocean_fcn_next);

    // steps_per_epoch is the number of timesteps in this epoch
    steps_per_epoch = E->ve_data_cont.stages_step[E->ve_data_cont.stage] - step_prev;

    // set incr_ocean_fcn
    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (j=1;j<=E->lmesh.nsf;j++) 
        E->slice_ve.incr_ocean_fcn[m][j] = ( E->slice_ve.ocean_fcn_next[m][j]
                                         -E->slice_ve.ocean_fcn_prev[m][j]) 
                                       /steps_per_epoch;

    // set incr_ocean_area
    E->ve_data_cont.incr_ocean_area=total_surface_integral(E,E->slice_ve.incr_ocean_fcn,1);

    // set incr_ice_volume (for the incremental non-dim volume of meltwater)

    ice_volume_new=total_surface_integral(E,E->slice_ve.ice_height_curr,1);
    ice_volume_old=total_surface_integral(E,E->slice_ve.ice_height_prev,1);

    E->ve_data_cont.incr_ice_volume = (ice_volume_new - ice_volume_old) 
                             /steps_per_epoch;

    if (E->parallel.me==0) fprintf(E->fp,"ICE_Volume %g %g\n",ice_volume_new,ice_volume_old);
                                 
    return;
}
