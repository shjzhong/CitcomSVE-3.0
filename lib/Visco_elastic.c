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

#include <math.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "parsing.h"

void process_new_velocity(struct All_variables *E,int ii)
{
    void output_velo_related();
    void get_STD_topo();
    void get_CBF_topo();
    void parallel_process_sync();
    void gravity_geoid();
    void interp_surf_velocity();
    void parallel_process_termination();

    void update_stress_strain();

    void print_surf_topo();
  void print_surf_topo_comp();
    void print_stress_visc();
    void deform_grid();
    void calculate_potential();
  void calculate_potential_comp();
    void calculate_Love_numbers();

    static int been_here=0;


    int m,i,it;
    FILE *fp;

    if(been_here==0) {
       E->monitor.length_scale = E->data.layer_km/E->mesh.layer[2]; /* km */
       E->monitor.time_scale = pow(E->data.layer_km*1000.0,2.0)/
            (E->data.therm_diff*3600.0*24.0*365.25*1.0e6);   /* Million years */
        been_here++;
    }
//printf("after first timestep: before update_stress_strain\n");
    update_stress_strain(E,ii);

         //deform grid, deal with CM apparent motion, update total_V 
//printf("after first timestep: beform deform_grid\n");
    deform_grid(E);
     
//>>>> modify by tao ---- update init_total_potential here, although geruo did it in deform_grid.
    if (E->ve_data_cont.SELFG)   {
      // calculate E->incr_potential[0/1] from slice.surf/botm[2]
      // (note: surf/botm[2] is set in deform_grid)
      /* comment by tao 
      Here I use **all_potential** to replace incr_potential[0/1];
      Here **all_potential** is incremental potential we get since solved the last timestep, which
        includes the effects from incremental displacement and dynamic ocean load.
        but not include the known incremental potential in last timestep, which is already included 
        before.
      *** all_potential here use the **all_stress** calculated in deform_grid.	
      Notice, in deform grid, we need to deal with CM motion, where we need to calculate potential.
      here again, we calculate it. Now the icre disp is shifted by r_cm.
    */
//printf("after first timestep: before calculate_potential_comp\n");
    if (E->ve_data_cont.compressible)
      {
      calculate_potential_comp(E, E->slice_ve.all_potential,0); 
      }
    else	/*comment by tao: incompressible*/
      {
      calculate_potential( E, E->slice_ve.surf[2], E->slice_ve.botm[2],
                       E->incr_potential[0], E->incr_potential[1], 1 );
      }
      

      // Add the incr_potential (due to deformation) into init_potential
      /* comment by tao: update init_total_potential with all_potential*/
    if (E->ve_data_cont.compressible)
      {
      for(m=1;m<=E->sphere.caps_per_proc;m++)
      for (i =1; i <=E->lmesh.nno; ++i)
        {
          E->slice_ve.init_total_potential[m][i] += E->slice_ve.all_potential[m][i];
        }
      }
    else /* comment by tao: incompressible*/
      {
      for (m=1;m<=E->sphere.caps_per_proc;m++)
            for (i=1;i<=E->lmesh.nsf;i++)   {
              E->init_potential[1][m][i] += E->incr_potential[1][m][i];
              E->init_potential[0][m][i] += E->incr_potential[0][m][i];
            }
      }


    //question by tao: should I also do this (below) in apply_new_loads?
      if (E->ve_data_cont.polar_wander) {
         E->ve_data_cont.PW[0] += E->ve_data_cont.PW_incr[0];
         E->ve_data_cont.PW[1] += E->ve_data_cont.PW_incr[1];
        if (E->parallel.me==0) fprintf(E->fp,"m0m1 in process_new_ve %g %g\n",E->ve_data_cont.PW_incr[0],E->ve_data_cont.PW_incr[1]);
         }
    }
  
      // Now init_potential has the potential due to the current (completed)
      // deformation and iceload.

//<<<< end tao ---- update init_potential

    print_stress_visc(E,ii);


  if (E->ve_data_cont.compressible)
      print_surf_topo_comp(E,ii);
  else
    print_surf_topo(E,ii);

    if (E->ve_data_cont.Heaviside==1) {  // single harmonic for benchmark
       calculate_Love_numbers(E,ii);
      }


    parallel_process_sync(E);
    return;
}

 //==========================================================================
 /* ========================================================================== */

void pickup_dt(struct All_variables *E)
{
    int i,j;
    void std2_timestep();

        std2_timestep(E);  // for viscoelastic models

    return;
}


/* ==========================================================================
 * deform_grid(E)
 * ==========================================================================
 * Calculates the new positions of each nodes by adding the displacements onto
 * the old coordinates for each nodes. And inject the new positions onto
 * different levels for multigrid.
 * Also, sets:  
 *       E->slice.surf/botm[1]    incr deformation (from E->U)
 *       E->slice.surf/botm[2]    incr stresses    (no iceload)
 *       E->slice.surf/botm[3]    cum  deformation (+= surf/botm[1])
 *       E->slice.load[0/2]       cum  stresses, w/ iceload (+= surf/botm[2])
 * ========================================================================== */
/* comment by tao:	
question about : shift_U_to_CM 
Here the shift includes the rm displacement calculated in apply_new_loads. 
this function actually shift the cap.V with r_cm, not U. although V is U modified to consistent with boundary condition

but, I think maybe we need modify X SX i.e. grid after this, i.e. with new incremental disp?
*/

void deform_grid(struct All_variables *E)
{
    void inject_grid();
    void mass_matrix();
    void remove_average();
    void get_cmb_topo();   
    void parallel_process_termination();
    void load_to_CM(), load_to_CM_grav(), shift_to_CM(), load_to_CM_grav_comp() ;
    void shift_U_to_CM() ;

    double sint,cost,sinf,cosf,ux,uy,uz,myatan();
    int node,count,i,j,lev,m;

    lev=E->mesh.levmax;

    for (m=1;m<=E->sphere.caps_per_proc;m++)  
      for (i=1;i<=E->lmesh.nno;i++)  {
         sint = E->SinCos[lev][m][0][i];
         sinf = E->SinCos[lev][m][1][i];
         cost = E->SinCos[lev][m][2][i];
         cosf = E->SinCos[lev][m][3][i];

         ux = E->sphere.cap[m].V[1][i]*cost*cosf
            - E->sphere.cap[m].V[2][i]*sinf
            + E->sphere.cap[m].V[3][i]*sint*cosf;
         uy = E->sphere.cap[m].V[1][i]*cost*sinf
            + E->sphere.cap[m].V[2][i]*cosf
            + E->sphere.cap[m].V[3][i]*sint*sinf;
         uz =-E->sphere.cap[m].V[1][i]*sint
            + E->sphere.cap[m].V[3][i]*cost;

         E->X[lev][m][1][i] += ux;
         E->X[lev][m][2][i] += uy;
         E->X[lev][m][3][i] += uz;

         E->SX[lev][m][3][i] += E->U[m][E->id[m][i].doff[3]];
         E->SX[lev][m][1][i] = acos(E->X[lev][m][3][i]/E->SX[lev][m][3][i]);
         E->SX[lev][m][2][i] = myatan(E->X[lev][m][2][i],E->X[lev][m][1][i]);
         
         E->SinCos[lev][m][0][i] = sin(E->SX[lev][m][1][i]);
         E->SinCos[lev][m][1][i] = sin(E->SX[lev][m][2][i]);
         E->SinCos[lev][m][2][i] = cos(E->SX[lev][m][1][i]);
         E->SinCos[lev][m][3][i] = cos(E->SX[lev][m][2][i]);
      }

  if (E->control.NMULTIGRID||E->control.EMULTIGRID)   {
      inject_grid(E); 
      }

    if (E->monitor.solution_cycles%E->ve_data_cont.KERNEL==0 || E->ve_data_cont.DIRECT)  {
        mass_matrix(E);
    }

    /* done with updating grids */

        /* We are at the end of the timestep, so do some post-processing.
         * Put the final incremental deformation into slice.surf/botm[1];
         * these fields will be used for writing out results and accumulating
         * into slice.surf/botm[3].
         * Also put incremental stresses from the current timestep into
         * slice.surf/botm[2] for calculation of CM motion (below).    */
         

    if(E->ve_data_cont.compressible){

      //>>>> add by tao: update the all_stress for compressible version. the compressible version of above.
      //debug: now I am not sure why I treat single-harmonic specially. clearly below j not given proper value.
//			 if (E->ve_data_cont.Heaviside==1)  {   /* single-harmonic */
//	                for (m=1;m<=E->sphere.caps_per_proc;m++)  
//	                for (i=1;i<=E->lmesh.nsf;i++)
//						node = j*E->lmesh.noz;
//	                    E->slice_ve.all_stress[m][node] = E->U[m][E->id[m][node].doff[3]];		// rho is 1? 	
//			 }
//			 else if (E->ve_data_cont.Heaviside==2) { 
//				get_boundary_forces(E,1); // this subroutine should fill all_stress with delta_rho * u and dynamic ocean load.
//			 }
      //<<<< end tao

      // comment by tao: below the code copy from get_boundary_forces without load_to_CM
      int m,ll,mm,i,j,p;
      double density_cmb,density_surf,total, length;
      double modified_plgndr_a(),con,t1;
      void remove_average();
      void parallel_process_termination();
      double total_surface_integral();
      double *TG[4],temp1,temp2;


      density_surf = 1.0;
      density_cmb = (E->data.density_below-E->data.density)/E->data.density;
  
            /* comment by tao: below is the compressible version.*/
            int e,n, node;
            double density_jump, stress_scaling;
            const int onedp = onedvpoints[E->mesh.nsd];
            const int vpts=vpoints[E->mesh.nsd];
            const int lev = E->mesh.levmax;
            stress_scaling =  E->data.density*E->data.grav_acc*E->sphere.dradius/E->ve_data_cont.shear_mod;

            for (m=1;m<=E->sphere.caps_per_proc;m++){
              for (e=1;e<=E->lmesh.nel;e++) {
                /* find delta-rho across the top surface of each element*/
                density_jump = 0.0;
                if ( e%E->lmesh.elz==0 /*&& E->parallel.me_loc[3]==E->parallel.nprocz-1*/ ) { /* at the surface */
                  //for (n=1;n<=vpts;n++){
                  if(E->parallel.me_loc[3]==E->parallel.nprocz-1)
                    density_jump += E->erho[lev][m][e];  //comment by tao: geruo define erho on vpt, here is on element
          else
          density_jump += E->erho[lev][m][e] - 1.0;
                  //}
                  //density_jump = density_jump/vpts;
                }
                else {
                  //for (n=1;n<=vpts;n++){
                    density_jump = density_jump + E->erho[lev][m][e] - E->erho[lev][m][e+1];
                  //}
                  //density_jump = density_jump/vpts;
                }
                /* compute all_stress at top nodes */
                for (j=1;j<=onedp;j++){
                  node = E->ien[m][e].node[j+onedp];
                  E->slice_ve.all_stress[m][node] = density_jump*E->U[m][E->id[m][node].doff[3]]; /* i do not multiply this by grav acc because I only need mass */
                }
              } // end e loop

              /* add dynamic_oceanload to the Earth's surface */
              if (E->ve_data_cont.SLE){
              if ( E->parallel.me_loc[3]==E->parallel.nprocz-1 ) {
                for (j=1;j<=E->lmesh.nsf;j++)    {
                  node = j*E->lmesh.noz;
                  E->slice_ve.all_stress[m][node] += E->slice_ve.dynamic_oceanload[m][j]/stress_scaling; /* note that the rho used in dynamic_oceanload is actually rho_ice */
                                              /* note that dynamic_oceanload is non-dim stress, so convert it to non-dim mass - rho_water * u */
                }
              }
              }
            }

            if (1 /*E->parallel.me_loc[3]==0*/)                              // CMB
              for (m=1;m<=E->sphere.caps_per_proc;m++) {
                for (j=1;j<=E->lmesh.SNEL[lev];j++)    {
                  e = (j-1)*E->lmesh.elz+1;

                  density_jump = 0.0;
                  //for (n=1;n<=vpts;n++)
                  if(E->parallel.me_loc[3]==0)
                    density_jump = density_jump + E->data.density_below/E->data.density -E->erho[lev][m][e];
          else
          density_jump = 1.0 -E->erho[lev][m][e];
                  //density_jump = density_jump/vpts;

                  for (n=1;n<=onedp;n++) {
                    node = E->ien[m][e].node[n];
                    E->slice_ve.all_stress[m][node] = density_jump*E->U[m][E->id[m][node].doff[3]];
                  }
                }
              }
    }
    else{
          // surface:
          if (E->parallel.me_loc[3]==E->parallel.nprocz-1)  {
              for (m=1;m<=E->sphere.caps_per_proc;m++)  
              for (i=1;i<=E->lmesh.nsf;i++)  {
                  node = i*E->lmesh.noz;
                  E->slice_ve.surf[1][m][i] = E->U[m][E->id[m][node].doff[3]];
              }
              remove_average(E,E->slice_ve.surf[1],1);

              if (E->ve_data_cont.Heaviside==1)  {   /* single-harmonic */
                  for (m=1;m<=E->sphere.caps_per_proc;m++)  
                  for (i=1;i<=E->lmesh.nsf;i++)
                      E->slice_ve.surf[2][m][i] = E->slice_ve.surf[1][m][i]
                                              *E->ve_data_cont.surf_scaling;
              }
              else if (E->ve_data_cont.Heaviside==2) {    /* ice-model  */
                  for (m=1;m<=E->sphere.caps_per_proc;m++)  
                  for (i=1;i<=E->lmesh.nsf;i++)
                      E->slice_ve.surf[2][m][i] = E->slice_ve.surf[1][m][i]
                                                *E->ve_data_cont.surf_scaling 
                                              + E->slice_ve.dynamic_oceanload[m][i];
              }
          }
          // cmb:
          if (E->parallel.me_loc[3]==0)  {
              for (m=1;m<=E->sphere.caps_per_proc;m++)  
              for (i=1;i<=E->lmesh.nsf;i++)  {
                  node = (i-1)*E->lmesh.noz+1;
                  E->slice_ve.botm[1][m][i] = E->U[m][E->id[m][node].doff[3]];
              }
              remove_average(E,E->slice_ve.botm[1],0);

              for (m=1;m<=E->sphere.caps_per_proc;m++)  
              for (i=1;i<=E->lmesh.nsf;i++)  {
                  E->slice_ve.botm[2][m][i] = E->slice_ve.botm[1][m][i]
                                          *E->ve_data_cont.botm_scaling;
              }
          }			

    }

        /* Now slice.surf/botm[2] are the dynamic loads due to deformation and 
         * dynamic ocean loads. 
         * Now change the loads in slice.surf/botm[2] so that their CM remains
         * at the origin with load_to_CM. This also adds to the current value
         * of data.CM_incr (that's what the 2 flag means).  Then remove the
         * incremental CM motion from the topography (with shift_to_CM).   */

        if (E->ve_data_cont.change_of_load==0) {
          count=1;   // no more new load or static ocean loads
          }
    else{
          count=2;
      }
//printf("after first timestep: before load_to_CM_grav_comp\n");
    if(E->ve_data_cont.compressible){
     // load_to_CM_grav_comp(E, 0, 1, 1 /*count*/); // Geruo's use 1 even change_of_load is 1, since the load effect is considered by shift_V_to_CM, although I am not sure whether its ok to implement like that.
    	load_to_CM_grav_comp(E, 0, 1, count); // Dec 1 2023 (tao), see also apply_new_loads. Above (Geruo's) has a problem: incremental disp doesn't include CM motion induced by new load
	}
    else{
      load_to_CM_grav( E, E->slice_ve.surf[2], E->slice_ve.botm[2], count );
    }
        

/* adjust topographies, displacement U, and total displacement total_V due to 
shift to CM. This is a shift of a constant vector of CM apparent motion */

    /* comment by tao:	Here the shift includes the rm displacement calculated in apply_new_loads. 
    this function actually shift the cap.V with r_cm, not U. although V is U modified to consistent with boundary condition
      do not shift total_VS, which is updated below

    but, I think maybe we need modify X SX i.e. grid after this, i.e. with new incremental disp?
    */
//printf("after first timestep: before shift_U_to_CM\n");
if(E->parallel.me==0){
  printf("debug deg-1 motion: now with shift_U_to_CM\n");
    printf("CM_incr: %e %e %e\n",E->ve_data_cont.CM_incr[0],E->ve_data_cont.CM_incr[1],E->ve_data_cont.CM_incr[2]);
}
        shift_U_to_CM(E, E->slice_ve.surf[1], E->slice_ve.botm[1]); //comment by tao: this actually shift E->cap.V

//        if (E->parallel.me_loc[3]==E->parallel.nprocz-1) 
//            shift_to_CM(E, E->slice_ve.surf[1] );
//        if (E->parallel.me_loc[3]==0) 
//            shift_to_CM(E, E->slice_ve.botm[1] );


        /* Accumulate the incremental topography into slice.surf/botm[3]
         * (which are the displacement fields written out in print_surf_topo,
         * and used for the dynamic_oceanload). Also accumulate
         * slice.surf/botm[2] into slice.load[0/2].    */


    if (E->ve_data_cont.compressible)
    {
      double stress_scaling = E->data.density*E->data.grav_acc*E->sphere.dradius/E->ve_data_cont.shear_mod;

      /* compute total displacement and total load */
      for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (node=1;node<=E->lmesh.nno;node++) {
          E->sphere.cap[m].total_V[1][node] += E->sphere.cap[m].V[1][node];
          E->sphere.cap[m].total_V[2][node] += E->sphere.cap[m].V[2][node];
          E->sphere.cap[m].total_V[3][node] += E->sphere.cap[m].V[3][node];
          E->slice_ve.total_load[m][node] += E->slice_ve.all_stress[m][node]*E->grav[E->mesh.levmax][m][node]*stress_scaling; // here we need non-dim stress
        }

              // surface:
          if (E->parallel.me_loc[3]==E->parallel.nprocz-1)  {
              for (m=1;m<=E->sphere.caps_per_proc;m++)  
              for (i=1;i<=E->lmesh.nsf;i++) {
                  node = i*E->lmesh.noz;
            //      E->slice.surf[3][m][i] += E->U[m][E->id[m][node].doff[3]];

          // comment by tao: total_VS is also updated in shift_U_to_CM above; here I just overwrite it with identical value.
          E->sphere.cap[m].total_VS[1][i]= E->sphere.cap[m].total_V[1][node];
          E->sphere.cap[m].total_VS[2][i]= E->sphere.cap[m].total_V[2][node];
          E->sphere.cap[m].total_VS[3][i]= E->sphere.cap[m].total_V[3][node];
          
                  E->slice_ve.surf[3][m][i] = E->sphere.cap[m].total_VS[3][i];
          E->slice_ve.surf[1][m][i] = E->sphere.cap[m].total_VS[1][i];
          E->slice_ve.surf[2][m][i] = E->sphere.cap[m].total_VS[2][i];
           //       E->slice.surf[1][m][i] = E->U[m][E->id[m][node].doff[3]];
           // I don't output the incr topo. Instead, I set surf[1] as the incr displacement in south direction,
           // and I set surf[2] as the incr displacement in east direction.
           //       E->slice.surf[1][m][i] = E->U[m][E->id[m][node].doff[1]];
           //       E->slice.surf[2][m][i] = E->U[m][E->id[m][node].doff[2]];

              }
          }
          // cmb:
          if (E->parallel.me_loc[3]==0)  {
              for (m=1;m<=E->sphere.caps_per_proc;m++)  
              for (i=1;i<=E->lmesh.nsf;i++)  {
                  node = (i-1)*E->lmesh.noz+1;
                  E->slice_ve.botm[3][m][i] += E->sphere.cap[m].V[3][node];
                  E->slice_ve.botm[1][m][i] = E->sphere.cap[m].V[3][node]; 
              }
          }

    }
    else{
      // comment by tao: for incompressible case, total_VS is updated in shift_U_to_CM called above. total_VS is cumulative surface disp (both hori & topo, i.e. surf[3]), used in Love Number calculation
          // surface:
          if (E->parallel.me_loc[3]==E->parallel.nprocz-1)  {
              for (m=1;m<=E->sphere.caps_per_proc;m++)  
              for (i=1;i<=E->lmesh.nsf;i++) {
                  E->slice_ve.surf[3][m][i] += E->slice_ve.surf[1][m][i];
                  E->slice_ve.load[0][m][i] += E->slice_ve.surf[2][m][i];
              }
          }
          // cmb:
          if (E->parallel.me_loc[3]==0)  {
              for (m=1;m<=E->sphere.caps_per_proc;m++)  
              for (i=1;i<=E->lmesh.nsf;i++)  {
                  E->slice_ve.botm[3][m][i] += E->slice_ve.botm[1][m][i];
                  E->slice_ve.load[2][m][i] += E->slice_ve.botm[2][m][i];
              }
          }

    }
    
    return;
}


//void deform_grid_mybadtao(struct All_variables *E)
//{
//    void inject_grid();
//    void mass_matrix();
//    void remove_average();
//    void get_cmb_topo();   
//    void parallel_process_termination();
//    void load_to_CM(), load_to_CM_grav(), shift_to_CM() ;
//    void shift_U_to_CM() ;
//
//    double sint,cost,sinf,cosf,ux,uy,uz,myatan();
//    int node,count,i,j,lev,m;
//
//    lev=E->mesh.levmax;
//
//    for (m=1;m<=E->sphere.caps_per_proc;m++)  
//      for (i=1;i<=E->lmesh.nno;i++)  {
//         sint = E->SinCos[lev][m][0][i];
//         sinf = E->SinCos[lev][m][1][i];
//         cost = E->SinCos[lev][m][2][i];
//         cosf = E->SinCos[lev][m][3][i];
//
//         ux = E->sphere.cap[m].V[1][i]*cost*cosf
//            - E->sphere.cap[m].V[2][i]*sinf
//            + E->sphere.cap[m].V[3][i]*sint*cosf;
//         uy = E->sphere.cap[m].V[1][i]*cost*sinf
//            + E->sphere.cap[m].V[2][i]*cosf
//            + E->sphere.cap[m].V[3][i]*sint*sinf;
//         uz =-E->sphere.cap[m].V[1][i]*sint
//            + E->sphere.cap[m].V[3][i]*cost;
//
//         E->X[lev][m][1][i] += ux;
//         E->X[lev][m][2][i] += uy;
//         E->X[lev][m][3][i] += uz;
//
//         E->SX[lev][m][3][i] += E->U[m][E->id[m][i].doff[3]];
//         E->SX[lev][m][1][i] = acos(E->X[lev][m][3][i]/E->SX[lev][m][3][i]);
//         E->SX[lev][m][2][i] = myatan(E->X[lev][m][2][i],E->X[lev][m][1][i]);
//         
//         E->SinCos[lev][m][0][i] = sin(E->SX[lev][m][1][i]);
//         E->SinCos[lev][m][1][i] = sin(E->SX[lev][m][2][i]);
//         E->SinCos[lev][m][2][i] = cos(E->SX[lev][m][1][i]);
//         E->SinCos[lev][m][3][i] = cos(E->SX[lev][m][2][i]);
//      }
//
//  if (E->control.NMULTIGRID||E->control.EMULTIGRID)   {
//      inject_grid(E); 
//      }
//
//    if (E->monitor.solution_cycles%E->ve_data_cont.KERNEL==0 || E->ve_data_cont.DIRECT)  {
//        mass_matrix(E);
//    }
//
//    /* done with updating grids */
//
//        /* We are at the end of the timestep, so do some post-processing.
//         * Put the final incremental deformation into slice.surf/botm[1];
//         * these fields will be used for writing out results and accumulating
//         * into slice.surf/botm[3].
//         * Also put incremental stresses from the current timestep into
//         * slice.surf/botm[2] for calculation of CM motion (below).    */
//        
//// from tao. 11/27 here I used Geruo's implementation: use E->slice.botm as 
//	// all horizontal slice.
//
//		// surface:
////        if (E->parallel.me_loc[3]==E->parallel.nprocz-1)  {
////            for (m=1;m<=E->sphere.caps_per_proc;m++)  
////            for (i=1;i<=E->lmesh.nsf;i++)  {
////                node = i*E->lmesh.noz;
////                E->slice_ve.surf[1][m][i] = E->U[m][E->id[m][node].doff[3]];
////            }
////            remove_average(E,E->slice_ve.surf[1],1);
////
////            if (E->ve_data_cont.Heaviside==1)  {   /* single-harmonic */
////                for (m=1;m<=E->sphere.caps_per_proc;m++)  
////                for (i=1;i<=E->lmesh.nsf;i++)
////                    E->slice_ve.surf[2][m][i] = E->slice_ve.surf[1][m][i]
////                                            *E->ve_data_cont.surf_scaling;
////            }
////            else if (E->ve_data_cont.Heaviside==2) {    /* ice-model  */
////                for (m=1;m<=E->sphere.caps_per_proc;m++)  
////                for (i=1;i<=E->lmesh.nsf;i++)
////                    E->slice_ve.surf[2][m][i] = E->slice_ve.surf[1][m][i]
////                                              *E->ve_data_cont.surf_scaling 
////                                            + E->slice_ve.dynamic_oceanload[m][i];
////            }
////        }
//        // cmb:
//
//	//<from tao: ++++---- here I follow geruo
////        if (E->parallel.me_loc[3]==0)  {  
//            for (m=1;m<=E->sphere.caps_per_proc;m++)
//			for (j=1;j<=E->lmesh.noz;j++){
//	            for (i=1;i<=E->lmesh.nsf;i++)  {
//	                node = (i-1)*E->lmesh.noz+j;
//	                E->slice_ve.botm[1][m][i] = E->U[m][E->id[m][node].doff[3]];
//	            }
//				if (j==1)
//	            	remove_average(E,E->slice_ve.botm[1],0);
//				else
//					remove_average(E,E->slice_ve.botm[1],1);
//				
//				for (i=1;i<=E->lmesh.nsf;i++)  {
//            		node = (i-1)*E->lmesh.noz+j;
//            		E->U[m][E->id[m][node].doff[3]] = E->slice_ve.botm[1][m][i];
//          		}
//			}
//	//>from tao: ++++----
//	
//// from tao ----
////	            for (m=1;m<=E->sphere.caps_per_proc;m++)  
////	            for (i=1;i<=E->lmesh.nsf;i++)  {
////	                E->slice_ve.botm[2][m][i] = E->slice_ve.botm[1][m][i]
////	                                        *E->ve_data_cont.botm_scaling;
//// from tao ----
////	            }
//
////        }
//
//if(0){
//	//< from tao ---- geruo did the CM calculation in calculate_potential_comp
//	////< from tao +++ below should use surf/botm[2]
//	//
//	//		// surface:
//	//        if (E->parallel.me_loc[3]==E->parallel.nprocz-1)  {
//	//            for (m=1;m<=E->sphere.caps_per_proc;m++)  
//	//            for (i=1;i<=E->lmesh.nsf;i++)  {
//	//                node = i*E->lmesh.noz;
//	//                E->slice_ve.surf[1][m][i] = E->U[m][E->id[m][node].doff[3]];
//	//            }
//	//            remove_average(E,E->slice_ve.surf[1],1);
//	//
//	//            if (E->ve_data_cont.Heaviside==1)  {   /* single-harmonic */
//	//                for (m=1;m<=E->sphere.caps_per_proc;m++)  
//	//                for (i=1;i<=E->lmesh.nsf;i++)
//	//                    E->slice_ve.surf[2][m][i] = E->slice_ve.surf[1][m][i]
//	//                                            *E->ve_data_cont.surf_scaling;
//	//            }
//	//            else if (E->ve_data_cont.Heaviside==2) {    /* ice-model  */
//	//                for (m=1;m<=E->sphere.caps_per_proc;m++)  
//	//                for (i=1;i<=E->lmesh.nsf;i++)
//	//                    E->slice_ve.surf[2][m][i] = E->slice_ve.surf[1][m][i]
//	//                                              *E->ve_data_cont.surf_scaling 
//	//                                            + E->slice_ve.dynamic_oceanload[m][i];
//	//            }
//	//        }
//	//		// cmb:
//	//        if (E->parallel.me_loc[3]==0)  { 
//	//            for (m=1;m<=E->sphere.caps_per_proc;m++)
//	//	            for (i=1;i<=E->lmesh.nsf;i++)  {
//	//	                node = (i-1)*E->lmesh.noz+1;
//	//	                E->slice_ve.botm[1][m][i] = E->U[m][E->id[m][node].doff[3]];
//	//	            }
//	//	            remove_average(E,E->slice_ve.botm[1],0);
//	//
//	//	            for (m=1;m<=E->sphere.caps_per_proc;m++)  
//	//	            for (i=1;i<=E->lmesh.nsf;i++)  {
//	//	                E->slice_ve.botm[2][m][i] = E->slice_ve.botm[1][m][i]
//	//	                                        *E->ve_data_cont.botm_scaling;
//	//	            }			
//	//		}
//	//
//	////> from tao +++ 
//	//
//	//
//	//        /* Now slice.surf/botm[2] are the dynamic loads due to deformation and 
//	//         * dynamic ocean loads. 
//	//         * Now change the loads in slice.surf/botm[2] so that their CM remains
//	//         * at the origin with load_to_CM. This also adds to the current value
//	//         * of data.CM_incr (that's what the 2 flag means).  Then remove the
//	//         * incremental CM motion from the topography (with shift_to_CM).   */
//	//
//	//        if (E->ve_data_cont.change_of_load==0) 
//	//          count=1;   // no more new load or static ocean loads
//	//        else
//	//          count=2;
//	//
//	//        load_to_CM_grav( E, E->slice_ve.surf[2], E->slice_ve.botm[2], count );
//	//
//	///* adjust topographies, displacement U, and total displacement total_V due to 
//	//shift to CM. This is a shift of a constant vector of CM apparent motion */
//	//
//	//        shift_U_to_CM(E, E->slice_ve.surf[1], E->slice_ve.botm[1]);
//	//
//	////        if (E->parallel.me_loc[3]==E->parallel.nprocz-1) 
//	////            shift_to_CM(E, E->slice_ve.surf[1] );
//	////        if (E->parallel.me_loc[3]==0) 
//	////            shift_to_CM(E, E->slice_ve.botm[1] );
//	//> from tao ---- geruo did the CM calculation in calculate_potential_comp
//}
//
////< from tao +++++  here geruo implemented CM degree-1 motion, then shift the incremental displacement U
//	// also, update the total potential with calculated incremental disp as post-process.
//		if (E->ve_data_cont.SELFG){
//	        get_boundary_forces(E,1);  // calculate all_stress
//
//	        for (i=1;i<=3;i++){
//	          calculate_potential_comp(E,E->slice_ve.all_potential,0,1);
//	          /* shift the U field according to the CM motion*/
//	          shift_U_to_CM(E);
//	        }
//
//	        get_potential_t_comp(E);
//        }
//
////> from tao +++++
//
//        /* Accumulate the incremental topography into slice.surf/botm[3]
//         * (which are the displacement fields written out in print_surf_topo,
//         * and used for the dynamic_oceanload). Also accumulate
//         * slice.surf/botm[2] into slice.load[0/2].    */
//        // surface:
//        if (E->parallel.me_loc[3]==E->parallel.nprocz-1)  {
//            for (m=1;m<=E->sphere.caps_per_proc;m++)  
//            for (i=1;i<=E->lmesh.nsf;i++) {
//                E->slice_ve.surf[3][m][i] += E->slice_ve.surf[1][m][i];
//                E->slice_ve.load[0][m][i] += E->slice_ve.surf[2][m][i];
//            }
//        }
//        // cmb:
//        if (E->parallel.me_loc[3]==0)  {
//            for (m=1;m<=E->sphere.caps_per_proc;m++)  
//            for (i=1;i<=E->lmesh.nsf;i++)  {
//                E->slice_ve.botm[3][m][i] += E->slice_ve.botm[1][m][i];
//                E->slice_ve.load[2][m][i] += E->slice_ve.botm[2][m][i];
//            }
//        }
//    
//    return;
//}

void grid_et_al (struct All_variables *E)
  {

  int lev,nno,nox,noy,noz,i,j,k,m;

  for (lev=E->mesh.levmin;lev<=E->mesh.levmax;lev++)  {
    nox = E->lmesh.NOX[lev];
    noy = E->lmesh.NOY[lev];
    noz = E->lmesh.NOZ[lev];
    nno = E->lmesh.NNO[lev];
    for (m=1;m<=E->sphere.caps_per_proc;m++)  
        for (i=1;i<=nno;i++)  {
          E->X[lev][m][1][i] = E->SX[lev][m][3][i]*sin(E->SX[lev][m][1][i])*cos(E->SX[lev][m][2][i]);
          E->X[lev][m][2][i] = E->SX[lev][m][3][i]*sin(E->SX[lev][m][1][i])*sin(E->SX[lev][m][2][i]);
          E->X[lev][m][3][i] = E->SX[lev][m][3][i]*cos(E->SX[lev][m][1][i]);
          } 

    for (m=1;m<=E->sphere.caps_per_proc;m++)  
      for (i=1;i<=nno;i++)  {
        E->SinCos[lev][m][0][i] = sin(E->SX[lev][m][1][i]);
        E->SinCos[lev][m][1][i] = sin(E->SX[lev][m][2][i]);
        E->SinCos[lev][m][2][i] = cos(E->SX[lev][m][1][i]);
        E->SinCos[lev][m][3][i] = cos(E->SX[lev][m][2][i]);
        }

    }


  return;
  }

// ===================================================================
// compute Love numbers h, k, and l for single harmonic cases 
// for benchmark purpose
// ===================================================================
void  calculate_Love_numbers(struct All_variables *E,int ii)
{
    char outfile[255];
    static FILE *fp1;
    int ll,mm,m,node,nsf,i,j;
    double total_surface_integral();
    void remove_average();
    void sphere_expansion_output();

    static int been_here=0;
    double *ve1[NCS],*ve2[NCS];
    double temp1,const1,const2,const3,h_love,k_love,l_love;
    double hdispersion_err,kdispersion_err;
 
    nsf  = E->lmesh.nsf;

    // save the initial load potential:
    if(E->monitor.solution_cycles==0)  // only save once
    {
    for (m=1;m<=E->sphere.caps_per_proc;m++)  
    for (j=1;j<=E->lmesh.nsf;j++)
        E->Xsurf[3][m][j] =  E->slice_ve.all_load_potential[m][j*E->lmesh.noz];
     // this is not good for deg-1, since deg-1 load is removed. 
    sphere_expansion_output(E,1,E->Xsurf[3],
            E->sphere.sphc_init_load,E->sphere.sphs_init_load,
            E->monitor.solution_cycles,"init_load_pttl");
    }

 if(E->ve_data_cont.apply_potential==0) {    // surface loading
    const1 = (2*E->convection.perturb_ll[0]+1)*E->data.grav_acc/(4.0*M_PI*E->data.grav_const*E->sphere.dradius*E->data.density);
    const2 = (2*E->convection.perturb_ll[0]+1);
    const3 = const1/sqrt(E->convection.perturb_ll[0]*(E->convection.perturb_ll[0]+1));
    }
 else if(E->ve_data_cont.apply_potential==1) {    // apply potential loading
    const1 = E->data.grav_acc/(4.0*M_PI*E->data.grav_const*E->sphere.dradius*E->data.density);
    const2 = 1;
    const3 = const1/sqrt(E->convection.perturb_ll[0]*(E->convection.perturb_ll[0]+1));
    }


 if ( been_here ==0 || (ii % (E->control.record_every) == 0) || ii==E->advection.max_timesteps ) {

   if (E->parallel.me==E->parallel.nprocz-1)  {  // only for one cpu
      hdispersion_err = -1e10;
      kdispersion_err = -1e10;
      for (ll=0;ll<=E->output.llmax;ll++)
      for(mm=0;mm<=ll;mm++)  {
        i = E->sphere.hindex[ll][mm];
        if (ll==E->convection.perturb_ll[0]&&mm==E->convection.perturb_mm[0]) {
          h_love = E->sphere.sphc[2][i];
          if (ll==1 || (ll==2 && mm==1))
            k_love = E->sphere.sphc[3][i];
          else
            k_love = E->sphere.sphc[3][i] - E->sphere.sphc_init_load[i];
          if(ll==1 || (ll==2 && mm==1))
            k_love = k_love*const2/E->convection.perturb_mag[0] - 1.0;
          else
            k_love = k_love*const2/E->convection.perturb_mag[0]; // - 1.0;
          // because the load is applied at cos term
          if (fabs(E->sphere.sphs[2][i]) > hdispersion_err) 
                       hdispersion_err = fabs(E->sphere.sphs[2][i]);
          if (fabs(E->sphere.sphs[3][i]) > kdispersion_err) 
                       kdispersion_err = fabs(E->sphere.sphs[3][i]);
          }
        else {
          if (fabs(E->sphere.sphc[2][i]) > hdispersion_err) 
                       hdispersion_err = fabs(E->sphere.sphc[2][i]);
          if (fabs(E->sphere.sphs[2][i]) > hdispersion_err) 
                       hdispersion_err = fabs(E->sphere.sphs[2][i]);
          if (fabs(E->sphere.sphc[3][i]) > kdispersion_err) 
                       kdispersion_err = fabs(E->sphere.sphc[3][i]);
          if (fabs(E->sphere.sphs[3][i]) > kdispersion_err) 
                       kdispersion_err = fabs(E->sphere.sphs[3][i]);
          }
        }

     h_love = h_love*const1/E->convection.perturb_mag[0];
     hdispersion_err = hdispersion_err*const1/E->convection.perturb_mag[0];
     hdispersion_err = hdispersion_err / h_love; // relative dispersion error
     kdispersion_err = kdispersion_err*const2/E->convection.perturb_mag[0];
     kdispersion_err = kdispersion_err / k_love; // should I use (k_love+1) or k_love here? or perhaps I should remove the sphc/s_init_load for every l,m
     }


    if (E->parallel.me_loc[3] == E->parallel.nprocz-1) {

      for (m=1;m<=E->sphere.caps_per_proc;m++)  {
        ve1[m] = (double *)malloc((nsf+2)*sizeof(double)); 
        ve2[m] = (double *)malloc((nsf+2)*sizeof(double)); 
        }

      for (m=1;m<=E->sphere.caps_per_proc;m++)
      for (i=1;i<=E->lmesh.nsf;i++)  {
        j = E->surf_node[m][i];
        ve1[m][i]=E->sphere.cap[m].total_VS[1][i];
        ve2[m][i]=E->sphere.cap[m].total_VS[2][i];
      }

//      remove_average(E,ve1,1);
//      remove_average(E,ve2,1);

      for (m=1;m<=E->sphere.caps_per_proc;m++)
      for (i=1;i<=E->lmesh.nsf;i++)  {
        ve1[m][i]=ve1[m][i]*ve1[m][i] + ve2[m][i]*ve2[m][i];
      }

    l_love = total_surface_integral(E,ve1,1);
    l_love = sqrt(l_love)*const3/E->convection.perturb_mag[0];

   if (E->parallel.me==E->parallel.nprocz-1)  {  // only for one cpu
     fprintf(E->fp_LN,"%d %.4e %.6e %.6e %.6e %.6e %.6e\n",ii,E->monitor.elapsed_time,h_love,hdispersion_err,k_love,kdispersion_err,l_love);
     fflush(E->fp_LN);
     }

    for (m=1;m<=E->sphere.caps_per_proc;m++) {
      free(ve1[m]);
      free(ve2[m]);
      }
   }
  }

  been_here = 1;

   return;
}
/* =============================================================================
 * print_stress_visc(E,ii)
 * =============================================================================
 * Write out stress and visc
 * ========================================================================== */
void print_stress_visc(struct All_variables *E,int ii)
{

    char outfile[255];
    FILE *fp;
    int m,node,i,j;
    static int been_here=0;
    double error, tot, amp, tau, tau1,tau2,tau3,ana1,ana3,ana2, 
           total,seconds_in_a_year,temp[601];
    int ll,mm;
    static float *stress[NCS],*visc_h,*stress_h;
    static float *stress_xx[NCS],*stress_yy[NCS],*stress_zz[NCS],
                 *stress_xy[NCS],*stress_xz[NCS],*stress_zy[NCS];  // nodal stress tensor for output

    void parallel_process_termination();
    void return_horiz_ave_f();
    void ele_to_nodes(struct All_variables *, float **, float **, int);
    void p_to_nodes(struct All_variables *, double **, float **, int); // should not use ele_to_nodes which is for float only

    if ( been_here==0) {
       for (m=1;m<=E->sphere.caps_per_proc;m++)  
         stress[m] = (float *) malloc ((E->lmesh.nno + 1)*sizeof(float));

       visc_h = (float *) malloc ((E->lmesh.noz + 1)*sizeof(float));
       stress_h = (float *) malloc ((E->lmesh.noz + 1)*sizeof(float));

        for (m=1;m<=E->sphere.caps_per_proc;m++){
          stress_xx[m] = (float *) malloc ((E->lmesh.nno + 1)*sizeof(float));
          stress_yy[m] = (float *) malloc ((E->lmesh.nno + 1)*sizeof(float));
          stress_zz[m] = (float *) malloc ((E->lmesh.nno + 1)*sizeof(float));
          stress_xy[m] = (float *) malloc ((E->lmesh.nno + 1)*sizeof(float));
          stress_xz[m] = (float *) malloc ((E->lmesh.nno + 1)*sizeof(float));
          stress_zy[m] = (float *) malloc ((E->lmesh.nno + 1)*sizeof(float));
        }

        // output coord z:
        if (E->parallel.me<E->parallel.nprocz) { // only for one column
            sprintf( outfile,"%s.coord_z.%d",
                     E->control.data_file,E->parallel.me);
            fp = fopen(outfile,"w");
            m=1;
            for (j=1;j<=E->lmesh.noz;j++)  { 
                fprintf(fp,"%.5e\n",E->sx[m][3][j]); 
            }
            fclose(fp);
        }
    } 

    if ( been_here ==0 || (ii % (E->control.record_every) == 0) || ii==E->advection.max_timesteps)    {

      been_here++;
      return_horiz_ave_f(E,E->Vi,visc_h);
      ele_to_nodes(E,E->S2inv,stress,E->mesh.levmax);
      return_horiz_ave_f(E,stress,stress_h);

      p_to_nodes(E, E->S2xx,stress_xx,E->mesh.levmax);
      p_to_nodes(E, E->S2yy,stress_yy,E->mesh.levmax);
      p_to_nodes(E, E->S2zz,stress_zz,E->mesh.levmax);
      p_to_nodes(E, E->S2xy,stress_xy,E->mesh.levmax);
      p_to_nodes(E, E->S2xz,stress_xz,E->mesh.levmax);
      p_to_nodes(E, E->S2zy,stress_zy,E->mesh.levmax);

      if (E->parallel.me<E->parallel.nprocz)  {
        sprintf(outfile,"%s.horiz_ave.%d.%d", 
                   E->control.data_file,E->parallel.me,ii);
        fp = fopen(outfile,"w");
        fprintf(fp,"%05d %d %.5e %.5e\n",
                    ii,  E->lmesh.noz,                                // timstep
                    E->ve_data_cont.tau*E->monitor.elapsed_time, // current time (years)
                    E->ve_data_cont.tau*E->advection.timestep);  // length of timestep
        m = 1;
        for (j=1;j<=E->lmesh.noz;j++)
          fprintf(fp,"%.5e %.5e %.5e\n",E->sx[m][3][j],visc_h[j],stress_h[j]);
        fclose(fp);
        }

      sprintf(outfile,"%s.visc.%d.%d", 
                   E->control.data_file,E->parallel.me,ii);
      fp = fopen(outfile,"w");
      fprintf(fp,"%05d %d %.5e %.5e\n",
                    ii,  E->lmesh.nno,                                // timstep
                    E->ve_data_cont.tau*E->monitor.elapsed_time, // current time (years)
                    E->ve_data_cont.tau*E->advection.timestep);  // length of timestep
      for (m=1;m<=E->sphere.caps_per_proc;m++)  
        for (j=1;j<=E->lmesh.nno;j++)
          fprintf(fp,"%.5e\n",E->Vi[m][j]);
      fclose(fp);

      // output stress tensor:
      sprintf(outfile,"%s.stress_tensor.%d.%d", 
                   E->control.data_file,E->parallel.me,ii);
      fp = fopen(outfile,"w");
      fprintf(fp,"%05d %d %.5e %.5e\n",
                    ii,  E->lmesh.nno,                                // timstep
                    E->ve_data_cont.tau*E->monitor.elapsed_time, // current time (years)
                    E->ve_data_cont.tau*E->advection.timestep);  // length of timestep
      for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nno;j++)
          fprintf(fp,"%.5e %.5e %.5e %.5e %.5e %.5e\n",
                  stress_xx[m][j],stress_yy[m][j],stress_zz[m][j],
                  stress_xy[m][j],stress_xz[m][j],stress_zy[m][j]);
      fclose(fp);
    }

   return;
   }


/* =============================================================================
 * print_surf_topo(E,ii)
 * =============================================================================
 * Write out topography and potential.
 * ========================================================================== */
void print_surf_topo(struct All_variables *E,int ii)
{

    char outfile[255];
    FILE *fp;
    int m,node,i,j;
    void sphere_expansion_output();
    static int been_here=0;
    double error, tot, amp, tau, tau1,tau2,tau3,ana1,ana3,ana2, 
           total,seconds_in_a_year,temp[601];
    int ll,mm;
    double modified_plgndr_a(),con,t1;
    void parallel_process_termination();

    if ( been_here == 0 )    {

        // write out mesh coordinates:
        if (E->parallel.me_loc[3]==E->parallel.nprocz-1)  {
            sprintf( outfile,"%s.coord_s.%d",
                     E->control.data_file,E->parallel.me);
            fp = fopen(outfile,"w");
            fprintf(fp,"%05d %.5e\n",ii,E->ve_data_cont.tau*E->monitor.elapsed_time);
            j = 0;
            for (m=1;m<=E->sphere.caps_per_proc;m++)  {
                for (j=1;j<=E->lmesh.nsf;j++)  { 
                    i=j*E->lmesh.noz;
                    E->Xsurf[1][m][j] = E->SX[E->mesh.levmax][m][1][i];
                    E->Xsurf[2][m][j] = E->SX[E->mesh.levmax][m][2][i];
                    fprintf(fp, "%.4e %.4e \n", 
                                E->Xsurf[1][m][j],E->Xsurf[2][m][j]);
                }
//                for (j=1;j<=E->lmesh.noz;j++)  { 
//                    fprintf(fp,"%.5e\n",E->SX[E->mesh.levmax][m][3][j]); 
//                }
              }
            fclose(fp);

            // the following is to observe the elastic behavior:
            if (1)  {

                sphere_expansion_output(E,1,E->slice_ve.surf[3],
                        E->sphere.sphc[2],E->sphere.sphs[2],
                        E->monitor.solution_cycles,"tps");

                // The following is incr_potential from both deformation and loads:
                for (m=1;m<=E->sphere.caps_per_proc;m++)  
                for (j=1;j<=E->lmesh.nsf;j++)
                    E->Xsurf[3][m][j] =  E->incr_potential[0][m][j]  
                                       + E->incr_potential[2][m][j] ;

                sphere_expansion_output(E,1,E->Xsurf[3],
                        E->sphere.sphc[3],E->sphere.sphs[3],
                        E->monitor.solution_cycles,"pttl");

                if (E->ve_data_cont.SLE) {
                    // get the oceanload:
                    for (m=1;m<=E->sphere.caps_per_proc;m++)  
                    for (j=1;j<=E->lmesh.nsf;j++)
                        E->Xsurf[3][m][j] =  E->slice_ve.static_oceanload[m][j]
                                           + E->slice_ve.dynamic_oceanload[m][j];
                    sphere_expansion_output(E,1,E->Xsurf[3],
                            E->sphere.sphc[0],E->sphere.sphs[0],
                            E->monitor.solution_cycles,"oceanload");
                }

                for (i=0;i<2;i++) {
                    fp = (i==0)? stderr : E->fp_out ;
                    fprintf(fp, "Wrote initial topo and potential files.\n");
                    fflush(fp);
                }
            }
        }
    }// end if not been_here

    been_here++;

    if ( ((ii % E->control.record_every) == 0) || (ii == E->advection.max_timesteps))    {

        if (E->parallel.me_loc[3]==E->parallel.nprocz-1)  {
            //sprintf(outfile,"%s.topo_s.%d.%d",
            //       E->control.data_file2,E->parallel.me,ii);
            sprintf(outfile,"%s.topo_s.%d.%d", 
                   E->control.data_file,E->parallel.me,ii);
            fp = fopen(outfile,"w");
            fprintf(fp,"%05d %.5e %.5e %.5e %.5e %.5e %.5e\n",
                    ii,                                  // timstep
                    E->ve_data_cont.tau*E->monitor.elapsed_time, // current time (years)
                    E->ve_data_cont.tau*E->advection.timestep,  // length of timestep
              E->ve_data_cont.PW[0],E->ve_data_cont.PW[1], // cumu. polar motion
              E->ve_data_cont.PW_incr[0],E->ve_data_cont.PW_incr[1]); // polar motion

            // The following is incr_potential from both deformation [0] and loads [2]:
            for (m=1;m<=E->sphere.caps_per_proc;m++)  
            for (j=1;j<=E->lmesh.nsf;j++)
                E->Xsurf[3][m][j] =  E->incr_potential[0][m][j]  
                                   + E->incr_potential[2][m][j] ;

            for (m=1;m<=E->sphere.caps_per_proc;m++)  
            for (j=1;j<=E->lmesh.nsf;j++)  { 
                i=j*E->lmesh.noz;
                /*// for writing out coordinates:
                E->Xsurf[1][m][j] = E->SX[E->mesh.levmax][m][1][i];
                E->Xsurf[2][m][j] = E->SX[E->mesh.levmax][m][2][i];
                E->Xsurf[3][m][j] = E->SX[E->mesh.levmax][m][3][i]; */
                /* output oceanloads for debugging
                   fprintf(fp,"%.4e %.4e %.4e %.4e %.4e %.4e\n",
                   E->slice_ve.surf[3][m][j],       // total topography
                   E->slice_ve.surf[1][m][j],       // incr topo
                   E->init_potential[0][m][j],   // total potential
                   E->incr_potential[0][m][j],   // incr potential
                   E->slice_ve.total_static_oceanload[m][j]  // oceanload
                   /E->ve_data_cont.ice_stress_scale,        // (nondim height)
                   E->slice_ve.total_dynamic_oceanload[m][j] //   
                   /E->ve_data_cont.ice_stress_scale         // 
                   ); */
                fprintf(fp,"%.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e\n",
                        E->slice_ve.surf[3][m][j],       // total topography
                        E->slice_ve.surf[1][m][j],       // incr topo
                        E->init_potential[0][m][j],   // total potential
                        E->Xsurf[3][m][j],           // incr potential
                        //E->incr_potential[0][m][j]);  // incr potential due to loads change only
                        E->sphere.cap[m].total_VS[1][j],        
                        E->sphere.cap[m].V[1][i],        
                        E->sphere.cap[m].total_VS[2][j],        
                        E->sphere.cap[m].V[2][i]);        
            }
            fclose(fp);

            sphere_expansion_output(E,1, E->slice_ve.surf[3],
                                         E->sphere.sphc[2], E->sphere.sphs[2],
                                         E->monitor.solution_cycles, "tps");
            sphere_expansion_output(E,1, E->slice_ve.surf[1],
                                         E->sphere.sphc[0],E->sphere.sphs[0],
                                         E->monitor.solution_cycles,"vtps");
            sphere_expansion_output(E,1, E->init_potential[0],
                                         E->sphere.sphc[3],E->sphere.sphs[3],
                                         E->monitor.solution_cycles,"pttl");
            sphere_expansion_output(E,1, E->Xsurf[3],
                                         E->sphere.sphc[0],E->sphere.sphs[0],
                                         E->monitor.solution_cycles,"pttldot");

            for (m=1;m<=E->sphere.caps_per_proc;m++)
                for (j=1;j<=E->lmesh.nsf;j++)  { 
                    i=j*E->lmesh.noz;
                    E->Xsurf[1][m][j] = E->sphere.cap[m].total_VS[1][j];
                    E->Xsurf[2][m][j] = E->sphere.cap[m].total_VS[2][j];
                    }
            sphere_expansion_output(E,1, E->Xsurf[1],
                                         E->sphere.sphc[4],E->sphere.sphs[4],
                                         E->monitor.solution_cycles,"utheta");
            sphere_expansion_output(E,1, E->Xsurf[2],
                                         E->sphere.sphc[5],E->sphere.sphs[5],
                                         E->monitor.solution_cycles,"ufi");

        }

        /*// write out CMB values:
        if (E->parallel.me_loc[3]==0)  {

            sprintf(outfile,"%s.topo_b.%d.%d",
                             E->control.data_file2,E->parallel.me,ii);
            //sprintf(outfile,"%s.topo_b.%d.%d",
            //                 E->control.data_file,E->parallel.me,ii);
            fp = fopen(outfile,"w");
            fprintf(fp,"%05d %.5e\n",ii,E->monitor.elapsed_time);
            for (m=1;m<=E->sphere.caps_per_proc;m++)  
                for (j=1;j<=E->lmesh.nsf;j++)  { 
                    i=(j-1)*E->lmesh.noz+1;
                    E->Xsurf[1][m][j] = E->SX[E->mesh.levmax][m][1][i];
                    E->Xsurf[2][m][j] = E->SX[E->mesh.levmax][m][2][i];
                    E->Xsurf[3][m][j] = E->SX[E->mesh.levmax][m][3][i];
                    fprintf(fp,"%.4e %.4e %.4e %.4e %.4e\n",
                               E->Xsurf[1][m][j],E->Xsurf[2][m][j],
                               E->slice_ve.botm[3][m][j],
                               E->U[m][E->id[m][i].doff[3]],E->potential[m][i]);
                }
            fclose(fp);

            sphere_expansion_output(E,0,E->slice_ve.botm[3],
                                        E->sphere.sphc[0],E->sphere.sphs[0],
                                        E->monitor.solution_cycles,"tpb");

            for (m=1;m<=E->sphere.caps_per_proc;m++)  
                for (j=1;j<=E->lmesh.nsf;j++)  { 
                    i=(j-1)*E->lmesh.noz+1;
                    E->Xsurf[3][m][j] = E->U[m][E->id[m][i].doff[3]];
                }

            sphere_expansion_output(E,0,E->Xsurf[3],
                                     E->sphere.sphc[0],E->sphere.sphs[0],
                                     E->monitor.solution_cycles,"vtpb");
        } */

    } // end if on a recording timestep

    return;
}


/*
 * note: for compressible case: slice_ve.surf[3] has been given values of total_VS[3], which is the cumulative topo.
 * 			but incr topo is not stored in any surf-like variables.
 *
 */
void print_surf_topo_comp(struct All_variables *E,int ii)
//    struct All_variables *E;
//    int ii;
{

    char outfile[255];
    FILE *fp;
    static FILE *fp1,*fp2;
    int m,node,i,j;
    void sphere_expansion_output();
    static int been_here=0;
    double error, tot, amp, tau, tau1,tau2,tau3,ana1,ana3,ana2, 
           total,seconds_in_a_year,temp[601];
    int ll,mm;
    double modified_plgndr_a(),con,t1;
    void parallel_process_termination();
    void output_strain();

    if ( been_here == 0 )    {

        // write out mesh coordinates:
        if (E->parallel.me_loc[3]==E->parallel.nprocz-1)  {
            sprintf( outfile,"%s.coord_s.%d",
                     E->control.data_file,E->parallel.me);
            fp = fopen(outfile,"w");
            fprintf(fp,"%05d %.5e\n",ii,E->ve_data_cont.tau*E->monitor.elapsed_time);
            j = 0;
            for (m=1;m<=E->sphere.caps_per_proc;m++)  
                for (j=1;j<=E->lmesh.nsf;j++)  { 
                    i=j*E->lmesh.noz;
                    E->Xsurf[1][m][j] = E->SX[E->mesh.levmax][m][1][i];
                    E->Xsurf[2][m][j] = E->SX[E->mesh.levmax][m][2][i];
                    fprintf(fp, "%.4e %.4e \n", 
                                E->Xsurf[1][m][j],E->Xsurf[2][m][j]);
                }
            fclose(fp);

            // the following is to observe the elastic behavior:
            if (1)  {

                sphere_expansion_output(E,1,E->slice_ve.surf[3],
                        E->sphere.sphc[2],E->sphere.sphs[2],
                        E->monitor.solution_cycles,"tps");

                // The following is incr_potential from both deformation and loads:
                for (m=1;m<=E->sphere.caps_per_proc;m++)  
                for (j=1;j<=E->lmesh.nsf;j++)
                  /* the original surface incr potential:
                    E->Xsurf[3][m][j] =  E->incr_potential[0][m][j]  
                                       + E->incr_potential[2][m][j] ; */
                    /* the compressible multiple layer version : */
                    E->Xsurf[3][m][j] =  E->slice_ve.all_potential[m][j*E->lmesh.noz]
                                       + E->slice_ve.all_load_potential[m][j*E->lmesh.noz];

                sphere_expansion_output(E,1,E->Xsurf[3],
                        E->sphere.sphc[3],E->sphere.sphs[3],
                        E->monitor.solution_cycles,"pttl");


                if (E->ve_data_cont.SLE) {
                    // get the oceanload:
                    for (m=1;m<=E->sphere.caps_per_proc;m++)  
                    for (j=1;j<=E->lmesh.nsf;j++)
                        E->Xsurf[3][m][j] =  E->slice_ve.static_oceanload[m][j]
                                           + E->slice_ve.dynamic_oceanload[m][j];
                    sphere_expansion_output(E,1,E->Xsurf[3],
                            E->sphere.sphc[0],E->sphere.sphs[0],
                            E->monitor.solution_cycles,"oceanload");
                }

                for (i=0;i<2;i++) {
                    fp = (i==0)? stderr : E->fp ;
                    fprintf(fp, "Wrote initial topo and potential files.\n");
                    fflush(fp);
                }
            }
        }
    }// end if not been_here

    been_here++;

    if ( ((ii % E->control.record_every) == 0)|| (ii == E->advection.max_timesteps))    {

        output_strain(E,ii);

        if (E->parallel.me_loc[3]==E->parallel.nprocz-1)  {
            sprintf(outfile,"%s.topo_s.%d.%d",
                    E->control.data_file,E->parallel.me,ii);
            //sprintf(outfile,"%s.topo_s.%d.%d", 
            //       E->control.data_file,E->parallel.me,ii);
            fp = fopen(outfile,"w");
 			fprintf(fp, "timestep, current year, dt, PW (x2), PW_incr (x2), eustatic_sea_level, barystatic_sea_level, surface net rotation (x3), CM_incr(x3), CM_incr_ice_static_ocean(x3) \n");
            fprintf(fp,"%05d %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n",
                    ii,                                  // timstep
                    E->ve_data_cont.tau*E->monitor.elapsed_time, // current time (years)
                    E->ve_data_cont.tau*E->advection.timestep,  // length of timestep
                    E->ve_data_cont.PW[0],E->ve_data_cont.PW[1], // cumu. polar motion
                    E->ve_data_cont.PW_incr[0],E->ve_data_cont.PW_incr[1], // polar motion
                    E->ve_data_cont.eustatic_sea_level, E->ve_data_cont.barystatic_sea_level, // RSL (m)
                    E->ve_data_cont.Omega_surface[0], E->ve_data_cont.Omega_surface[1], E->ve_data_cont.Omega_surface[2], // Surface Net rotation rate, deg/yr
					E->ve_data_cont.CM_incr[0],E->ve_data_cont.CM_incr[1],E->ve_data_cont.CM_incr[2],	// CM_incr
					E->ve_data_cont.CM_incr_ice_static_ocean[0], E->ve_data_cont.CM_incr_ice_static_ocean[1],E->ve_data_cont.CM_incr_ice_static_ocean[2]);

            // The following is incr_potential from both deformation and loads:
            for (m=1;m<=E->sphere.caps_per_proc;m++)  
            for (j=1;j<=E->lmesh.nsf;j++)
                E->Xsurf[3][m][j] =  E->slice_ve.all_potential[m][j*E->lmesh.noz]
                                   + E->slice_ve.all_load_potential[m][j*E->lmesh.noz];



            for (m=1;m<=E->sphere.caps_per_proc;m++)  
            for (j=1;j<=E->lmesh.nsf;j++)  { 
            /* the compressible multiple layer version */
                fprintf(fp,"%.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e\n",
                        E->slice_ve.surf[3][m][j],                             // total topography
                        E->sphere.cap[m].V[3][j*E->lmesh.noz],                             // incr topo 
                        E->slice_ve.init_total_potential[m][j*E->lmesh.noz],      // total potential
                        E->Xsurf[3][m][j],                                 // incr potential
                        E->sphere.cap[m].total_VS[1][j],        
                        E->sphere.cap[m].V[1][j*E->lmesh.noz],        
                        E->sphere.cap[m].total_VS[2][j],        
                        E->sphere.cap[m].V[2][j*E->lmesh.noz],
                        // the following are the strain tensor components:
                        E->StrainTensorNodeXX[m][j*E->lmesh.noz],
                        E->StrainTensorNodeYY[m][j*E->lmesh.noz],
                        E->StrainTensorNodeZZ[m][j*E->lmesh.noz],
                        E->StrainTensorNodeXY[m][j*E->lmesh.noz],
                        E->StrainTensorNodeXZ[m][j*E->lmesh.noz],
                        E->StrainTensorNodeZY[m][j*E->lmesh.noz]);
            }

      
            fclose(fp);


            for (m=1;m<=E->sphere.caps_per_proc;m++)
                     for (j=1;j<=E->lmesh.nsf;j++)  {
                          i=j*E->lmesh.noz;
                          E->Xsurf[1][m][j] = E->sphere.cap[m].V[3][i]; // incr topo
                          E->Xsurf[2][m][j] = E->slice_ve.init_total_potential[m][i]; // total potential
                     }

            sphere_expansion_output(E,1, E->slice_ve.surf[3],
                                         E->sphere.sphc[2], E->sphere.sphs[2],
                                         E->monitor.solution_cycles, "tps");

            sphere_expansion_output(E,1, E->Xsurf[1],  
                                         E->sphere.sphc[0],E->sphere.sphs[0],
                                         E->monitor.solution_cycles,"vtps");

            sphere_expansion_output(E,1, E->Xsurf[2],
                                         E->sphere.sphc[3],E->sphere.sphs[3],
                                         E->monitor.solution_cycles,"pttl");

            sphere_expansion_output(E,1, E->Xsurf[3],
                                         E->sphere.sphc[0],E->sphere.sphs[0],
                                         E->monitor.solution_cycles,"pttldot");


            for (m=1;m<=E->sphere.caps_per_proc;m++)
                     for (j=1;j<=E->lmesh.nsf;j++)  {
                          // i=j*E->lmesh.noz;
                          E->Xsurf[1][m][j] = E->slice_ve.iceload[0][m][j]; // iceload
                          E->Xsurf[2][m][j] = E->slice_ve.static_oceanload[m][j]; // static ocean load
                          E->Xsurf[3][m][j] = E->slice_ve.dynamic_oceanload[m][j]; // dynamic ocean load
                     }
            
            sphere_expansion_output(E,1, E->Xsurf[1],
                                         E->sphere.sphc[0],E->sphere.sphs[0],
                                         E->monitor.solution_cycles,"iceload");
            
            sphere_expansion_output(E,1, E->Xsurf[2],
                                          E->sphere.sphc[0],E->sphere.sphs[0],
                                          E->monitor.solution_cycles,"static_oceanload");
            
            sphere_expansion_output(E,1, E->Xsurf[3],
                                          E->sphere.sphc[0],E->sphere.sphs[0],
                                          E->monitor.solution_cycles,"dynamic_oceanload");


            for (m=1;m<=E->sphere.caps_per_proc;m++)
                  for (j=1;j<=E->lmesh.nsf;j++)  { 
                    i=j*E->lmesh.noz;
                    E->Xsurf[1][m][j] = E->sphere.cap[m].total_VS[1][j];
                    E->Xsurf[2][m][j] = E->sphere.cap[m].total_VS[2][j];
                  }

            sphere_expansion_output(E,1, E->Xsurf[1],
                              E->sphere.sphc[4],E->sphere.sphs[4],
                              E->monitor.solution_cycles,"utheta");
            sphere_expansion_output(E,1, E->Xsurf[2],
                              E->sphere.sphc[5],E->sphere.sphs[5],
                              E->monitor.solution_cycles,"ufi");

              }

        /*// write out CMB values:
          if (E->parallel.me_loc[3]==0)  {

              sprintf(outfile,"%s.topo_b.%d.%d",
                              E->control.data_file2,E->parallel.me,ii);
              //sprintf(outfile,"%s.topo_b.%d.%d",
              //                 E->control.data_file,E->parallel.me,ii);
              fp = fopen(outfile,"w");
              fprintf(fp,"%05d %.5e\n",ii,E->monitor.elapsed_time);
              for (m=1;m<=E->sphere.caps_per_proc;m++)  
                  for (j=1;j<=E->lmesh.nsf;j++)  { 
                      i=(j-1)*E->lmesh.noz+1;
                      E->Xsurf[1][m][j] = E->SX[E->mesh.levmax][1][m][i];
                      E->Xsurf[2][m][j] = E->SX[E->mesh.levmax][2][m][i];
                      E->Xsurf[3][m][j] = E->SX[E->mesh.levmax][3][m][i];
                      fprintf(fp,"%.4e %.4e %.4e %.4e %.4e\n",
                                E->Xsurf[1][m][j],E->Xsurf[2][m][j],
                                E->slice.botm[3][m][j],
                                E->U[m][E->id[m][i].doff[3]],E->potential[m][i]);
                  }
              fclose(fp);

              sphere_expansion_output(E,0,E->slice.botm[3],
                                          E->sphere.sphc[0],E->sphere.sphs[0],
                                          E->monitor.solution_cycles,"tpb");

              for (m=1;m<=E->sphere.caps_per_proc;m++)  
                  for (j=1;j<=E->lmesh.nsf;j++)  { 
                      i=(j-1)*E->lmesh.noz+1;
                      E->Xsurf[3][m][j] = E->U[m][E->id[m][i].doff[3]];
                  }

              sphere_expansion_output(E,0,E->Xsurf[3],
                                      E->sphere.sphc[0],E->sphere.sphs[0],
                                      E->monitor.solution_cycles,"vtpb");
          } 
        */

    } // end if on a recording timestep

    return;
}


/*========================================================= */
/*========================================================= */
void get_temperature_half_space(struct All_variables *E)
 {

 double timea;

 timea = E->monitor.elapsed_time;

 return;
 }

/*========================================================= */
/*========================================================= */

  void add_viscoelasticity(E, evisc)
  struct All_variables *E;
  float **evisc;
  {

  double visc, alpha;
  const int vpts = vpoints[E->mesh.nsd];
  int m,i,j;

for (m=1;m<=E->sphere.caps_per_proc;m++)
  for (i=1;i<=E->lmesh.nel;i++)   {
    visc = 0.0;
    for(j=1;j<=vpts;j++) 
       visc += evisc[m][(i-1)*vpts+j];
    visc = visc/vpts;

    alpha = E->Maxwelltime[m][i]/(2.0*visc)*E->advection.timestep;
    for(j=1;j<=vpts;j++) 
      evisc[m][(i-1)*vpts+j]=E->Maxwelltime[m][i]/(1.0+alpha);

    E->Maxwelltime[m][i] = (1.0-alpha)/(1.0+alpha);
    }


                           /* Maxwelltime is now alpha/(dt+alpha) */
                           /* evisc is now evisc/(dt+alpha) */

  return;
  }

  /*========================================================= */
/*========================================================= */
// add by tao: copy from geruo
// note: geruo define those at vpts, but here we only define at element
/*  add viscoelasticity for compressible rheology
    now the equation reads sigma_ij(n+1) =    Maxwell*sigma_ij(n)+esmu*sigma_kk(n)*delta_ij
                                           +  2*EVI*de_ij
                                           +  elambda*de_kk*delta_ij
*/

  void add_viscoelasticity_comp(struct All_variables *E)
  //struct All_variables *E;
  {

  double visc, alpha, lambda, smu;
  const int vpts = vpoints[E->mesh.nsd];
  int m,i,j;
  static int been = 0;
  FILE *fp;
//printf("add_viscoelasticity_comp is called \n");
  for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (i=1;i<=E->lmesh.nel;i++) {
      visc = 0.0;
      for (j=1;j<=vpts;j++)
  visc +=E->EVI[E->mesh.levmax][m][(i-1)*vpts+j];
      visc = visc/vpts;

      alpha = E->Maxwelltime[m][i]/(2.0*visc)* E->advection.timestep;
      for (j=1;j<=vpts;j++){
        E->EVI[E->mesh.levmax][m][(i-1)*vpts+j] = E->esmu[E->mesh.levmax][m][i]/(1.0+alpha);
      }
        lambda = E->elambda[E->mesh.levmax][m][i];
        smu = E->esmu[E->mesh.levmax][m][i];
        E->elambda[E->mesh.levmax][m][i] = (lambda+(lambda+(2.0/3.0)*smu)*alpha)/(1.0+alpha);
        E->esmu[E->mesh.levmax][m][i] = 2.0/3.0*alpha/(1.0+alpha);
//        E->esmu[E->mesh.levmax][m][(i-1)*vpts+j] = (lambda+(2.0/3.0)*smu)*alpha/(1.0+alpha);
     // }
      E->Maxwelltime[m][i] = (1.0-alpha)/(1.0+alpha);
    }

//  	if (E->control.DIRECT) {
//	  fp = fopen("/home/geruo/VE_CitcomS/test_comp/chk_par","w");
//	  for (m=1;m<=E->sphere.caps_per_proc;m++){
//	    for (i=1;i<=E->lmesh.elz;i++){
//	        j=1;
//	        fprintf(fp,"%d %d %.5e %.5e %.5e %.5e\n",i,j,E->EVI[E->mesh.levmax][m][(i-1)*vpts+j], E->elambda[E->mesh.levmax][m][(i-1)*vpts+j], E->esmu[E->mesh.levmax][m][(i-1)*vpts+j], E->Maxwelltime[m][i]);
//	        
//	    }
//  	}
//  	fclose(fp);
//  	fprintf(stderr,"wrote chk_par\n");
//  		been++;
//  		}
  return;
  }


// change end

void update_stress_strain(struct All_variables *E, int ii)
//    struct All_variables *E;
//    int ii;
{

    void get_rtf_at_vpts();
    void velo_from_element_d();
    void construct_c3x3matrix_el();
    int i,j,k,e,node,snode,m,nel2;
    
    double *SXX[NCS],*SYY[NCS],*SXY[NCS],*SXZ[NCS],*SZY[NCS],*SZZ[NCS];
    double VV[4][9],Szz,Sxx,Syy,Sxy,Sxz,Szy;
    double Vxyz1,Vxyz2,Vxyz3,Vxyz4,Vxyz5,Vxyz6;
    double Sxyz1,Sxyz2,Sxyz3,Sxyz4,Sxyz5,Sxyz6;
    double sinaa,cosaa,ct,el_volume,Visc,a,b,rtf[4][9];
    struct Shape_function GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;
    static struct CC Cc;
    static struct CCX Ccx;

    
    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
    const int vpts=vpoints[dims];
    const int ppts=ppoints[dims];
    const int ends=enodes[dims];
    const int nno=E->lmesh.nno;
    const int lev=E->mesh.levmax;
    const int sphere_key=1;

    double *gnxx,*gnda;
    double old_Skk;

  for(m=1;m<=E->sphere.caps_per_proc;m++)      {

    for(e=1;e<=E->lmesh.nel;e++)  {

      E->StrainTensorXX[m][e] = E->StrainTensorYY[m][e] = E->StrainTensorZZ[m][e] = 0.0;
      E->StrainTensorXY[m][e] = E->StrainTensorXZ[m][e] = E->StrainTensorZY[m][e] = 0.0;

      gnxx = E->gNX[m][e].vpt;

      get_rtf_at_vpts(E,m,E->mesh.levmax,e,rtf);

      if ((e-1)%E->lmesh.elz==0)
        construct_c3x3matrix_el(E,e,&Cc,&Ccx,E->mesh.levmax,m,0);

      velo_from_element_d(E,VV,m,e,sphere_key);
  
      for(j=1;j<=vpts;j++)   {
        Sxyz1 = Sxyz2 = Sxyz3 = Sxyz4 = Sxyz5 = Sxyz6 = 0.0;
        sinaa = sin(rtf[1][j]);
        cosaa = cos(rtf[1][j]);
        ct = cosaa/sinaa;
        Visc = E->EVi[m][(e-1)*vpts+j];
        for(i=1;i<=ends;i++)   {
          for(k=1;k<=dims;k++)   {
          
          Sxyz1 += VV[k][i]*rtf[3][j]*
                 (gnxx[GNVXINDEX(0,i,j)]*Cc.vpt[BVINDEX(1,k,i,j)]
                 +E->N.vpt[GNVINDEX(i,j)]*Ccx.vpt[BVXINDEX(1,k,1,i,j)]
                 +E->N.vpt[GNVINDEX(i,j)]*Cc.vpt[BVINDEX(3,k,i,j)]);

          Sxyz2 += VV[k][i]*rtf[3][j]*
                (E->N.vpt[GNVINDEX(i,j)]*Cc.vpt[BVINDEX(1,k,i,j)]*ct
                +E->N.vpt[GNVINDEX(i,j)]*Cc.vpt[BVINDEX(3,k,i,j)]
                +(gnxx[GNVXINDEX(1,i,j)]*Cc.vpt[BVINDEX(2,k,i,j)]
                 +E->N.vpt[GNVINDEX(i,j)]*Ccx.vpt[BVXINDEX(2,k,2,i,j)])/sinaa);

          Sxyz3 += VV[k][i]*gnxx[GNVXINDEX(2,i,j)]*Cc.vpt[BVINDEX(3,k,i,j)];

          Sxyz4 += VV[k][i]*rtf[3][j]*
                (gnxx[GNVXINDEX(0,i,j)]*Cc.vpt[BVINDEX(2,k,i,j)]
                +E->N.vpt[GNVINDEX(i,j)]*Ccx.vpt[BVXINDEX(2,k,1,i,j)]
                -ct*Cc.vpt[BVINDEX(2,k,i,j)]*E->N.vpt[GNVINDEX(i,j)]
                +(E->N.vpt[GNVINDEX(i,j)]*Ccx.vpt[BVXINDEX(1,k,2,i,j)]
                 +gnxx[GNVXINDEX(1,i,j)]*Cc.vpt[BVINDEX(1,k,i,j)])/sinaa);

          Sxyz5 += VV[k][i]*rtf[3][j]*
                (gnxx[GNVXINDEX(2,i,j)]*Cc.vpt[BVINDEX(1,k,i,j)]/rtf[3][j]
                +(E->N.vpt[GNVINDEX(i,j)]*Ccx.vpt[BVXINDEX(3,k,1,i,j)]
                 +gnxx[GNVXINDEX(0,i,j)]*Cc.vpt[BVINDEX(3,k,i,j)]
                 -E->N.vpt[GNVINDEX(i,j)]*Cc.vpt[BVINDEX(1,k,i,j)]));

          Sxyz6 += VV[k][i]*
                (gnxx[GNVXINDEX(2,i,j)]*Cc.vpt[BVINDEX(2,k,i,j)]
                -rtf[3][j]*E->N.vpt[GNVINDEX(i,j)]*Cc.vpt[BVINDEX(2,k,i,j)]
                +rtf[3][j]/sinaa*(E->N.vpt[GNVINDEX(i,j)]*Ccx.vpt[BVXINDEX(3,k,2,i,j)]+gnxx[GNVXINDEX(1,i,j)]*Cc.vpt[BVINDEX(3,k,i,j)]));
          }
        }

        E->StrainTensorXX[m][e] += Sxyz1;
        E->StrainTensorYY[m][e] += Sxyz2;
        E->StrainTensorZZ[m][e] += Sxyz3;
        E->StrainTensorXY[m][e] += Sxyz4;
        E->StrainTensorXZ[m][e] += Sxyz5;
        E->StrainTensorZY[m][e] += Sxyz6;

        //>>>> add by tao:
        if(E->ve_data_cont.compressible){
              old_Skk = E->Sxx[m][(e-1)*vpts+j] + E->Syy[m][(e-1)*vpts+j] + E->Szz[m][(e-1)*vpts+j];
                E->Sxx[m][(e-1)*vpts+j] = 2.0*Visc*Sxyz1 + E->elambda[lev][m][e]*(Sxyz1+Sxyz2+Sxyz3)
                                        + E->Maxwelltime[m][e]*E->Sxx[m][(e-1)*vpts+j] + E->esmu[lev][m][e]*old_Skk;
                E->Syy[m][(e-1)*vpts+j] = 2.0*Visc*Sxyz2 + E->elambda[lev][m][e]*(Sxyz1+Sxyz2+Sxyz3)
                                        + E->Maxwelltime[m][e]*E->Syy[m][(e-1)*vpts+j] + E->esmu[lev][m][e]*old_Skk;
                E->Szz[m][(e-1)*vpts+j] = 2.0*Visc*Sxyz3 + E->elambda[lev][m][e]*(Sxyz1+Sxyz2+Sxyz3)
                                        + E->Maxwelltime[m][e]*E->Szz[m][(e-1)*vpts+j] + E->esmu[lev][m][e]*old_Skk;
                E->Sxy[m][(e-1)*vpts+j] =     Visc*Sxyz4 + E->Maxwelltime[m][e]*E->Sxy[m][(e-1)*vpts+j];
                E->Sxz[m][(e-1)*vpts+j] =     Visc*Sxyz5 + E->Maxwelltime[m][e]*E->Sxz[m][(e-1)*vpts+j];
                E->Szy[m][(e-1)*vpts+j] =     Visc*Sxyz6 + E->Maxwelltime[m][e]*E->Szy[m][(e-1)*vpts+j];

                E->Ekk[m][(e-1)*vpts+j] = E->Ekk[m][(e-1)*vpts+j] + Sxyz1+Sxyz2+Sxyz3;
        }
        else{
                E->Sxx[m][(e-1)*vpts+j]=2.0*Visc*Sxyz1 + E->Maxwelltime[m][e]*E->Sxx[m][(e-1)*vpts+j];
                E->Syy[m][(e-1)*vpts+j]=2.0*Visc*Sxyz2 + E->Maxwelltime[m][e]*E->Syy[m][(e-1)*vpts+j];
                E->Szz[m][(e-1)*vpts+j]=2.0*Visc*Sxyz3 + E->Maxwelltime[m][e]*E->Szz[m][(e-1)*vpts+j];
                E->Sxy[m][(e-1)*vpts+j]=    Visc*Sxyz4 + E->Maxwelltime[m][e]*E->Sxy[m][(e-1)*vpts+j];
                E->Sxz[m][(e-1)*vpts+j]=    Visc*Sxyz5 + E->Maxwelltime[m][e]*E->Sxz[m][(e-1)*vpts+j];
                E->Szy[m][(e-1)*vpts+j]=    Visc*Sxyz6 + E->Maxwelltime[m][e]*E->Szy[m][(e-1)*vpts+j];
        }
        //<<<< end tao

      } // end for j (vpts)

      E->StrainTensorXX[m][e] = E->StrainTensorXX[m][e]/vpts;
      E->StrainTensorYY[m][e] = E->StrainTensorYY[m][e]/vpts;
      E->StrainTensorZZ[m][e] = E->StrainTensorZZ[m][e]/vpts;
      E->StrainTensorXY[m][e] = E->StrainTensorXY[m][e]/(2.0*vpts);  // 
      E->StrainTensorXZ[m][e] = E->StrainTensorXZ[m][e]/(2.0*vpts);
      E->StrainTensorZY[m][e] = E->StrainTensorZY[m][e]/(2.0*vpts);


      E->S2xx[m][e] = E->S2yy[m][e] = E->S2zz[m][e] = E->S2xy[m][e] = E->S2xz[m][e] = E->S2zy[m][e] = 0.0;
      for(j=1;j<=vpts;j++)   {
         E->S2xx[m][e] += E->Sxx[m][(e-1)*vpts+j];
         E->S2yy[m][e] += E->Syy[m][(e-1)*vpts+j];
         E->S2zz[m][e] += E->Szz[m][(e-1)*vpts+j];
         E->S2xy[m][e] += E->Sxy[m][(e-1)*vpts+j];
         E->S2xz[m][e] += E->Sxz[m][(e-1)*vpts+j];
         E->S2zy[m][e] += E->Szy[m][(e-1)*vpts+j];
         }
      E->S2xx[m][e] = E->S2xx[m][e]/vpts;
      E->S2yy[m][e] = E->S2yy[m][e]/vpts;
      E->S2zz[m][e] = E->S2zz[m][e]/vpts;
      E->S2xy[m][e] = E->S2xy[m][e]/vpts;
      E->S2xz[m][e] = E->S2xz[m][e]/vpts;
      E->S2zy[m][e] = E->S2zy[m][e]/vpts;

      if (!E->viscosity.SDEPV)  {
      // compute the 2nd invariant of stress
      E->S2inv[m][e]=E->S2xx[m][e]*E->S2xx[m][e]+E->S2xy[m][e]*E->S2xy[m][e]*2.0
                    +E->S2yy[m][e]*E->S2yy[m][e]+E->S2zy[m][e]*E->S2zy[m][e]*2.0
                    +E->S2zz[m][e]*E->S2zz[m][e]+E->S2xz[m][e]*E->S2xz[m][e]*2.0;
      E->S2inv[m][e] = sqrt(0.5*E->S2inv[m][e]);   
      }
    }    /* end for el */

    // project E->StrainTensor to E->StrainTensorNode
    void ele_to_nodes();
    ele_to_nodes(E,E->StrainTensorXX,E->StrainTensorNodeXX,lev);
    ele_to_nodes(E,E->StrainTensorYY,E->StrainTensorNodeYY,lev);
    ele_to_nodes(E,E->StrainTensorZZ,E->StrainTensorNodeZZ,lev);
    ele_to_nodes(E,E->StrainTensorXY,E->StrainTensorNodeXY,lev);
    ele_to_nodes(E,E->StrainTensorXZ,E->StrainTensorNodeXZ,lev);
    ele_to_nodes(E,E->StrainTensorZY,E->StrainTensorNodeZY,lev);

  }     /* end for m */

    return; 
}

void output_strain(struct All_variables *E, int ii)
{
    void get_rtf_at_vpts();
    void velo_from_element_d();
    void construct_c3x3matrix_el();
    void sphere_expansion_output();
    int i,j,k,e,node,snode,m,n,nel2;
    
    double *SXX[NCS],*SYY[NCS],*SXY[NCS],*SXZ[NCS],*SZY[NCS],*SZZ[NCS];
    double VV[4][9],Szz,Sxx,Syy,Sxy,Sxz,Szy;
    double Vxyz1,Vxyz2,Vxyz3,Vxyz4,Vxyz5,Vxyz6;

    // temp value for each element
    double Sxyz1,Sxyz2,Sxyz3,Sxyz4,Sxyz5,Sxyz6;
    double S2Dxy1, S2Dxy2, S2Dxy3; // for 2D strain
    double S2Dxy1_V2, S2Dxy2_V2, S2Dxy3_V2; // use only surface velocity
    double div, div2D; // for divergence
    double vor; // for vorticity
    double vor_V2; // use only surface velocity
    double div2D_V2; // use only surface velocity

    double sinaa,cosaa,ct,el_volume,Visc,a,b,rtf[4][9];
    struct Shape_function GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;
    static struct CC Cc;
    static struct CCX Ccx;

    // define variables for output
    static double *Sxx_out[NCS],*Syy_out[NCS],*Szz_out[NCS],*Sxy_out[NCS],*Sxz_out[NCS],*Szy_out[NCS]; // 3D strain

    static double *S2Dxx_out[NCS],*S2Dyy_out[NCS] ,*S2Dxy_out[NCS] ; // 2D strain in theta phi (surface)
    static double *S2Dxx_out_V2[NCS],*S2Dyy_out_V2[NCS] ,*S2Dxy_out_V2[NCS] ; // use only surface velocity
    static double *div_out[NCS]; // divergence
    static double *div2D_out[NCS]; // divergence in theta phi (surface)
    static double *div2D_out_V2[NCS]; // use only surface velocity
    static double *vor_out[NCS]; // vorticity
    static double *vor_out_V2[NCS]; // use only surface velocity
    static double *surf_var[NCS]; // surface variable

    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
    const int vpts=vpoints[dims];
    const int ppts=ppoints[dims];
    const int ends=enodes[dims];
    const int nno=E->lmesh.nno;
    const int lev=E->mesh.levmax;
    const int sphere_key=1;

    double *gnxx,*gnda;
    double old_Skk;
    char outfile[255];
    FILE *fp;


    static int been_here=0;

    if(been_here == 0){
      for (j=1;j<=E->sphere.caps_per_proc;j++)  {
        Sxx_out[j] = (double *)malloc((E->lmesh.nno+2)*sizeof(double));
        Syy_out[j] = (double *)malloc((E->lmesh.nno+2)*sizeof(double));
        Szz_out[j] = (double *)malloc((E->lmesh.nno+2)*sizeof(double));
        Sxy_out[j] = (double *)malloc((E->lmesh.nno+2)*sizeof(double));
        Sxz_out[j] = (double *)malloc((E->lmesh.nno+2)*sizeof(double));
        Szy_out[j] = (double *)malloc((E->lmesh.nno+2)*sizeof(double));
        S2Dxx_out[j] = (double *)malloc((E->lmesh.nno+2)*sizeof(double));
        S2Dyy_out[j] = (double *)malloc((E->lmesh.nno+2)*sizeof(double));
        S2Dxy_out[j] = (double *)malloc((E->lmesh.nno+2)*sizeof(double));
        S2Dxx_out_V2[j] = (double *)malloc((E->lmesh.nno+2)*sizeof(double));
        S2Dyy_out_V2[j] = (double *)malloc((E->lmesh.nno+2)*sizeof(double));
        S2Dxy_out_V2[j] = (double *)malloc((E->lmesh.nno+2)*sizeof(double));
        div_out[j] = (double *)malloc((E->lmesh.nno+2)*sizeof(double));
        div2D_out[j] = (double *)malloc((E->lmesh.nno+2)*sizeof(double));
        div2D_out_V2[j] = (double *)malloc((E->lmesh.nno+2)*sizeof(double));
        vor_out[j] = (double *)malloc((E->lmesh.nno+2)*sizeof(double));
        vor_out_V2[j] = (double *)malloc((E->lmesh.nno+2)*sizeof(double));
        surf_var[j] = (double *)malloc((E->lmesh.nsf+2)*sizeof(double));

      }
      been_here = 1;
    }
    
    for(m=1;m<=E->sphere.caps_per_proc;m++)      {
      for(n=1;n<=E->lmesh.NNO[lev];n++){
      Sxx_out[m][n] = 0.0; Syy_out[m][n] = 0.0; Szz_out[m][n] = 0.0;
      Sxy_out[m][n] = 0.0; Sxz_out[m][n] = 0.0; Szy_out[m][n] = 0.0;
      S2Dxx_out[m][n] = 0.0; S2Dyy_out[m][n] = 0.0; S2Dxy_out[m][n] = 0.0;
      div_out[m][n] = 0.0; div2D_out[m][n] = 0.0; vor_out[m][n] = 0.0;
      S2Dxx_out_V2[m][n] = 0.0; S2Dyy_out_V2[m][n] = 0.0; S2Dxy_out_V2[m][n] = 0.0;
      div2D_out_V2[m][n] = 0.0; vor_out_V2[m][n] = 0.0;
      }
    }

  for(m=1;m<=E->sphere.caps_per_proc;m++)      {

    for(e=1;e<=E->lmesh.nel;e++)  {

      if(e%E->lmesh.elz!=0){
        continue;  // only calculate top layer
      }

      gnxx = E->gNX[m][e].vpt;

      get_rtf_at_vpts(E,m,E->mesh.levmax,e,rtf);

      // if ((e-1)%E->lmesh.elz==0)
      construct_c3x3matrix_el(E,e,&Cc,&Ccx,E->mesh.levmax,m,0);

      velo_from_element_d(E,VV,m,e,sphere_key);
  
      Sxyz1 = Sxyz2 = Sxyz3 = Sxyz4 = Sxyz5 = Sxyz6 = 0.0;  // keep record of elemental value
      S2Dxy1 = S2Dxy2 = S2Dxy3 = 0.0;
      div = div2D = vor = 0.0;

      for(j=1;j<=vpts;j++)   {
        sinaa = sin(rtf[1][j]);
        cosaa = cos(rtf[1][j]);
        ct = cosaa/sinaa;
        Visc = E->EVi[m][(e-1)*vpts+j];
        for(i=1;i<=ends;i++)   {
          for(k=1;k<=dims;k++)   {
          
          Sxyz1 += VV[k][i]*rtf[3][j]*
                 (gnxx[GNVXINDEX(0,i,j)]*Cc.vpt[BVINDEX(1,k,i,j)]
                 +E->N.vpt[GNVINDEX(i,j)]*Ccx.vpt[BVXINDEX(1,k,1,i,j)]
                 +E->N.vpt[GNVINDEX(i,j)]*Cc.vpt[BVINDEX(3,k,i,j)]);

          S2Dxy1 += VV[k][i]*rtf[3][j]*
                 (gnxx[GNVXINDEX(0,i,j)]*Cc.vpt[BVINDEX(1,k,i,j)]
                 +E->N.vpt[GNVINDEX(i,j)]*Ccx.vpt[BVXINDEX(1,k,1,i,j)]
                 ); 
                 // epsion_theta_theta for (u_theta,u_phi)
                 // the Sxyz1 is epsion_theta_theta for (u_r,u_theta,u_phi)
                 
          Sxyz2 += VV[k][i]*rtf[3][j]*
                (E->N.vpt[GNVINDEX(i,j)]*Cc.vpt[BVINDEX(1,k,i,j)]*ct
                +E->N.vpt[GNVINDEX(i,j)]*Cc.vpt[BVINDEX(3,k,i,j)]
                +(gnxx[GNVXINDEX(1,i,j)]*Cc.vpt[BVINDEX(2,k,i,j)]
                 +E->N.vpt[GNVINDEX(i,j)]*Ccx.vpt[BVXINDEX(2,k,2,i,j)])/sinaa);
          
          S2Dxy2 += VV[k][i]*rtf[3][j]*
                (E->N.vpt[GNVINDEX(i,j)]*Cc.vpt[BVINDEX(1,k,i,j)]*ct
               // +E->N.vpt[GNVINDEX(i,j)]*Cc.vpt[BVINDEX(3,k,i,j)]
                +(gnxx[GNVXINDEX(1,i,j)]*Cc.vpt[BVINDEX(2,k,i,j)]
                 +E->N.vpt[GNVINDEX(i,j)]*Ccx.vpt[BVXINDEX(2,k,2,i,j)])/sinaa);
                 // epsion_phi_phi for (u_theta,u_phi)

          Sxyz3 += VV[k][i]*gnxx[GNVXINDEX(2,i,j)]*Cc.vpt[BVINDEX(3,k,i,j)];

          Sxyz4 += VV[k][i]*rtf[3][j]*
                (gnxx[GNVXINDEX(0,i,j)]*Cc.vpt[BVINDEX(2,k,i,j)]
                +E->N.vpt[GNVINDEX(i,j)]*Ccx.vpt[BVXINDEX(2,k,1,i,j)]
                -ct*Cc.vpt[BVINDEX(2,k,i,j)]*E->N.vpt[GNVINDEX(i,j)]
                +(E->N.vpt[GNVINDEX(i,j)]*Ccx.vpt[BVXINDEX(1,k,2,i,j)]
                 +gnxx[GNVXINDEX(1,i,j)]*Cc.vpt[BVINDEX(1,k,i,j)])/sinaa);
          
          S2Dxy3 = Sxyz4; // epsion_theta_phi for (u_theta,u_phi)

          vor += VV[k][i]*rtf[3][j]*
                (gnxx[GNVXINDEX(0,i,j)]*Cc.vpt[BVINDEX(2,k,i,j)]
                +E->N.vpt[GNVINDEX(i,j)]*Ccx.vpt[BVXINDEX(2,k,1,i,j)]    // first 2 terms: 1/r* d(V_phi)/dtheta
                +ct*Cc.vpt[BVINDEX(2,k,i,j)]*E->N.vpt[GNVINDEX(i,j)]      // third term: cot(theta)*V_phi/r
                -(E->N.vpt[GNVINDEX(i,j)]*Ccx.vpt[BVXINDEX(1,k,2,i,j)]    
                 +gnxx[GNVXINDEX(1,i,j)]*Cc.vpt[BVINDEX(1,k,i,j)])/sinaa);  // last 2 terms: [1/r*sin(theta)]*d(V_theta)/dphi

          Sxyz5 += VV[k][i]*rtf[3][j]*
                (gnxx[GNVXINDEX(2,i,j)]*Cc.vpt[BVINDEX(1,k,i,j)]/rtf[3][j]
                +(E->N.vpt[GNVINDEX(i,j)]*Ccx.vpt[BVXINDEX(3,k,1,i,j)]
                 +gnxx[GNVXINDEX(0,i,j)]*Cc.vpt[BVINDEX(3,k,i,j)]
                 -E->N.vpt[GNVINDEX(i,j)]*Cc.vpt[BVINDEX(1,k,i,j)]));

          Sxyz6 += VV[k][i]*
                (gnxx[GNVXINDEX(2,i,j)]*Cc.vpt[BVINDEX(2,k,i,j)]
                -rtf[3][j]*E->N.vpt[GNVINDEX(i,j)]*Cc.vpt[BVINDEX(2,k,i,j)]
                +rtf[3][j]/sinaa*(E->N.vpt[GNVINDEX(i,j)]*Ccx.vpt[BVXINDEX(3,k,2,i,j)]+gnxx[GNVXINDEX(1,i,j)]*Cc.vpt[BVINDEX(3,k,i,j)]));
          }
        }
      } // end for j (vpts)
       
       Sxyz1 /= vpts;
       Sxyz2 /= vpts;
       Sxyz3 /= vpts;
       Sxyz4 /= 2.0*vpts;
       Sxyz5 /= 2.0*vpts;
       Sxyz6 /= 2.0*vpts;
 
       S2Dxy1 /= vpts;
       S2Dxy2 /= vpts;
       S2Dxy3 /= 2.0*vpts;

       vor /= vpts;
       
       div2D = S2Dxy1 + S2Dxy2;
       div = Sxyz1 + Sxyz2 + Sxyz3;

      /////////////////
      //   re-calculate for surface strain (V2)
      /////////////////
      if(1){
        // overwrite those values
        S2Dxy1_V2 = S2Dxy2_V2 = S2Dxy3_V2 = 0.0;
        div2D_V2 = vor_V2 = 0.0;
        
        assert(onedvpoints[dims]==4);
        for(j=1;j<=onedvpoints[dims];j++){  // ends 
          for(i=1; i<=dims; i++){
            VV[i][j] = VV[i][j+4];  // set the bottom nodes's velocity to be the same as the top nodes
          }
        }

        for(j=1;j<=vpts;j++)   {
          sinaa = sin(rtf[1][j]);
          cosaa = cos(rtf[1][j]);
          ct = cosaa/sinaa;
          for(i=1;i<=ends;i++)   {
            for(k=1;k<=dims;k++)   {

            S2Dxy1_V2 += VV[k][i]*rtf[3][j]*
                  (gnxx[GNVXINDEX(0,i,j)]*Cc.vpt[BVINDEX(1,k,i,j)]
                  +E->N.vpt[GNVINDEX(i,j)]*Ccx.vpt[BVXINDEX(1,k,1,i,j)]
                  ); 
                  // epsion_theta_theta for (u_theta,u_phi)
                  // the Sxyz1 is epsion_theta_theta for (u_r,u_theta,u_phi)
            
            
            S2Dxy2_V2 += VV[k][i]*rtf[3][j]*
                  (E->N.vpt[GNVINDEX(i,j)]*Cc.vpt[BVINDEX(1,k,i,j)]*ct
                // +E->N.vpt[GNVINDEX(i,j)]*Cc.vpt[BVINDEX(3,k,i,j)]
                  +(gnxx[GNVXINDEX(1,i,j)]*Cc.vpt[BVINDEX(2,k,i,j)]
                  +E->N.vpt[GNVINDEX(i,j)]*Ccx.vpt[BVXINDEX(2,k,2,i,j)])/sinaa);
                  // epsion_phi_phi for (u_theta,u_phi)

            S2Dxy3_V2 += VV[k][i]*rtf[3][j]*
                  (gnxx[GNVXINDEX(0,i,j)]*Cc.vpt[BVINDEX(2,k,i,j)]
                  +E->N.vpt[GNVINDEX(i,j)]*Ccx.vpt[BVXINDEX(2,k,1,i,j)]
                  -ct*Cc.vpt[BVINDEX(2,k,i,j)]*E->N.vpt[GNVINDEX(i,j)]
                  +(E->N.vpt[GNVINDEX(i,j)]*Ccx.vpt[BVXINDEX(1,k,2,i,j)]
                  +gnxx[GNVXINDEX(1,i,j)]*Cc.vpt[BVINDEX(1,k,i,j)])/sinaa);

            vor_V2 += VV[k][i]*rtf[3][j]*
                  (gnxx[GNVXINDEX(0,i,j)]*Cc.vpt[BVINDEX(2,k,i,j)]
                  +E->N.vpt[GNVINDEX(i,j)]*Ccx.vpt[BVXINDEX(2,k,1,i,j)]    // first 2 terms: 1/r* d(V_phi)/dtheta
                  +ct*Cc.vpt[BVINDEX(2,k,i,j)]*E->N.vpt[GNVINDEX(i,j)]      // third term: cot(theta)*V_phi/r
                  -(E->N.vpt[GNVINDEX(i,j)]*Ccx.vpt[BVXINDEX(1,k,2,i,j)]    
                  +gnxx[GNVXINDEX(1,i,j)]*Cc.vpt[BVINDEX(1,k,i,j)])/sinaa);  // last 2 terms: [1/r*sin(theta)]*d(V_theta)/dphi
            }
          }
        } // end for j (vpts)
      }
 
       S2Dxy1_V2 /= vpts;
       S2Dxy2_V2 /= vpts;
       S2Dxy3_V2 /= 2.0*vpts;

       vor_V2 /= vpts;
       
       div2D_V2 = S2Dxy1_V2 + S2Dxy2_V2;

       

      // from element to nodes, but don't forget to multiply by the mass matrix
      for(j=1;j<=ends;j++) {
        n = E->IEN[lev][m][e].node[j];
        Sxx_out[m][n] += E->TWW[lev][m][e].node[j] * Sxyz1;
        Syy_out[m][n] += E->TWW[lev][m][e].node[j] * Sxyz2;
        Szz_out[m][n] += E->TWW[lev][m][e].node[j] * Sxyz3;
        Sxy_out[m][n] += E->TWW[lev][m][e].node[j] * Sxyz4;
        Sxz_out[m][n] += E->TWW[lev][m][e].node[j] * Sxyz5;
        Szy_out[m][n] += E->TWW[lev][m][e].node[j] * Sxyz6;

        S2Dxx_out[m][n] += E->TWW[lev][m][e].node[j] * S2Dxy1;
        S2Dyy_out[m][n] += E->TWW[lev][m][e].node[j] * S2Dxy2;
        S2Dxy_out[m][n] += E->TWW[lev][m][e].node[j] * S2Dxy3;

        S2Dxx_out_V2[m][n] += E->TWW[lev][m][e].node[j] * S2Dxy1_V2;
        S2Dyy_out_V2[m][n] += E->TWW[lev][m][e].node[j] * S2Dxy2_V2;
        S2Dxy_out_V2[m][n] += E->TWW[lev][m][e].node[j] * S2Dxy3_V2;

        div_out[m][n] += E->TWW[lev][m][e].node[j] * div;
        div2D_out[m][n] += E->TWW[lev][m][e].node[j] * div2D;
        vor_out[m][n] += E->TWW[lev][m][e].node[j] * vor;

        div2D_out_V2[m][n] += E->TWW[lev][m][e].node[j] * div2D_V2;
        vor_out_V2[m][n] += E->TWW[lev][m][e].node[j] * vor_V2;
      }

    }    /* end for el */

    (E->exchange_node_d)(E,Sxx_out,lev);
    (E->exchange_node_d)(E,Syy_out,lev);
    (E->exchange_node_d)(E,Szz_out,lev);
    (E->exchange_node_d)(E,Sxy_out,lev);
    (E->exchange_node_d)(E,Sxz_out,lev);
    (E->exchange_node_d)(E,Szy_out,lev);
    (E->exchange_node_d)(E,S2Dxx_out,lev);
    (E->exchange_node_d)(E,S2Dyy_out,lev);
    (E->exchange_node_d)(E,S2Dxy_out,lev);
    (E->exchange_node_d)(E,div_out,lev);
    (E->exchange_node_d)(E,div2D_out,lev);
    (E->exchange_node_d)(E,vor_out,lev);
    (E->exchange_node_d)(E,S2Dxx_out_V2,lev);
    (E->exchange_node_d)(E,S2Dyy_out_V2,lev);
    (E->exchange_node_d)(E,S2Dxy_out_V2,lev);
    (E->exchange_node_d)(E,div2D_out_V2,lev);
    (E->exchange_node_d)(E,vor_out_V2,lev);
  }     /* end for m */

  for(m=1;m<=E->sphere.caps_per_proc;m++)
    for(n=1;n<=E->lmesh.NNO[lev];n++){
      Sxx_out[m][n] *= E->MASS[lev][m][n];
      Syy_out[m][n] *= E->MASS[lev][m][n];
      Szz_out[m][n] *= E->MASS[lev][m][n];
      Sxy_out[m][n] *= E->MASS[lev][m][n];
      Sxz_out[m][n] *= E->MASS[lev][m][n];
      Szy_out[m][n] *= E->MASS[lev][m][n];
      S2Dxx_out[m][n] *= E->MASS[lev][m][n];
      S2Dyy_out[m][n] *= E->MASS[lev][m][n];
      S2Dxy_out[m][n] *= E->MASS[lev][m][n];
      div_out[m][n] *= E->MASS[lev][m][n];
      div2D_out[m][n] *= E->MASS[lev][m][n];
      vor_out[m][n] *= E->MASS[lev][m][n];
      S2Dxx_out_V2[m][n] *= E->MASS[lev][m][n];
      S2Dyy_out_V2[m][n] *= E->MASS[lev][m][n];
      S2Dxy_out_V2[m][n] *= E->MASS[lev][m][n];
      div2D_out_V2[m][n] *= E->MASS[lev][m][n];
      vor_out_V2[m][n] *= E->MASS[lev][m][n];
    }

  // save Sxx Syy Szz Sxy Sxz Szy to file
    if ( ((ii % E->control.record_every) == 0) || (ii == E->advection.max_timesteps))    {

      if (E->parallel.me_loc[3]==E->parallel.nprocz-1)  {
          sprintf(outfile,"%s.surf_strain.%d.%d",
                  E->control.data_file,E->parallel.me,ii);
          fp = fopen(outfile,"w");
          fprintf(fp,"%05d %.5e %.5e\n",
                  ii,                                  // timstep
                  E->ve_data_cont.tau*E->monitor.elapsed_time, // current time (years)
                  E->ve_data_cont.tau*E->advection.timestep  // length of timestep
          );
          for (m=1;m<=E->sphere.caps_per_proc;m++)  
          for (j=1;j<=E->lmesh.nsf;j++)  { 
              fprintf(fp,"%.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e\n",
                  //first: 6 components of strain tensor
                      Sxx_out[m][j*E->lmesh.noz],
                      Syy_out[m][j*E->lmesh.noz],
                      Szz_out[m][j*E->lmesh.noz],
                      Sxy_out[m][j*E->lmesh.noz],
                      Sxz_out[m][j*E->lmesh.noz],
                      Szy_out[m][j*E->lmesh.noz],
                  //second: 3 components of 2D strain tensor
                      S2Dxx_out[m][j*E->lmesh.noz],
                      S2Dyy_out[m][j*E->lmesh.noz],
                      S2Dxy_out[m][j*E->lmesh.noz],
                  //third: divergence
                      div_out[m][j*E->lmesh.noz],
                      div2D_out[m][j*E->lmesh.noz],
                  //fourth: vorticity
                      vor_out[m][j*E->lmesh.noz]


              );
          }
          fclose(fp);
      }

      if (E->parallel.me_loc[3]==E->parallel.nprocz-1)  {
          sprintf(outfile,"%s.surf_strain_V2.%d.%d",
                  E->control.data_file,E->parallel.me,ii);
          fp = fopen(outfile,"w");
          fprintf(fp,"%05d %.5e %.5e\n",
                  ii,                                  // timstep
                  E->ve_data_cont.tau*E->monitor.elapsed_time, // current time (years)
                  E->ve_data_cont.tau*E->advection.timestep  // length of timestep
          );
          for (m=1;m<=E->sphere.caps_per_proc;m++)  
          for (j=1;j<=E->lmesh.nsf;j++)  { 
              fprintf(fp,"%.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e\n",
                  //first: 6 components of strain tensor
                      Sxx_out[m][j*E->lmesh.noz],
                      Syy_out[m][j*E->lmesh.noz],
                      Szz_out[m][j*E->lmesh.noz],
                      Sxy_out[m][j*E->lmesh.noz],
                      Sxz_out[m][j*E->lmesh.noz],
                      Szy_out[m][j*E->lmesh.noz],
                  //second: 3 components of 2D strain tensor
                      S2Dxx_out_V2[m][j*E->lmesh.noz],
                      S2Dyy_out_V2[m][j*E->lmesh.noz],
                      S2Dxy_out_V2[m][j*E->lmesh.noz],
                  //third: divergence
                      div_out[m][j*E->lmesh.noz],
                      div2D_out_V2[m][j*E->lmesh.noz],
                  //fourth: vorticity
                      vor_out_V2[m][j*E->lmesh.noz]

                      
              );
          }
          fclose(fp);
      }

      for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)  {
              i=j*E->lmesh.noz;
              surf_var[m][j] =  div2D_out[m][i];
        }

      sphere_expansion_output(E,1, surf_var,  
                              E->sphere.sphc[0],E->sphere.sphs[0],
                              E->monitor.solution_cycles,"div2D");

      for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)  {
              i=j*E->lmesh.noz;
              surf_var[m][j] =  vor_out[m][i];
        }

      sphere_expansion_output(E,1, surf_var,  
                              E->sphere.sphc[0],E->sphere.sphs[0],
                              E->monitor.solution_cycles,"vor");



      // V2
      for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)  {
              i=j*E->lmesh.noz;
              surf_var[m][j] =  div2D_out_V2[m][i];
        }


      sphere_expansion_output(E,1, surf_var,  
                              E->sphere.sphc[0],E->sphere.sphs[0],
                              E->monitor.solution_cycles,"div2D_V2");

      for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)  {
              i=j*E->lmesh.noz;
              surf_var[m][j] =  vor_out_V2[m][i];
        }

      sphere_expansion_output(E,1, surf_var,  
                              E->sphere.sphc[0],E->sphere.sphs[0],
                              E->monitor.solution_cycles,"vor_V2");
        
    }
    return; 
}

void std2_timestep(E)
     struct All_variables *E;

{
    static int been = 0;
    static int stages = 0;
    static higher_precision init_u,diff_timestep,root3,root2;
    int i,d,n,nel,el,node;

    higher_precision adv_timestep, ts,uc1,uc2,uc3,uc,size,step,VV[4][9];

    const int dims=E->mesh.nsd;
    const int dofs=E->mesh.dof;
    const int nno=E->mesh.nno;
    const int lev=E->mesh.levmax;
    const int ends=enodes[dims];

   been ++;
    if(E->advection.fixed_timestep != 0.0) {
      E->advection.timestep = E->advection.fixed_timestep;
      return;
    }

  E->ve_data_cont.DIRECT=0;    // 1: means updating viscosity/stiffness matrix

 for (i=stages;i<E->ve_data_cont.stages;i++) {
    E->ve_data_cont.stage = i;
    E->advection.timestep = E->ve_data_cont.stages_timestep[i];

    if (i == 0)  
       been = 1;
    else 
       been = E->ve_data_cont.stages_step[i-1]+1;

    if ( E->monitor.solution_cycles < E->ve_data_cont.stages_step[i] )  {
       // neither 1st nor the last step od the stage 
       E->advection.next_timestep = E->ve_data_cont.stages_timestep[i];
       if (E->monitor.solution_cycles==been) {  // the 1st step of the stage
             E->ve_data_cont.DIRECT=1;
             }
       break;
    }
    else if ( E->monitor.solution_cycles == E->ve_data_cont.stages_step[i] )  {
       // before last timestep of current stage
       if ( i == E->ve_data_cont.stages-1 ) // last stage
           E->advection.next_timestep = E->advection.timestep ;
       else 
           E->advection.next_timestep = E->ve_data_cont.stages_timestep[i+1];
       break;
    }
  }
  stages = E->ve_data_cont.stage;

    return;
  }
