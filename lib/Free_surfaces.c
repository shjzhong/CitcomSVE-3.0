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
#include "parsing.h"
#include <math.h>
#include <sys/types.h>

#include <stdio.h>

/* ========================================================================== 
 * set_free_surfaces(E)
 * Called only once in read_instructions (in Instructions.c)
 * ========================================================================== */
void set_free_surfaces(E)
    struct All_variables *E;
{

    void parallel_process_termination();
    void  get_surface_loads();
    int lev,i,node1,node2,m;

    /* mark the vertical d.o.f. for nodes on density boundaries */

    if (E->parallel.me_loc[3] == 0)                      { 
        for (lev=E->mesh.gridmax;lev>=E->mesh.gridmin;lev--)
        for (m=1;m<=E->sphere.caps_per_proc;m++)  
        for (i=1;i<=E->lmesh.NEL[lev];i++)    
        if ((i-1)%E->lmesh.ELZ[lev]==0)   {
            E->ELEMENT[lev][m][i] = E->ELEMENT[lev][m][i] | BOUNDARIES;
            E->ELEMENT[lev][m][i] = E->ELEMENT[lev][m][i] | BOUNDARY2;    /*CMB*/
            E->NODE[lev][m][E->IEN[lev][m][i].node[1]] = E->NODE[lev][m][E->IEN[lev][m][i].node[1]] | RESTORE;    /*CMB*/
            E->NODE[lev][m][E->IEN[lev][m][i].node[2]] = E->NODE[lev][m][E->IEN[lev][m][i].node[2]] | RESTORE;    /*CMB*/
            E->NODE[lev][m][E->IEN[lev][m][i].node[3]] = E->NODE[lev][m][E->IEN[lev][m][i].node[3]] | RESTORE;    /*CMB*/
            E->NODE[lev][m][E->IEN[lev][m][i].node[4]] = E->NODE[lev][m][E->IEN[lev][m][i].node[4]] | RESTORE;    /*CMB*/
        }
    }
    if (E->parallel.me_loc[3] == E->parallel.nprocz-1)   {
        for (lev=E->mesh.gridmax;lev>=E->mesh.gridmin;lev--)
        for (m=1;m<=E->sphere.caps_per_proc;m++)  
        for (i=1;i<=E->lmesh.NEL[lev];i++)    
        if (i%E->lmesh.ELZ[lev]==0)   {
            E->ELEMENT[lev][m][i] = E->ELEMENT[lev][m][i] | BOUNDARIES;
            E->ELEMENT[lev][m][i] = E->ELEMENT[lev][m][i] | BOUNDARY1;    /*Surface*/
            E->NODE[lev][m][E->IEN[lev][m][i].node[5]] = E->NODE[lev][m][E->IEN[lev][m][i].node[5]] | RESTORE;
            E->NODE[lev][m][E->IEN[lev][m][i].node[6]] = E->NODE[lev][m][E->IEN[lev][m][i].node[6]] | RESTORE;
            E->NODE[lev][m][E->IEN[lev][m][i].node[7]] = E->NODE[lev][m][E->IEN[lev][m][i].node[7]] | RESTORE;
            E->NODE[lev][m][E->IEN[lev][m][i].node[8]] = E->NODE[lev][m][E->IEN[lev][m][i].node[8]] | RESTORE;
        }
    }

    void get_system_viscosity();//add by tao: need E->grav in get_surface_loads
	get_system_viscosity(E,1,E->EVI[E->mesh.levmax],E->VI[E->mesh.levmax]);
    
	// get the load (single harmonic, iceModel, or ice caps):
    /* get the initial rho*g*h for nodes on density boundaries */
    get_surface_loads(E);

    return;
}


/* ======================================================= 
  before assembling a force vector, determine topography
  at each density interface and convert it into stress 
  drho*g*h
 =======================================================  */

void get_boundary_forces(struct All_variables *E, int count) 
//    struct All_variables *E;
//    int count;
{
    int m,ll,mm,i,j,p;
    double density_cmb,density_surf,total, length;
    double modified_plgndr_a(),con,t1;
    void remove_average();
    void parallel_process_termination();
    double total_surface_integral();
    double *TG[4],temp1,temp2;

    void load_to_CM();
    void load_to_CM_grav();
	void load_to_CM_grav_comp();

    density_surf = 1.0;
    density_cmb = (E->data.density_below-E->data.density)/E->data.density;

        /* If this is the second (or later) self-gravity iteration for this
         * timestep, put the deformation from the previous iteration into
         * slice.surf/botm[2], then zero the CM for these loads (these loads
         * will be used in calculate_potential).
         * Then add them to the initial slice.load[0/2] to get the loads for 
         * this timestep, slice.load[1/3] (although I believe these are just
         * over-written before they're used).      */

    if (count==0)   {   
        /* slice.load[0] and slice.load[2] are the surf/cmb loads
         * from the end of the previous timestep, so return. */
      return;
      }
    else if (count)   {  //>>>> modify by tao:
		
		if(E->ve_data_cont.compressible == 0){ /* comment by tao: use surf/botm when incompress */
        // surface:
	        if (E->parallel.me_loc[3]==E->parallel.nprocz-1)   {
	            for (m=1;m<=E->sphere.caps_per_proc;m++)
	            for (j=1;j<=E->lmesh.nsf;j++)    {
	                i = j*E->lmesh.noz;
	                E->slice_ve.surf[2][m][j] = E->U[m][E->id[m][i].doff[3]]
	                                               * E->ve_data_cont.surf_scaling
	                                        + E->slice_ve.dynamic_oceanload[m][j];
	            }
	            remove_average(E,E->slice_ve.surf[2],1);
	        }
	        // cmb:
	        if (E->parallel.me_loc[3]==0)   {
	            for (m=1;m<=E->sphere.caps_per_proc;m++)
	            for (j=1;j<=E->lmesh.nsf;j++)    {
	                i = (j-1)*E->lmesh.noz + 1;
	                E->slice_ve.botm[2][m][j] = E->U[m][E->id[m][i].doff[3]];
	                E->slice_ve.botm[2][m][j] *= E->ve_data_cont.botm_scaling;
	            }
	            remove_average(E,E->slice_ve.botm[2],0);
	        }
		}
		else{
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
				if ( e%E->lmesh.elz==0 /* && E->parallel.me_loc[3]==E->parallel.nprocz-1*/ ) { /* at the surface */
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
				for (j=1;j<=E->lmesh.nsf;j++)	 {
				  node = j*E->lmesh.noz;
				  E->slice_ve.all_stress[m][node] += E->slice_ve.dynamic_oceanload[m][j]/stress_scaling; /* note that the rho used in dynamic_oceanload is actually rho_ice */
											  /* note that dynamic_oceanload is non-dim stress, so convert it to non-dim mass - rho_water * u */			
				}
			  }
			  }    
			}
			
			if (1 /*E->parallel.me_loc[3]==0*/)							   // CMB
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

		if(E->ve_data_cont.compressible){
			load_to_CM_grav_comp(E, 0, 1, 0);
		}
		else{
			load_to_CM_grav( E, E->slice_ve.surf[2], E->slice_ve.botm[2], 0 ); 
		}

        // todo_tao: 

    }// end if count
    	//<<<< end tao

    return;
}



/* ======================================================= 
   for each element el on boundaries, determine contributions
   from topographic deflection to elt_f
 =======================================================  */

void load_boundary_forces(struct All_variables *E, int el, double elt_f[24],int m)
//    struct All_variables *E;
 //   int el,m;
//    double elt_f[24];
{

    void get_global_1d_shape_fn_2();
//    void construct_c3x3matrix_el();
    void get_rtf_at_vpts();

    double radius2,temp,x[3][3],force[9],force_at_gs[9];
    int ic,nn[5],e, i, k, p, p2, lnode[5],node[9];
    double sinaa[9],ct[9],rtf[4][9];
    double div_w[4][9][9];
    double *gnxx,*gnda;

    const int lev = E->mesh.levmax;
    const int onedp = onedvpoints[E->mesh.nsd];
    const int vpts=vpoints[E->mesh.nsd];
    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
    const int ends=enodes[dims];

    struct Shape_function1 GM;
    struct Shape_function1_dA dGammax;

    static struct CC Cc;
    static struct CCX Ccx;

	//>>>> add by tao:
	double con;
	int a,j;
    void construct_c3x3matrix_el();
	
	//<<<< end tao

     // for incompressible
  if (!E->ve_data_cont.compressible)  {

    if (E->parallel.me_loc[3]==0 || E->parallel.me_loc[3]==E->parallel.nprocz-1)  {

        if (el%E->lmesh.elz==0 && E->parallel.me_loc[3]==E->parallel.nprocz-1)   {
            e = el/E->lmesh.elz;
            for(k=1;k<=onedvpoints[E->mesh.nsd];k++)   {
                lnode[k] = E->sien[m][e].node[k];
                force[k] = E->slice_ve.load[1][m][ lnode[k] ];
            }
            ic = 1;     /*      top  */
        }
        else if ((el-1)%E->lmesh.elz==0 && E->parallel.me_loc[3]==0)  {
            e = (el-1)/E->lmesh.elz+1;
            for(k=1;k<=onedvpoints[E->mesh.nsd];k++)   {
                lnode[k] = E->sien[m][e].node[k];
                force[k] = E->slice_ve.load[3][m][ lnode[k] ];
            }
            ic = 0;     /*      bottom */
        }
        else {
            return;
        }

        for(k=1;k<=onedvpoints[E->mesh.nsd];k++)   
            nn[k] = k+ic*4;

        for (i=1;i<=onedvpoints[E->mesh.nsd];i++)   {
            force_at_gs[i] = 0.0;
            for(k=1;k<=onedvpoints[E->mesh.nsd];k++)
                force_at_gs[i] += force[k] * E->M.vpt[GMVINDEX(k,i)];
        }

        for(k=1;k<=onedp;k++) {
            for (i=1;i<=onedp;i++)  
                elt_f[3*nn[k]-1] += E->B_R[lev][m][nn[k]][(e-1)*onedp+i] * force_at_gs[i];
            }
    }
  }
  else if (E->ve_data_cont.compressible) {
    con = 0.0;
    if (E->ve_data_cont.SELFG)
       con = E->ve_data_cont.Rsg;

    //e = (el-1)/E->lmesh.elz+1;

//if(E->parallel.me==0)
//printf("element %d \n",el);  //debug 20221213
    for(k=1;k<=onedvpoints[E->mesh.nsd];k++)   {
      lnode[k] = E->ien[m][el].node[k+onedp];
      force[k] = E->slice_ve.all_load[m][lnode[k]];
//if(E->parallel.me==0)
//printf("\t node: %d; all_load: %e \n", lnode[k], force[k]);
    }
   
    radius2 = E->sx[m][3][lnode[1]]*E->sx[m][3][lnode[1]];

    for (i=1;i<=onedvpoints[E->mesh.nsd];i++)   {
      force_at_gs[i] = 0.0;
      for(k=1;k<=onedvpoints[E->mesh.nsd];k++)
        force_at_gs[i] += force[k] * E->M.vpt[GMVINDEX(k,i)];
    }

    /* elt_f has been set to zero before this subroutine is called, grep assemble_forces */
    /* and only the r components of elt_f are changed */
    for(k=1;k<=onedp;k++)
      for (i=1;i<=onedp;i++)
        elt_f[3*(k+onedp)-1] += /*radius2* */E->B_RR[lev][m][k+onedp][(el-1)*onedp+i] * force_at_gs[i]; // here is el, not e; e is 0 ~ #of surface
// modified by tao: before shijie used radius2*E->B_R, which will be same as E->B_RR.

    /* also pick the bottom surface if it is at CMB */
    if ((el-1)%E->lmesh.elz==0 /* && E->parallel.me_loc[3]==0*/)  {
      for(k=1;k<=onedvpoints[E->mesh.nsd];k++)   {
        lnode[k] = E->ien[m][el].node[k];
        force[k] = E->slice_ve.all_load[m][lnode[k]];
      }
      for (i=1;i<=onedvpoints[E->mesh.nsd];i++)   {
        force_at_gs[i] = 0.0;
        for(k=1;k<=onedvpoints[E->mesh.nsd];k++)
          force_at_gs[i] += force[k] * E->M.vpt[GMVINDEX(k,i)];

      }
      for(k=1;k<=onedp;k++)
        for (i=1;i<=onedp;i++)
          elt_f[3*k-1] += E->B_RR[lev][m][k][(el-1)*onedp+i] * force_at_gs[i];
    }   // no need for multiplying radius^2 here for the CMB, as it was already, 

	// also add in the div_w*rho*phi*dV term */
	//  compute the rho*phi at gaussian points */

	    for (a=1;a<=ends;a++) {
	      node[a] = E->ien[m][el].node[a];
	      force[a] = con*E->slice_ve.total_potential[m][node[a]];
	    }
	    for (i=1;i<=vpts;i++) {
	      force_at_gs[i] = 0.0;
	      for (k=1;k<=ends;k++)
	        force_at_gs[i] += E->erho[lev][m][el]*force[k]*E->N.vpt[GNVINDEX(k,i)];
	    }
	// compute div_w at each dimension 
	    gnxx = E->gNX[m][el].vpt;
	    gnda = E->gDA[m][el].vpt;
	    get_rtf_at_vpts(E,m,E->mesh.levmax,el,rtf);

	    if ((el-1)%E->lmesh.elz==0)
	      construct_c3x3matrix_el(E,el,&Cc,&Ccx,E->mesh.levmax,m,0);

	    for (j=1;j<=vpts;j++) {
	      sinaa[j] = sin(rtf[1][j]);
	      ct[j] = cos(rtf[1][j])/sinaa[j];
	      for (i=1;i<=ends;i++) {
	        for (k=1;k<=dims;k++){
	          div_w[k][i][j] = rtf[3][j]*
	                 (gnxx[GNVXINDEX(0,i,j)]*Cc.vpt[BVINDEX(1,k,i,j)]
	                 +E->N.vpt[GNVINDEX(i,j)]*Ccx.vpt[BVXINDEX(1,k,1,i,j)]
	                 +E->N.vpt[GNVINDEX(i,j)]*Cc.vpt[BVINDEX(3,k,i,j)])
	                          +rtf[3][j]*
	                (E->N.vpt[GNVINDEX(i,j)]*Cc.vpt[BVINDEX(1,k,i,j)]*ct[j]
	                +E->N.vpt[GNVINDEX(i,j)]*Cc.vpt[BVINDEX(3,k,i,j)]
	                +(gnxx[GNVXINDEX(1,i,j)]*Cc.vpt[BVINDEX(2,k,i,j)]
	                 +E->N.vpt[GNVINDEX(i,j)]*Ccx.vpt[BVXINDEX(2,k,2,i,j)])/sinaa[j])
	                          +gnxx[GNVXINDEX(2,i,j)]*Cc.vpt[BVINDEX(3,k,i,j)];
	        }
	      }
    }
	//      add to elf_f 
	    for(a=1;a<=ends;a++) {
	      p = dims*(a-1);
	      for(j=1;j<=vpts;j++) {
	        elt_f[p  ] -= div_w[1][a][j]*force_at_gs[j]*gnda[j]*g_point[j].weight[dims-1];
	        elt_f[p+1] -= div_w[2][a][j]*force_at_gs[j]*gnda[j]*g_point[j].weight[dims-1];
	        elt_f[p+2] -= div_w[3][a][j]*force_at_gs[j]*gnda[j]*g_point[j].weight[dims-1];
	      }
	    }

  }


    return;
}


/* ========================================================================== 
 * get_surface_loads(E)
 * Called only once in set_free_surfaces (above).
 * Gets the ice load (single harmonic, iceModel, or ice caps).
 * ========================================================================== */
void get_surface_loads(struct All_variables *E) 
 //   struct All_variables *E;
{
    FILE *fp;
    int oldstep,m,ll,mm,p,i,j,n;
    double tau,total, length, stress_scale;
    double dist,distance(),modified_plgndr_a(),con,t1,f1;
    double x1,y1,z1,density_ice,r,t,f,ri,ro,icet,icef,icer,oldampt;
    int m1 = E->parallel.me;

    double itemp[100], temp[100];

    char outfile[255];
    void sphere_expansion_output();
    void remove_average();
    void get_iceModel();  
    void parallel_process_termination();
    void apply_new_loads(); 

    void get_static_oceanload();
    void load_to_CM() ;
    void load_to_CM_grav() ;
	void load_to_CM_grav_comp();

    E->ve_data_cont.potential_scaling = 4.0*M_PI*E->data.grav_const*E->data.density
        *E->sphere.dradius*E->sphere.dradius;
    E->ve_data_cont.rotation_rate = E->ve_data_cont.rotation_rate*2.0*M_PI/(24.0*3600.0);

        E->ve_data_cont.Rsg = 4.0*M_PI*E->data.grav_const*E->data.density
            *E->data.density*E->sphere.dradius*E->sphere.dradius
            /E->ve_data_cont.shear_mod;
        // the following are used to convert nondim height -> nondim stress
        E->ve_data_cont.botm_scaling = (E->data.density_below-E->data.density)
                              *E->data.grav_acc
                              *E->sphere.dradius/E->ve_data_cont.shear_mod;
        E->ve_data_cont.surf_scaling = E->data.density*E->data.grav_acc
                              *E->sphere.dradius /E->ve_data_cont.shear_mod;

    density_ice = 917.4;
    E->ve_data_cont.ice_stress_scale = E->ve_data_cont.surf_scaling*density_ice/E->data.density;

    if (E->parallel.me==0)  {
        fprintf(E->fp,"scale %.5e %.5e %.5e %.5e\n",
                      E->ve_data_cont.botm_scaling,
                      E->ve_data_cont.surf_scaling,
                      E->ve_data_cont.ice_stress_scale,E->ve_data_cont.tau);
        fflush(E->fp); 
    }

// the following three lines are moved to Instructions.c
//    input_int("stages",&(E->ve_data_cont.stages),"1",m1);
//    input_int_vector("step",E->ve_data_cont.stages,(E->ve_data_cont.stages_step),m1);
//    input_double_vector("timestep",E->ve_data_cont.stages,(E->ve_data_cont.stages_time),m1);

    oldstep = 0;
    for (i=0;i<E->ve_data_cont.stages;i++)  {
        if (E->ve_data_cont.Heaviside==1)
            E->ve_data_cont.stages_timestep[i] = E->ve_data_cont.stages_time[i]
                /((E->ve_data_cont.stages_step[i]-oldstep));
        else 
            E->ve_data_cont.stages_timestep[i] = E->ve_data_cont.stages_time[i]
                /(E->ve_data_cont.tau*(E->ve_data_cont.stages_step[i]-oldstep));
        oldstep = E->ve_data_cont.stages_step[i];
    }

    if (E->ve_data_cont.Heaviside==1)   { 
        E->monitor.elapsed_time=0.0;
        E->ve_data_cont.DIRECT=1;
        }
    else  {
        // set the current timstep to 1, and the current elapsed_time to dt:
        E->monitor.elapsed_time=E->ve_data_cont.stages_timestep[0];
        E->advection.timestep=E->ve_data_cont.stages_timestep[0];
        E->monitor.solution_cycles = 1;
        E->ve_data_cont.DIRECT=1;
    }

    // Set the initial loads:
    
    E->ve_data_cont.change_of_load = 1;  // used by apply_new_loads

    if (E->ve_data_cont.Heaviside==1)    {
        /* loading at a single harmonic */
        mm = E->convection.perturb_mm[0];
        ll = E->convection.perturb_ll[0];
        con = E->convection.perturb_mag[0];

        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)    {
            i = j*E->lmesh.noz;
            t1 = E->SX[E->mesh.levmax][m][1][i];
            f1 = E->SX[E->mesh.levmax][m][2][i];
            if (E->ve_data_cont.apply_potential==0)
                E->Xsurf[3][m][j] = con*modified_plgndr_a(ll,mm,t1)*cos(mm*f1);
            else if (E->ve_data_cont.apply_potential==1) { 
                E->Xsurf[3][m][j] = 0.0;
                E->init_potential[0][m][j] = con*E->sphere.ro*E->sphere.ro*modified_plgndr_a(ll,mm,t1)*cos(mm*f1);
                E->init_potential[1][m][j] = con*E->sphere.ri*E->sphere.ri*modified_plgndr_a(ll,mm,t1)*cos(mm*f1);
                }
        }
        remove_average(E,E->Xsurf[3],1);

        sphere_expansion_output(E,1,E->Xsurf[3],
                                E->sphere.sphc[0],E->sphere.sphs[0],
                                E->monitor.solution_cycles,"init_surf");


		if(E->ve_data_cont.compressible){  // add by tao: compressible version
			int node;
			
	        for (m=1;m<=E->sphere.caps_per_proc;m++)
	          for (i=1;i<=E->lmesh.nno;i++)
	            E->slice_ve.all_stress[m][i] = 0.0;

			if(E->parallel.me_loc[3]==E->parallel.nprocz-1)
			for (m=1;m<=E->sphere.caps_per_proc;m++)
	        for (j=1;j<=E->lmesh.nsf;j++)    {
					node = j*E->lmesh.noz;
					E->slice_ve.all_stress[m][node] = E->Xsurf[3][m][j]; // comment by tao: here I follow geruo, all_stress is rho*u, or dynamic_oceanload / stress_scaling.
																			// so do not *E->ve_data_cont.surf_scaling.
//if(E->parallel.me==0) //debug 20221213 and below
//	printf("node: %d; all_stress: %e \n", node, E->slice_ve.all_stress[m][node]);
	        	}

/*        if(0||E->parallel.me==0 || E->parallel.me==1){
            printf("---debug E->all_stress (count = %d) -----\n", E->monitor.count);
            char debug_outfile[255];
            sprintf(debug_outfile, "debug_all_stress_Mine.%d.%d.proc%d",E->monitor.solution_cycles,E->monitor.count, E->parallel.me);
            FILE * debug_fp;
            debug_fp = fopen(debug_outfile, "w");
            int debug_i;
            for(debug_i=1; debug_i<=E->lmesh.nno; debug_i++){
                fprintf(debug_fp,"%e\n", E->slice_ve.all_stress[1][debug_i]);
            }
            fclose(debug_fp);
        }

        double global_tdot_d_noskip();
        double _temp1;
        _temp1 =  global_tdot_d_noskip(E,E->slice_ve.all_stress,E->slice_ve.all_stress, E->mesh.levmax);
        printf("all_stress(before shiftV): %e; from core %d \n", _temp1, E->parallel.me);
*/
			for(i=1; i<=1; i++){
				load_to_CM_grav_comp(E, 1, 1, 1);
				if(E->parallel.me==0){
    				printf("debug deg-1 motion: now with shift_V_to_CM\n");
        			printf("CM_incr: %e %e %e\n",E->ve_data_cont.CM_incr[0],E->ve_data_cont.CM_incr[1],E->ve_data_cont.CM_incr[2]);
				}
    			{// add by tao: follow geruo: note for ice loading case, this has not implemented
        		void shift_V_to_CM();
				if(E->parallel.me==0)
    				printf("shift_V_to_CM called\n");
        		shift_V_to_CM(E);
    			}

				//if(E->parallel.me==0)
				/*{
					printf("---debug E->all_stress (after shift_V_to_CM) -----\n");
					char debug_outfile[255];
					sprintf(debug_outfile, "debug_all_stress_Mine.%d.%d",E->monitor.solution_cycles,E->parallel.me);
					FILE * debug_fp;
					debug_fp = fopen(debug_outfile, "w");
					int debug_i;
					for(debug_i=1; debug_i<=E->lmesh.nno; debug_i++){
						fprintf(debug_fp,"%e\n", E->slice_ve.all_stress[1][debug_i]);
					}
					fclose(debug_fp);
				}*/
			}
			for (m=1;m<=E->sphere.caps_per_proc;m++)
				for (i=1;i<=E->lmesh.nno;i++)
					E->slice_ve.total_load[m][i] = 0.0;
				
			double stress_scaling = E->data.density*E->data.grav_acc*E->sphere.dradius/E->ve_data_cont.shear_mod;
			for (m=1;m<=E->sphere.caps_per_proc;m++)
			for (j=1;j<=E->lmesh.nno;j++)	 {
					node = j;
					E->slice_ve.total_load[m][node] = E->slice_ve.all_stress[m][node] * E->grav[E->mesh.levmax][m][node]*stress_scaling;
	
/*			
if(E->parallel.me==0)
{
	printf("node: %d; grav: %e \n", node, E->grav[E->mesh.levmax][m][node]);	
	printf("node: %d; total_load: %e \n", node, E->slice_ve.total_load[m][node]);
}*/
			}

			
			}
		else{
	        for (m=1;m<=E->sphere.caps_per_proc;m++)
	        for (j=1;j<=E->lmesh.nsf;j++)    {
	            E->slice_ve.surf[2][m][j] = E->Xsurf[3][m][j]*E->ve_data_cont.surf_scaling;
	            E->slice_ve.botm[2][m][j] = 0.0 ;
	        }

			load_to_CM_grav( E, E->slice_ve.surf[2], E->slice_ve.botm[2], 1 );
			for (m=1;m<=E->sphere.caps_per_proc;m++)
			for (j=1;j<=E->lmesh.nsf;j++)	 {
				E->slice_ve.load[0][m][j] = E->slice_ve.surf[2][m][j] ;
				E->slice_ve.load[2][m][j] = E->slice_ve.botm[2][m][j] ;
			}
		}
		
      
    }
    else if (E->ve_data_cont.Heaviside==2)  {   /* ice-model */
        // read 1st epoch of iceModel and (if appropriate) the oceanload:
        get_iceModel(E,E->slice_ve.iceload[0]);  
        if (E->ve_data_cont.SLE) get_static_oceanload(E);

    } // end ice loading

/*        double global_tdot_d_noskip();
        double _temp;
        _temp =  global_tdot_d_noskip(E,E->slice_ve.all_stress,E->slice_ve.all_stress, E->mesh.levmax);
        printf("all_stress: %e; from core %d \n", _temp, E->parallel.me);
*/
    // Put the initial loads into slice.load[0/2], and set init_potential:
    apply_new_loads(E);

    if (E->ve_data_cont.Heaviside==1)  
        E->ve_data_cont.change_of_load = 0; // no more load changes in this case
/*
if(E->parallel.me==0){
    printf("---debug E->init_total_potential (count = 0) -----\n");
    char debug_outfile[255];
    sprintf(debug_outfile, "debug_init_total_potential.%d",E->monitor.solution_cycles);
    FILE * debug_fp;
    debug_fp = fopen(debug_outfile, "w");
    int debug_i;
    for(debug_i=1; debug_i<=E->lmesh.nno; debug_i++){
        //if(debug_i%3==0) fprintf(debug_fp,"\n");
          fprintf(debug_fp,"%e \n", E->slice_ve.init_total_potential[1][debug_i]);
          }
    fclose(debug_fp);
	}
*/
    return;
}

/* =============================================================================
 * double distance(x0,y0,z0,t1,f1)
 *    dist = 2.0*arcsin(dist/diameter) = 2.0*arcsin(0.5*dist/radius) 
 * ===========================================================================*/
double distance(x0,y0,z0,t1,f1)
    double x0,y0,z0,t1,f1;
{

    double dist,x2,y2,z2;

    x2 = sin(t1)*cos(f1); 
    y2 = sin(t1)*sin(f1); 
    z2 = cos(t1);
    dist = sqrt( (x2-x0)*(x2-x0) + (y2-y0)*(y2-y0) + (z2-z0)*(z2-z0) );
    dist = 0.5*dist;
    dist = 2.0*asin(dist);
    /* dist = 2.0*arcsin(dist/diameter) = 2.0*arcsin(0.5*dist/radius) */
    return (dist);
}


/* ======================================================= 
  for self-G cases, the relation between stress srr and boundary
  topography needs to add the effects of potential for CMB.
 =======================================================  */
void  append_boundary_forces(E)
    struct All_variables *E;
{
    int m,is,ib,j;
    double con,con1;
    int debug_write_load = 0;
    void sphere_expansion_output();

    if (E->ve_data_cont.SELFG)  {// for conversion of potential -> geoid 
        con1 = E->ve_data_cont.Rsg; // surface
        con =  E->ve_data_cont.Rsg  // CMB
               * (E->data.density_below-E->data.density) / E->data.density;
    }
    else  con = con1 = 0.0 ;

    if (E->parallel.me_loc[3]==E->parallel.nprocz-1)          // surface
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)    {
            E->slice_ve.load[1][m][j] = -(  E->slice_ve.load[0][m][j] 
                                       + E->slice_ve.dynamic_oceanload[m][j] )
                                     + con1*E->potential[0][m][j];
            /* E->slice_ve.load[1][m][j] = -E->slice_ve.load[1][m][j] 
                                     + con1*E->potential[m][is];*/
        }
    if (E->parallel.me_loc[3]==0)                              // CMB
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)    {
            /* the negative sign out of ( ... ) is from the surface
               normal at the bottom is negative or points to the
               -r direction */
            E->slice_ve.load[3][m][j] = - ( E->slice_ve.load[2][m][j] 
                                         - con*E->potential[1][m][j] ) ; 
            /* E->slice_ve.load[3][m][j] = - ( E->slice_ve.load[3][m][j] 
                                        - con*E->potential[m][ib]) ; */
        }

    if (debug_write_load) 
        sphere_expansion_output(E,1,E->slice_ve.load[1],
                                E->sphere.sphc[0],E->sphere.sphs[0],
                                E->monitor.solution_cycles,"surfaceload");

    return;
}

void  append_multi_boundary_forces(struct All_variables * E)
//    struct All_variables *E;
{
    int m,is,ib,j,i,e;
    double con,con1;
    int debug_write_load = 0;
    void sphere_expansion_output();

    int n,node;
    double density_jump, con2;
    const int onedp = onedvpoints[E->mesh.nsd];
    const int vpts=vpoints[E->mesh.nsd];
    const int lev = E->mesh.levmax;

    if (E->ve_data_cont.SELFG)  { 
        con = E->ve_data_cont.Rsg; 
    }
    else  con = 0.0 ;

    for (m=1;m<=E->sphere.caps_per_proc;m++) {
      for (e=1;e<=E->lmesh.nel;e++){
        /* find delta-rho across the top surface of each element*/
        density_jump = 0.0;
        if ( e%E->lmesh.elz==0 /*&& E->parallel.me_loc[3]==E->parallel.nprocz-1 */) { /* at the surface */
          //for (n=1;n<=vpts;n++)
          if(E->parallel.me_loc[3]==E->parallel.nprocz-1)
            density_jump += E->erho[lev][m][e];
		  else
			density_jump += E->erho[lev][m][e] - 1.0;
          //density_jump = density_jump/vpts;
        }
        else {
          //for (n=1;n<=vpts;n++)
            density_jump = density_jump + E->erho[lev][m][e] - E->erho[lev][m][e+1];
          //density_jump = density_jump/vpts;
        }
       /* compute all_load at top nodes */
       for (j=1;j<=onedp;j++){
         node = E->ien[m][e].node[j+onedp];
         E->slice_ve.all_load[m][node] = -E->slice_ve.total_load[m][node] + con*density_jump*E->slice_ve.total_potential[m][node];
       }
      } // end e loop

      /* add dynamic_oceanload to the Earth's surface */
      if ( E->parallel.me_loc[3]==E->parallel.nprocz-1 ) {
        for (j=1;j<=E->lmesh.nsf;j++)    {
          node = j*E->lmesh.noz;
          E->slice_ve.all_load[m][node] -= E->slice_ve.dynamic_oceanload[m][j];
          // ideally, should mupltiply by non-dim surface g, but it is 1.0
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
              E->slice_ve.all_load[m][node] = -E->slice_ve.total_load[m][node] + con*density_jump*E->slice_ve.total_potential[m][node];
            }
          }
        }


    return;
}


/*==============================================================================
 * calculate_potential (E,X_surf,X_cmb,potential_surf,potential_cmb)
 * =============================================================================
 * Takes the surface masses given in X_surf/cmb, computes their gravitational
 * potential at the surface and CMB, and puts it in potential_surf/cmb. 
 * If E->control.polar_wander is true, it further changes the potential arrays
 * (in spherical harmonic l,m=2,1) to via polar_wander_effects.
 * Notes:
 *     X_surf/cmb must be nondimensional stresses (such as E->iceload).
 *     This routine will not change the values of X_surf and X_cmb.
 *=========================================================================== */
void calculate_potential(E,X_surf,X_cmb,potential_surf,potential_cmb,icon)
    struct All_variables *E;
    double **X_surf;
    double **X_cmb;
    double **potential_surf;
    double **potential_cmb;
    int icon;
{
    int ib,is,m,ll,ll1,ll2,mm,p,i,j,k,n;
    double con,r,t,f,rir,ri,ro,modified_plgndr_a();
    double density_surf,density_cmb;

    static int been=0;
    static double brll1[65],brll2[65];
    static double srll1[65],srll2[65];
    void  sphere_expansion_VE();
    void  polar_wander_effects();
    void  exchange_sphcs();

    ri = E->sphere.ri;
    ro = E->sphere.ro;
    density_surf = 1.0;
    density_cmb = (E->data.density_below-E->data.density)/E->data.density;

    if (been==0)  {
        been = 1;
        rir = 1.0;
        for (ll=1;ll<=E->output.llmax;ll++)    {
            brll1[ll] = pow(ri,(double)(ll))/(2.0*ll+1.0); 
            brll2[ll] = ri*pow(rir,(double)(ll+1))/(2.0*ll+1.0); 
            srll1[ll] = pow(ro,(double)(ll))/(2.0*ll+1.0); 
            srll2[ll] = ri*pow(ri,(double)(ll+1))/(2.0*ll+1.0); 
        }
    }

   for (m=1;m<=E->sphere.caps_per_proc;m++)
   for (j=1;j<=E->lmesh.nsf;j++) {
        potential_surf[m][j] = 0.0; 
        potential_cmb[m][j] = 0.0; 
        E->Xsurf[1][m][j] = 0.0;
        E->Xsurf[2][m][j] = 0.0;
        }

    // Ylm expansion of surface mass:
    if (E->parallel.me_loc[3]==E->parallel.nprocz-1)   {
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)   {
            // scale nondim stress -> nondim height of rock
            E->Xsurf[1][m][j] = X_surf[m][j]/E->ve_data_cont.surf_scaling;
        }
        sphere_expansion_VE(E,1,E->Xsurf[1],
                         E->sphere.sphc[0],E->sphere.sphs[0],E->output.llmax);
        for (ll=0;ll<=E->output.llmax;ll++)
        for (mm=0;mm<=ll;mm++)   {
            p = E->sphere.hindex[ll][mm];
            E->sphere.sphc[0][p] *= density_surf;
            E->sphere.sphs[0][p] *= density_surf;
        }
    }

    // Ylm expansion of CMB mass:
    if (E->parallel.me_loc[3]==0)  {
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)   {
            // scale nondim stress -> nondim height of cmb boundary
            E->Xsurf[2][m][j] = X_cmb[m][j]/E->ve_data_cont.botm_scaling;
        }
        sphere_expansion_VE(E,0,E->Xsurf[2],
                          E->sphere.sphc[1],E->sphere.sphs[1],E->output.llmax);
        for (ll=0;ll<=E->output.llmax;ll++)
        for (mm=0;mm<=ll;mm++)   {
            p = E->sphere.hindex[ll][mm];
            E->sphere.sphc[1][p] *= density_cmb;
            E->sphere.sphs[1][p] *= density_cmb;
        }
    }

    if (E->parallel.nprocz>1)  // if >1 cap in z-direction (probably not)
        exchange_sphcs(E,E->sphere.sphc[0],E->sphere.sphs[0],
                         E->sphere.sphc[1],E->sphere.sphs[1]);

    // Now calculate potential:
        // cmb:
    if (E->parallel.me_loc[3]==0)  {
      for (m=1;m<=E->sphere.caps_per_proc;m++)
      for (j=1;j<=E->lmesh.nsf;j++) {
        potential_cmb[m][j] = 0.0;
        for (ll=1;ll<=E->output.llmax;ll++)  
        for (mm=0;mm<=ll;mm++)   { 
            p = E->sphere.hindex[ll][mm]; 
            potential_cmb[m][j] += ( E->Tbl_cs[m][mm][j]*
              (brll1[ll]*E->sphere.sphc[0][p] + brll2[ll]*E->sphere.sphc[1][p])
               + E->Tbl_sn[m][mm][j]*
              (brll1[ll]*E->sphere.sphs[0][p] + brll2[ll]*E->sphere.sphs[1][p]))
                *E->Tbl_lm[m][p][j];
        }
      }
    }
        // surface:
    if (E->parallel.me_loc[3]==E->parallel.nprocz-1)   {
      for (m=1;m<=E->sphere.caps_per_proc;m++)
      for (j=1;j<=E->lmesh.nsf;j++) {
        potential_surf[m][j] = 0.0; 
        for (ll=1;ll<=E->output.llmax;ll++)   
        for (mm=0;mm<=ll;mm++)   { 
            p = E->sphere.hindex[ll][mm]; 
            potential_surf[m][j] += ( E->Tbl_cs[m][mm][j]*
             (srll1[ll]*E->sphere.sphc[0][p] + srll2[ll]*E->sphere.sphc[1][p])
                    + E->Tbl_sn[m][mm][j]*
             (srll1[ll]*E->sphere.sphs[0][p] + srll2[ll]*E->sphere.sphs[1][p]) )
                *E->Tbl_lm[m][p][j];
        }
      }
    }

    // alter the Y_21 potential, if polar_wander is on:
    if (E->ve_data_cont.polar_wander && icon==1)  
        // Xsurf[1/2] have already been converted to nondim height of boundary
//        polar_wander_effects( E, E->Xsurf[1], E->Xsurf[2], 
//                                 potential_surf, potential_cmb );
        polar_wander_effects( E, X_surf, X_cmb, 
                                 potential_surf, potential_cmb );

    return;
}




//>>>> add by tao: only calculate the phi, the CM effect should be corrected before separately
// todo_tao: not implemented.

void exchange_sphcs_to_all(struct All_variables *E, double **sphc, int i) // i means the sending layer, index of sphs, global z index.
// struct All_variables *E;
// double *sphc , *sphs ;
{

 void parallel_process_termination();
 int p,target_proc,msginfo[8];


 static int been_here = 0;
 static int idb,sizeofk;
 static double *RV, *RV_positive, *RV_negative;

 MPI_Status status[100];
 MPI_Status status1;
 MPI_Request request[100];

 if (been_here ==0 )   {

   idb =  E->sphere.hindice;
   sizeofk = (idb+1)*sizeof(double);
	RV = (double *)malloc( sizeofk );
//   RV_positive=(double *)malloc( sizeofk );
//   RV_negative=(double *)malloc( sizeofk );
   been_here ++;
   }

 //   MPI_Barrier(E->parallel.vertical_comm);
    MPI_Allreduce(sphc[i], RV, idb, MPI_DOUBLE, MPI_SUM , E->parallel.vertical_comm);
    
//    MPI_Allreduce(sphc[i], RV_positive, idb, MPI_DOUBLE, MPI_MAX , E->parallel.vertical_comm);
//    MPI_Allreduce(sphc[i], RV_negative, idb, MPI_DOUBLE, MPI_MIN , E->parallel.vertical_comm);
    // Synchronize again before obtaining final time
//    MPI_Barrier(E->parallel.vertical_comm);

    for (p=0;p<E->sphere.hindice;p++)  {
      sphc[i][p] = RV[p]; // _positive[p] + RV_negative[p];
      }

	
//    MPI_Barrier(E->parallel.vertical_comm);
//	
////    MPI_Allreduce(sphs[i], RV, idb, MPI_DOUBLE, MPI_SUM , E->parallel.vertical_comm);
//    MPI_Allreduce(sphs[i], RV_positive, idb, MPI_DOUBLE, MPI_MAX , E->parallel.vertical_comm);
//    MPI_Allreduce(sphs[i], RV_negative, idb, MPI_DOUBLE, MPI_MIN , E->parallel.vertical_comm);
//    // Synchronize again before obtaining final time
//    MPI_Barrier(E->parallel.vertical_comm);
//    for (p=0;p<E->sphere.hindice;p++)  {
//      sphs[p] = RV_positive[p] + RV_negative[p];
//      }
	
	return;
}

void exchange_depth_to_all(struct All_variables *E, double *depth) 
// double *sphc , *sphs ;
{

 void parallel_process_termination();
 int p,target_proc,msginfo[8];


 static int been_here = 0;
 static int idb,sizeofk;
// static double *RV_positive; //, *RV_negative;
 static double *RV;



 MPI_Status status[100];
 MPI_Status status1;
 MPI_Request request[100];

 if (been_here ==0 )   {

   idb =  E->mesh.noz + 1;
   sizeofk = idb*sizeof(double);

   RV =(double *)malloc( sizeofk );
 //  RV_negative=(double *)malloc( sizeofk );
   been_here ++;
   }

 //   MPI_Barrier(E->parallel.vertical_comm);
    MPI_Allreduce(depth, RV, idb, MPI_DOUBLE, MPI_SUM , E->parallel.vertical_comm);
    
  //  MPI_Allreduce(depth, RV_positive, idb, MPI_DOUBLE, MPI_MAX , E->parallel.vertical_comm);
 //   MPI_Allreduce(depth, RV_negative, idb, MPI_DOUBLE, MPI_MIN , E->parallel.vertical_comm);
    // Synchronize again before obtaining final time
 //   MPI_Barrier(E->parallel.vertical_comm);

    for (p=1;p<=E->mesh.noz;p++)  {
//      depth[p] = RV_positive[p] + RV_negative[p];
		depth[p] = RV[p];

      }

	
//    MPI_Barrier(E->parallel.vertical_comm);
//	
////    MPI_Allreduce(sphs[i], RV, idb, MPI_DOUBLE, MPI_SUM , E->parallel.vertical_comm);
//    MPI_Allreduce(sphs[i], RV_positive, idb, MPI_DOUBLE, MPI_MAX , E->parallel.vertical_comm);
//    MPI_Allreduce(sphs[i], RV_negative, idb, MPI_DOUBLE, MPI_MIN , E->parallel.vertical_comm);
//    // Synchronize again before obtaining final time
//    MPI_Barrier(E->parallel.vertical_comm);
//    for (p=0;p<E->sphere.hindice;p++)  {
//      sphs[p] = RV_positive[p] + RV_negative[p];
//      }
	
	return;
}


void calculate_potential_comp(struct All_variables *E, double ** all_potential, int load_only)
{
	void calculate_potential_comp_OR_deg1_2();

	calculate_potential_comp_OR_deg1_2(E,all_potential,load_only,1,0,0);
	return;
}

// add by tao: also used in load_to_CM_grav_comp, that's where only_deg1_2 comes. add_incr_CM only works with only_deg1_2
void calculate_potential_comp_OR_deg1_2(struct All_variables *E, double ** all_potential, int load_only, int icon, int only_deg1_2, int add_incr_CM)
//    struct All_variables *E;
//    double **X_surf;
//    double **X_cmb;
//    double **potential_surf;
//    double **potential_cmb;
//    int icon;
{
    int ib,is,m,ll,ll1,ll2,mm,p,i,j,k,n;
    double con,r,t,f,rir,ri,ro,modified_plgndr_a();
    double density_surf,density_cmb;

    static int been=0;
//    static double brll1[65],brll2[65];
//    static double srll1[65],srll2[65];
    void  sphere_expansion_VE();
    void  polar_wander_effects();
    void  exchange_sphcs();
	void remove_average();
    void construct_c3x3matrix_el();

	int global_z_id;
	int node;
//	static double local_r[400];
	static double global_r[400]; 

	static double * rho_div_u[NCS];
	static double M_conv, C_conv, omega;


    ri = E->sphere.ri;
    ro = E->sphere.ro;
    density_surf = 1.0;
    density_cmb = (E->data.density_below-E->data.density)/E->data.density;

	double time, time0;
	double CPU_time0();
        if(E->parallel.me==0 && E->advection.timesteps%50==0){
                time0 = CPU_time0() ;
                fprintf(E->fp,"DEBUG (In calculate_potential_comp_OR_deg1) : begin  \n");
        }

    if (been==0)  {
        been = 1;
        rir = 1.0;

		for(m=1; m<=E->sphere.caps_per_proc;m++){
			rho_div_u[m] = (double *) malloc((E->lmesh.nno+2)*vpoints[E->mesh.nsd]*sizeof(double));
		}

		omega = 2.0*M_PI/(24.0*3600.0);  // rotation_rate
        M_conv = (3.0*E->data.grav_const*E->data.density)/(omega*omega*E->ve_data_cont.kf)*5.0*sqrt(4.0*M_PI/15.0); // 3G*rho/(w^2*kf)*5sqrt(4pi/15)
          // convert non-dim phi into non-dim M
        C_conv = -omega*omega/(4.0*M_PI*E->data.grav_const*E->data.density);
          // convert non-dim r^2(m1*cosf+m2*sinf)*sintcost into non-dim centrifugal potential, which can be added into the potential terms;
          // to dimensionalize the potential, multiply by 4*pi*G*rho*R^2.
		
		
		//>>>> add by tao: exchange global radius

		for (i = 1; i <= E->mesh.noz; ++i)
			{
			global_r[i] = 0.0;
			}
		
		for(i=2; i<=E->lmesh.noz;i++){
			for(m=1;m<=E->sphere.caps_per_proc;m++){
				node = i; // first column
				global_z_id = E->lmesh.nzs + i -1;
				global_r[global_z_id] = E->SX[E->mesh.levmax][m][3][node];
			}
		}
		
		if(E->parallel.me_loc[3]==0){
			global_r[1] = E->SX[E->mesh.levmax][1][3][1];
		}
		
		if(E->parallel.nprocz>1){
			exchange_depth_to_all( E, global_r);
		}
		/*		
		if(E->sphere.capid[1] == 1)
		{
			printf("debug global_r \n");
			printf("nzs: %d\n", E->lmesh.nzs);
			for(i=1; i<=E->mesh.noz;i++){
				printf("global[%d] = %e from proc %d \n", i, global_r[i], E->parallel.me);	
			}
		}
		*/
		//<<<< end tao
    }

//   for (m=1;m<=E->sphere.caps_per_proc;m++)
//   for (j=1;j<=E->lmesh.nsf;j++) {
////        potential_surf[m][j] = 0.0; 
////        potential_cmb[m][j] = 0.0; 
//        E->Xsurf[1][m][j] = 0.0;
//        E->Xsurf[2][m][j] = 0.0;
//        }


	//comment by tao: should set to 0 every time called, since only layer in this processor is modified, others will be lefted.
	for (i = 1; i <= E->mesh.noz; ++i) // commented by tao: only worked for compressible, otherwise out of range
	{
		for (p = 0; p < E->sphere.hindice+2; ++p)
		{
			E->sphere.sphc_c[i][p] = 0.;
			E->sphere.sphs_c[i][p] = 0.;
			E->sphere.sphc_v[i][p] = 0.;
			E->sphere.sphs_v[i][p] = 0.;
		}
	}


assert(E->mesh.noz<300); // add by tao: the sphc_c was defined haing a max layer of 300

int llmax;

if(!only_deg1_2)
	llmax=E->output.llmax;
else
	llmax=2; //2;

	// Ylm expansion of the mass perturbation at each boundary - do the Ylm expansion of all_stress:
        if(E->parallel.me==0 && E->advection.timesteps%50==0){
                time = CPU_time0() -time0;
                fprintf(E->fp,"DEBUG (In calculate_potential_comp_OR_deg1_2: !only_deg1_2) : before sphere_expansion_VE , elapsed time = %f \n",time);
               // time = CPU_time0();
        }	
	for (i=1;i<=E->lmesh.noz;i++) {
      // save the ith layer into E->Xsurf[1]
      		
      //>>>> modify by tao: below from geruo
	      for (m=1;m<=E->sphere.caps_per_proc;m++) {
	        for (j=1;j<=E->lmesh.nsf;j++) {
	          node = (j-1)*E->lmesh.noz+i;
	          E->Xsurf[1][m][j] = E->slice_ve.all_stress[m][node];
	        }
	      }
	      // if I(Geruo) wasn't mistaken, it will produce identical results for ic=1 or 0
			/*      if (i!=1)
	        raw_sphere_expansion(E,1,E->Xsurf[1],
	                             E->sphere.sc[i],E->sphere.ss[i]);
	      	else*/

		  // comment by tao: here should use global id on z dir
		  global_z_id = E->lmesh.nzs + i -1;
	        sphere_expansion_VE(E,0,E->Xsurf[1],
	                             E->sphere.sphc_c[global_z_id],E->sphere.sphs_c[global_z_id],llmax);  //comment by tao: here should be global_id
		//<<<< end tao: above from geruo
		//
    }
        if(E->parallel.me==0 && E->advection.timesteps%50==0){
                time = CPU_time0() -time0;
                fprintf(E->fp,"DEBUG (In calculate_potential_comp_OR_deg1_2: !only_deg1_2) : after sphere_expansion_VE , delta time = %f \n",time);
                //time = CPU_time0();
        }	
/*
        if(E->parallel.me_loc[1]==0 && E->parallel.me_loc[2]==0){
            printf("debug calculate potential, hindice=%d\n", E->sphere.hindice);
            int debug_k;
  //          for(debug_k=1; debug_k<=E->mesh.noz; debug_k++)
//                printf("global_r[%d]= %e , from proc %d\n",debug_k, global_r[debug_k], E->parallel.me);

            char debug_outfile[255];
            sprintf(debug_outfile, "debug_sphere_expan_Before_MPI.%d.%d.%d",E->monitor.solution_cycles,E->monitor.count,E->parallel.me);
            FILE * debug_fp;
            debug_fp = fopen(debug_outfile, "w");
            int debug_i, debug_j, debug_global_z;
            for(debug_j=1; debug_j<=E->mesh.noz; debug_j++)
            {
                for(debug_i=0; debug_i<E->sphere.hindice; debug_i++)
                {
                    fprintf(debug_fp,"%e %e %e %e\n", E->sphere.sphc_c[debug_j][debug_i],
                                                      E->sphere.sphs_c[debug_j][debug_i],
                                                      E->sphere.sphc_v[debug_j][debug_i],
                                                      E->sphere.sphs_v[debug_j][debug_i]);
                }
            }
            fclose(debug_fp);

        }
void parallel_process_termination();
parallel_process_termination();
*/
	//>>>> add by tao: here not sure whether this works.
	if (E->parallel.nprocz>1)  // if >1 cap in z-direction (probably not)
		for(i=1; i<=E->mesh.noz; i++) //comment by tao: global z;
			{
			exchange_sphcs_to_all(E,E->sphere.sphc_c, i);
			exchange_sphcs_to_all(E,E->sphere.sphs_c, i);
			}
	

/*	
        if(E->parallel.me_loc[1]==0 && E->parallel.me_loc[2]==0){
            printf("debug calculate potential, hindice=%d\n", E->sphere.hindice);
            int debug_k;
  //          for(debug_k=1; debug_k<=E->mesh.noz; debug_k++)
//                printf("global_r[%d]= %e , from proc %d\n",debug_k, global_r[debug_k], E->parallel.me);

            char debug_outfile[255];
            sprintf(debug_outfile, "debug_sphere_expan_Before_MPI.%d.%d.%d",E->monitor.solution_cycles,E->monitor.count,E->parallel.me);
            FILE * debug_fp;
            debug_fp = fopen(debug_outfile, "w");
            int debug_i, debug_j, debug_global_z;
            for(debug_j=1; debug_j<=E->mesh.noz; debug_j++)
            {
                for(debug_i=0; debug_i<E->sphere.hindice; debug_i++)
                {
                    fprintf(debug_fp,"%e %e %e %e\n", E->sphere.sphc_c[debug_j][debug_i],
                                                      E->sphere.sphs_c[debug_j][debug_i],
                                                      E->sphere.sphc_v[debug_j][debug_i],
                                                      E->sphere.sphs_v[debug_j][debug_i]);
                }
            }
            fclose(debug_fp);

        }
void parallel_process_termination();
parallel_process_termination();
*/		
	//<<<< end tao: above not sure.





	// !!! compute rho*div_u*dV at each gaussian points, then spherical expansion
	
	double *gnxx, *gnda;
	double rtf[4][9],sinaa[9],ct[9],div_u_gp[9],div_u[4][9][9];
	static struct CC Cc;
    static struct CCX Ccx;
	double temp;
	int a, el;
	const int dims=E->mesh.nsd,dofs=E->mesh.dof;
    const int vpts=vpoints[E->mesh.nsd];
    const int ends=enodes[E->mesh.nsd];
	
	void get_rtf_at_vpts();
	// here will get sphere.sphc_v and sphs_v as geruo use sc_v ss_v.

    if (!load_only){
		  for (m=1;m<=E->sphere.caps_per_proc;m++)
            for(i=1;i<=E->lmesh.nno;i++)
              rho_div_u[m][i] = 0.0;

		//if(E->parallel.me==0)
		//{
		/*
		if(E->parallel.me==0)
			printf("debug div_u sphere expan. Part 1: define outputfile\n");
			char debug_outfile1[255];
			sprintf(debug_outfile1, "debug_rho_div_u_ele_Mine.%d.%d.%d",E->parallel.me,E->monitor.solution_cycles,E->monitor.count);
			char debug_outfile2[255];
			sprintf(debug_outfile2, "debug_TWW_Mine.%d.%d.%d",E->parallel.me,E->monitor.solution_cycles,E->monitor.count);
			FILE * debug_fp1, *debug_fp2;
			debug_fp1 = fopen(debug_outfile1, "w");
			debug_fp2 = fopen(debug_outfile2, "w");
		*/
		//}

          for (m=1;m<=E->sphere.caps_per_proc;m++)
          for (el=1;el<=E->lmesh.nel;el++) {
            gnxx = E->gNX[m][el].vpt;
            gnda = E->gDA[m][el].vpt;
            get_rtf_at_vpts(E,m,E->mesh.levmax,el,rtf); // get r,t,f at gaussian points
            if ((el-1)%E->lmesh.elz==0) {
              construct_c3x3matrix_el(E,el,&Cc,&Ccx,E->mesh.levmax,m,0);
            }
            temp = 0.;
            for (j=1;j<=vpts;j++) {
              // compute div_u at gaussian points
              sinaa[j] = sin(rtf[1][j]);
              ct[j] = cos(rtf[1][j])/sinaa[j];
              div_u_gp[j] = 0.0;
              for (a=1;a<=ends;a++) {
                node = E->ien[m][el].node[a];
                for (k=1;k<=dims;k++) {
                  div_u[k][a][j] = rtf[3][j]*
                           (gnxx[GNVXINDEX(0,a,j)]*Cc.vpt[BVINDEX(1,k,a,j)]
                           +E->N.vpt[GNVINDEX(a,j)]*Ccx.vpt[BVXINDEX(1,k,1,a,j)]
                           +E->N.vpt[GNVINDEX(a,j)]*Cc.vpt[BVINDEX(3,k,a,j)])
                                  +rtf[3][j]*
                           (E->N.vpt[GNVINDEX(a,j)]*Cc.vpt[BVINDEX(1,k,a,j)]*ct[j]
                           +E->N.vpt[GNVINDEX(a,j)]*Cc.vpt[BVINDEX(3,k,a,j)]
                           +(gnxx[GNVXINDEX(1,a,j)]*Cc.vpt[BVINDEX(2,k,a,j)]
                           +E->N.vpt[GNVINDEX(a,j)]*Ccx.vpt[BVXINDEX(2,k,2,a,j)])/sinaa[j])
                                  +gnxx[GNVXINDEX(2,a,j)]*Cc.vpt[BVINDEX(3,k,a,j)];
                  div_u_gp[j] +=div_u[k][a][j]*E->U[m][E->id[m][node].doff[k]];
                }  
              }

              temp += div_u_gp[j]*E->erho[E->mesh.levmax][m][el];

/*              rho_div_u_gp[m][(el-1)*vpts+j] = div_u_gp[j]
                                           *E->erho[E->mesh.levmax][m][(el-1)*vpts+j]*gnda[j]*g_point[j].weight[dims-1];
                                           // rho*div_u*dV
              t_gp[m][(el-1)*vpts+j] = rtf[1][j];
              f_gp[m][(el-1)*vpts+j] = rtf[2][j];
              r_gp[m][(el-1)*vpts+j] = 1./rtf[3][j]; */
            } // end j loop
            temp = temp/vpts;
            for(j=1;j<=ends;j++)                {
              node = E->ien[m][el].node[j];
              rho_div_u[m][node] += E->TWW[E->mesh.levmax][m][el].node[j]*temp;
            }
				/*
				//if(E->parallel.me==0)
				//{
				//    printf("debug div_u sphere expan. Part 2: output to file\n");
				//	fprintf(debug_fp1, "%e \n", temp);
					for(j=1;j<=ends;j++)
					{
						node = E->ien[m][el].node[j];
						fprintf(debug_fp2, "%e ",  E->TWW[E->mesh.levmax][m][el].node[j]);
					}
					fprintf(debug_fp2, "\n");

				//}
				*/

          } // end el loop
				/*
				//if(E->parallel.me==0)
				//{
				if(E->parallel.me==0)
				printf("debug div_u sphere expan. Part 3: close files\n");
				fclose(debug_fp1);
				fclose(debug_fp2);
				//}
				*/
/*          // compute Ylm at the guassian points
          for (m=1;m<=E->sphere.caps_per_proc;m++){
          for (ll=1;ll<=E->sphere.output_llmax;ll++)
          for (mm=0;mm<=ll;mm++) {
            for (el=1;el<=E->lmesh.nel;el++){
              for (j=1;j<=vpts;j++){
                ylm_cs[m][ll][mm][(el-1)*vpts+j] = cos(mm*f_gp[m][(el-1)*vpts+j])*modified_plgndr_a(ll,mm,t_gp[m][(el-1)*vpts+j]);
                ylm_sn[m][ll][mm][(el-1)*vpts+j] = sin(mm*f_gp[m][(el-1)*vpts+j])*modified_plgndr_a(ll,mm,t_gp[m][(el-1)*vpts+j]);
              }
            }
          }
          }*/

          (E->exchange_node_d) (E,rho_div_u,E->mesh.levmax);

//add by tao: the TWW definition is different between Geruo's and Shijie's
   		  for(m=1;m<=E->sphere.caps_per_proc;m++)
     		for(node=1;node<=E->lmesh.nno;node++)
        		rho_div_u[m][node] *= E->MASS[E->mesh.levmax][m][node];

	      for (i=1;i<E->lmesh.noz;i++) {
	        // save the ith layer into E->Xsurf[1]
	        for (m=1;m<=E->sphere.caps_per_proc;m++) {
	          for (j=1;j<=E->lmesh.nsf;j++) {
	            node = (j-1)*E->lmesh.noz+i;
	            E->Xsurf[1][m][j] = (rho_div_u[m][node]+rho_div_u[m][node+1])/2.0;
	          }
	        }
			global_z_id = E->lmesh.nzs + i -1;
	       // raw_sphere_expansion(E,0,E->Xsurf[1],E->sphere.sphc_v[i],E->sphere.sphs_v[i]);
		   sphere_expansion_VE(E,0,E->Xsurf[1],
		                             E->sphere.sphc_v[global_z_id],E->sphere.sphs_v[global_z_id], llmax);  
	      }  

		  if (E->parallel.nprocz>1)  // if >1 cap in z-direction (probably not)
		    for(i=1; i<E->mesh.noz; i++) //comment by tao: global z; i < noz
			{
				exchange_sphcs_to_all(E,E->sphere.sphc_v, i);
				exchange_sphcs_to_all(E,E->sphere.sphs_v, i);
			}
    
    	}

		if(0)	
		if(E->sphere.capid==0){
			printf("debug calculate potential, hindice=%d\n", E->sphere.hindice);
			int debug_k;
//			for(debug_k=1; debug_k<=E->mesh.noz; debug_k++)
//				printf("global_r[%d]= %e , from proc %d\n",debug_k, global_r[debug_k], E->parallel.me);

			char debug_outfile[255];
			sprintf(debug_outfile, "debug_sphere_expan_Mine.%d.%d.%d",E->monitor.solution_cycles,E->monitor.count,E->parallel.me);
			FILE * debug_fp;
			debug_fp = fopen(debug_outfile, "w");
			int debug_i, debug_j, debug_global_z;
			for(debug_j=1; debug_j<=E->mesh.noz; debug_j++)
			{	
				for(debug_i=0; debug_i<E->sphere.hindice; debug_i++)
				{
					fprintf(debug_fp,"%e %e %e %e\n", E->sphere.sphc_c[debug_j][debug_i],
													  E->sphere.sphs_c[debug_j][debug_i],
													  E->sphere.sphc_v[debug_j][debug_i],
													  E->sphere.sphs_v[debug_j][debug_i]);
				}
			}
			fclose(debug_fp);	

		}
		
        if(E->parallel.me==0 && E->advection.timesteps%50==0){
                time = CPU_time0() -time0;
                fprintf(E->fp,"DEBUG (In calculate_potential_comp_OR_deg1_2: !only_deg1_2) : after sphere_expansion_VE (div v) , elapsed time = %f \n",time);
               // time = CPU_time0();
        }

	double phi_lm_cs, phi_lm_sn, tmp, dlayer;
	double M[3],centri_pot;

	if(!only_deg1_2){
		
			for (m=1;m<=E->sphere.caps_per_proc;m++) {
		      for (i=1;i<=E->lmesh.nno;i++) {
	    	    all_potential[m][i] = 0.0;
	      	  }
	    	}
	    // Now calculate potential:
	    //>>>> add by tao: all layer:

			for (ll=1;ll<=E->output.llmax;ll++) // degree 0 term is supposed to be zero at surface, but here we are excluding it for each layer
			for (mm=0;mm<=ll;mm++) {
			  p = E->sphere.hindex[ll][mm];
			  
			  for (m=1;m<=E->sphere.caps_per_proc;m++) {
				for (i=1;i<=E->lmesh.noz;i++) {
					  // compute phi_lm(r)
					  phi_lm_cs = 0.0;
					  phi_lm_sn = 0.0;
					  r = E->SX[E->mesh.levmax][m][3][i];
					  global_z_id = E->lmesh.nzs + i -1;

					  // 1. the boundary terms
					  for (k=1;k <= global_z_id;k++) {  //comment by tao: inside layers
						ri = global_r[k]; // the smaller radius // todo_tao: get from global_k
						ro = r; 							 // the larger radius
						tmp = ri*pow(ri/ro,(double)(ll+1))/(2.0*ll+1.0);
						phi_lm_cs += tmp * E->sphere.sphc_c[k][p];
						phi_lm_sn += tmp * E->sphere.sphs_c[k][p];
					  }
					  if (global_z_id < E->mesh.noz) { // comment by tao: outside layers
						for (k=global_z_id+1;k<=E->mesh.noz;k++) {
						  ri = r;
						  ro = global_r[k];
						  tmp = ri*pow(ri/ro,(double)(ll-1))/(2.0*ll+1.0);
						  phi_lm_cs += tmp*  E->sphere.sphc_c[k][p];
						  phi_lm_sn += tmp*  E->sphere.sphs_c[k][p];
						}
					  }

					  // 2. volume div_v   //done tao

					  if(!load_only){
					  //if (i>1)
						for(k=1;k<global_z_id;k++) {
						  ri = (global_r[k] + global_r[k+1])/2.0;
						  ro = r;
						  dlayer = global_r[k+1] - global_r[k];
						  tmp = ri*pow(ri/ro,(double)(ll+1))/(2.0*ll+1.0);
						  phi_lm_cs -= tmp*E->sphere.sphc_v[k][p]*dlayer;
						  phi_lm_sn -= tmp*E->sphere.sphs_v[k][p]*dlayer;
						}
					  //if (i<E->lmesh.noz)
						for(k=global_z_id;k<E->mesh.noz;k++) {
						  ri = r;
						  ro = (global_r[k] + global_r[k+1])/2.0;
						  dlayer = global_r[k+1] - global_r[k];
						  tmp = ri*pow(ri/ro,(double)(ll-1))/(2.0*ll+1.0);
						  phi_lm_cs -= tmp*E->sphere.sphc_v[k][p]*dlayer;
						  phi_lm_sn -= tmp*E->sphere.sphs_v[k][p]*dlayer;
						}
					  }
						//end of div_v

					// comment by tao: polar_wander, copy from geruo
					  if (E->ve_data_cont.polar_wander && E->parallel.me_loc[3] == E->parallel.nprocz-1 && i==E->lmesh.noz && ll==2 && mm==1){
			            M[0] = M_conv * phi_lm_cs * pow(E->SX[E->mesh.levmax][m][3][i],3.0);
			            M[1] = M_conv * phi_lm_sn * pow(E->SX[E->mesh.levmax][m][3][i],3.0);
//						M[0] = E->ve_data_cont.damping * M[0] + (1-E->ve_data_cont.damping) * E->ve_data_cont.PW_incr[0];
//						M[1] = E->ve_data_cont.damping * M[1] + (1-E->ve_data_cont.damping) * E->ve_data_cont.PW_incr[1];
			            E->ve_data_cont.PW_incr[0] = M[0];
			            E->ve_data_cont.PW_incr[1] = M[1];
	                 }
					  
					  // sum 1+2
					   for (j=1;j<=E->lmesh.nsf;j++) {
		            	node = (j-1)*E->lmesh.noz+i;
		            		all_potential[m][node] += (phi_lm_cs*E->Tbl_cs[m][mm][j] + phi_lm_sn*E->Tbl_sn[m][mm][j])
		                                         *E->Tbl_lm[m][p][j];          
		          		}
					 }
				
			  	}
			}

            // !!! only top processors have PW_incr now. should send to all other processors below.
            if(E->parallel.nprocz>1){
                //printf("MPI: before pass M, M[0]= %e, M[1]= %e, proc %d\n", M[0], M[1], E->parallel.me);
                if(E->parallel.me_loc[3] == E->parallel.nprocz-1){
                    //send to E->parallel.me + procz (1 to E->parallel.nprocz)
                    int procz;
                    for(procz=0; procz<E->parallel.nprocz-1; procz++)
                        MPI_Send(M, 2, MPI_DOUBLE, procz, 0, E->parallel.vertical_comm);
                }
                if(E->parallel.me_loc[3] != E->parallel.nprocz-1){
                    MPI_Recv(M,2, MPI_DOUBLE, E->parallel.nprocz-1, 0, E->parallel.vertical_comm, MPI_STATUS_IGNORE);
                    //receive from E->parallel.me - E->parallel.me_loc[3];
                }
                //printf("MPI: finish pass M, M[0]= %e, M[1]=%e, proc %d\n", M[0], M[1], E->parallel.me);
				E->ve_data_cont.PW_incr[0] = M[0];
                E->ve_data_cont.PW_incr[1] = M[1];
            }
			
	        if(E->parallel.me==0 && E->advection.timesteps%50==0){
                	time = CPU_time0() -time0;
                	fprintf(E->fp,"DEBUG (In calculate_potential_comp_OR_deg1_2: !only_deg1_2) : after sum all spherical harms to get nodal all_potential, elapsed time = %f \n",time);
               		// time = CPU_time0();
        	}	
		    if (E->ve_data_cont.polar_wander && icon==1){
		      // compute: -r^2*omega^2*sin(t)*cos(t)*[M[0]cos(f) + M[1]sin(f)]
		      for(m=1;m<=E->sphere.caps_per_proc;m++)
		      for(node=1;node<=E->lmesh.nno;node++){
		        t = E->SX[E->mesh.levmax][m][1][node];
		        f = E->SX[E->mesh.levmax][m][2][node];
		        r = E->SX[E->mesh.levmax][m][3][node];
		        centri_pot = r*r*sin(t)*cos(t)*(M[0]*cos(f)+M[1]*sin(f))*C_conv;
		        all_potential[m][node] += centri_pot;
		      }
		    }


			
			//calculate E->ve_data_cont.potential_vary_PW . what is icon?
			if (E->ve_data_cont.polar_wander && icon==1){

				static int been_here = 0;
				double  potential[3],prod[3];
	    		static double *oldpotl_surf[NCS], *oldpotl_cmb[NCS];
				

				    if (!been_here) {
				        been_here = 1;
				        
				       for (m=1;m<=E->sphere.caps_per_proc;m++)   {
				          oldpotl_surf[m] =(double*)malloc((E->lmesh.nsf+2)*sizeof(double));
				          oldpotl_cmb[m] = (double*)malloc((E->lmesh.nsf+2)*sizeof(double));
				          for (n=1;n<=E->lmesh.nsf;n++)  { 
				             oldpotl_surf[m][n] = 0.0;
				             oldpotl_cmb[m][n] = 0.0;
				        }

				        }
				    }

				    potential[0]=potential[1]=0.0;
				    prod[0]=prod[1]=0.0;
				    for (m=1;m<=E->sphere.caps_per_proc;m++)
				    for (n=1;n<=E->lmesh.nsf;n++)   {
				        t = E->SX[E->mesh.levmax][m][1][n*E->lmesh.NOZ[E->mesh.levmax]];
				        f = E->SX[E->mesh.levmax][m][2][n*E->lmesh.NOZ[E->mesh.levmax]];
						
				        if (E->parallel.me_loc[3]==E->parallel.nprocz-1)   {
								node = (n-1)*E->lmesh.noz + E->lmesh.noz;
				                //potential_surf[m][n] = all_potential[m][node];
				                potential[0] += all_potential[m][node]*all_potential[m][node];
				                temp = all_potential[m][node]-oldpotl_surf[m][n];
				                potential[1] += temp*temp;
								oldpotl_surf[m][n] = all_potential[m][node];
				                }
				        if (E->parallel.me_loc[3]==0)   {
				                //potential_cmb [m][n] = ;
				                node = (n-1)*E->lmesh.noz + 1;
				                potential[0] += all_potential[m][node]*all_potential[m][node];
				                temp = all_potential[m][node]-oldpotl_cmb[m][n];
				                potential[1] += temp*temp;
								oldpotl_cmb[m][n] = all_potential[m][node];
				                }
	     			 }

					MPI_Allreduce(potential,prod,3,MPI_DOUBLE,MPI_SUM,E->parallel.world);
					E->ve_data_cont.potential_vary_PW = sqrt(prod[1]/prod[0]);
				}// end calculate E->ve_data_cont.potential_vary_PW
				
		}// end if !only_deg1_2



	//comment by tao: deal with degree-1 motion. why only one surface point (for each processor) is used in calculate CM? might cause CM different for each processor.
		if (only_deg1_2 ){ 
//		if(0){			
//if(E->parallel.me==0)
//	printf("deg 1 correction\n");
			//// variables related to CM calculation, from geruo
					double cm[3];
					static double Earth_mass;
					double xy_mag, r0, t0, f0, density_jump, t1, f1, rel_cos;
					
					for (i=0;i<3;i++) cm[i] = 0.0;
					// compute Earth's mass
					if (add_incr_CM != 0 || been==1){
					  Earth_mass = 0.0;
					  tmp = 0.0;
	
					  for (m=1;m<=E->sphere.caps_per_proc;m++)
						for (el=1;el<=E->lmesh.nel;el++) {
						  for (j=1;j<=vpts;j++) {
							gnda = E->gDA[m][el].vpt;
							tmp += E->erho[E->mesh.levmax][m][el]*gnda[j]*g_point[j].weight[dims-1];
						  } // end j loop
						} // end el loop		  
					  MPI_Allreduce(&tmp,&Earth_mass,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
					  Earth_mass += 4./3.*M_PI*E->data.density_below/E->data.density*pow(E->sphere.ri,3);
					
					  if(E->parallel.me==0)
							printf("Geruo's Earth mass: %g \n", Earth_mass);
				    /*  
					double surface_g = 9.82;
    	Earth_mass = (surface_g* pow(E->sphere.dradius,2) / E->data.grav_const)/(E->data.density * pow(E->sphere.dradius,3));
					   if(E->parallel.me==0)
							printf("Use mine Earth mass: %g \n", Earth_mass); 
					*/
					}
			//////////

			if(E->parallel.me_loc[3] == E->parallel.nprocz-1){
			// compute deg-1 motion from ll=1 and i=E->lmesh.noz, and then update all_stress
			ll=1;
			for (mm=0;mm<=ll;mm++){
			  p = E->sphere.hindex[ll][mm];
			  for(m=1;m<=E->sphere.caps_per_proc;m++){
				i=E->lmesh.noz; // surface. only one point as this grid. first column of this grid.
				phi_lm_cs = 0.0;
				phi_lm_sn = 0.0;
				r = E->SX[E->mesh.levmax][m][3][i];
				global_z_id = E->lmesh.nzs + i -1;
				for (k=1;k<=global_z_id;k++) {
				  ri = global_r[k]; // the smaller radius
				  ro = r;							   // the larger radius
				  tmp = ri*pow(ri/ro,(double)(ll+1))/(2.0*ll+1.0);
				  phi_lm_cs += tmp*E->sphere.sphc_c[k][p];
				  phi_lm_sn += tmp*E->sphere.sphs_c[k][p];
				}
				if (!load_only){
				  for(k=1;k<global_z_id;k++) {
					  ri = (global_r[k] + global_r[k+1])/2.0;
					  ro = r;
					  dlayer = global_r[k+1] - global_r[k];
					  tmp = ri*pow(ri/ro,(double)(ll+1))/(2.0*ll+1.0);
					  phi_lm_cs -= tmp*E->sphere.sphc_v[k][p]*dlayer;
					  phi_lm_sn -= tmp*E->sphere.sphs_v[k][p]*dlayer;
				  }
				}

				

				if (mm==1) {
				  cm[0] = -phi_lm_cs*sqrt(12.*M_PI)/Earth_mass*pow(E->SX[E->mesh.levmax][m][3][i],2);
				  cm[1] = -phi_lm_sn*sqrt(12.*M_PI)/Earth_mass*pow(E->SX[E->mesh.levmax][m][3][i],2);
				}
				else {
				  cm[2] = phi_lm_cs*sqrt(12.*M_PI)/Earth_mass;
						//			fprintf(stderr,"m = 0: phi_cs = %g, phi_sn = %g\n",phi_lm_cs, phi_lm_sn);
								}
						  }
						}
		//		printf("in proc %d \n", E->parallel.me);
		//		  printf("cm[0] = %g, cm[1] = %g, cm[2] = %g \n", cm[0],cm[1],cm[2]);
		//		  printf("m = 1: phi_cs = %g, phi_sn = %g\n",phi_lm_cs,phi_lm_sn);
				}


			// !!! only top processors have cm now. should send cm to all other vertical processors.
			if(E->parallel.nprocz>1){
				//printf("MPI: before pass cm, %d\n", E->parallel.me);
				if(E->parallel.me_loc[3] == E->parallel.nprocz-1){
					//send to E->parallel.me + procz (1 to E->parallel.nprocz)
					int procz;
					for(procz=0; procz<E->parallel.nprocz-1; procz++)
						MPI_Send(cm, 3, MPI_DOUBLE, procz, 0, E->parallel.vertical_comm);
				}
				if(E->parallel.me_loc[3] != E->parallel.nprocz-1){
					MPI_Recv(cm,3, MPI_DOUBLE, E->parallel.nprocz-1, 0, E->parallel.vertical_comm, MPI_STATUS_IGNORE);
					//receive from E->parallel.me - E->parallel.me_loc[3];
				}
				//printf("MPI: finish pass cm, %d\n", E->parallel.me);
			}
			
			if (add_incr_CM == 1)
			  for (j=0;j<3;j++) E->ve_data_cont.CM_incr[j] = cm[j];
			else if (add_incr_CM == 2)
			  for (j=0;j<3;j++) E->ve_data_cont.CM_incr[j] += cm[j];

if( E->parallel.me==0){
    printf("debug deg-1 motion: in calculate_potential_comp\n");
	printf("cm: %e %e %e\n",cm[0],cm[1],cm[2]);
}
////

///
			
			if (1){
			// convert CM coordinates to spherical (r0,t0,f0):
			xy_mag = sqrt( cm[0]*cm[0] + cm[1]*cm[1] );
			r0 = sqrt( cm[2]*cm[2] + xy_mag*xy_mag );
			if (r0) {
				t0 = acos( cm[2]/r0 );
				f0 = acos( cm[0]/xy_mag );
				if ( cm[1]<0 ) f0 = 2.0*M_PI - f0 ;
			}
			else t0 = f0 = 0.0 ;
		//	  fprintf(stderr,"t,f,r = %f, %f, %f \n",t0*180./M_PI,f0*180./M_PI,r0);
//	if(E->parallel.me==0)
//		printf("-- deg-1 correction: here only work for one core on z direction.\n");	
			// apply additional loads to all_stress
			for (m=1;m<=E->sphere.caps_per_proc;m++)
			for (i=1;i<=E->lmesh.noz;i++) {
			  density_jump = 0.0;
			  if (i==1){ // at the cmb
				//for (j=1;j<=vpts;j++)
				if(E->parallel.me_loc[3]==0)
				  density_jump +=  E->data.density_below/E->data.density-E->erho[E->mesh.levmax][m][i];
				else
				  density_jump = 1. -E->erho[E->mesh.levmax][m][i];
			  }
			  else if (i==E->lmesh.noz){ // at the surface
				//for (j=1;j<=vpts;j++)
				if(E->parallel.me_loc[3]== E->parallel.nprocz-1)
				  density_jump += E->erho[E->mesh.levmax][m][i-1];
				else
					density_jump += E->erho[E->mesh.levmax][m][i-1] -1.0;
			  }
			  else
				//for (j=1;j<=vpts;j++)
				  density_jump += E->erho[E->mesh.levmax][m][i-1] - E->erho[E->mesh.levmax][m][i];
			  //density_jump = density_jump/vpts;
		 // 	if (i==1) density_jump += E->data.density_below/E->data.density;
		
			  for (k=1;k<=E->lmesh.nsf;k++) {
				node = (k-1)*E->lmesh.noz+i;
				t1 = E->SX[E->mesh.levmax][m][1][node];
				f1 = E->SX[E->mesh.levmax][m][2][node];
				rel_cos = -1.0 * (cos(t1)*cos(t0) + sin(t1)*sin(t0)*cos(f0-f1)) ;
				E->Xsurf[1][m][k] = density_jump*r0*rel_cos;
		/*		  E->all_stress[m][node] += density_jump*r0*rel_cos;
				E->Xsurf[1][m][k] = E->all_stress[m][node];*/
			  }
			  if (i!=1)
				remove_average(E,E->Xsurf[1],1);
			  else
				remove_average(E,E->Xsurf[1],0);
		
			  for (k=1;k<=E->lmesh.nsf;k++) {
				node = (k-1)*E->lmesh.noz+i;
				E->slice_ve.all_stress[m][node] += E->Xsurf[1][m][k];
			  
				// remove below line by tao: dont need this because I will do it later when call calculate_potential_comp to calculate potential
				//			E->Xsurf[1][m][k] = E->slice_ve.all_stress[m][node];
			  }
					/*		if (i!=1)
							raw_sphere_expansion(E,1,E->Xsurf[1],E->sphere.sc[i],E->sphere.ss[i]);
						  else*/
				//remove below line by tao: I think here is no need to 
				//	raw_sphere_expansion(E,0,E->Xsurf[1],E->sphere.sc[i],E->sphere.ss[i]);
		
			}
			}
/*        double global_tdot_d_noskip();
        double _temp;
        _temp =  global_tdot_d_noskip(E,E->slice_ve.all_stress,E->slice_ve.all_stress, E->mesh.levmax);
        printf("all_stress(just after shift V): %e; from core %d \n", _temp, E->parallel.me);
		
    
        if(1||E->parallel.me==0 || E->parallel.me==1){
            printf("---debug E->all_stress (just after shift V) (count = %d) -----\n", E->monitor.count);
            char debug_outfile[255];
            sprintf(debug_outfile, "debug_all_stress_Mine.%d.%d.proc%d",E->monitor.solution_cycles,E->monitor.count, E->parallel.me);
            FILE * debug_fp;
            debug_fp = fopen(debug_outfile, "w");
            int debug_i;
            for(debug_i=1; debug_i<=E->lmesh.nno; debug_i++){
                fprintf(debug_fp,"%e\n", E->slice_ve.all_stress[1][debug_i]);
            }
            fclose(debug_fp);
        }
*/
    
	}
		/* remove by tao:
		    // alter the Y_21 potential, if polar_wander is on:
		    if (E->ve_data_cont.polar_wander && icon==1)  
		        // Xsurf[1/2] have already been converted to nondim height of boundary
		//        polar_wander_effects( E, E->Xsurf[1], E->Xsurf[2], 
		//                                 potential_surf, potential_cmb );
		        polar_wander_effects( E, X_surf, X_cmb, 
		                                 potential_surf, potential_cmb );
		*/
			//<<<< end tao
        if(E->parallel.me==0 && E->advection.timesteps%50==0){
                time = CPU_time0() -time0;
                fprintf(E->fp,"DEBUG (In calculate_potential_comp_OR_deg1_2: !only_deg1_2) : End of this subroutine, elapsed time = %f \n",time);
               // time = CPU_time0();
        }
    return;
}

////<<<< end tao: only calculate the phi, the CM effect should be corrected before separately
///*
//void calculate_potential_comp_Geruo(E,all_potential,load_only,add_incr_CM)
//    struct All_variables *E;
//    higher_precision **all_potential;
//    int load_only,add_incr_CM;
///* I probably should add an option for apply_new_loads, if we can simply process degree-1 motion from the degree-1 potential terms.
//   Since in that way, when it is called in apply_new loads, we will only have a surface load as the input for calculate_potential,
//   and the degree-1 motion can be determined after we obtain the potential. */
//{
//    int ib,is,m,ll,ll1,ll2,mm,p,i,j,k,n,node,el,a,m1;
//    double con,r,t,f,rir,ri,ro,modified_plgndr_a(),ra,phi_lm_cs, phi_lm_sn,div_u[4][9][9];
//    double div_u_gp[9], cs_gp[9], sn_gp[9], rtf[4][9], sinaa[9], ct[9];
//    static double Earth_mass;
//    static double *rho_div_u_gp[NCS], *r_gp[NCS], *t_gp[NCS], *f_gp[NCS], *ylm_cs[NCS][100][100], *ylm_sn[NCS][100][100],*rho_div_u[NCS];
//    double density_surf,density_cmb,density_jump;
//    double tmp,v_cs,v_sn,v_total_cs,v_total_sn,dlayer,temp,cm[3],r0,t0,f0,xy_mag,t1,f1,rel_cos;
//    double M[3],centri_pot;
//
//    static int been=0;
//    static double M_conv, C_conv, omega;
//    void  raw_sphere_expansion();
//    void  polar_wander_effects();
//    void  exchange_sphcs();
//    void construct_c3x3matrix_el();
//    void get_rtf();
//    void visc_from_gint_to_nodes();
//    void remove_average();
//
//    static struct CC Cc;
//    static struct CCX Ccx;
//    float *gnxx,*gnda;
//
//
//    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
//    const int vpts=vpoints[E->mesh.nsd];
//    const int ends=enodes[E->mesh.nsd];
//
//   if (been == 0 ) { 
//    for(m=1;m<=E->sphere.caps_per_proc;m++){
//    rho_div_u[m] = (double *) malloc((E->lmesh.nno+2)*vpoints[E->mesh.nsd]*sizeof(double));
///*    rho_div_u_gp[m] = (double *) malloc((E->lmesh.nel+2)*vpoints[E->mesh.nsd]*sizeof(double));
//    r_gp[m]         = (double *) malloc((E->lmesh.nel+2)*vpoints[E->mesh.nsd]*sizeof(double));
//    t_gp[m]         = (double *) malloc((E->lmesh.nel+2)*vpoints[E->mesh.nsd]*sizeof(double));
//    f_gp[m]         = (double *) malloc((E->lmesh.nel+2)*vpoints[E->mesh.nsd]*sizeof(double));
//
//    for (ll=1;ll<=E->sphere.output_llmax;ll++)
//    for (mm=0;mm<=ll;mm++) {
//      ylm_cs[m][ll][mm] = (double *) malloc((E->lmesh.nel+2)*vpoints[E->mesh.nsd]*sizeof(double));
//      ylm_sn[m][ll][mm] = (double *) malloc((E->lmesh.nel+2)*vpoints[E->mesh.nsd]*sizeof(double));
//    }*/
//    }
//
//    omega = 2.0*M_PI/(24.0*3600.0);  // rotation_rate
//    M_conv = (3.0*E->data.grav_const*E->data.density)/(omega*omega*E->control.kf)*5.0*sqrt(4.0*M_PI/15.0); // 3G*rho/(w^2*kf)*5sqrt(4pi/15)
//          // convert non-dim phi into non-dim M
//    C_conv = -omega*omega/(4.0*M_PI*E->data.grav_const*E->data.density);
//          // convert non-dim r^2(m1*cosf+m2*sinf)*sintcost into non-dim centrifugal potential, which can be added into the potential terms;
//          // to dimensionalize the potential, multiply by 4*pi*G*rho*R^2.
//
//    been++;
//    }
//    for (i=0;i<3;i++) cm[i] = 0.0;
//    // compute Earth's mass
//    if (add_incr_CM != 0 || been==1){
//      Earth_mass = 0.0;
//      tmp = 0.0;
//
//      for (m=1;m<=E->sphere.caps_per_proc;m++)
//        for (el=1;el<=E->lmesh.nel;el++) {
//          for (j=1;j<=vpts;j++) {
//            gnda = E->gDA[m][el].vpt;
//            tmp += E->erho[E->mesh.levmax][m][(el-1)*vpts+j]*gnda[j]*g_point[j].weight[dims-1];
//          } // end j loop
//        } // end el loop          
//      MPI_Allreduce(&tmp,&Earth_mass,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
//      Earth_mass += 4./3.*M_PI*E->data.density_below/E->data.density*pow(E->sphere.ri,3);
//    }
//    
//
//    for (m=1;m<=E->sphere.caps_per_proc;m++) {
//      for (i=1;i<=E->lmesh.nno;i++) {
//        all_potential[m][i] = 0.0;
//      }
//    }
//
//
//    // Ylm expansion of the mass perturbation at each boundary - do the Ylm expansion of all_stress:
//    for (i=1;i<=E->lmesh.noz;i++) {
//      // save the ith layer into E->Xsurf[1]
//      for (m=1;m<=E->sphere.caps_per_proc;m++) {
//        for (j=1;j<=E->lmesh.nsf;j++) {
//          node = (j-1)*E->lmesh.noz+i;
//          E->Xsurf[1][m][j] = E->all_stress[m][node];
//        }
//      }
//      // if I wasn't mistaken, it will produce identical results for ic=1 or 0
///*      if (i!=1)
//        raw_sphere_expansion(E,1,E->Xsurf[1],
//                             E->sphere.sc[i],E->sphere.ss[i]);
//      else*/
//        raw_sphere_expansion(E,0,E->Xsurf[1],
//                             E->sphere.sc[i],E->sphere.ss[i]);
//
//    }
//
//    // compute rho*div_u*dV at each gaussian points
//    if (!load_only){
//          for (m=1;m<=E->sphere.caps_per_proc;m++)
//            for(i=1;i<=E->lmesh.nno;i++)
//              rho_div_u[m][i] = 0.0;
//
//          for (m=1;m<=E->sphere.caps_per_proc;m++)
//          for (el=1;el<=E->lmesh.nel;el++) {
//            gnxx = E->gNX[m][el].vpt;
//            gnda = E->gDA[m][el].vpt;
//            get_rtf(E,el,0,rtf,E->mesh.levmax,m); // get r,t,f at gaussian points
//            if ((el-1)%E->lmesh.elz==0) {
//              construct_c3x3matrix_el(E,el,&Cc,&Ccx,E->mesh.levmax,m,0);
//            }
//            temp = 0.;
//            for (j=1;j<=vpts;j++) {
//              // compute div_u at gaussian points
//              sinaa[j] = sin(rtf[1][j]);
//              ct[j] = cos(rtf[1][j])/sinaa[j];
//              div_u_gp[j] = 0.0;
//              for (a=1;a<=ends;a++) {
//                node = E->ien[m][el].node[a];
//                for (k=1;k<=dims;k++) {
//                  div_u[k][a][j] = rtf[3][j]*
//                           (gnxx[GNVXINDEX(0,a,j)]*Cc.vpt[BVINDEX(1,k,a,j)]
//                           +E->N.vpt[GNVINDEX(a,j)]*Ccx.vpt[BVXINDEX(1,k,1,a,j)]
//                           +E->N.vpt[GNVINDEX(a,j)]*Cc.vpt[BVINDEX(3,k,a,j)])
//                                  +rtf[3][j]*
//                           (E->N.vpt[GNVINDEX(a,j)]*Cc.vpt[BVINDEX(1,k,a,j)]*ct[j]
//                           +E->N.vpt[GNVINDEX(a,j)]*Cc.vpt[BVINDEX(3,k,a,j)]
//                           +(gnxx[GNVXINDEX(1,a,j)]*Cc.vpt[BVINDEX(2,k,a,j)]
//                           +E->N.vpt[GNVINDEX(a,j)]*Ccx.vpt[BVXINDEX(2,k,2,a,j)])/sinaa[j])
//                                  +gnxx[GNVXINDEX(2,a,j)]*Cc.vpt[BVINDEX(3,k,a,j)];
//                  div_u_gp[j] +=div_u[k][a][j]*E->U[m][E->id[m][node].doff[k]];
//                }  
//              }
//
//              temp += div_u_gp[j]*E->erho[E->mesh.levmax][m][(el-1)*vpts+j];
//
///*              rho_div_u_gp[m][(el-1)*vpts+j] = div_u_gp[j]
//                                           *E->erho[E->mesh.levmax][m][(el-1)*vpts+j]*gnda[j]*g_point[j].weight[dims-1];
//                                           // rho*div_u*dV
//              t_gp[m][(el-1)*vpts+j] = rtf[1][j];
//              f_gp[m][(el-1)*vpts+j] = rtf[2][j];
//              r_gp[m][(el-1)*vpts+j] = 1./rtf[3][j]; */
//            } // end j loop
//            temp = temp/vpts;
//            for(j=1;j<=ends;j++)                {
//              node = E->ien[m][el].node[j];
//              rho_div_u[m][node] += E->TWW[E->mesh.levmax][m][el].node[j]*temp;
//            }
//
//          } // end el loop
//
///*          // compute Ylm at the guassian points
//          for (m=1;m<=E->sphere.caps_per_proc;m++){
//          for (ll=1;ll<=E->sphere.output_llmax;ll++)
//          for (mm=0;mm<=ll;mm++) {
//            for (el=1;el<=E->lmesh.nel;el++){
//              for (j=1;j<=vpts;j++){
//                ylm_cs[m][ll][mm][(el-1)*vpts+j] = cos(mm*f_gp[m][(el-1)*vpts+j])*modified_plgndr_a(ll,mm,t_gp[m][(el-1)*vpts+j]);
//                ylm_sn[m][ll][mm][(el-1)*vpts+j] = sin(mm*f_gp[m][(el-1)*vpts+j])*modified_plgndr_a(ll,mm,t_gp[m][(el-1)*vpts+j]);
//              }
//            }
//          }
//          }*/
//
//          exchange_node_d (E,rho_div_u,E->mesh.levmax);
//
//
//      for (i=1;i<E->lmesh.noz;i++) {
//        // save the ith layer into E->Xsurf[1]
//        for (m=1;m<=E->sphere.caps_per_proc;m++) {
//          for (j=1;j<=E->lmesh.nsf;j++) {
//            node = (j-1)*E->lmesh.noz+i;
//            E->Xsurf[1][m][j] = (rho_div_u[m][node]+rho_div_u[m][node+1])/2.0;
//          }
//        }
//       raw_sphere_expansion(E,0,E->Xsurf[1],E->sphere.sc_v[i],E->sphere.ss_v[i]);
//      }  
//    }
//
//    // if >1 cap in z-direction (probably not), additional change is needed
//
//    // compute deg-1 motion from ll=1 and i=E->lmesh.noz, and then update all_stress
//    ll=1;
//    for (mm=0;mm<=ll;mm++){
//      p = E->sphere.hindex[ll][mm];
//      for(m=1;m<=E->sphere.caps_per_proc;m++){
//        i=E->lmesh.noz;
//        phi_lm_cs = 0.0;
//        phi_lm_sn = 0.0;
//        r = E->SX[E->mesh.levmax][3][m][i];
//        for (k=1;k<=i;k++) {
//          ri = E->SX[E->mesh.levmax][3][m][k]; // the smaller radius
//          ro = r;                              // the larger radius
//          tmp = ri*pow(ri/ro,(double)(ll+1))/(2.0*ll+1.0);
//          phi_lm_cs += tmp*E->sphere.sc[k][p];
//          phi_lm_sn += tmp*E->sphere.ss[k][p];
//        }
//        if (!load_only){
//          for(k=1;k<i;k++) {
//              ri = (E->SX[E->mesh.levmax][3][m][k] + E->SX[E->mesh.levmax][3][m][k+1])/2.0;
//              ro = r;
//              dlayer = E->SX[E->mesh.levmax][3][m][k+1] - E->SX[E->mesh.levmax][3][m][k];
//              tmp = ri*pow(ri/ro,(double)(ll+1))/(2.0*ll+1.0);
//              phi_lm_cs -= tmp*E->sphere.sc_v[k][p]*dlayer;
//              phi_lm_sn -= tmp*E->sphere.ss_v[k][p]*dlayer;
//          }
//        }
//        if (mm==1) {
//          cm[0] = -phi_lm_cs*sqrt(12.*M_PI)/Earth_mass*pow(E->SX[E->mesh.levmax][3][m][i],2);
//          cm[1] = -phi_lm_sn*sqrt(12.*M_PI)/Earth_mass*pow(E->SX[E->mesh.levmax][3][m][i],2);
//        }
//        else {
//          cm[2] = phi_lm_cs*sqrt(12.*M_PI)/Earth_mass;
////          fprintf(stderr,"m = 0: phi_cs = %g, phi_sn = %g\n",phi_lm_cs, phi_lm_sn);
//        }
//      }
//    }
////    fprintf(stderr,"cm[0] = %g, cm[1] = %g, cm[2] = %g \n", cm[0],cm[1],cm[2]);
////    fprintf(stderr,"m = 1: phi_cs = %g, phi_sn = %g\n",phi_lm_cs,phi_lm_sn);
//
//    if (add_incr_CM == 1)
//      for (j=0;j<3;j++) E->data.CM_incr[j] = cm[j];
//    else if (add_incr_CM == 2)
//      for (j=0;j<3;j++) E->data.CM_incr[j] += cm[j];
//
//    if (1){
//    // convert CM coordinates to spherical (r0,t0,f0):
//    xy_mag = sqrt( cm[0]*cm[0] + cm[1]*cm[1] );
//    r0 = sqrt( cm[2]*cm[2] + xy_mag*xy_mag );
//    if (r0) {
//        t0 = acos( cm[2]/r0 );
//        f0 = acos( cm[0]/xy_mag );
//        if ( cm[1]<0 ) f0 = 2.0*M_PI - f0 ;
//    }
//    else t0 = f0 = 0.0 ;
////    fprintf(stderr,"t,f,r = %f, %f, %f \n",t0*180./M_PI,f0*180./M_PI,r0);
//
//    // apply additional loads to all_stress
//    for (m=1;m<=E->sphere.caps_per_proc;m++)
//    for (i=1;i<=E->lmesh.noz;i++) {
//      density_jump = 0.0;
//      if (i==1) // at the cmb
//        for (j=1;j<=vpts;j++)
//          density_jump +=  E->data.density_below/E->data.density-E->erho[E->mesh.levmax][m][(i-1)*vpts+j];
//      else if (i==E->lmesh.noz) // at the surface
//        for (j=1;j<=vpts;j++)
//          density_jump += E->erho[E->mesh.levmax][m][(i-2)*vpts+j];
//      else
//        for (j=1;j<=vpts;j++)
//          density_jump += E->erho[E->mesh.levmax][m][(i-2)*vpts+j] - E->erho[E->mesh.levmax][m][(i-1)*vpts+j];
//      density_jump = density_jump/vpts;
// //     if (i==1) density_jump += E->data.density_below/E->data.density;
//
//      for (k=1;k<=E->lmesh.nsf;k++) {
//        node = (k-1)*E->lmesh.noz+i;
//        t1 = E->SX[E->mesh.levmax][1][m][node];
//        f1 = E->SX[E->mesh.levmax][2][m][node];
//        rel_cos = -1.0 * (cos(t1)*cos(t0) + sin(t1)*sin(t0)*cos(f0-f1)) ;
//        E->Xsurf[1][m][k] = density_jump*r0*rel_cos;
///*        E->all_stress[m][node] += density_jump*r0*rel_cos;
//        E->Xsurf[1][m][k] = E->all_stress[m][node];*/
//      }
//      if (i!=1)
//        remove_average(E,E->Xsurf[1],1);
//      else
//        remove_average(E,E->Xsurf[1],0);
//
//      for (k=1;k<=E->lmesh.nsf;k++) {
//        node = (k-1)*E->lmesh.noz+i;
//        E->all_stress[m][node] += E->Xsurf[1][m][k];
//        E->Xsurf[1][m][k] = E->all_stress[m][node];
//      }
///*      if (i!=1)
//        raw_sphere_expansion(E,1,E->Xsurf[1],E->sphere.sc[i],E->sphere.ss[i]);
//      else*/
//        raw_sphere_expansion(E,0,E->Xsurf[1],E->sphere.sc[i],E->sphere.ss[i]);
//
//    }
//    }
//
//    // calculate potential
////    for (ll=1;ll<=10;ll++) //VEdeg10.mpi
//    for (ll=1;ll<=E->sphere.output_llmax;ll++) // degree 0 term is supposed to be zero at surface, but here we are excluding it for each layer
//    for (mm=0;mm<=ll;mm++) {
//      p = E->sphere.hindex[ll][mm];
//      for (m=1;m<=E->sphere.caps_per_proc;m++) {
//        for (i=1;i<=E->lmesh.noz;i++) {
//          // compute phi_lm(r)
//          phi_lm_cs = 0.0;
//          phi_lm_sn = 0.0;
//          v_cs      = 0.0;
//          v_sn      = 0.0;
//          r = E->SX[E->mesh.levmax][3][m][i];
//          // 1. the boundary terms
//          for (k=1;k<=i;k++) {
//            ri = E->SX[E->mesh.levmax][3][m][k]; // the smaller radius
//            ro = r;                              // the larger radius
//            tmp = ri*pow(ri/ro,(double)(ll+1))/(2.0*ll+1.0);
//            phi_lm_cs += tmp*E->sphere.sc[k][p];
//            phi_lm_sn += tmp*E->sphere.ss[k][p];
//          }
//          if (i<E->lmesh.noz) {
//            for (k=i+1;k<=E->lmesh.noz;k++) {
//              ri = r;
//              ro = E->SX[E->mesh.levmax][3][m][k];
//              tmp = ri*pow(ri/ro,(double)(ll-1))/(2.0*ll+1.0);
//              phi_lm_cs += tmp*E->sphere.sc[k][p];
//              phi_lm_sn += tmp*E->sphere.ss[k][p];
//            }
//          }
//          // 2. the volumn integral term
//          if (!load_only) {
///* my old way to do this, it is extremely slow !!*/
///*          for (m1=1;m1<=E->sphere.caps_per_proc;m1++)
//          for (el=1;el<=E->lmesh.nel;el++) {
//            for (j=1;j<=vpts;j++) {
//              // compute (r<)^l/(r>)^(l+1) * Ylm  at gaussian points
//              if ( r_gp[m1][(el-1)*vpts+j] <= r) {
//                ri = r_gp[m1][(el-1)*vpts+j];
//                ro = r;
//              }
//              else {
//                ri = r;
//                ro = r_gp[m1][(el-1)*vpts+j];
//              }
//              tmp = pow(ri/ro,(double)(ll))/ro/(2.0*ll+1.0);
//              cs_gp[j] = tmp*ylm_cs[m1][ll][mm][(el-1)*vpts+j];
//              sn_gp[j] = tmp*ylm_sn[m1][ll][mm][(el-1)*vpts+j];
//              // multiply by div_u at gaussian points              
//              v_cs -= cs_gp[j]*rho_div_u_gp[m1][(el-1)*vpts+j];
//              v_sn -= sn_gp[j]*rho_div_u_gp[m1][(el-1)*vpts+j];            
//            } // end j loop
//          } // end el loop          
//
//          MPI_Allreduce(&v_cs,&v_total_cs,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
//          MPI_Allreduce(&v_sn,&v_total_sn,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
//          phi_lm_cs += v_total_cs;
//          phi_lm_sn += v_total_sn;*/
//
//          if (i>1)
//            for(k=1;k<i;k++) {
//              ri = (E->SX[E->mesh.levmax][3][m][k] + E->SX[E->mesh.levmax][3][m][k+1])/2.0;
//              ro = r;
//              dlayer = E->SX[E->mesh.levmax][3][m][k+1] - E->SX[E->mesh.levmax][3][m][k];
//              tmp = ri*pow(ri/ro,(double)(ll+1))/(2.0*ll+1.0);
//              phi_lm_cs -= tmp*E->sphere.sc_v[k][p]*dlayer;
//              phi_lm_sn -= tmp*E->sphere.ss_v[k][p]*dlayer;
//            }
//          if (i<E->lmesh.noz)
//            for(k=i;k<E->lmesh.noz;k++) {
//              ri = r;
//              ro = (E->SX[E->mesh.levmax][3][m][k] + E->SX[E->mesh.levmax][3][m][k+1])/2.0;
//              dlayer = E->SX[E->mesh.levmax][3][m][k+1] - E->SX[E->mesh.levmax][3][m][k];
//              tmp = ri*pow(ri/ro,(double)(ll-1))/(2.0*ll+1.0);
//              phi_lm_cs -= tmp*E->sphere.sc_v[k][p]*dlayer;
//              phi_lm_sn -= tmp*E->sphere.ss_v[k][p]*dlayer;
//            }
//          }
//          // compute the change in rotation axis: m0, m1
//          if (E->control.polar_wander && i==E->lmesh.noz && ll==2 && mm==1){
//            M[0] = M_conv * phi_lm_cs * pow(E->SX[E->mesh.levmax][3][m][i],3.0);
//            M[1] = M_conv * phi_lm_sn * pow(E->SX[E->mesh.levmax][3][m][i],3.0);
//            E->data.rot[0] = M[0];
//            E->data.rot[1] = M[1];
//          }
//          // check the net mass increase
//          if (ll==0 && i==E->lmesh.noz){
//            tmp = phi_lm_cs*2*E->SX[E->mesh.levmax][3][m][i]*sqrt(M_PI);
////            fprintf(stderr,"@@ earth_mass = %g, incr mass = %g, ratio = %g \n",Earth_mass, tmp, tmp/Earth_mass);
//          }
//
//          // sum over ll, mm for each nodes within the ith layer
//          for (j=1;j<=E->lmesh.nsf;j++) {
//            node = (j-1)*E->lmesh.noz+i;
//            all_potential[m][node] += (phi_lm_cs*E->Tbl_cs[m][mm][j] + phi_lm_sn*E->Tbl_sn[m][mm][j])
//                                         *E->Tbl_lm[m][p][j];          
//          }
//        } // end the ith layer loop
//      }
//    }
//    
//
//    // alter the Y_21 potential, if polar_wander is on:
//    if (E->control.polar_wander){
//      // compute: -r^2*omega^2*sin(t)*cos(t)*[M[0]cos(f) + M[1]sin(f)]
//      for(m=1;m<=E->sphere.caps_per_proc;m++)
//      for(node=1;node<=E->lmesh.nno;node++){
//        t = E->SX[E->mesh.levmax][1][m][node];
//        f = E->SX[E->mesh.levmax][2][m][node];
//        r = E->SX[E->mesh.levmax][3][m][node];
//        centri_pot = r*r*sin(t)*cos(t)*(M[0]*cos(f)+M[1]*sin(f))*C_conv;
//        all_potential[m][node] += centri_pot;
//      }
//    }
//
//    return;
//}
//
//*/

void calculate_potential_deg1_2(E,X_surf,X_cmb)
    struct All_variables *E;
    double **X_surf;
    double **X_cmb;
{
    int llmax,ib,is,m,ll,ll1,ll2,mm,p,i,j,k,n;
    double con,r,t,f,rir,ri,ro,modified_plgndr_a();
    double density_surf,density_cmb;

    static int been=0;
    static double brll1[65],brll2[65];
    static double srll1[65],srll2[65];
    void  sphere_expansion_VE();
    void  exchange_sphcs();

    ri = E->sphere.ri;
    ro = E->sphere.ro;
    density_surf = 1.0;
    density_cmb = (E->data.density_below-E->data.density)/E->data.density;

    llmax = 2;
    if (been==0)  {
        been = 1;
        for (ll=1;ll<=llmax;ll++)    {
            srll1[ll] = pow(ro,(double)(ll))/(2.0*ll+1.0); 
            srll2[ll] = ri*pow(ri,(double)(ll+1))/(2.0*ll+1.0); 
        }
    }

    // Ylm expansion of surface mass:
    if (E->parallel.me_loc[3]==E->parallel.nprocz-1)   {
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)   {
            // scale nondim stress -> nondim height of rock
            E->Xsurf[1][m][j] = X_surf[m][j]/E->ve_data_cont.surf_scaling;
        }
        sphere_expansion_VE(E,1,E->Xsurf[1],
                             E->sphere.sphc[0],E->sphere.sphs[0],llmax);
        for (ll=0;ll<=llmax;ll++)
        for (mm=0;mm<=ll;mm++)   {
            p = E->sphere.hindex[ll][mm];
            E->sphere.sphc[0][p] *= density_surf;
            E->sphere.sphs[0][p] *= density_surf;
        }
    }

    // Ylm expansion of CMB mass:
    if (E->parallel.me_loc[3]==0)  {
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)   {
            // scale nondim stress -> nondim height of cmb boundary
            E->Xsurf[2][m][j] = X_cmb[m][j]/E->ve_data_cont.botm_scaling;
        }
        sphere_expansion_VE(E,0,E->Xsurf[2],
                             E->sphere.sphc[1],E->sphere.sphs[1],llmax);
        for (ll=0;ll<=llmax;ll++)
        for (mm=0;mm<=ll;mm++)   {
            p = E->sphere.hindex[ll][mm];
            E->sphere.sphc[1][p] *= density_cmb;
            E->sphere.sphs[1][p] *= density_cmb;
        }
    }

    if (E->parallel.nprocz>1)  // if >1 cap in z-direction (probably not)
        exchange_sphcs(E,E->sphere.sphc[0],E->sphere.sphs[0],
                         E->sphere.sphc[1],E->sphere.sphs[1]);

    // Now calculate potential:
	//comment by tao: seems here sphere.sphc[0][1] is only implemented in top and bottom processors? todo_tao
    ll=1; mm=0;  
    p = E->sphere.hindex[ll][mm]; 
    E->ve_data_cont.CM_pot[0]= srll1[ll]*E->sphere.sphc[0][p]+srll2[ll]*E->sphere.sphc[1][p];

    ll=1; mm=1;  
    p = E->sphere.hindex[ll][mm]; 
    E->ve_data_cont.CM_pot[1]= srll1[ll]*E->sphere.sphc[0][p]+srll2[ll]*E->sphere.sphc[1][p];
    E->ve_data_cont.CM_pot[2]= srll1[ll]*E->sphere.sphs[0][p]+srll2[ll]*E->sphere.sphs[1][p];

    ll=2; mm=1;  
    p = E->sphere.hindex[ll][mm]; 
    E->ve_data_cont.PW_pot[0]= srll1[ll]*E->sphere.sphc[0][p]+srll2[ll]*E->sphere.sphc[1][p];
    E->ve_data_cont.PW_pot[1]= srll1[ll]*E->sphere.sphs[0][p]+srll2[ll]*E->sphere.sphs[1][p];

  return;

  }

/*===========================================================================
 * get_potential(E,count)
 *===========================================================================
 * Sets E->potential[0/1], the nondim grav'l potential at the surface/cmb.
 *               |  pred      (if there is no self-grav iteration)
 *   potential = |  init      (for the first step in a self-grav iteration)
 *               |  init+incr (if in the middle of self-grav interation)
 *=========================================================================== */
void get_potential(struct All_variables *E, int count)
//    struct All_variables *E;
 //   int count;
{
    int ib,is,m,ll,ll1,ll2,mm,p,i,j,k,n;
    void calculate_potential();
    void calculate_potential_comp();
    void parallel_process_termination();
    double total_surface_integral();

/*    double temp1,*TG[4];

  for (m=1;m<=E->sphere.caps_per_proc;m++)
    TG[m] = (double *)malloc((E->lmesh.nsf+1)*sizeof(double));
*/
    if (count==0)   {  // => first (perhaps only) visit here for this timestep
            // We are performing self-grav iteration. So initialize E->potential
            // with the potential at the beginning of the timestep.

			//>>>> modify by tao:
            if (E->ve_data_cont.compressible)
				{
					for(m=1;m<=E->sphere.caps_per_proc;m++)
						for (i=1; i<=E->lmesh.nno; ++i)
						{
							E->slice_ve.all_potential[m][i] = 0.0;
							E->slice_ve.total_potential[m][i] = E->slice_ve.init_total_potential[m][i];
						}
            	}
			else
				{
				for (m=1;m<=E->sphere.caps_per_proc;m++)
            		for (j=1;j<=E->lmesh.nsf;j++)   {
                	E->potential[1][m][j] = E->init_potential[1][m][j];
                	E->potential[0][m][j] = E->init_potential[0][m][j];
            		}
				}			
			//<<<< end tao

    }

    else  {   // second time or later in self-grav iteration

		//>>>> add by tao:
		if (E->ve_data_cont.compressible)
			{
				calculate_potential_comp(E, E->slice_ve.all_potential,0); 

				for(m=1;m<=E->sphere.caps_per_proc;m++)
			  	for (i =1; i <=E->lmesh.nno; ++i)
				  {
					  E->slice_ve.total_potential[m][i] = E->slice_ve.init_total_potential[m][i]
					  				+ E->slice_ve.all_potential[m][i];
				  }

			}		
		//<<<< end tao
		
		else  // comment by tao: here is original incompressible version.
			{
		        // calculate E->incr_potential[0/1] from E->slice.surf/botm[2]:
		        calculate_potential( E, E->slice_ve.surf[2], E->slice_ve.botm[2], 
		                             E->incr_potential[0], E->incr_potential[1], 1);

		        for (m=1;m<=E->sphere.caps_per_proc;m++)
		        for (j=1;j<=E->lmesh.nsf;j++)   {
		            E->potential[1][m][j] =  E->init_potential[1][m][j] 
		                                   + E->incr_potential[1][m][j];  // cmb
		            E->potential[0][m][j] =  E->init_potential[0][m][j] 
		                                   + E->incr_potential[0][m][j];  // surface
		        }
			}
    }

	if(E->ve_data_cont.compressible){
		    // save the surface potential as in the original version; used in dynamic ocean calculation.
		    for (m=1;m<=E->sphere.caps_per_proc;m++)
		      for (j=1;j<=E->lmesh.nsf;j++)   {
		        E->potential[0][m][j] = E->slice_ve.total_potential[m][j*E->lmesh.noz];
		      }
		}

/*
for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (j=1;j<=E->lmesh.nsf;j++)  { 
         TG[m][j] = fabs(E->potential[0][m][j]);
            }
temp1=total_surface_integral(E,TG,1);
if (E->parallel.me==0) fprintf(stderr,"in get_potential %g count %d\n",temp1,count);
    parallel_process_termination();

  free ((void *)TG);
*/

    return;
}

/* ======================================================= 
   for each element el on boundaries, determine contributions
   from incremental displacement. Only the matrix is determined
   here.
 =======================================================  */
void  construct_B_R(E)
  struct All_variables *E;
 {

  void get_global_1d_shape_fn_2();
  void  exchange_id_d();

  const int onedp=onedvpoints[E->mesh.nsd];

  struct Shape_function1 GM;
  struct Shape_function1_dA dGammax;

  double con,slope;
  double temp,x[3][3],force[9],force_at_gs[9];
  int lev,m,ic,nn[5],e,el, i, k, p, p2, lnode[5];

 if (E->parallel.me_loc[3]==E->parallel.nprocz-1 || E->parallel.me_loc[3]==0)  {
 for (lev=E->mesh.levmax;lev>=E->mesh.levmin;lev--)     
 for (m=1;m<=E->sphere.caps_per_proc;m++)   {

   if (E->parallel.me_loc[3]==E->parallel.nprocz-1)  
     for (e=1;e<=E->lmesh.SNEL[lev];e++)    {
       ic = 1;     /*top  */
       el = e*E->lmesh.ELZ[lev];

       get_global_1d_shape_fn_2(E,el,&GM,&dGammax,ic,lev,m);

       for(k=1;k<=onedp;k++)
          nn[k] = k+ic*onedp;

       for(k=1;k<=onedp;k++)
         for (i=1;i<=onedp;i++)
           E->B_R[lev][m][nn[k]][(e-1)*onedp+i] =
                     E->M.vpt[GMVINDEX(k,i)]
                     * dGammax.vpt[GMVGAMMA(0,i)];
       }
   if (E->parallel.me_loc[3]==0)  
     for (e=1;e<=E->lmesh.SNEL[lev];e++)     {
       ic = 0;     /*bottom  */
       el = (e-1)*E->lmesh.ELZ[lev]+1;

       get_global_1d_shape_fn_2(E,el,&GM,&dGammax,ic,lev,m);

       for(k=1;k<=onedp;k++)
          nn[k] = k+ic*onedp;

       for(k=1;k<=onedp;k++)
         for (i=1;i<=onedp;i++)
           E->B_R[lev][m][nn[k]][(e-1)*onedp+i] =
                     E->M.vpt[GMVINDEX(k,i)]
                     * dGammax.vpt[GMVGAMMA(0,i)];
       }
    }
  }

 return;
 }

void  construct_B_RR(struct All_variables *E)
 // struct All_variables *E;
 {

  void get_global_1d_shape_fn_2();
  void  exchange_id_d();

  const int onedp=onedvpoints[E->mesh.nsd];

  struct Shape_function1 GM;
  struct Shape_function1_dA dGammax;

  double con,slope;
  double temp,x[3][3],force[9],force_at_gs[9];
  int lev,m,ic,nn[5],e,el, i, k, p, p2, lnode[5], j;

  if (1)  {
   for (lev=E->mesh.levmax;lev>=E->mesh.levmin;lev--)
     for (m=1;m<=E->sphere.caps_per_proc;m++)   {
       for (e=1;e<=E->lmesh.SNEL[lev];e++)    {
         for (j=1;j<=E->lmesh.ELZ[lev];j++) {
            ic = 1;
            el = (e-1)*E->lmesh.ELZ[lev] + j;
            get_global_1d_shape_fn_2(E,el,&GM,&dGammax,ic,lev,m);
            
            for(k=1;k<=onedp;k++)
              nn[k] = k+ic*onedp;
            for(k=1;k<=onedp;k++)
              for (i=1;i<=onedp;i++)
                E->B_RR[lev][m][nn[k]][(el-1)*onedp+i] =
                     E->M.vpt[GMVINDEX(k,i)]
                     * dGammax.vpt[GMVGAMMA(0,i)];
         }
       }
       if (1 /*E->parallel.me_loc[3]==0*/) {
         for (e=1;e<=E->lmesh.SNEL[lev];e++)     {
           ic = 0;     /*bottom  */
           el = (e-1)*E->lmesh.ELZ[lev]+1;
           get_global_1d_shape_fn_2(E,el,&GM,&dGammax,ic,lev,m);

           for(k=1;k<=onedp;k++)
             nn[k] = k+ic*onedp;
           for(k=1;k<=onedp;k++)
             for (i=1;i<=onedp;i++)
               E->B_RR[lev][m][nn[k]][(el-1)*onedp+i] =
                     E->M.vpt[GMVINDEX(k,i)]
                     * dGammax.vpt[GMVGAMMA(0,i)];
         }
       }
     }
  }
  return;
}

/*===========================================================================
 * add_restoring(E,lev,m,Au,u)
 *===========================================================================*/
void add_restoring(E,lev,m,Au,u)
    struct All_variables *E;
    double **Au,**u;
    int lev,m;
{

    const int onedp=onedvpoints[E->mesh.nsd];

    double temp,x[3][3],force[9],force_at_gs[9];
    int ic,nn[5],el,e, i, k, p, p2, lnode[5];

    if (E->parallel.me_loc[3]==E->parallel.nprocz-1)  {

        for (e=1;e<=E->lmesh.SNEL[lev];e++)   {
            el=e*E->lmesh.ELZ[lev];

            for(k=1;k<=onedp;k++)   {
                nn[k] = E->ID[lev][m][E->IEN[lev][m][el].node[k+onedp]].doff[3];
                force[k] = E->ve_data_cont.surf_scaling*u[m][nn[k]];
            }

            for (i=1;i<=onedp;i++)   {
                force_at_gs[i] = 0.0;
                for(k=1;k<=onedp;k++)
                    force_at_gs[i] += force[k] * E->M.vpt[GMVINDEX(k,i)];
            }

            for (k=1;k<=onedp;k++)
            for (i=1;i<=onedp;i++)
                Au[m][nn[k]] += E->B_R[lev][m][k+onedp][(e-1)*onedp+i] * force_at_gs[i];

        }
    }

    if (E->parallel.me_loc[3]==0)  {

        for (e=1;e<=E->lmesh.SNEL[lev];e++)   {
            el=(e-1)*E->lmesh.ELZ[lev]+1;

            for(k=1;k<=onedp;k++)   {
                nn[k] = E->ID[lev][m][E->IEN[lev][m][el].node[k]].doff[3];
                force[k] = E->ve_data_cont.botm_scaling*u[m][nn[k]];
            }

            for (i=1;i<=onedp;i++)   {
                force_at_gs[i] = 0.0;
                for(k=1;k<=onedp;k++)
                    force_at_gs[i] += force[k] * E->M.vpt[GMVINDEX(k,i)];
            }

            for(k=1;k<=onedp;k++)
            for (i=1;i<=onedp;i++)
                Au[m][nn[k]] += E->B_R[lev][m][k][(e-1)*onedp+i] * force_at_gs[i];

        }
    }

    return;
}

//add by tao: from geruo
void add_restoring_comp(struct All_variables *E, int lev, int m, double **Au, double** u)
//    struct All_variables *E;
//    double **Au,**u;
//    int lev,m;
{

    void get_global_1d_shape_fn_2();
    const int onedp=onedvpoints[E->mesh.nsd];
    const int vpts=vpoints[E->mesh.nsd];

    double temp,x[3][3],force[9],force_at_gs[9];
    int ic,nn[5],el,e, i, k, p, p2, lnode[5];
    int j,n;
    double density_jump, grav, stress_scaling;

    struct Shape_function1 GM;
    struct Shape_function1_dA dGammax;
    static int been=0;


/*    if (E->parallel.me_loc[3]==E->parallel.nprocz-1)  {


        for (e=1;e<=E->lmesh.SNEL[lev];e++)   {
            el=e*E->lmesh.ELZ[lev];
            for(k=1;k<=onedp;k++)   {
                nn[k] = E->ID[lev][m][E->IEN[lev][m][el].node[k+onedp]].doff[3];
                force[k] = E->data.surf_scaling*u[m][nn[k]];
            }

            for (i=1;i<=onedp;i++)   {
                force_at_gs[i] = 0.0;
                for(k=1;k<=onedp;k++)
                    force_at_gs[i] += force[k] * E->M.vpt[GMVINDEX(k,i)];
            }

            for (k=1;k<=onedp;k++)
            for (i=1;i<=onedp;i++)
                Au[m][nn[k]] += E->B_R[lev][m][k+onedp][(e-1)*onedp+i] * force_at_gs[i];
          }

        }
    }
*/
/* change start - add the restoring force at each density discontinuity */

    stress_scaling = E->data.density*E->data.grav_acc*E->sphere.dradius/E->ve_data_cont.shear_mod;

    if (1) { /* I probably need to add an option to identify the proc_loc */
        for (e=1;e<=E->lmesh.SNEL[lev];e++)  {
          for (j=1;j<=E->lmesh.ELZ[lev];j++) {
            el = (e-1)*E->lmesh.ELZ[lev] + j;

            /* find the delta-rho */
            density_jump = 0.0;
            if (j == E->lmesh.ELZ[lev]) { /* at the surface */
//              for (n=1;n<=vpts;n++)
				if(E->parallel.me_loc[3]==E->parallel.nprocz-1)
                	density_jump += E->erho[lev][m][el];
				else
					density_jump += E->erho[lev][m][el] -1.0;
//              density_jump = density_jump/vpts;
            }
            else {
 //             for (n=1;n<=vpts;n++)
                density_jump = density_jump + E->erho[lev][m][el] - E->erho[lev][m][el+1];
 //             density_jump = density_jump/vpts;
            }

//            if (e==10 && E->parallel.me==10 && been==0) fprintf(stderr,"%d  %.8e\n",j+1,density_jump*E->data.density);


            /* find the forces acting on the top nodes of the specific element */
            for (k=1;k<=onedp;k++) {
              nn[k] = E->ID[lev][m][E->IEN[lev][m][el].node[k+onedp]].doff[3];
              grav = E->grav[lev][m][E->IEN[lev][m][el].node[k+onedp]];
            /* always pick the top surface of each element so: k+onedp; doff[3] means the r-direction component*/
              force[k] = density_jump*grav*u[m][nn[k]]*stress_scaling;
            }

            /* find the forces at each gaussian point */
            for (i=1;i<=onedp;i++)   {
                force_at_gs[i] = 0.0;
                for(k=1;k<=onedp;k++)
                    force_at_gs[i] += force[k] * E->M.vpt[GMVINDEX(k,i)];
            }

            /* find the 2d jacobian, I probably should write a subroutine similar to construct_B_R for this purpose
            ic = 1;     
            get_global_1d_shape_fn_2(E,el,&GM,&dGammax,ic,lev,m);
            */

            /* add the forcing terms ( drho*g*u_r, here it only changes the nn(k) terms in Au) to the Au vector */
            for (k=1;k<=onedp;k++)
            for (i=1;i<=onedp;i++)
                  Au[m][nn[k]] += E->B_RR[lev][m][k+onedp][(el-1)*onedp+i] * force_at_gs[i];
//                Au[m][nn[k]] += E->M.vpt[GMVINDEX(k,i)] * dGammax.vpt[GMVGAMMA(0,i)] * force_at_gs[i];

          }

        }
    }

    if ( 1 /*E->parallel.me_loc[3]==0*/)  { /* we also need to consider the restoring forces at the CMB */

        for (e=1;e<=E->lmesh.SNEL[lev];e++)   {
            el=(e-1)*E->lmesh.ELZ[lev]+1;

            /* find the delta-rho */
            density_jump = 0.0;
			if(E->parallel.me_loc[3]==0) {
            	//for (n=1;n<=vpts;n++)
              density_jump = density_jump + E->data.density_below/E->data.density-E->erho[lev][m][el];
            	//density_jump = density_jump/vpts ;
			}
			else{
				density_jump = 1.0 -E->erho[lev][m][el];  // this will be compensated by proc below when MPI sum.
			}
//            if (e==10 && E->parallel.me==10 && been==0) {fprintf(stderr,"%d  %.8e\n",1,density_jump*E->data.density);been++;}

            
            for(k=1;k<=onedp;k++)   {
                nn[k] = E->ID[lev][m][E->IEN[lev][m][el].node[k]].doff[3];
                grav = E->grav[lev][m][E->IEN[lev][m][el].node[k]];
//                force[k] = E->data.botm_scaling*u[m][nn[k]];
                force[k] = density_jump*grav*u[m][nn[k]]*stress_scaling;
            }

            for (i=1;i<=onedp;i++)   {
                force_at_gs[i] = 0.0;
                for(k=1;k<=onedp;k++)
                    force_at_gs[i] += force[k] * E->M.vpt[GMVINDEX(k,i)];
            }


            for(k=1;k<=onedp;k++)
            for (i=1;i<=onedp;i++)
                Au[m][nn[k]] += E->B_RR[lev][m][k][(el-1)*onedp+i] * force_at_gs[i];
//                Au[m][nn[k]] += E->B_R[lev][m][k][(e-1)*onedp+i] * force_at_gs[i];

        }
    }

    return;
}


/* =============================================================================
 * apply_new_loads (E)
 * =============================================================================
 * At the beginning of each timestep, just after setting the current value for
 * iceload (in get_iceModel, for example), this routine will prepare the following
 * for the stokes solver:
 *     slice.load[0/2]
 *     init_potential[1/0]
 * Also, the degree-1 iceload is changed via load_to_CM.
 * About the flag E->control.change_of_load:
 *     0 => incr loads are same as last timestep (eg, within an epoch)
 *     1 => incr loads are changing from prev timestep
 *     2 => incr loads are all zero (eg, the last 4kyrs, post-glacial)
 * ===========================================================================*/

void apply_new_loads(struct All_variables * E)  
 //   struct All_variables *E;
{
    int m,i,node;
    void calculate_potential();
    void calculate_potential_comp();
    double total_surface_integral();
    double temp1,temp2;
    void parallel_process_termination();
    void load_to_CM();
    void load_to_CM_grav();
	void load_to_CM_grav_comp();

    if (E->ve_data_cont.change_of_load==0) {  // no more ice or static ocean load
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (i=1;i<=E->lmesh.nsf;i++)   {
            E->incr_potential[2][m][i] = 0.0 ;
            E->incr_potential[3][m][i] = 0.0 ;
        }
		//>>>> add by tao:
		if (E->ve_data_cont.compressible){
			for(m=1;m<=E->sphere.caps_per_proc;m++)
			for (i=1; i<=E->lmesh.nno; ++i)
				{
					E->slice_ve.all_load_potential[m][i] = 0.0;
				}
		}
		//<<<< end tao
		

        return;
    }

    if (E->ve_data_cont.Heaviside==2) {  // only for ice-model type loading

        // Put loads into slice.surf/botm[2]:
        if(E->ve_data_cont.compressible){
	        // reset all_stress to zero
	        double stress_scaling = E->data.density*E->data.grav_acc*E->sphere.dradius/E->ve_data_cont.shear_mod;
	        for (m=1;m<=E->sphere.caps_per_proc;m++)
	          for (i=1;i<=E->lmesh.nno;i++)
	            E->slice_ve.all_stress[m][i] = 0.0;

	        // Put loads into the surface layer of all_stress:
	        if (E->parallel.me_loc[3]==E->parallel.nprocz-1)
	            for (m=1;m<=E->sphere.caps_per_proc;m++)
	            for (i=1;i<=E->lmesh.nsf;i++){
	              node = i*E->lmesh.noz;
	              E->slice_ve.all_stress[m][node] = (E->slice_ve.iceload[0][m][i]
	                                        + E->slice_ve.static_oceanload[m][i])/stress_scaling; /* convert to non-dim rho*u from non-dim stress */
	            }
		}
		else{
	        if (E->parallel.me_loc[3]==E->parallel.nprocz-1)
	            for (m=1;m<=E->sphere.caps_per_proc;m++)  
	            for (i=1;i<=E->lmesh.nsf;i++)
	                E->slice_ve.surf[2][m][i] =  E->slice_ve.iceload[0][m][i]
	                                        + E->slice_ve.static_oceanload[m][i];
	        if (E->parallel.me_loc[3]==0)
	            for (m=1;m<=E->sphere.caps_per_proc;m++)  
	            for (i=1;i<=E->lmesh.nsf;i++)
	                E->slice_ve.botm[2][m][i] = 0.0;
			}
		
        // Adjust loads to prevent CM motion (the 1 flag will overwrite the
        // current data.CM_incr, starting fresh for this timestep):
        if(E->ve_data_cont.compressible){
			load_to_CM_grav_comp(E, 1, 1, 1); // load_only=1; icon=1; add_incr_CM=1, which means overwrite data.CM_incr
    		//{// add by tao: follow geruo
				// Dec 1 2023 (tao): don't shift total_V here. Instead, deform U (incremental) in deform_grid, where we not add this to the new CM_incr when change_of_load==1.
        		//	void shift_V_to_CM();
				// if(E->parallel.me==0)
    			//		printf("shift_V_to_CM called\n");
        		//	shift_V_to_CM(E);
    		//}
			for(i=0;i<3;i++)
				E->ve_data_cont.CM_incr_ice_static_ocean[i] = E->ve_data_cont.CM_incr[i]; // record ice and static ocean induced CM motion.			
		}
		else{
			load_to_CM_grav( E, E->slice_ve.surf[2], E->slice_ve.botm[2], 1 );
		}

		// comment by tao: update total_load, which is delta_rho * g * Ur + ice and ocean loads. the delta_rho*g*Ur is add in deform_grid
		if(E->ve_data_cont.compressible){
			double stress_scaling = E->data.density*E->data.grav_acc*E->sphere.dradius/E->ve_data_cont.shear_mod;
		
		   if (E->ve_data_cont.compressible)
			for (m=1;m<=E->sphere.caps_per_proc;m++)
			  for (i=1;i<=E->lmesh.nno;i++){
				  E->slice_ve.total_load[m][i] +=  E->slice_ve.all_stress[m][i]*E->grav[E->mesh.levmax][m][i]*stress_scaling; // here we need a non-dim stress
			  }

			}
		else{
	        // Add these loads into slice.load[0/2] (which is used in
	        // assemble_forces to create force vector):
	        if (E->parallel.me_loc[3]==E->parallel.nprocz-1)
	            for (m=1;m<=E->sphere.caps_per_proc;m++)  
	            for (i=1;i<=E->lmesh.nsf;i++)
	                E->slice_ve.load[0][m][i] +=  E->slice_ve.surf[2][m][i];
	        if (E->parallel.me_loc[3]==0)
	            for (m=1;m<=E->sphere.caps_per_proc;m++)  
	            for (i=1;i<=E->lmesh.nsf;i++)
	                E->slice_ve.load[2][m][i] += E->slice_ve.botm[2][m][i];
			}

    } // end if ice-model

// add by tao: follow geruo, I think here should include the load-induced CM correction
/*
	if(E->ve_data_cont.compressible)
	{
		void shift_V_to_CM();
if(E->parallel.me==0)
	printf("shift_V_to_CM called\n");	
		shift_V_to_CM(E);
	}
*/
// end add by tao
	
	//>>>> modify by tao: calculate the load induced potential change,
		// should let the variable for incre_potential be zero before using those vars.
		// var: compressible version: all_potential.
    if (E->ve_data_cont.SELFG)   {

        // if the loads have changed, find their new incremental potential due
        // to incremental loads from ice and static ocean loads:
        if (E->ve_data_cont.change_of_load==1 || E->ve_data_cont.SLE)
            // Note: if we're using the SLE, the static_oceanload does not
            // change in a uniform way (due to the time-dependant ocean
            // function).  So we have to recalculate the incr_potential every
            // timestep.


		if (E->ve_data_cont.compressible)
		  {
		  calculate_potential_comp(E, E->slice_ve.all_load_potential,1); /* */
		  }
	  	else	/*comment by tao: incompressible*/
		  {
		  calculate_potential( E, E->slice_ve.surf[2], E->slice_ve.botm[2],
                       E->incr_potential[2], E->incr_potential[3], 1 );
		  }

/*        double global_tdot_d();
        double _temp;
        _temp =  global_tdot_d(E,E->slice_ve.all_load_potential,E->slice_ve.all_load_potential, E->mesh.levmax);
        printf("all_load_potential: %e; from core %d \n", _temp, E->parallel.me);
*/		
/*
if(E->parallel.me==0){
    printf("----- debug potential (in apply_new_loads) ------\n");
    int m, i_debug, node_debug, eqn_debug,level_debug, ele;
    for(m=1;m<=E->sphere.caps_per_proc;m++)
    for(i_debug=1; i_debug<=E->lmesh.noz;i_debug++){
    ele=i_debug;
    printf("node: %d, all_load_potential: %e \n", i_debug, E->slice_ve.all_load_potential[m][i_debug]);
    //printf("node: %d, all_stress: %e \n", i_debug, E->slice_ve.all_stress[m][i_debug]);
    }
}
*/
        // Add load's incr_potenial to init_ and pred_:

	      /* comment by tao: update init_total_potential with all_potential*/
	    if (E->ve_data_cont.compressible)
		  {
		  for(m=1;m<=E->sphere.caps_per_proc;m++)
			for (i =1; i <=E->lmesh.nno; ++i)
				{
					E->slice_ve.init_total_potential[m][i] += E->slice_ve.all_load_potential[m][i];
				}
				/*
				if(E->parallel.me==0){
				int node, debug_i;
				printf("---- debug init_total_potential ----\n");
				for(node=1; node<=E->lmesh.noz; node++){
					printf("node: %d, all_stress: %e \n", node, E->slice_ve.all_stress[1][node]);
					printf("node: %d, init_total_potential: %e\n",node,E->slice_ve.init_total_potential[1][node]);
				}
				}
				*/

		  }
	    else /* comment by tao: incompressible*/
		  {
			for (m=1;m<=E->sphere.caps_per_proc;m++)
			for (i=1;i<=E->lmesh.nsf;i++)	{
				E->init_potential[0][m][i] += E->incr_potential[2][m][i];
				E->init_potential[1][m][i] += E->incr_potential[3][m][i];
			}

		  }
    }     // end if self-gravitation

//>>>> add by tao: I think here I also need update the PW
    if(E->parallel.me==0) fprintf(E->fp,"Question in apply_new_loads: should I update PW[0/1] with PW_incr here?\n");
	//if (E->ve_data_cont.polar_wander) {
	  // E->ve_data_cont.PW[0] += E->ve_data_cont.PW_incr[0];
	  // E->ve_data_cont.PW[1] += E->ve_data_cont.PW_incr[1];
	//	  if (E->parallel.me==0) fprintf(E->fp,"m0m1 in process_new_ve %g %g\n",E->ve_data_cont.PW_incr[0],E->ve_data_cont.PW_incr[1]);
	  // }
//<<<< end polar wander

    //<<<< end tao

    return;
}


/* =============================================================================
 * get_iceModel(E,iceload)  
 * =============================================================================
 * Specifically for iceModel loading.
 * Puts new load into iceload[m][n] (n = surface node). Specifically, 
 * iceload is the incremental (ie, per timestep) nondimensional stress, 
 * after removing its average.
 * Only sets (passed in) array 'iceload'.
 * About the flag E->control.change_of_load:
 *     0 => incr loads are same as last timestep (eg, within an epoch)
 *     1 => incr loads are changing from prev timestep
 *     2 => incr loads are all zero (eg, the last 4kyrs, post-glacial)
 * ========================================================================== */

void get_iceModel(E,iceload)  
    struct All_variables *E;
    double **iceload;
{
    FILE *fp;
    int m,i,j,n;
    static int ifile;  // this is the iceModel 'epoch'
    static int been_here=0;
    static int step_prev=0;
    void sphere_expansion_output();
    void parallel_process_termination();
    void read_reg_grids();
    void remove_average();
    double total,temp1, temp2;
    char outfile[255],input_s[200];

	void get_iceModel_withInitialLoad(struct All_variables * ,double ** );
	
	if(E->ve_data_cont.NonZero_Init_ICELoad == 1){
		get_iceModel_withInitialLoad(E,iceload);
		return;
	}
	
    if (been_here==0)  {
        // Allocate and initialize some arrays:
        for (m=1;m<=E->sphere.caps_per_proc;m++)   {
            E->slice_ve.ice_height_prev[m] = (double *)
                              malloc((E->lmesh.nsf+2)*sizeof(double));
            E->slice_ve.ice_height_curr[m] = (double *)
                              malloc((E->lmesh.nsf+2)*sizeof(double));
            for (n=1;n<=E->lmesh.nsf;n++)  {  
                E->slice_ve.ice_height_prev[m][n] = 0.0;
                E->slice_ve.ice_height_curr[m][n] = 0.0;
            }
        } // end allocate & initialize

    }

   if (E->ve_data_cont.DIRECT == 0 && been_here!=0) 
        return;

    been_here ++;

        // Read in the iceload for each epoch. 

    ifile = E->ve_data_cont.stage + 1;

    // The following flag will force a re-calculation of load potential in
    // apply_new_loads:
    E->ve_data_cont.change_of_load = 1;

    if (ifile==1)  {       // for the first epoch 
        sprintf(outfile,"%s%dx%d.%d",
            E->ve_data_cont.ice_file, E->sphere.elx, E->sphere.ely, ifile-1);
        read_reg_grids(E,outfile,E->slice_ve.ice_height_prev);
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)  {
           E->slice_ve.ice_height_prev[m][j] = E->slice_ve.ice_height_prev[m][j]/E->sphere.dradius; // non-dimensionalized by the Earth's radius
           E->slice_ve.ice_height_prev[m][j] = 0.0; // reset to be ice free
           }
        step_prev = 0;
        }
    else if (ifile>1)  {       // if we're beyond the first epoch, curr->prev
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)  
            E->slice_ve.ice_height_prev[m][j] = E->slice_ve.ice_height_curr[m][j];   
        step_prev = E->ve_data_cont.stages_step[E->ve_data_cont.stage-1];
        }

    // read new epoch data into slice.ice_height_curr
    //
    sprintf(outfile,"%s%dx%d.%d",
                E->ve_data_cont.ice_file, E->sphere.elx, E->sphere.ely, ifile);
    read_reg_grids(E,outfile,E->slice_ve.ice_height_curr);

    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (j=1;j<=E->lmesh.nsf;j++)  {
        E->slice_ve.ice_height_curr[m][j] = E->slice_ve.ice_height_curr[m][j]/E->sphere.dradius; // non-dimensionalized by the Earth's radius
    }
    
    // temp1 is height->stress scaling / number of timesteps in this epoch
    temp1 = E->ve_data_cont.ice_stress_scale/(E->ve_data_cont.stages_step[E->ve_data_cont.stage] - step_prev);

    // iceload now becomes incremental (ie, per timestep) nondimensional stress
    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (j=1;j<=E->lmesh.nsf;j++) 
        iceload[m][j] = ( E->slice_ve.ice_height_curr[m][j]
                         -E->slice_ve.ice_height_prev[m][j])*temp1;

    // remove average and save the spherical harmonic coefficients:
    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (n=1;n<=E->lmesh.nsf;n++)    
        E->Xsurf[3][m][n] = iceload[m][n];
    remove_average(E,E->Xsurf[3],1);
    sprintf(outfile,"iceload_stage%d",ifile);

    sphere_expansion_output( E, 1, E->Xsurf[3], 
                             E->sphere.sphc[0], E->sphere.sphs[0],
                             E->monitor.solution_cycles,outfile );
    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (n=1;n<=E->lmesh.nsf;n++)
        iceload[m][n] = E->Xsurf[3][m][n];

    return;
}

// version of get_iceModel: have initial non-zero loading.
void get_iceModel_withInitialLoad(struct All_variables * E,double ** iceload)  
//    struct All_variables *E;
//    double **iceload;
{
    FILE *fp;
    int m,i,j,n;
    static int ifile;  // this is the iceModel 'epoch'
    static int been_here=0;
    static int step_prev=0;
    void sphere_expansion_output();
    void parallel_process_termination();
    void read_reg_grids();
    void remove_average();
    double total,temp1, temp2;
    char outfile[255],input_s[200];

    if (been_here==0)  {
        // Allocate and initialize some arrays:
        for (m=1;m<=E->sphere.caps_per_proc;m++)   {
            E->slice_ve.ice_height_prev[m] = (double *)
                              malloc((E->lmesh.nsf+2)*sizeof(double));
            E->slice_ve.ice_height_curr[m] = (double *)
                              malloc((E->lmesh.nsf+2)*sizeof(double));
            for (n=1;n<=E->lmesh.nsf;n++)  {  
                E->slice_ve.ice_height_prev[m][n] = 0.0;
                E->slice_ve.ice_height_curr[m][n] = 0.0;
            }
        } // end allocate & initialize

		ifile = 0;
		sprintf(outfile,"%s%dx%d.%d",
            E->ve_data_cont.ice_file, E->sphere.elx, E->sphere.ely, ifile);
        read_reg_grids(E,outfile,E->slice_ve.ice_height_prev);
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)  {
           E->slice_ve.ice_height_curr[m][j] = E->slice_ve.ice_height_prev[m][j]/E->sphere.dradius; // non-dimensionalized by the Earth's radius
           E->slice_ve.ice_height_prev[m][j] = 0.0; 
           }
        //step_prev = 0;
		E->ve_data_cont.change_of_load = 1;

		// set elapsed time and solution_cycle back to step 0, since the 1st epoch has non-zero load
        E->monitor.elapsed_time=0;
        E->monitor.solution_cycles = 0;
		
        E->advection.timestep=E->ve_data_cont.stages_timestep[0];


		for (m=1;m<=E->sphere.caps_per_proc;m++)
	    for (j=1;j<=E->lmesh.nsf;j++) 
	        iceload[m][j] = E->slice_ve.ice_height_curr[m][j];

	    // remove average and save the spherical harmonic coefficients:
	    for (m=1;m<=E->sphere.caps_per_proc;m++)
	    for (n=1;n<=E->lmesh.nsf;n++)    
	        E->Xsurf[3][m][n] = iceload[m][n];
	    remove_average(E,E->Xsurf[3],1);
	    sprintf(outfile,"iceload_stage%d",ifile);

	    sphere_expansion_output( E, 1, E->Xsurf[3], 
	                             E->sphere.sphc[0], E->sphere.sphs[0],
	                             E->monitor.solution_cycles,outfile );
	    for (m=1;m<=E->sphere.caps_per_proc;m++)
	    for (n=1;n<=E->lmesh.nsf;n++)
	        iceload[m][n] = E->Xsurf[3][m][n];
		
		been_here ++;

	    return;
    }

   if (E->ve_data_cont.DIRECT == 0 && been_here!=0) 
        return;

    

        // Read in the iceload for each epoch. 

    ifile = E->ve_data_cont.stage + 1;

    // The following flag will force a re-calculation of load potential in
    // apply_new_loads:
    E->ve_data_cont.change_of_load = 1;

    if (ifile==1)  {       // for the first epoch 
        sprintf(outfile,"%s%dx%d.%d",
            E->ve_data_cont.ice_file, E->sphere.elx, E->sphere.ely, ifile-1);
        read_reg_grids(E,outfile,E->slice_ve.ice_height_prev);
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)  {
           E->slice_ve.ice_height_prev[m][j] = E->slice_ve.ice_height_prev[m][j]/E->sphere.dradius; // non-dimensionalized by the Earth's radius
           E->slice_ve.ice_height_prev[m][j] = 0.0; // reset to be ice free
           }
        step_prev = 0;
        }
    else if (ifile>1)  {       // if we're beyond the first epoch, curr->prev
        for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (j=1;j<=E->lmesh.nsf;j++)  
            E->slice_ve.ice_height_prev[m][j] = E->slice_ve.ice_height_curr[m][j];   
        step_prev = E->ve_data_cont.stages_step[E->ve_data_cont.stage-1];
        }

    // read new epoch data into slice.ice_height_curr
    //
    sprintf(outfile,"%s%dx%d.%d",
                E->ve_data_cont.ice_file, E->sphere.elx, E->sphere.ely, ifile);
    read_reg_grids(E,outfile,E->slice_ve.ice_height_curr);

    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (j=1;j<=E->lmesh.nsf;j++)  {
        E->slice_ve.ice_height_curr[m][j] = E->slice_ve.ice_height_curr[m][j]/E->sphere.dradius; // non-dimensionalized by the Earth's radius
    }
    
    // temp1 is height->stress scaling / number of timesteps in this epoch
    temp1 = E->ve_data_cont.ice_stress_scale/(E->ve_data_cont.stages_step[E->ve_data_cont.stage] - step_prev);

    // iceload now becomes incremental (ie, per timestep) nondimensional stress
    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (j=1;j<=E->lmesh.nsf;j++) 
        iceload[m][j] = ( E->slice_ve.ice_height_curr[m][j]
                         -E->slice_ve.ice_height_prev[m][j])*temp1;

    // remove average and save the spherical harmonic coefficients:
    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (n=1;n<=E->lmesh.nsf;n++)    
        E->Xsurf[3][m][n] = iceload[m][n];
    remove_average(E,E->Xsurf[3],1);
    sprintf(outfile,"iceload_stage%d",ifile);

    sphere_expansion_output( E, 1, E->Xsurf[3], 
                             E->sphere.sphc[0], E->sphere.sphs[0],
                             E->monitor.solution_cycles,outfile );
    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (n=1;n<=E->lmesh.nsf;n++)
        iceload[m][n] = E->Xsurf[3][m][n];

    return;
}


////////////////////////////////////////////////////////
 void read_reg_grids_old(E,outfile,field)
   struct All_variables *E;
   char *outfile;
   double **field;   
  {

char input_s[1000];
FILE *fp2;

void  parallel_process_termination();
int numtheta,numphi,nregnodes;
int kk,j;
int ntheta,nphi,iregnode,iregnode2;
int node;
int c1,c2,c3,c4;
int isurface_element,lev;

double del_degree;
double del_theta;
double del_phi;
double lat_min,lat_max;
double theta_min,theta_max;
double phi_min,phi_max;
double vphi;
double rlat,rlong;
double theta,phi;
double phiprime,thetaprime;
double xi,eta;
double shape1,shape2,shape3,shape4;

static double *regular_f;
static int been_here=0;

/* INPUT FILE SPECIFIC STUFF */
/* format longitude, latitude, vphi, vtheta (vtheta is opposite in convention) */
/* file sequence given in specific order (see below):                          */
/* inputfile specific parameters:                                              */

 if (E->sphere.nox==181) {
   del_degree=1.0;           /* increments given in input file */
   lat_min=-89.5;            /* minimum lattitude */
   lat_max=89.5;             /* maximum lattitude */
   phi_min=0.5;             /* minimum longitude */
   phi_max=359.5;           /* maximum longitude */
   numtheta=180;             /* number of theta increments */
   numphi=360;               /* number of phi increments */
  }
 else if (E->sphere.nox==361) {
   del_degree=0.5;           /* increments given in input file */
   lat_min=-89.75;            /* minimum lattitude */
   lat_max=89.75;             /* maximum lattitude */
   phi_min=0.25;             /* minimum longitude */
   phi_max=359.75;           /* maximum longitude */
   numtheta=360;             /* number of theta increments */
   numphi=720;               /* number of phi increments */
  }
 
 
   del_theta=del_degree*M_PI/180.0;
   del_phi=del_degree*M_PI/180.0;
   theta_min=-1.0*(lat_max-90.0)*M_PI/180.0;
   theta_max=-1.0*(lat_min-90.0)*M_PI/180.0;
   phi_min=phi_min*M_PI/180.0;
   phi_max=phi_max*M_PI/180.0;

   nregnodes=numtheta*numphi;


/* adjust parameters for ghosting - extra regular elements are added on all sides */
/* This is required to interpolate nodes between thetamin and thetamax, and phimin and phimax */
/* Domain is padded by extra +2 numphi and +2 numtheta */

   numtheta=numtheta+2;
   numphi=numphi+2;
   theta_min=theta_min-del_theta;
   theta_max=theta_max+del_theta;
   phi_min=phi_min-del_phi;
   phi_max=phi_max+del_phi;

   nregnodes=numtheta*numphi;

 if (been_here==0)  {
   regular_f = (double *)malloc((nregnodes+1)*sizeof(double));
   been_here = 1;
   }


/* read from file */

   if ( (fp2=fopen(outfile,"r"))==NULL)   {
      fprintf(stderr,"ERROR(read_reg_grids)- file %s not found\n",outfile);
      fflush(stderr);
      parallel_process_termination();
   }

   for (kk=1;kk<=nregnodes;kk++)
      regular_f[kk]=-999999;

/* input file specific order of reading  right here */
/* (note - leaving room for ghosting) */

//   for (ntheta=(numtheta-1);ntheta>1;ntheta--)
   for (ntheta=2;ntheta<=(numtheta-1);ntheta++)
     for (nphi=2;nphi<=(numphi-1);nphi++)     {

      iregnode=ntheta+(nphi-1)*numtheta;

      if ((iregnode>nregnodes)||(iregnode<1))   {
        fprintf(stderr,"ERROR(plate velocities)-wrong iregnode %d %d\n",iregnode,nregnodes);
        fflush(stderr);
        exit(10);
      }

      fgets(input_s,1000,fp2);
      sscanf(input_s,"%lf %lf %le",&rlong,&rlat,&vphi);

     regular_f[iregnode]=vphi;

   }

   fclose(fp2);

/* ghost sides of regular domain */

   nphi=1;
   for (ntheta=2;ntheta<=(numtheta-1);ntheta++)  {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=ntheta+((numphi-1)-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }

   nphi=numphi;
   for (ntheta=2;ntheta<=(numtheta-1);ntheta++)   {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=ntheta+(2-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }

   ntheta=1;
   for (nphi=1;nphi<=(numphi/2);nphi++)   {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=(ntheta+1)+(nphi+numphi/2-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }
   for (nphi=(numphi/2)+1;nphi<=numphi;nphi++)   {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=(ntheta+1)+(nphi-(numphi/2)-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }

   ntheta=(numtheta);
   for (nphi=1;nphi<=(numphi/2);nphi++)   {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=(ntheta-1)+(nphi+numphi/2-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }
   for (nphi=(numphi/2)+1;nphi<=numphi;nphi++)  {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=(ntheta-1)+(nphi-(numphi/2)-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }

/* interpolate real mesh */

   isurface_element=0;
   for (j=1;j<=E->sphere.caps_per_proc;j++)
   for (node=1;node<=E->lmesh.nno;node++)   {
     theta=E->sx[j][1][node];
     phi  =E->sx[j][2][node];
     if (node%E->lmesh.noz==0)  {
        isurface_element++;
/* find position on regular mesh */

        ntheta=((theta-theta_min)/del_theta)+1;
        nphi=((phi-phi_min)/del_phi)+1;

/* iregnode is the lower left node in the regular element */

        iregnode=ntheta+(nphi-1)*numtheta;

/* find theta and phi distance from node */

        thetaprime=theta-(theta_min+(ntheta-1)*del_theta);
        phiprime=phi-(phi_min+(nphi-1)*del_phi);

        if (thetaprime<0.0 || thetaprime>del_theta)  {
           fprintf(stderr,"ERROR(get plate velocities)-wierd thetaprime %f %f\n",thetaprime,del_theta);
           fflush(stderr);
           exit(10);
        }
        if (phiprime<0.0 || phiprime>del_phi)   {
           fprintf(stderr,"ERROR(get plate velocities)-wierd phiprime %f %f\n",phiprime,del_phi);
           fflush(stderr);
           exit(10);
        }

/* determine shape functions */

       xi=2.0/del_theta*thetaprime-1.0;
       eta=2.0/del_phi*phiprime-1.0;

       if (xi<(-1-1e-6) || xi > (1+1e-6) || eta < (-1-1e-6) || eta > (1+1e-6))  {
          fprintf(stderr,"bad local variables %f %f (%f %f) (%f %f)\n",xi,eta,thetaprime,phiprime,del_theta,del_phi);
          fflush(stderr);
          exit(10);
       }
       shape1=0.25*(1.0-xi)*(1.0-eta);
       shape2=0.25*(1.0+xi)*(1.0-eta);
       shape3=0.25*(1.0+xi)*(1.0+eta);
       shape4=0.25*(1.0-xi)*(1.0+eta);

       c1=iregnode;
       c2=iregnode+1;
       c3=iregnode+1+numtheta;
       c4=iregnode+numtheta;

       if ( (c3>(numtheta*numphi)) || (c1<1) )  {
         fprintf(stderr,"ERROR(get plate...) - %d %d %d %d\n",c1,c2,c3,c4);
         fflush(stderr);
         exit(10);
       }

       vphi = regular_f[c1]*shape1+regular_f[c2]*shape2
            + regular_f[c3]*shape3+regular_f[c4]*shape4;

       field[j][isurface_element]=vphi;

     }
   }

return;
}

/////////////////////////////////////////////////////
//
 void read_reg_grids(E,outfile,field)
   struct All_variables *E;
   char *outfile;
   double **field;   
  {

char input_s[1000];
FILE *fp2;

void  parallel_process_termination();
int numtheta,numphi,nregnodes;
int kk,j;
int ntheta,nphi,iregnode,iregnode2;
int node;
int c1,c2,c3,c4;
int isurface_element,lev;

double del_degree;
double del_theta;
double del_phi;
double lat_min,lat_max;
double theta_min,theta_max;
double phi_min,phi_max;
double vphi;
double rlat,rlong;
double theta,phi;
double phiprime,thetaprime;
double xi,eta;
double shape1,shape2,shape3,shape4;

static double *regular_f;
static int been_here=0;

/* INPUT FILE SPECIFIC STUFF */
/* format longitude, latitude, vphi, vtheta (vtheta is opposite in convention) */
/* file sequence given in specific order (see below):                          */
/* inputfile specific parameters:                                              */

   del_degree=180.0/(E->sphere.nox-1); /* increments given in input file */
   lat_min=-(90.0-del_degree/2.0);         /* minimum lattitude */
   lat_max=+(90.0-del_degree/2.0);         /* maximum lattitude */
   phi_min=del_degree/2.0;             /* minimum longitude */
   phi_max=360.0-del_degree/2.0;           /* maximum longitude */
   numtheta=E->sphere.nox-1;       /* number of theta increments */
   numphi=2*numtheta;               /* number of phi increments */
 
   del_theta=del_degree*M_PI/180.0;
   del_phi=del_degree*M_PI/180.0;
   theta_min=-1.0*(lat_max-90.0)*M_PI/180.0;
   theta_max=-1.0*(lat_min-90.0)*M_PI/180.0;
   phi_min=phi_min*M_PI/180.0;
   phi_max=phi_max*M_PI/180.0;

   nregnodes=numtheta*numphi;


/* adjust parameters for ghosting - extra regular elements are added on all sides */
/* This is required to interpolate nodes between thetamin and thetamax, and phimin and phimax */
/* Domain is padded by extra +2 numphi and +2 numtheta */

   numtheta=numtheta+2;
   numphi=numphi+2;
   theta_min=theta_min-del_theta;
   theta_max=theta_max+del_theta;
   phi_min=phi_min-del_phi;
   phi_max=phi_max+del_phi;

   nregnodes=numtheta*numphi;

 if (been_here==0)  {
   regular_f = (double *)malloc((nregnodes+1)*sizeof(double));
   been_here = 1;
   }


/* read from file */

   if ( (fp2=fopen(outfile,"r"))==NULL)   {
      fprintf(stderr,"ERROR(read_reg_grids)- file %s not found\n",outfile);
      fflush(stderr);
      parallel_process_termination();
   }

   for (kk=1;kk<=nregnodes;kk++)
      regular_f[kk]=-999999;

/* input file specific order of reading  right here */
/* (note - leaving room for ghosting) */

//   for (ntheta=(numtheta-1);ntheta>1;ntheta--)
   for (ntheta=2;ntheta<=(numtheta-1);ntheta++)
     for (nphi=2;nphi<=(numphi-1);nphi++)     {

      iregnode=ntheta+(nphi-1)*numtheta;

      if ((iregnode>nregnodes)||(iregnode<1))   {
        fprintf(stderr,"ERROR(plate velocities)-wrong iregnode %d %d\n",iregnode,nregnodes);
        fflush(stderr);
        exit(10);
      }

      fgets(input_s,1000,fp2);
      sscanf(input_s,"%lf %lf %le",&rlong,&rlat,&vphi);

     regular_f[iregnode]=vphi;

   }

   fclose(fp2);

/* ghost sides of regular domain */

   nphi=1;
   for (ntheta=2;ntheta<=(numtheta-1);ntheta++)  {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=ntheta+((numphi-1)-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }

   nphi=numphi;
   for (ntheta=2;ntheta<=(numtheta-1);ntheta++)   {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=ntheta+(2-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }

   ntheta=1;
   for (nphi=1;nphi<=(numphi/2);nphi++)   {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=(ntheta+1)+(nphi+numphi/2-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }
   for (nphi=(numphi/2)+1;nphi<=numphi;nphi++)   {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=(ntheta+1)+(nphi-(numphi/2)-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }

   ntheta=(numtheta);
   for (nphi=1;nphi<=(numphi/2);nphi++)   {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=(ntheta-1)+(nphi+numphi/2-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }
   for (nphi=(numphi/2)+1;nphi<=numphi;nphi++)  {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=(ntheta-1)+(nphi-(numphi/2)-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }

/* interpolate real mesh */

   isurface_element=0;
   for (j=1;j<=E->sphere.caps_per_proc;j++)
   for (node=1;node<=E->lmesh.nno;node++)   {
     theta=E->sx[j][1][node];
     phi  =E->sx[j][2][node];
     if (node%E->lmesh.noz==0)  {
        isurface_element++;
/* find position on regular mesh */

        ntheta=((theta-theta_min)/del_theta)+1;
        nphi=((phi-phi_min)/del_phi)+1;

/* iregnode is the lower left node in the regular element */

        iregnode=ntheta+(nphi-1)*numtheta;

/* find theta and phi distance from node */

        thetaprime=theta-(theta_min+(ntheta-1)*del_theta);
        phiprime=phi-(phi_min+(nphi-1)*del_phi);

        if (thetaprime<0.0 || thetaprime>del_theta)  {
           fprintf(stderr,"ERROR(get plate velocities)-wierd thetaprime %f %f\n",thetaprime,del_theta);
           fflush(stderr);
           exit(10);
        }
        if (phiprime<0.0 || phiprime>del_phi)   {
           fprintf(stderr,"ERROR(get plate velocities)-wierd phiprime %f %f\n",phiprime,del_phi);
           fflush(stderr);
           exit(10);
        }

/* determine shape functions */

       xi=2.0/del_theta*thetaprime-1.0;
       eta=2.0/del_phi*phiprime-1.0;

       if (xi<(-1-1e-6) || xi > (1+1e-6) || eta < (-1-1e-6) || eta > (1+1e-6))  {
          fprintf(stderr,"bad local variables %f %f (%f %f) (%f %f)\n",xi,eta,thetaprime,phiprime,del_theta,del_phi);
          fflush(stderr);
          exit(10);
       }
       shape1=0.25*(1.0-xi)*(1.0-eta);
       shape2=0.25*(1.0+xi)*(1.0-eta);
       shape3=0.25*(1.0+xi)*(1.0+eta);
       shape4=0.25*(1.0-xi)*(1.0+eta);

       c1=iregnode;
       c2=iregnode+1;
       c3=iregnode+1+numtheta;
       c4=iregnode+numtheta;

       if ( (c3>(numtheta*numphi)) || (c1<1) )  {
         fprintf(stderr,"ERROR(get plate...) - %d %d %d %d\n",c1,c2,c3,c4);
         fflush(stderr);
         exit(10);
       }

       vphi = regular_f[c1]*shape1+regular_f[c2]*shape2
            + regular_f[c3]*shape3+regular_f[c4]*shape4;

       field[j][isurface_element]=vphi;

     }
   }

return;
}


void interpolate_to_FE_grid(E,field,regular_f,numtheta,numphi)
  struct All_variables *E;
  double *field,*regular_f;
  int numtheta, numphi;
{

void  parallel_process_termination();
int nregnodes;
int kk,j;
int ntheta,nphi,iregnode,iregnode2;
int node;
int c1,c2,c3,c4;
int isurface_element,lev;

double del_degree;
double del_theta;
double del_phi;
double lat_min,lat_max;
double theta_min,theta_max;
double phi_min,phi_max;
double vphi;
double rlat,rlong;
double theta,phi;
double phiprime,thetaprime;
double xi,eta;
double shape1,shape2,shape3,shape4;

   del_degree=180.0/(numtheta-2); /* increments given in input file */
   lat_min=-(90.0-del_degree/2.0);         /* minimum lattitude */
   lat_max=+(90.0-del_degree/2.0);         /* maximum lattitude */
   phi_min=del_degree/2.0;             /* minimum longitude */
   phi_max=360.0-del_degree/2.0;           /* maximum longitude */
 
   del_theta=del_degree*M_PI/180.0;
   del_phi=del_degree*M_PI/180.0;
   theta_min=-1.0*(lat_max-90.0)*M_PI/180.0;
   theta_max=-1.0*(lat_min-90.0)*M_PI/180.0;
   phi_min=phi_min*M_PI/180.0;
   phi_max=phi_max*M_PI/180.0;

   nregnodes=numtheta*numphi;

   theta_min=theta_min-del_theta;
   theta_max=theta_max+del_theta;
   phi_min=phi_min-del_phi;
   phi_max=phi_max+del_phi;

/* ghost sides of regular domain */

   nphi=1;
   for (ntheta=2;ntheta<=(numtheta-1);ntheta++)  {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=ntheta+((numphi-1)-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }

   nphi=numphi;
   for (ntheta=2;ntheta<=(numtheta-1);ntheta++)   {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=ntheta+(2-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }

   ntheta=1;
   for (nphi=1;nphi<=(numphi/2);nphi++)   {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=(ntheta+1)+(nphi+numphi/2-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }
   for (nphi=(numphi/2)+1;nphi<=numphi;nphi++)   {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=(ntheta+1)+(nphi-(numphi/2)-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }

   ntheta=(numtheta);
   for (nphi=1;nphi<=(numphi/2);nphi++)   {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=(ntheta-1)+(nphi+numphi/2-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }
   for (nphi=(numphi/2)+1;nphi<=numphi;nphi++)  {
      iregnode=ntheta+(nphi-1)*numtheta;
      iregnode2=(ntheta-1)+(nphi-(numphi/2)-1)*numtheta;
      regular_f[iregnode]=regular_f[iregnode2];
   }

/* interpolate real mesh */

   isurface_element=0;
   for (j=1;j<=E->sphere.caps_per_proc;j++)
   for (node=1;node<=E->lmesh.nno;node++)   {
     theta=E->sx[j][1][node];
     phi  =E->sx[j][2][node];
     if (node%E->lmesh.noz==0)  {
        isurface_element++;
/* find position on regular mesh */

        ntheta=((theta-theta_min)/del_theta)+1;
        nphi=((phi-phi_min)/del_phi)+1;

/* iregnode is the lower left node in the regular element */

        iregnode=ntheta+(nphi-1)*numtheta;

/* find theta and phi distance from node */

        thetaprime=theta-(theta_min+(ntheta-1)*del_theta);
        phiprime=phi-(phi_min+(nphi-1)*del_phi);

        if (thetaprime<0.0 || thetaprime>del_theta)  {
           fprintf(stderr,"ERROR(get plate velocities)-wierd thetaprime %f %f\n",thetaprime,del_theta);
           fflush(stderr);
           exit(10);
        }
        if (phiprime<0.0 || phiprime>del_phi)   {
           fprintf(stderr,"ERROR(get plate velocities)-wierd phiprime %f %f\n",phiprime,del_phi);
           fflush(stderr);
           exit(10);
        }

/* determine shape functions */

       xi=2.0/del_theta*thetaprime-1.0;
       eta=2.0/del_phi*phiprime-1.0;

       if (xi<(-1-1e-6) || xi > (1+1e-6) || eta < (-1-1e-6) || eta > (1+1e-6))  {
          fprintf(stderr,"bad local variables %f %f (%f %f) (%f %f)\n",xi,eta,thetaprime,phiprime,del_theta,del_phi);
          fflush(stderr);
          exit(10);
       }
       shape1=0.25*(1.0-xi)*(1.0-eta);
       shape2=0.25*(1.0+xi)*(1.0-eta);
       shape3=0.25*(1.0+xi)*(1.0+eta);
       shape4=0.25*(1.0-xi)*(1.0+eta);

       c1=iregnode;
       c2=iregnode+1;
       c3=iregnode+1+numtheta;
       c4=iregnode+numtheta;

       if ( (c3>(numtheta*numphi)) || (c1<1) )  {
         fprintf(stderr,"ERROR(get plate...) - %d %d %d %d\n",c1,c2,c3,c4);
         fflush(stderr);
         exit(10);
       }

       vphi = regular_f[c1]*shape1+regular_f[c2]*shape2
            + regular_f[c3]*shape3+regular_f[c4]*shape4;

       field[isurface_element]=vphi;

     }
   }

return;
}
