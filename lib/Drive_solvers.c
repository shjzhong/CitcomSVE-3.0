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



/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *<LicenseText>
 *
 * CitcomS by Louis Moresi, Shijie Zhong, Lijie Han, Eh Tan,
 * Clint Conrad, Michael Gurnis, and Eun-seo Choi.
 * Copyright (C) 1994-2005, California Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *</LicenseText>
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
#include <math.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "drive_solvers.h"

double global_vdot();
double global_vdot_e();
double vnorm_nonnewt();
int need_visc_update(struct All_variables *);
void parallel_process_termination();



/************************************************************/

void general_equations_of_motions_solver_setup(struct All_variables *E)
{
  int i, m;
  void construct_node_maps();

  if (E->control.NMULTIGRID || E->control.NASSEMBLE)
    construct_node_maps(E);
  else
    for (i=E->mesh.gridmin;i<=E->mesh.gridmax;i++)
      for (m=1;m<=E->sphere.caps_per_proc;m++)
	E->elt_k[i][m]=(struct EK *)malloc((E->lmesh.NEL[i]+1)*sizeof(struct EK));

  if (E->parallel.me==0) fprintf(stderr,"done for equations_of_motions_solver_setup\n");

  return;
}




void general_equations_of_motions_solver(struct All_variables *E)
{
  void solve_constrained_flow_iterative();
  void construct_stiffness_B_matrix();
  void velocities_conform_bcs();
  void assemble_forces();
  void get_system_viscosity();
  void remove_rigid_rot();

  float vmag;

  double Udot_mag, dUdot_mag,omega[3];
  double Udot_mag1, dUdot_mag1,omega1[3];
  double PW_mag,dPW_mag;
  int m,count,count1,i,j,k;

  static double *oldU[NCS],*oldU1[NCS],pwold[3],time0;
  static int visits=0;

  double *delta_U[NCS],*delta_U1[NCS];
  double CPU_time0(),time;

  const int nno = E->lmesh.nno;
  const int nel = E->lmesh.nel;
  const int neq = E->lmesh.neq;
  const int vpts = vpoints[E->mesh.nsd];
  const int dims = E->mesh.nsd;
  const int addi_dof = additional_dof[dims];

  if (visits==0) {
    for (m=1;m<=E->sphere.caps_per_proc;m++)  {
      oldU[m] = (double *)malloc(neq*sizeof(double));
      oldU1[m] = (double *)malloc(neq*sizeof(double));
      for(i=0;i<neq;i++) {
	oldU[m][i]=0.0;
	oldU1[m][i]=0.0;
        }
    }
   time0 = CPU_time0();
   visits++;
   } 

  for (m=1;m<=E->sphere.caps_per_proc;m++) 
      delta_U[m] = (double *)malloc(neq*sizeof(double));

  if(E->parallel.me==0)  time=CPU_time0();

  E->viscosity.iterate=0;

  velocities_conform_bcs(E,E->U);

            // for Newtonian rheology (i.e., non-stress-dependent)
if (  !E->viscosity.SDEPV  )   {



    if(need_visc_update(E)){
      get_system_viscosity(E,1,E->EVI[E->mesh.levmax],E->VI[E->mesh.levmax]);

      construct_stiffness_B_matrix(E);
    }

    count=0;
    do {
	  E->monitor.count=count; // add by tao;
	  assemble_forces(E,count);
		
          solve_constrained_flow_iterative(E);

    for (m=1;m<=E->sphere.caps_per_proc;m++)
	for (i=0;i<neq;i++) {
	  delta_U[m][i] = E->U[m][i] - oldU[m][i];
	  oldU[m][i] = E->U[m][i];
	}

      dPW_mag = PW_mag = 0.0;
      if (E->ve_data_cont.polar_wander) {
        PW_mag = sqrt(pwold[0]*pwold[0] + pwold[1]*pwold[1]); 
        dPW_mag = sqrt( ((E->ve_data_cont.PW_incr[0]-pwold[0])*(E->ve_data_cont.PW_incr[0]-pwold[0]))
                      + ((E->ve_data_cont.PW_incr[1]-pwold[1])*(E->ve_data_cont.PW_incr[1]-pwold[1])) );
        dPW_mag = dPW_mag/PW_mag;
        pwold[0] = E->ve_data_cont.PW_incr[0];
        pwold[1] = E->ve_data_cont.PW_incr[1];
        }

      Udot_mag  = sqrt(global_vdot(E,oldU,oldU,E->mesh.levmax));
      dUdot_mag = sqrt(global_vdot(E,delta_U,delta_U,E->mesh.levmax));
      dUdot_mag = dUdot_mag/Udot_mag;

      if(E->parallel.me==0) {
	fprintf(stderr,"dU %.4e (%.4e) dPW %.4e (%.4e) for iteration %d\n",
		dUdot_mag,Udot_mag,dPW_mag,PW_mag, count);
        time=CPU_time0()-time0;
	fprintf(E->fp,"!!! iteration=%d relative change %g  %g m0 m1 %g %g potl_d %g timestep %d elapse time %g\n",
		count,dUdot_mag,dPW_mag,pwold[0],pwold[1],E->ve_data_cont.potential_vary_PW,E->monitor.solution_cycles,time);
        fflush(E->fp);
      }

      count++;

//     } while ( (count<50) && ((dUdot_mag > E->viscosity.sdepv_misfit) || (dPW_mag > (E->viscosity.sdepv_misfit*3))) && E->ve_data_cont.SELFG);     // for BM
//     } while ( (count<50) && ((dUdot_mag > E->viscosity.sdepv_misfit) || (dPW_mag > 0.03)) && E->ve_data_cont.SELFG);
//     } while ( (count<50) && ((dUdot_mag > E->viscosity.sdepv_misfit) || (E->ve_data_cont.potential_vary_PW > 3*E->viscosity.sdepv_misfit)) && E->ve_data_cont.SELFG);
     } while ( (count<50) && (dUdot_mag > E->viscosity.sdepv_misfit) && E->ve_data_cont.SELFG);

  } /*end if !SDEPV */

else if (E->viscosity.SDEPV ) {

   count1 = 0;
   E->viscosity.iterate=0;

//   do {
      
     count=0;
     do {
       get_system_viscosity(E,1,E->EVI[E->mesh.levmax],E->VI[E->mesh.levmax]);
       construct_stiffness_B_matrix(E);

       assemble_forces(E,count);
       solve_constrained_flow_iterative(E);

       for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (i=0;i<neq;i++) {
          delta_U[m][i] = E->U[m][i] - oldU[m][i];
          oldU[m][i] = E->U[m][i];
        }

       dPW_mag = PW_mag = 0.0;
       if (E->ve_data_cont.polar_wander) {
         PW_mag = sqrt(pwold[0]*pwold[0] + pwold[1]*pwold[1]); 
         dPW_mag = sqrt( ((E->ve_data_cont.PW_incr[0]-pwold[0])*(E->ve_data_cont.PW_incr[0]-pwold[0]))
                       + ((E->ve_data_cont.PW_incr[1]-pwold[1])*(E->ve_data_cont.PW_incr[1]-pwold[1])) );
         dPW_mag = dPW_mag/PW_mag;
         pwold[0] = E->ve_data_cont.PW_incr[0];
         pwold[1] = E->ve_data_cont.PW_incr[1];
         }

       Udot_mag  = sqrt(global_vdot(E,oldU,oldU,E->mesh.levmax));
       dUdot_mag = sqrt(global_vdot(E,delta_U,delta_U,E->mesh.levmax));
       dUdot_mag = dUdot_mag/Udot_mag;

       Udot_mag1  = sqrt(global_vdot_e(E,oldU,oldU,E->mesh.levmax));
       dUdot_mag1 = sqrt(global_vdot_e(E,delta_U,delta_U,E->mesh.levmax));
       dUdot_mag1 = dUdot_mag1/Udot_mag1;

       if (E->parallel.me == 0) {
         fprintf(stderr,
                 "dU %.4e (%.4e) %.4e dPW %.4e (%.4e) for iteration %d\n",
                 dUdot_mag, Udot_mag, dUdot_mag1, dPW_mag, PW_mag, count);
         time = CPU_time0() - time0;
         fprintf(E->fp,
                 "!!! iteration=%d relative change %g %g %g m0 m1 %g %g "
                 "timestep %d elapse time %g\n",
                 count, dUdot_mag, dUdot_mag1, dPW_mag, pwold[0], pwold[1],
                 E->monitor.solution_cycles, time);
         fflush(E->fp);
       }

      count++;
      E->viscosity.iterate++;

    // } while ( (count<50) && ((dUdot_mag1 > E->viscosity.sdepv_misfit) || (dPW_mag >0.03)) && E->ve_data_cont.SELFG);
      } while ( (count<50) && ((dUdot_mag1 > E->viscosity.sdepv_misfit) ) && E->ve_data_cont.SELFG);

  } // end for SDEPV 




  /* remove the rigid rotation component from the velocity solution */
  if((E->sphere.caps == 12) &&
     (E->control.remove_rigid_rotation || E->control.remove_angular_momentum)) {
      remove_rigid_rot(E);
  }


    for (m=1;m<=E->sphere.caps_per_proc;m++)  {
      free((void *) delta_U[m]);
    }


  return;
}

int need_visc_update(struct All_variables *E)
{
  if(E->ve_data_cont.DIRECT){
    /* always update */
    return 1;
  }else{
    /* deal with first time called */
    if(E->control.restart){	
      /* restart step - when this function is called, the cycle has
	 already been incremented */
      if(E->monitor.solution_cycles ==  E->monitor.solution_cycles_init + 1)
	return 1;
      else
	return 0;
    }else{
      /* regular step */
      if(E->monitor.solution_cycles == 0)
	return 1;
      else
	return 0;
    }
  }
}

void general_equations_of_motions_solver_pseudo_surf(struct All_variables *E)
{
  void solve_constrained_flow_iterative_pseudo_surf();
  void construct_stiffness_B_matrix();
  void velocities_conform_bcs();
  void assemble_forces_pseudo_surf();
  void get_system_viscosity();
  void std_timestep();
  void remove_rigid_rot();
  void get_STD_freesurf(struct All_variables *, float**);

  float vmag;

  double Udot_mag, dUdot_mag;
  int m,count,i,j,k,topo_loop;

  double *oldU[NCS], *delta_U[NCS];

  const int nno = E->lmesh.nno;
  const int nel = E->lmesh.nel;
  const int neq = E->lmesh.neq;
  const int vpts = vpoints[E->mesh.nsd];
  const int dims = E->mesh.nsd;
  const int addi_dof = additional_dof[dims];

  velocities_conform_bcs(E,E->U);

  E->monitor.stop_topo_loop = 0;
  E->monitor.topo_loop = 0;
  if(E->monitor.solution_cycles==0) std_timestep(E);
  while(E->monitor.stop_topo_loop == 0) {
	  assemble_forces_pseudo_surf(E,0);
	  if(need_visc_update(E)){
	    get_system_viscosity(E,1,E->EVI[E->mesh.levmax],E->VI[E->mesh.levmax]);
	    construct_stiffness_B_matrix(E);
	  }
	  solve_constrained_flow_iterative_pseudo_surf(E);

	  if (E->viscosity.SDEPV || E->viscosity.PDEPV) {

		  for (m=1;m<=E->sphere.caps_per_proc;m++)  {
			  delta_U[m] = (double *)malloc(neq*sizeof(double));
			  oldU[m] = (double *)malloc(neq*sizeof(double));
			  for(i=0;i<neq;i++)
				  oldU[m][i]=0.0;
		  }

		  Udot_mag=dUdot_mag=0.0;
		  count=1;

		  while (1) {

			  for (m=1;m<=E->sphere.caps_per_proc;m++)
				  for (i=0;i<neq;i++) {
					  delta_U[m][i] = E->U[m][i] - oldU[m][i];
					  oldU[m][i] = E->U[m][i];
				  }

			  Udot_mag  = sqrt(global_vdot(E,oldU,oldU,E->mesh.levmax));
			  dUdot_mag = vnorm_nonnewt(E,delta_U,oldU,E->mesh.levmax);

			  if(E->parallel.me==0){
				  fprintf(stderr,"Stress dependent viscosity: DUdot = %.4e (%.4e) for iteration %d\n",dUdot_mag,Udot_mag,count);
				  fprintf(E->fp,"Stress dependent viscosity: DUdot = %.4e (%.4e) for iteration %d\n",dUdot_mag,Udot_mag,count);
				  fflush(E->fp);
			  }

			  if (count>50 || dUdot_mag<E->viscosity.sdepv_misfit)
				  break;

			  get_system_viscosity(E,1,E->EVI[E->mesh.levmax],E->VI[E->mesh.levmax]);
			  construct_stiffness_B_matrix(E);
			  solve_constrained_flow_iterative_pseudo_surf(E);

			  count++;

		  } /*end while */
		  for (m=1;m<=E->sphere.caps_per_proc;m++)  {
			  free((void *) oldU[m]);
			  free((void *) delta_U[m]);
		  }

	  } /*end if SDEPV or PDEPV */
	  E->monitor.topo_loop++;
  }

  /* remove the rigid rotation component from the velocity solution */
  if((E->sphere.caps == 12) &&
     (E->control.remove_rigid_rotation || E->control.remove_angular_momentum)) {
      remove_rigid_rot(E);
  }

  get_STD_freesurf(E,E->slice.freesurf);

  return;
}
