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
/* Functions relating to the determination of viscosity field either
   as a function of the run, as an initial condition or as specified from
   a previous file */


#include <math.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "parsing.h"
#include <stdio.h>


void myerror(struct All_variables *,char *);

static void apply_low_visc_wedge_channel(struct All_variables *E, float **evisc);
static void low_viscosity_channel_factor(struct All_variables *E, float *F);
static void low_viscosity_wedge_factor(struct All_variables *E, float *F);
void parallel_process_termination();


void elastic_viscosity_input(struct All_variables *E)
{
    int m=E->parallel.me;
    int i;

    /* default values .... */
    for(i=0;i < CITCOM_MAX_VISC_LAYER;i++) {
        E->viscosity.N0[i]=1.0;
        E->viscosity.T[i] = 0.0;
        E->viscosity.Z[i] = 0.0;
        E->viscosity.E[i] = 0.0;

	E->viscosity.pdepv_a[i] = 1.e20; /* \sigma_y = min(a + b * (1-r),y) */
	E->viscosity.pdepv_b[i] = 0.0;
	E->viscosity.pdepv_y[i] = 1.e20;


    }
    for(i=0;i<10;i++)
      E->viscosity.cdepv_ff[i] = 1.0; /* flavor factors for CDEPV */


    /* read in information */
    input_boolean("VISC_UPDATE",&(E->viscosity.update_allowed),"on",m);
    input_int("rheol",&(E->viscosity.RHEOL),"3",m);

    E->viscosity.FROM_FILE = 0;
    input_int("3d_visc_from_file",&(E->viscosity.FROM_FILE),"0",m);
    if (E->viscosity.FROM_FILE) {
      input_string("3d_visc_datafile",E->ve_data_cont.visc_file,"initialize",m);
      }

    E->ve_data_cont.ELASTIC_FROM_FILE = 0;
    input_int("3d_elastic_from_file", &(E->ve_data_cont.ELASTIC_FROM_FILE), "0",m);
    if (E->ve_data_cont.ELASTIC_FROM_FILE){
 	 input_string("3d_elastic_datafile", E->ve_data_cont.elastic_file, "initialize",m);
	}

    input_int("compressible",&(E->ve_data_cont.compressible),"0",m);

    E->ve_data_cont.OneDmodel_read=1;
    E->ve_data_cont.OneDmodel_read_g = 0;
//    input_int("1Dmodel_from_file",&(E->ve_data_cont.OneDmodel_read),"0",m);
    if(E->ve_data_cont.OneDmodel_read) {
      input_string("1Dmodel_datafile",E->ve_data_cont.OneDmodel_file,"initialize",m);
      input_int("1Dmodel_g_from_file",&(E->ve_data_cont.OneDmodel_read_g),"0",m);
      } 


    input_float_vector("visc0",E->viscosity.num_mat,(E->viscosity.N0),m);
    input_float_vector("shearModulus",E->viscosity.num_mat,(E->viscosity.G),m);

    input_boolean("TDEPV",&(E->viscosity.TDEPV),"on",m);
    if (E->viscosity.TDEPV) {
        input_float_vector("viscT",E->viscosity.num_mat,(E->viscosity.T),m);
        input_float_vector("viscE",E->viscosity.num_mat,(E->viscosity.E),m);
        input_float_vector("viscZ",E->viscosity.num_mat,(E->viscosity.Z),m);
	/* for viscosity 8 */
        input_float("T_sol0",&(E->viscosity.T_sol0),"0.6",m);
        input_float("ET_red",&(E->viscosity.ET_red),"0.1",m);
    }


    E->viscosity.beta=1.0;

    input_float("polar_wander_relax",&(E->viscosity.beta),"1.000",m);

    E->viscosity.sdepv_misfit = 1.0;
    input_boolean("SDEPV",&(E->viscosity.SDEPV),"off",m);
    if (E->viscosity.SDEPV) {
      E->viscosity.sdepv_visited = 0;
      input_float_vector("sdepv_expt",E->viscosity.num_mat,(E->viscosity.sdepv_expt),m);
      input_float_vector("sdepv_trns",E->viscosity.num_mat,(E->viscosity.sdepv_trns),m);
      input_float_vector("sdepv_bg",E->viscosity.num_mat,(E->viscosity.sdepv_bg),m);
      input_float("sdepv_relax_alpha",&(E->viscosity.alpha),"1.000",m);
      for(i=0;i<E->viscosity.num_mat;i++) {
        E->viscosity.sdepv_trns[i]=E->viscosity.sdepv_trns[i]/E->ve_data_cont.shear_mod;
        E->viscosity.sdepv_bg[i]  =E->viscosity.sdepv_bg[i]  /E->ve_data_cont.shear_mod;
        }
    }


    input_boolean("PDEPV",&(E->viscosity.PDEPV),"off",m); /* plasticity addition by TWB */
    if (E->viscosity.PDEPV) {
      E->viscosity.pdepv_visited = 0;

      input_boolean("pdepv_eff",&(E->viscosity.pdepv_eff),"on",m);
      input_float_vector("pdepv_a",E->viscosity.num_mat,(E->viscosity.pdepv_a),m);
      input_float_vector("pdepv_b",E->viscosity.num_mat,(E->viscosity.pdepv_b),m);
      input_float_vector("pdepv_y",E->viscosity.num_mat,(E->viscosity.pdepv_y),m);

      input_float("pdepv_offset",&(E->viscosity.pdepv_offset),"0.0",m);
    }
      input_float("iteration_tolerance",&(E->viscosity.sdepv_misfit),"0.001",m);


    input_boolean("CDEPV",&(E->viscosity.CDEPV),"off",m);
    if(E->viscosity.CDEPV){
      /* compositional viscosity */
      if(E->control.tracer < 1){
	fprintf(stderr,"error: CDEPV requires tracers, but tracer is off\n");
	parallel_process_termination();
      }
      if(E->trace.nflavors > 10)
	myerror(E,"error: too many flavors for CDEPV");
      /* read in flavor factors */
      input_float_vector("cdepv_ff",E->trace.nflavors,
			 (E->viscosity.cdepv_ff),m);
      /* and take the log because we're using a geometric avg */
      for(i=0;i<E->trace.nflavors;i++)
	E->viscosity.cdepv_ff[i] = log(E->viscosity.cdepv_ff[i]);
    }


    input_boolean("low_visc_channel",&(E->viscosity.channel),"off",m);
    input_boolean("low_visc_wedge",&(E->viscosity.wedge),"off",m);

    input_float("lv_min_radius",&(E->viscosity.lv_min_radius),"0.9764",m);
    input_float("lv_max_radius",&(E->viscosity.lv_max_radius),"0.9921",m);
    input_float("lv_channel_thickness",&(E->viscosity.lv_channel_thickness),"0.0047",m);
    input_float("lv_reduction",&(E->viscosity.lv_reduction),"0.5",m);

    input_boolean("VMAX",&(E->viscosity.MAX),"off",m);
    if (E->viscosity.MAX)
        input_float("visc_max",&(E->viscosity.max_value),"1e22,1,nomax",m);

    input_boolean("VMIN",&(E->viscosity.MIN),"off",m);
    if (E->viscosity.MIN)
        input_float("visc_min",&(E->viscosity.min_value),"1e20",m);

    input_string("Viscosity",E->viscosity.STRUCTURE,"system",m);
    input_int ("visc_smooth_method",&(E->viscosity.smooth_cycles),"0",m);

    E->viscosity.FROM_SYSTEM = 1;

    return;
}


void set_mech_model_from_files(struct All_variables *E)
{
    void visc_from_file();
    void mech_model_from_file();
    void elastic_model_from_file_3D();

    mech_model_from_file(E);

    if (E->viscosity.FROM_FILE) // read 3D viscosity
        visc_from_file(E);

    if (E->ve_data_cont.ELASTIC_FROM_FILE)
        elastic_model_from_file_3D(E, E->elambda_o[E->mesh.levmax],E->esmu_o[E->mesh.levmax],E->ve_data_cont.elastic_file,E->ve_data_cont.shear_mod);


    return;
}


/* ============================================ */
/* get_system_viscosity
    always called as get_system_viscosity(E,1,E->EVI[E->mesh.levmax],E->VI[E->mesh.levmax]);
    However, does not rely on the initial values of E->EVI and E->VI. 
    It constructs them.

Construct those variables:
E->EVI: was viscosity on the gauss points before the last step (add_viscoelasticity), 
          afterward it is mu_hat: mu/(1+alpha) (Geruo A 2012)
E->VI: viscosity on nodes. (constructed from E->EVI before add_viscoelasticity).

E->Maxwelltime[m][i]: was non-dimensionlized mu before the last step (add_viscoelasticity), 
          afterward it is (1.0-alpha)/(1.0+alpha);

*/
void get_system_viscosity(E,propogate,evisc,visc)
     struct All_variables *E;
     int propogate;
     float **evisc,**visc;
{
    void visc_from_mat();
    void visc_from_T();
    void visc_from_S();

    void visc_from_P();
    void visc_from_C();

    void apply_viscosity_smoother();
    void visc_from_gint_to_nodes();
    void add_viscoelasticity();
    void add_viscoelasticity_comp();
    void parallel_process_termination();

    static int been_here=0;
    int i,j,m,e;
    float temp1,temp2,*vvvis;
    double *TG;

    const int vpts = vpoints[E->mesh.nsd];

    if (E->viscosity.SDEPV)  {
      if (been_here==0) {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
          for(e=1;e<=E->lmesh.nel;e++) {
            E->EVolder[m][e] = 1.0;
            E->EVold[m][e] = 1.0;
           }
        been_here ++;
        }
      else  {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
          for(e=1;e<=E->lmesh.nel;e++) {
            E->EVolder[m][e] = E->EVold[m][e];
            E->EVold[m][e] = 0.125*
              (visc[m][E->ien[m][e].node[1]] + visc[m][E->ien[m][e].node[2]]
              +visc[m][E->ien[m][e].node[3]] + visc[m][E->ien[m][e].node[4]]
              +visc[m][E->ien[m][e].node[5]] + visc[m][E->ien[m][e].node[6]]
              +visc[m][E->ien[m][e].node[7]] + visc[m][E->ien[m][e].node[8]]);
           }
         }
    }

//       visc_from_mat(E,evisc);
//

//    if(E->viscosity.TDEPV)
//        visc_from_T(E,evisc,propogate);

//    if(E->viscosity.CDEPV)	/* compositional prefactor */
//      visc_from_C(E,evisc);
//
    for(m=1;m<=E->sphere.caps_per_proc;m++)
       for(i=1;i<=E->lmesh.nel;i++)
            for(j=1;j<=vpts;j++)
                  evisc[m][(i-1)*vpts + j] = E->evi_o[E->mesh.levmax][m][i];

    if(E->viscosity.SDEPV)
      visc_from_S(E,evisc,propogate);

    if(E->viscosity.PDEPV)	/* "plasticity" */
      visc_from_P(E,evisc);

    /* min/max cut-off */

    if(E->viscosity.MAX) {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(i=1;i<=E->lmesh.nel;i++)
                for(j=1;j<=vpts;j++)
                    if(evisc[m][(i-1)*vpts + j] > E->viscosity.max_value)
                        evisc[m][(i-1)*vpts + j] = E->viscosity.max_value;
    }

    if(E->viscosity.MIN) {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(i=1;i<=E->lmesh.nel;i++)
                for(j=1;j<=vpts;j++)
                    if(evisc[m][(i-1)*vpts + j] < E->viscosity.min_value)
                        evisc[m][(i-1)*vpts + j] = E->viscosity.min_value;
    }

    if (E->control.verbose)  {
      fprintf(E->fp_out,"output_evisc0 %d\n",E->monitor.solution_cycles);
      for(m=1;m<=E->sphere.caps_per_proc;m++) {
        fprintf(E->fp_out,"output_evisc for cap %d\n",E->sphere.capid[m]);
      for(i=1;i<=E->lmesh.elz;i++)
          fprintf(E->fp_out,"%d %d %f\n",i,E->mat[m][i],evisc[m][(i-1)*vpts+1]);
      }
      fflush(E->fp_out);
    }
    /* interpolate from gauss quadrature points to node points for output */
    visc_from_gint_to_nodes(E,evisc,visc,E->mesh.levmax);

    if (!E->ve_data_cont.compressible)
        add_viscoelasticity(E,evisc);
    else if (E->ve_data_cont.compressible)
        add_viscoelasticity_comp(E,evisc);

    return;
}



void initial_viscosity(struct All_variables *E)
{
    void report(struct All_variables*, char*);

    report(E,"Initialize viscosity field");

    return;
}

void mech_model_from_file(struct All_variables *E)
{
    int i,d,ii,m,jj,j,k,nn;
    char output_file[255],input_s[200];
    float tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
    float rho,vp,vs,g,visc;
    void parallel_process_termination();
    void ele_to_nodes();
    FILE *fp1;

    static int been_here=0;
    int ee,el,ndepth,i1,i2,idepth,mgroup;
    static double *radius, *model[6];
    static int *matgroup;
    double rr;
    const int vpts = vpoints[E->mesh.nsd];
    const int lev = E->mesh.levmax;

    double temp1, temp2;

    // variables for compute g
    void compute_g();
    //	double r1, g1, con, r2, g2, rho2;
    //	con = 3.1415926536*6.673e-11;

    sprintf(output_file,"%s",E->ve_data_cont.OneDmodel_file);
    fp1 = fopen(output_file,"r");
    fgets(input_s,200,fp1);   // read description
    fgets(input_s,200,fp1);
    sscanf(input_s,"%d",&ndepth);


    if (been_here==0)  {
      matgroup = (int *)malloc((ndepth+1)*sizeof(int));
      radius = (double *)malloc((ndepth+1)*sizeof(double));
      for (i=1;i<=5;i++)
        model[i] = (double *)malloc((ndepth+1)*sizeof(double));
      been_here = 1;
      }


    for (i=1;i<=ndepth;i++)  {
      fgets(input_s,200,fp1);
      sscanf(input_s,"%g %g %g %g %g %g %d",&tmp1,&rho,&vp,&vs,&g,&visc,&mgroup);
      radius[i] = tmp1/E->sphere.dradius;
      model[1][i] = rho;
      model[2][i] = rho*vp*vp - 2.0*rho*vs*vs;
      model[3][i] = rho*vs*vs; 
      model[4][i] = g;
      model[5][i] = visc;
      matgroup[i] = mgroup;
      if (E->parallel.me==0) fprintf(stderr,"%g %g %g %g %g %g %d\n",radius[i],rho,model[3][i],model[2][i],g,visc,mgroup);
      }

    fclose(fp1);

  for(m=1;m<=E->sphere.caps_per_proc;m++)
    for (jj=1;jj<=E->lmesh.elz;jj++)  {
      rr = (E->sx[m][3][jj]+E->sx[m][3][jj+1])/2.0;

      if (rr>=radius[ndepth]) {
        i1 = ndepth;
        i2 = ndepth;
        }
      else if (rr<=radius[1])  {
        i1=1;
        i2=1;
        }
      else {
        for (i=1;i<ndepth;i++)  {
          if (rr>=radius[i] && rr<=radius[i+1]) {
            i1=i;
            i2=i+1;
            break;
            } 
          }
        }

      if (i1!=i2) {     // read one more layer and interpolate radially
          temp1 = (radius[i2] - rr)/(radius[i2] - radius[i1]); 
          temp2 = (rr - radius[i1])/(radius[i2] - radius[i1]); 
          tmp1 = model[1][i1]*temp1 + model[1][i2]*temp2; 
          tmp2 = model[2][i1]*temp1 + model[2][i2]*temp2; 
          tmp3 = model[3][i1]*temp1 + model[3][i2]*temp2; 
          tmp4 = model[4][i1]*temp1 + model[4][i2]*temp2; 
          tmp5 = model[5][i1]*temp1 + model[5][i2]*temp2; 
          tmp6 = matgroup[i1]*temp1 + matgroup[i2]*temp2; 
	  mgroup = tmp6;
          }
            
       for (i = 1; i <= E->lmesh.ely; i++)
          for (j = 1; j <= E->lmesh.elx; j++) {
              el = jj + (j - 1) * E->lmesh.elz +
                   (i - 1) * E->lmesh.elz * E->lmesh.elx;
              E->erho[lev][m][el] = tmp1 / E->data.density;
              E->elambda_o[lev][m][el] = tmp2 / E->ve_data_cont.shear_mod;
              E->esmu_o[lev][m][el] = tmp3 / E->ve_data_cont.shear_mod;
              E->egrav[lev][m][el] = tmp4 / E->data.grav_acc;
              E->evi_o[lev][m][el] = tmp5 / E->data.ref_viscosity;
              E->mat[m][el] = mgroup;
            }
    } // end for jj of lmesh.elz
  // end for m for-loop
 
  if (E->ve_data_cont.OneDmodel_read_g == 1)
    ele_to_nodes(E, E->egrav[lev], E->grav[lev], lev); // add by tao
  else
    compute_g(E, lev);

	// note: for lower levels, those needed variables are filled by project_viscosity.

  if (been_here == 1 && E->parallel.me < E->parallel.nprocz) {
    fprintf(E->fp, "debug earth model, elemental, from proc %d\n",
            E->parallel.me);
    int element;
    int m = 1;
    // print header
    fprintf(E->fp, "z(element-center) dz(element-size) rho lambda mu g visc \n");
    for (element = 1; element <= E->lmesh.elz; element++) {

      fprintf(E->fp, "%e %e %e %e %e %e %e %d",
              E->ECO[E->mesh.levmax][m][element].centre[3] * E->sphere.dradius,
              E->ECO[E->mesh.levmax][m][element].size[3] * E->sphere.dradius,
              E->erho[E->mesh.levmax][m][element] * E->data.density,
              E->elambda_o[E->mesh.levmax][m][element] * E->ve_data_cont.shear_mod,
              E->esmu_o[E->mesh.levmax][m][element] * E->ve_data_cont.shear_mod,
              E->egrav[E->mesh.levmax][m][element] * E->data.grav_acc,
              E->evi_o[E->mesh.levmax][m][element] * E->data.ref_viscosity,
              E->mat[m][element]);
      fprintf(E->fp, "\n");
    }
    been_here = 2;
  }

  return;
 }

void compute_g(struct All_variables *E, int lev)
{
	//void get_rtf();
	void parallel_process_termination();
	int m,i,e,k,j,node;
	double r1,r2,con,g1,g2,rho2,tmp[9];
	const int vpts = vpoints[E->mesh.nsd];
	

// calculate Earth_mass_below: total mass below this processor
	double Earth_mass_below = 0.0;
	double Earth_mass_hori = 0.0;
	double tmp_mass = 0.0;
	int el;
	double* gnda;
	const int dims=E->mesh.nsd;

	con = 3.1415926536*6.673e-11;

	for (m=1;m<=E->sphere.caps_per_proc;m++)
		for (el=1;el<=E->lmesh.NEL[lev];el++) {
			for (j=1;j<=vpts;j++) {
				gnda = E->GDA[lev][m][el].vpt;
				tmp_mass += E->erho[lev][m][el]*gnda[j]*g_point[j].weight[dims-1];
			 } // end j loop
		} // end el loop		  
	MPI_Allreduce(&tmp_mass, &Earth_mass_hori,1,MPI_DOUBLE,MPI_SUM, E->parallel.horizontal_comm);
	
	MPI_Scan(&Earth_mass_hori, &Earth_mass_below, 1, MPI_DOUBLE,MPI_SUM, E->parallel.vertical_comm);
	Earth_mass_below -= Earth_mass_hori;
//	if(E->parallel.me_loc[3] ==0)
//		Earth_mass_below = 0;
	Earth_mass_below += 4./3.*M_PI*E->data.density_below/E->data.density*pow(E->sphere.ri,3);
	Earth_mass_below *= E->data.density * pow(E->sphere.dradius,3);

//MPI_Barrier(MPI_COMM_WORLD);	
	// element 
	for(m=1;m<=E->sphere.caps_per_proc;m++)
	  for(i=1;i<=E->lmesh.ELZ[lev];i++) {
		
		if (i == 1) {
			r1 = E->SX[lev][m][3][i] * E->sphere.dradius; //E->sx[m][3][jj]
			g1 = 6.673e-11* Earth_mass_below /(r1 * r1);
			r2 = (E->SX[lev][m][3][i]+E->SX[lev][m][3][i+1])/2.0 *E->sphere.dradius;
			rho2 = E->erho[lev][m][i]*E->data.density;
			g2 = g1*r1*r1/r2/r2 + 4./3.*con*rho2*(r2-pow(r1,3.)/r2/r2);
			E->egrav[lev][m][i] = g2/E->data.grav_acc;
			//r2 = E->SX[lev][m][3][i+1]*E->sphere.dradius;
			//g1 = g1*r1*r1/r2/r2 + 4./3.*con*rho2*(r2-pow(r1,3.)/r2/r2);			
		}
		else {
			r1 = (E->SX[lev][m][3][i]+E->SX[lev][m][3][i-1])/2.0*E->sphere.dradius;
			r2 = (E->SX[lev][m][3][i])*E->sphere.dradius;
			g1 = E->egrav[lev][m][i-1]*E->data.grav_acc;
			rho2 = E->erho[lev][m][i-1]*E->data.density;
			g2 = g1*r1*r1/r2/r2 + 4./3.*con*rho2*(r2-pow(r1,3.)/r2/r2); // 

			r1 = (E->SX[lev][m][3][i])*E->sphere.dradius;
			g1 = g2;
			r2 = (E->SX[lev][m][3][i]+E->SX[lev][m][3][i+1])/2.0*E->sphere.dradius;
			rho2 = E->erho[lev][m][i]*E->data.density;
			g2 = g1*r1*r1/r2/r2 + 4./3.*con*rho2*(r2-pow(r1,3.)/r2/r2); 
			E->egrav[lev][m][i] = g2/E->data.grav_acc;

		}
		for (j=1;j<=E->lmesh.SNEL[lev];j++){
		  e = (j-1)*E->lmesh.ELZ[lev]+i;
		  E->egrav[lev][m][e] = E->egrav[lev][m][i];
		}
	 }
	// nodes
	for(m=1;m<=E->sphere.caps_per_proc;m++)
	  for(i=1;i<=E->lmesh.NOZ[lev];i++) {
		if (i == 1) {
			r1 = E->SX[lev][m][3][i]*E->sphere.dradius;
			g1 = 6.673e-11 * Earth_mass_below /(r1 * r1);
			E->grav[lev][m][i] = g1/E->data.grav_acc;
		}
		else {
		  g1 = E->grav[lev][m][i-1]*E->data.grav_acc;
		  r1 = E->SX[lev][m][3][i-1]*E->sphere.dradius;
		  r2 = E->SX[lev][m][3][i]*E->sphere.dradius;
		  rho2 = E->erho[lev][m][i-1]*E->data.density;
		  g2 = g1*r1*r1/r2/r2 + 4./3.*con*rho2*(r2-pow(r1,3.)/r2/r2);
		  E->grav[lev][m][i] = g2/E->data.grav_acc;
		}
		for (j=1;j<=E->lmesh.NOX[lev]*E->lmesh.NOY[lev];j++){
		  node = (j-1)*E->lmesh.NOZ[lev]+i;
		  E->grav[lev][m][node] = E->grav[lev][m][i];
		}
	 }
//parallel_process_termination();	
	return;


}


void visc_from_file(E)
     struct All_variables *E;
{

    int i,d,ii,m,jj,j,k,nn;
    const int vpts = vpoints[E->mesh.nsd];
    char output_file[255],input_s[200];
    float temp1,temp2,temp3;
    void ele_to_nodes();
    void return_horiz_ave_f();
    void parallel_process_termination();
    void interpolate_to_FE_grid();
    FILE *fp1;

    static int been_here=0;
    int ee,el,ndepth,nfi,nth,nnfi,nnth,i1,i2,iskip,idepth;
    static double *radius_visc,*visc_reg,*visc_FE[2];
    static float *ave_visc;
    double rr;

    sprintf(output_file,"%s",E->ve_data_cont.visc_file);
    fp1 = fopen(output_file,"r");
    fgets(input_s,200,fp1);   // read description
    fgets(input_s,200,fp1);
    sscanf(input_s,"%d %d %d",&ndepth,&nfi,&nth);
    fgets(input_s,200,fp1);   // read description
    nn = nfi*nth;
    nnth = nth +2;    // add 1 row and 1 column to each side 
    nnfi = nfi +2;

    if (been_here==0)  {
      visc_reg = (double *)malloc(((nfi+2)*(nth+2)+1)*sizeof(double));
      radius_visc = (double *)malloc((ndepth+1)*sizeof(double));
      ave_visc = (float *)malloc((E->lmesh.noz+1)*sizeof(float));
      for (m=1;m<=E->sphere.caps_per_proc;m++)
         visc_FE[m] = (double *)malloc((E->lmesh.nsf+1)*sizeof(double));
      been_here = 1;
      }

    for (i=1;i<=ndepth;i++)  {
      fgets(input_s,200,fp1);
      sscanf(input_s,"%g",&temp1);
      radius_visc[i] = E->sphere.ro-temp1*1000.0/E->sphere.dradius;
//if(E->parallel.me==0) fprintf(E->fp_out,"visc model radius info %d %g\n",i,radius_visc[i]);
      }

    fclose(fp1);


  for(m=1;m<=E->sphere.caps_per_proc;m++)
    for (jj=1;jj<=E->lmesh.elz;jj++)  {
      rr = (E->sx[m][3][jj]+E->sx[m][3][jj+1])/2.0;

      if (rr>=radius_visc[ndepth]) {
       //  continue;
        i1 = ndepth;
        i2 = ndepth;
        }
      else if (rr<=radius_visc[1])  {
        //continue;
         i1=1;
         i2=1;
        }
      else {
        for (i=1;i<ndepth;i++)  {
          if (rr>=radius_visc[i] && rr<=radius_visc[i+1]) {
            i1=i;
            i2=i+1;
            break;
            } 
          }
        }

      sprintf(output_file,"%s",E->ve_data_cont.visc_file);
      fp1 = fopen(output_file,"r");
      iskip = 3+ndepth+1+nfi+1+nth + (i1-1)*(nn+1); 
      for (i=1;i<=iskip;i++)  
        fgets(input_s,200,fp1);

      fgets(input_s,200,fp1);
      for (i=2;i<=nnth-1;i++)  
      for (j=2;j<=nnfi-1;j++)  {
        k = i + (j-1)*nnth;
        fgets(input_s,200,fp1);
        sscanf(input_s,"%g",&temp1);
        // visc_reg[k] = log10f(temp1);       // another possibility is visc_reg[k] = log10f(temp1);
//        visc_reg[k] = 1/temp1;       // use Harmonic Mean
        visc_reg[k] = temp1;      
        }

      if (i1!=i2) {     // read one more layer and interpolate radially
        temp1 = (radius_visc[i2] - rr)/(radius_visc[i2] - radius_visc[i1]); 
        temp2 = (rr - radius_visc[i1])/(radius_visc[i2] - radius_visc[i1]); 
//if(E->parallel.me==0) fprintf(stderr,"%d %d %d rr= %g %g %g\n",jj,i1,i2,rr,temp1,temp2);
        fgets(input_s,200,fp1);
        for (i=2;i<=nnth-1;i++)  
        for (j=2;j<=nnfi-1;j++)  {
          k = i + (j-1)*nnth;
          fgets(input_s,200,fp1);
          sscanf(input_s,"%g",&temp3);
          // visc_reg[k] = visc_reg[k]*temp1 + log10f(temp3)*temp2;  // also use log10f here if used above
//          visc_reg[k] = visc_reg[k]*temp1 + 1/temp3*temp2;  // use Harmonic Mean
          visc_reg[k] = visc_reg[k]*temp1 + temp3*temp2;  // use Harmonic Mean
          }
        }
      fclose(fp1);

      interpolate_to_FE_grid(E,visc_FE[m],visc_reg,nnth,nnfi);

      for (i=1;i<=E->lmesh.ely;i++)  
      for (j=1;j<=E->lmesh.elx;j++)  {
         ee = j + (i-1)*E->lmesh.elx;      
         temp1 =( visc_FE[m][E->sien[m][ee].node[1]] 
                + visc_FE[m][E->sien[m][ee].node[2]] 
                + visc_FE[m][E->sien[m][ee].node[3]] 
                + visc_FE[m][E->sien[m][ee].node[4]])/4.0; 
         el = jj + (j-1)*E->lmesh.elz + (i-1)*E->lmesh.elz*E->lmesh.elx;
	 E->evi_o[E->mesh.levmax][m][el] = temp1/E->data.ref_viscosity;
//	 E->evi_o[E->mesh.levmax][m][el] = 1.0/temp1/E->data.ref_viscosity;
            // E->evi_o[E->mesh.levmax][m][el] = pow(10,temp1)/E->data.ref_viscosity;  // use power if log10 used above
         }

      }   // end for jj of lmesh.elz



  ele_to_nodes(E,E->evi_o[E->mesh.levmax],E->VI[E->mesh.levmax],E->mesh.levmax);

  return_horiz_ave_f(E,E->VI[E->mesh.levmax],ave_visc);


   if (E->parallel.me==0)  {
     fprintf(E->fp_out,"average visc vs radius after the interpolation to FE grid\n"); 
     for(k=1;k<=E->lmesh.noz;k++)   {
       fprintf(E->fp_out,"%d %.5e %.5e\n",k,E->sx[1][3][k],ave_visc[k]*E->data.ref_viscosity); 
       }
     }

//  parallel_process_termination();

  return;
  }



/*
read model (elambda, esmu, etc) from file (multiple layers of data on regular theta, phi grid).
Computed at element center, then assign to element (element_var).

This can be called after a inital model has been set up. It will overwrite the region that 
the input 3D model file covers.

So the 3D model file doesn't need to cover the whole depth range.
*/
void elastic_model_from_file_3D(E,element_var1, element_var2,filename, normalization)
     struct All_variables *E;
     float **element_var1,**element_var2; // 
     char * filename;
	 double normalization;  // value used for normalization
{

    int i,d,ii,m,jj,j,k,nn;
    const int vpts = vpoints[E->mesh.nsd];
    char output_file[255],input_s[200];
    float temp1,temp2,temp3,temp4;
    // void visc_from_gint_to_nodes();
	void ele_to_nodes();
    void return_horiz_ave_f();
    void parallel_process_termination();
    void interpolate_to_FE_grid();
    FILE *fp1;

    static int been_here=0;
    int ee,el,ndepth,nfi,nth,nnfi,nnth,i1,i2,iskip,idepth;
    static double *radius_visc,*visc_reg1,*visc_reg2,*visc_FE[2];
    static float *ave_visc;
    double rr;

    sprintf(output_file,"%s",filename);
    fp1 = fopen(output_file,"r");
    fgets(input_s,200,fp1);   // read description
    fgets(input_s,200,fp1);
    sscanf(input_s,"%d %d %d",&ndepth,&nfi,&nth);
    fgets(input_s,200,fp1);   // read description
    nn = nfi*nth;
    nnth = nth +2;    // add 1 row and 1 column to each side 
    nnfi = nfi +2;

    if (been_here==0)  {
      visc_reg1 = (double *)malloc(((nfi+2)*(nth+2)+1)*sizeof(double));
      visc_reg2 = (double *)malloc(((nfi+2)*(nth+2)+1)*sizeof(double));
      radius_visc = (double *)malloc((ndepth+1)*sizeof(double));
      ave_visc = (float *)malloc((E->lmesh.noz+1)*sizeof(float));
      for (m=1;m<=E->sphere.caps_per_proc;m++)
         visc_FE[m] = (double *)malloc((E->lmesh.nsf+1)*sizeof(double));
      been_here = 1;
      }

    for (i=1;i<=ndepth;i++)  {
      fgets(input_s,200,fp1);
      sscanf(input_s,"%g",&temp1);
      radius_visc[i] = E->sphere.ro-temp1*1000.0/E->sphere.dradius;
//if(E->parallel.me==0) fprintf(E->fp_out,"visc model radius info %d %g\n",i,radius_visc[i]);
      }

    fclose(fp1);


  for(m=1;m<=E->sphere.caps_per_proc;m++)
    for (jj=1;jj<=E->lmesh.elz;jj++)  {
      rr = (E->sx[m][3][jj]+E->sx[m][3][jj+1])/2.0;

      if (rr>=radius_visc[ndepth]) {
//          continue;
        i1 = ndepth;
        i2 = ndepth;
        }
      else if (rr<=radius_visc[1])  {
//        continue;
         i1=1;
         i2=1;
        }
      else {
        for (i=1;i<ndepth;i++)  {
          if (rr>=radius_visc[i] && rr<=radius_visc[i+1]) {
            i1=i;
            i2=i+1;
            break;
            } 
          }
        }

      sprintf(output_file,"%s",filename);
      fp1 = fopen(output_file,"r");
      iskip = 3+ndepth+1+nfi+1+nth + (i1-1)*(nn+1); 
      for (i=1;i<=iskip;i++)  
        fgets(input_s,200,fp1);

      fgets(input_s,200,fp1);
      for (i=2;i<=nnth-1;i++)  
      for (j=2;j<=nnfi-1;j++)  {
        k = i + (j-1)*nnth;
        fgets(input_s,200,fp1);
        sscanf(input_s,"%g %g",&temp1,&temp2);
        visc_reg1[k] = temp1;       // another possibility is visc_reg[k] = log10f(temp1);
        visc_reg2[k] = temp2;       // another possibility is visc_reg[k] = log10f(temp1);
        }

      if (i1!=i2) {     // read one more layer and interpolate radially
        temp1 = (radius_visc[i2] - rr)/(radius_visc[i2] - radius_visc[i1]); 
        temp2 = (rr - radius_visc[i1])/(radius_visc[i2] - radius_visc[i1]); 
//if(E->parallel.me==0) fprintf(stderr,"%d %d %d rr= %g %g %g\n",jj,i1,i2,rr,temp1,temp2);
        fgets(input_s,200,fp1);
        for (i=2;i<=nnth-1;i++)  
        for (j=2;j<=nnfi-1;j++)  {
          k = i + (j-1)*nnth;
          fgets(input_s,200,fp1);
          sscanf(input_s,"%g %g",&temp3,&temp4);
          visc_reg1[k] = visc_reg1[k]*temp1 + temp3*temp2;  // also use log10f here if used above
          visc_reg2[k] = visc_reg2[k]*temp1 + temp4*temp2;  // also use log10f here if used above
          }
        }
      fclose(fp1);

      interpolate_to_FE_grid(E,visc_FE[m],visc_reg1,nnth,nnfi);

      for (i=1;i<=E->lmesh.ely;i++)  
      for (j=1;j<=E->lmesh.elx;j++)  {
         ee = j + (i-1)*E->lmesh.elx;      
         temp1 =( visc_FE[m][E->sien[m][ee].node[1]] 
                + visc_FE[m][E->sien[m][ee].node[2]] 
                + visc_FE[m][E->sien[m][ee].node[3]] 
                + visc_FE[m][E->sien[m][ee].node[4]])/4.0; 
         el = jj + (j-1)*E->lmesh.elz + (i-1)*E->lmesh.elz*E->lmesh.elx;
         element_var1[m][el] = temp1 / normalization;
         }

      interpolate_to_FE_grid(E,visc_FE[m],visc_reg2,nnth,nnfi);

      for (i=1;i<=E->lmesh.ely;i++)  
      for (j=1;j<=E->lmesh.elx;j++)  {
         ee = j + (i-1)*E->lmesh.elx;      
         temp2 =( visc_FE[m][E->sien[m][ee].node[1]] 
                + visc_FE[m][E->sien[m][ee].node[2]] 
                + visc_FE[m][E->sien[m][ee].node[3]] 
                + visc_FE[m][E->sien[m][ee].node[4]])/4.0; 
         el = jj + (j-1)*E->lmesh.elz + (i-1)*E->lmesh.elz*E->lmesh.elx;
         element_var2[m][el] = temp2 / normalization;
         }

      }   // end for jj of lmesh.elz

  return;
  }


void visc_from_mat(E,EEta)
     struct All_variables *E;
     float **EEta;
{

    int i,m,jj;
    const int vpts = vpoints[E->mesh.nsd];

    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(i=1;i<=E->lmesh.nel;i++)
            for(jj=1;jj<=vpoints[E->mesh.nsd];jj++)
                EEta[m][ (i-1)*vpts+jj ] = E->viscosity.N0[E->mat[m][i]-1];

   /*
      fprintf(E->fp_out,"output_evisc0 %d\n",E->monitor.solution_cycles);
      for(m=1;m<=E->sphere.caps_per_proc;m++) {
        fprintf(E->fp_out,"output_evisc0 for cap %d\n",E->sphere.capid[m]);
      for(i=1;i<=E->lmesh.elz;i++)
          fprintf(E->fp_out,"%d %d %f %f\n",i,E->mat[m][i],EEta[m][(i-1)*8+1],EEta[m][(i-1)*8+7]);
      }
      fflush(E->fp_out);
*/

    return;
}

void visc_from_T(E,EEta,propogate)
     struct All_variables *E;
     float **EEta;
     int propogate;
{
    int m,i,k,l,z,jj,kk;
    float zero,one,eta0,temp,tempa,TT[9];
    float zzz,zz[9],dr;
    float visc1, visc2;
    const int vpts = vpoints[E->mesh.nsd];
    const int ends = enodes[E->mesh.nsd];
    const int nel = E->lmesh.nel;

    one = 1.0;
    zero = 0.0;

    /* consistent handling : l is (material number - 1) to allow
       addressing viscosity arrays, which are all 0...n-1  */
    switch (E->viscosity.RHEOL)   {
    case 1:
        /* eta = N_0 exp( E * (T_0 - T))  */
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(i=1;i<=nel;i++)   {
                l = E->mat[m][i] - 1;

                if(E->control.mat_control==0)
                    tempa = E->viscosity.N0[l];
                else if(E->control.mat_control==1)
                    tempa = E->viscosity.N0[l]*E->VIP[m][i];

                for(kk=1;kk<=ends;kk++) {
                    TT[kk] = E->T[m][E->ien[m][i].node[kk]];
                }

                for(jj=1;jj<=vpts;jj++) {
                    temp=0.0;
                    for(kk=1;kk<=ends;kk++)   {
                        temp += TT[kk] * E->N.vpt[GNVINDEX(kk,jj)];
                    }

                    EEta[m][ (i-1)*vpts + jj ] = tempa*
                        exp( E->viscosity.E[l] * (E->viscosity.T[l] - temp));

                }
            }
        break;

    case 2:
        /* eta = N_0 exp(-T/T_0) */
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(i=1;i<=nel;i++)   {
                l = E->mat[m][i] - 1;

                if(E->control.mat_control==0)
                    tempa = E->viscosity.N0[l];
                else if(E->control.mat_control==1)
                    tempa = E->viscosity.N0[l]*E->VIP[m][i];

                for(kk=1;kk<=ends;kk++) {
                    TT[kk] = E->T[m][E->ien[m][i].node[kk]];
                }

                for(jj=1;jj<=vpts;jj++) {
                    temp=0.0;
                    for(kk=1;kk<=ends;kk++)   {
                        temp += TT[kk] * E->N.vpt[GNVINDEX(kk,jj)];
                    }

                    EEta[m][ (i-1)*vpts + jj ] = tempa*
                        exp( -temp / E->viscosity.T[l]);

                }
            }
        break;

    case 3:
        /* eta = N_0 exp(E/(T+T_0) - E/(1+T_0)) 

	   where T is normalized to be within 0...1

	 */
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(i=1;i<=nel;i++)   {
                l = E->mat[m][i] - 1;
		if(E->control.mat_control) /* switch moved up here TWB */
		  tempa = E->viscosity.N0[l] * E->VIP[m][i];
		else
		  tempa = E->viscosity.N0[l];

                for(kk=1;kk<=ends;kk++) {
		  TT[kk] = E->T[m][E->ien[m][i].node[kk]];
                }

                for(jj=1;jj<=vpts;jj++) {
                    temp=0.0;
                    for(kk=1;kk<=ends;kk++)   {	/* took out
						   computation of
						   depth, not needed
						   TWB */
		      TT[kk]=max(TT[kk],zero);
		      temp += min(TT[kk],one) * E->N.vpt[GNVINDEX(kk,jj)];
                    }
		    EEta[m][ (i-1)*vpts + jj ] = tempa*
		      exp( E->viscosity.E[l]/(temp+E->viscosity.T[l])
			   - E->viscosity.E[l]/(one +E->viscosity.T[l]) );
                }
            }
        break;

    case 4:
        /* eta = N_0 exp( (E + (1-z)Z_0) / (T+T_0) ) */
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(i=1;i<=nel;i++)   {
                l = E->mat[m][i] - 1;
		if(E->control.mat_control) /* moved this up here TWB */
		  tempa = E->viscosity.N0[l] * E->VIP[m][i];
		else
		  tempa = E->viscosity.N0[l];

                for(kk=1;kk<=ends;kk++) {
                    TT[kk] = E->T[m][E->ien[m][i].node[kk]];
                    zz[kk] = (1.-E->sx[m][3][E->ien[m][i].node[kk]]);
                }

                for(jj=1;jj<=vpts;jj++) {
                    temp=0.0;
                    zzz=0.0;
                    for(kk=1;kk<=ends;kk++)   {
                        TT[kk]=max(TT[kk],zero);
                        temp += min(TT[kk],one) * E->N.vpt[GNVINDEX(kk,jj)];
                        zzz += zz[kk] * E->N.vpt[GNVINDEX(kk,jj)];
                    }


		    EEta[m][ (i-1)*vpts + jj ] = tempa*
		      exp( (E->viscosity.E[l] +  E->viscosity.Z[l]*zzz )
			   / (E->viscosity.T[l]+temp) );

                }
            }
        break;

    case 100:
        /* user-defined viscosity law goes here */
        fprintf(stderr, "Need user definition for viscosity law: 'rheol=%d'\n",
                E->viscosity.RHEOL);
        parallel_process_termination();
        break;

    default:
        /* unknown option */
        fprintf(stderr, "Invalid value of 'rheol=%d'\n", E->viscosity.RHEOL);

        parallel_process_termination();
        break;
    }


    return;
}


void visc_from_S(E,EEta,propogate)
     struct All_variables *E;
     float **EEta;
     int propogate;
{
    double one,two,scale,temp2,temp1,scale1,depth,exponent1;

    void stress_2_inv();
    int m,e,l,z,jj,kk;
    static int been_here=0;

    const int vpts = vpoints[E->mesh.nsd];
    const int nel = E->lmesh.nel;

    one = 1.0;
    two = 2.0;

    if (been_here==0)  {
      for(m=1;m<=E->sphere.caps_per_proc;m++)  
        for(e=1;e<=nel;e++)/* initialize with unity if no velocities around */
          E->S2inv[m][e] = 1.0;
      been_here = 1;
      }
    else 
      stress_2_inv(E,E->S2inv,1,E->viscosity.iterate);

    for(m=1;m<=E->sphere.caps_per_proc;m++)  
        for(e=1;e<=nel;e++)   {
           l = E->mat[m][e] - 1;
           temp1 = (E->S2inv[m][e] + E->viscosity.sdepv_bg[l])/E->viscosity.sdepv_trns[l];
           exponent1= E->viscosity.sdepv_expt[l]-one;
           scale=one/(one+pow(temp1,exponent1));

           temp2 = E->viscosity.sdepv_bg[l]/E->viscosity.sdepv_trns[l];
           scale1=(one+pow(temp2,exponent1));
           scale = scale*scale1;

           for(jj=1;jj<=vpts;jj++)
                EEta[m][(e-1)*vpts + jj] = scale*EEta[m][(e-1)*vpts + jj];
        }

    return;
}

void visc_from_P(E,EEta) /* "plasticity" implementation

			 viscosity will be limited by a yield stress

			 \sigma_y  = min(a + b * (1-r), y)

			 where a,b,y are parameters input via pdepv_a,b,y

			 and

			 \eta_y = \sigma_y / (2 \eps_II)

			 where \eps_II is the second invariant. Then

			 \eta_eff = (\eta_0 \eta_y)/(\eta_0 + \eta_y)

			 for pdepv_eff = 1

			 or

			 \eta_eff = min(\eta_0,\eta_y)

			 for pdepv_eff = 0

			 where \eta_0 is the regular viscosity


			 TWB

			 */
     struct All_variables *E;
     float **EEta;
{
    float *eedot,zz[9],zzz,tau,eta_p,eta_new;
    int m,e,l,z,jj,kk;

    const int vpts = vpoints[E->mesh.nsd];
    const int nel = E->lmesh.nel;
    const int ends = enodes[E->mesh.nsd];

    void strain_rate_2_inv();

    eedot = (float *) malloc((2+nel)*sizeof(float));

    for(m=1;m<=E->sphere.caps_per_proc;m++)  {

      if(E->viscosity.pdepv_visited){

        strain_rate_2_inv(E,m,eedot,1);	/* get second invariant for all elements */

      }else{
	for(e=1;e<=nel;e++)	/* initialize with unity if no velocities around */
	  eedot[e] = 1.0;
	if(m == E->sphere.caps_per_proc)
	  E->viscosity.pdepv_visited = 1;
	if((E->parallel.me == 0)&&(E->control.verbose)){
	  for(e=0;e < E->viscosity.num_mat;e++)
	    fprintf(stderr,"num mat: %i a: %g b: %g y: %g\n",
		    e,E->viscosity.pdepv_a[e],E->viscosity.pdepv_b[e],E->viscosity.pdepv_y[e]);
	}
      }

      for(e=1;e <= nel;e++)   {	/* loop through all elements */

	l = E->mat[m][e] -1 ;	/* material of this element */

	for(kk=1;kk <= ends;kk++) /* nodal depths */
	  zz[kk] = (1.0 - E->sx[m][3][E->ien[m][e].node[kk]]); /* for depth, zz = 1 - r */

	for(jj=1;jj <= vpts;jj++){ /* loop through integration points */

	  zzz = 0.0;		/* get mean depth of integration point */
	  for(kk=1;kk<=ends;kk++)
	    zzz += zz[kk] * E->N.vpt[GNVINDEX(kk,jj)];

	  /* depth dependent yield stress */
	  tau = E->viscosity.pdepv_a[l] + zzz * E->viscosity.pdepv_b[l];

	  /* min of depth dep. and constant yield stress */
	  tau = min(tau,  E->viscosity.pdepv_y[l]);

	  /* yield viscosity */
	  eta_p = tau/(2.0 * eedot[e] + 1e-7) + E->viscosity.pdepv_offset;


	  if(E->viscosity.pdepv_eff){
	    /* two dashpots in series */
	    eta_new  = 1.0/(1.0/EEta[m][ (e-1)*vpts + jj ] + 1.0/eta_p);
	  }else{
	    /* min viscosities*/
	    eta_new  = min(EEta[m][ (e-1)*vpts + jj ], eta_p);
	  }
	  //fprintf(stderr,"z: %11g mat: %i a: %11g b: %11g y: %11g ee: %11g tau: %11g eta_p: %11g eta_new: %11g eta_old: %11g\n",
	  //zzz,l,E->viscosity.pdepv_a[l], E->viscosity.pdepv_b[l],E->viscosity.pdepv_y[l],
	  //eedot[e],tau,eta_p,eta_new,EEta[m][(e-1)*vpts + jj]);
	  EEta[m][(e-1)*vpts + jj] = eta_new;
        } /* end integration point loop */
      }	/* end element loop */

    } /* end caps loop */
    free ((void *)eedot);
    return;
}

/*

multiply with compositional factor which is determined by a geometric
mean average from the tracer composition, assuming two flavors and
compositions between zero and unity

*/
void visc_from_C( E, EEta)
     struct All_variables *E;
     float **EEta;
{
  double vmean,cc_loc[10],CC[10][9],cbackground;
  int m,l,z,jj,kk,i,p,q;


  const int vpts = vpoints[E->mesh.nsd];
  const int nel = E->lmesh.nel;
  const int ends = enodes[E->mesh.nsd];

  for(m=1;m <= E->sphere.caps_per_proc;m++)  {
    for(i = 1; i <= nel; i++){
      /* determine composition of each of the nodes of the
	 element */
        for(p=0; p<E->composition.ncomp; p++) {
            for(kk = 1; kk <= ends; kk++){
                CC[p][kk] = E->composition.comp_node[m][p][E->ien[m][i].node[kk]];
                if(CC[p][kk] < 0)CC[p][kk]=0.0;
                if(CC[p][kk] > 1)CC[p][kk]=1.0;
            }
        }
        for(jj = 1; jj <= vpts; jj++) {
            /* concentration of background material */
            cbackground = 1;
            for(p=0; p<E->composition.ncomp; p++) {
                /* compute mean composition  */
                cc_loc[p] = 0.0;
                for(kk = 1; kk <= ends; kk++) {
                    cc_loc[p] += CC[p][kk] * E->N.vpt[GNVINDEX(kk, jj)];
                }
                cbackground -= cc_loc[p];
            }

            /* geometric mean of viscosity */
            vmean = cbackground * E->viscosity.cdepv_ff[0];
            for(p=0; p<E->composition.ncomp; p++) {
                vmean += cc_loc[p] * E->viscosity.cdepv_ff[p+1];
            }
            vmean = exp(vmean);

            /* multiply the viscosity with this prefactor */
            EEta[m][ (i-1)*vpts + jj ] *= vmean;

        } /* end jj loop */
    } /* end el loop */
  } /* end cap */
}

void strain_rate_2_inv(E,m,EEDOT,SQRT)
     struct All_variables *E;
     float *EEDOT;
     int m,SQRT;
{
    void get_rtf_at_ppts();
    void velo_from_element();
    void construct_c3x3matrix_el();
    void get_ba_p();

    struct Shape_function_dx *GNx;

    double edot[4][4], rtf[4][9];
    double theta;
    double ba[9][9][4][7];
    float VV[4][9], Vxyz[7][9], dilation[9];
    
    int e, i, j, p, q, n;

    const int nel = E->lmesh.nel;
    const int dims = E->mesh.nsd;
    const int ends = enodes[dims];
    const int lev = E->mesh.levmax;
    const int ppts = ppoints[dims];
    const int sphere_key = 1;

    for(e=1; e<=nel; e++) {

        get_rtf_at_ppts(E, m, lev, e, rtf);
        velo_from_element(E, VV, m, e, sphere_key);
        GNx = &(E->gNX[m][e]);

        theta = rtf[1][1];


        /* Vxyz is the strain rate vector, whose relationship with
         * the strain rate tensor (e) is that:
         *    Vxyz[1] = e11
         *    Vxyz[2] = e22
         *    Vxyz[3] = e33
         *    Vxyz[4] = 2*e12
         *    Vxyz[5] = 2*e13
         *    Vxyz[6] = 2*e23
         * where 1 is theta, 2 is phi, and 3 is r
         */
        for(j=1; j<=ppts; j++) {
            Vxyz[1][j] = 0.0;
            Vxyz[2][j] = 0.0;
            Vxyz[3][j] = 0.0;
            Vxyz[4][j] = 0.0;
            Vxyz[5][j] = 0.0;
            Vxyz[6][j] = 0.0;
            dilation[j] = 0.0;
        }

        if ((E->control.precise_strain_rate) || (theta < 0.09) || (theta > 3.05)) {
            /* When the element is close to the poles, use a more
             * precise method to compute the strain rate. 
	     
	     if precise_strain_rate=on, will always choose this option

	    */

            if ((e-1)%E->lmesh.elz==0) {
                construct_c3x3matrix_el(E,e,&E->element_Cc,&E->element_Ccx,lev,m,1);
            }

            get_ba_p(&(E->N), GNx, &E->element_Cc, &E->element_Ccx,
                     rtf, E->mesh.nsd, ba);

            for(j=1;j<=ppts;j++)
                for(p=1;p<=6;p++)
                    for(i=1;i<=ends;i++)
                        for(q=1;q<=dims;q++) {
                            Vxyz[p][j] += ba[i][j][q][p] * VV[q][i];
                        }

        }
        else {
            for(j=1; j<=ppts; j++) {
                for(i=1; i<=ends; i++) {
                    Vxyz[1][j] += (VV[1][i] * GNx->ppt[GNPXINDEX(0, i, j)]
                                   + VV[3][i] * E->N.ppt[GNPINDEX(i, j)])
                        * rtf[3][j];
                    Vxyz[2][j] += ((VV[2][i] * GNx->ppt[GNPXINDEX(1, i, j)]
                                    + VV[1][i] * E->N.ppt[GNPINDEX(i, j)]
                                    * cos(rtf[1][j])) / sin(rtf[1][j])
                                   + VV[3][i] * E->N.ppt[GNPINDEX(i, j)])
                        * rtf[3][j];
                    Vxyz[3][j] += VV[3][i] * GNx->ppt[GNPXINDEX(2, i, j)];

                    Vxyz[4][j] += ((VV[1][i] * GNx->ppt[GNPXINDEX(1, i, j)]
                                    - VV[2][i] * E->N.ppt[GNPINDEX(i, j)]
                                    * cos(rtf[1][j])) / sin(rtf[1][j])
                                   + VV[2][i] * GNx->ppt[GNPXINDEX(0, i, j)])
                        * rtf[3][j];
                    Vxyz[5][j] += VV[1][i] * GNx->ppt[GNPXINDEX(2, i, j)]
                        + rtf[3][j] * (VV[3][i] * GNx->ppt[GNPXINDEX(0, i, j)]
                                       - VV[1][i] * E->N.ppt[GNPINDEX(i, j)]);
                    Vxyz[6][j] += VV[2][i] * GNx->ppt[GNPXINDEX(2, i, j)]
                        + rtf[3][j] * (VV[3][i]
                                       * GNx->ppt[GNPXINDEX(1, i, j)]
                                       / sin(rtf[1][j])
                                       - VV[2][i] * E->N.ppt[GNPINDEX(i, j)]);
                }
            }
        } /* end of else */

        if(E->control.inv_gruneisen != 0) {
            for(j=1; j<=ppts; j++)
                dilation[j] = (Vxyz[1][j] + Vxyz[2][j] + Vxyz[3][j]) / 3.0;
        }

        edot[1][1] = edot[2][2] = edot[3][3] = 0;
        edot[1][2] = edot[1][3] = edot[2][3] = 0;

        /* edot is 2 * (the deviatoric strain rate tensor) */
        for(j=1; j<=ppts; j++) {
            edot[1][1] += 2.0 * (Vxyz[1][j] - dilation[j]);
            edot[2][2] += 2.0 * (Vxyz[2][j] - dilation[j]);
            edot[3][3] += 2.0 * (Vxyz[3][j] - dilation[j]);
            edot[1][2] += Vxyz[4][j];
            edot[1][3] += Vxyz[5][j];
            edot[2][3] += Vxyz[6][j];
        }

        EEDOT[e] = edot[1][1] * edot[1][1]
            + edot[1][2] * edot[1][2] * 2.0
            + edot[2][2] * edot[2][2]
            + edot[2][3] * edot[2][3] * 2.0
            + edot[3][3] * edot[3][3]
            + edot[1][3] * edot[1][3] * 2.0;
    }

    if(SQRT)
	for(e=1;e<=nel;e++)
	    EEDOT[e] =  sqrt(0.5 *EEDOT[e]);
    else
	for(e=1;e<=nel;e++)
	    EEDOT[e] *=  0.5;

    return;
}


static void apply_low_visc_wedge_channel(struct All_variables *E, float **evisc)
{
    void parallel_process_termination();

    int i,j,m;
    const int vpts = vpoints[E->mesh.nsd];
    float *F;

    /* low viscosity channel/wedge require tracers to work */
    if(E->control.tracer == 0) {
        if(E->parallel.me == 0) {
            fprintf(stderr, "Error: low viscosity channel/wedge is turned on, "
                   "but tracer is off!\n");
            fprintf(E->fp, "Error: low viscosity channel/wedge is turned on, "
                   "but tracer is off!\n");
            fflush(E->fp);
        }
        parallel_process_termination();
    }


    F = (float *)malloc((E->lmesh.nel+1)*sizeof(float));
    for(i=1 ; i<=E->lmesh.nel ; i++)
        F[i] = 0.0;

    /* if low viscosity channel ... */
    if(E->viscosity.channel)
        low_viscosity_channel_factor(E, F);


    /* if low viscosity wedge ... */
    if(E->viscosity.wedge)
        low_viscosity_wedge_factor(E, F);


    for(i=1 ; i<=E->lmesh.nel ; i++) {
        if (F[i] != 0.0)
            for(m = 1 ; m <= E->sphere.caps_per_proc ; m++) {
                for(j=1;j<=vpts;j++) {
                    evisc[m][(i-1)*vpts + j] = F[i];
            }
        }
    }


    free(F);

    return;
}




static void low_viscosity_channel_factor(struct All_variables *E, float *F)
{
    int i, ii, k, m, e, ee;
    int nz_min[NCS], nz_max[NCS];
    const int flavor = 0;
    double rad_mean, rr;

    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        /* find index of radius corresponding to lv_min_radius */
        for(e=1; e<=E->lmesh.elz; e++) {
            rad_mean = 0.5 * (E->sx[m][3][E->ien[m][e].node[1]] +
                              E->sx[m][3][E->ien[m][e].node[8]]);
            if(rad_mean >= E->viscosity.lv_min_radius) break;
        }
        nz_min[m] = e;

        /* find index of radius corresponding to lv_max_radius */
        for(e=E->lmesh.elz; e>=1; e--) {
            rad_mean = 0.5 * (E->sx[m][3][E->ien[m][e].node[1]] +
                              E->sx[m][3][E->ien[m][e].node[8]]);
            if(rad_mean <= E->viscosity.lv_max_radius) break;
        }
        nz_max[m] = e;
    }



    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        for(k=1; k<=E->lmesh.elx*E->lmesh.ely; k++) {
            for(i=nz_min[m]; i<=nz_max[m]; i++) {
                e = (k-1)*E->lmesh.elz + i;

                rad_mean = 0.5 * (E->sx[m][3][E->ien[m][e].node[1]] +
                                  E->sx[m][3][E->ien[m][e].node[8]]);

                /* loop over elements below e */
                for(ii=i; ii>=nz_min[m]; ii--) {
                    ee = (k-1)*E->lmesh.elz + ii;

                    rr = 0.5 * (E->sx[m][3][E->ien[m][ee].node[1]] +
                                E->sx[m][3][E->ien[m][ee].node[8]]);

                    /* if ee has tracers in it and is within the channel */
                    if((E->trace.ntracer_flavor[m][flavor][ee] > 0) &&
                       (rad_mean <= rr + E->viscosity.lv_channel_thickness)) {
                           F[e] = E->viscosity.lv_reduction;
                           break;
                       }
                }
            }
        }
    }


    /** debug **
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(e=1; e<=E->lmesh.nel; e++)
            fprintf(stderr, "lv_reduction: %d %e\n", e, F[e]);
    */

    return;
}


static void low_viscosity_wedge_factor(struct All_variables *E, float *F)
{
    int i, ii, k, m, e, ee;
    int nz_min[NCS], nz_max[NCS];
    const int flavor = 0;
    double rad_mean, rr;

    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        /* find index of radius corresponding to lv_min_radius */
        for(e=1; e<=E->lmesh.elz; e++) {
            rad_mean = 0.5 * (E->sx[m][3][E->ien[m][e].node[1]] +
                              E->sx[m][3][E->ien[m][e].node[8]]);
            if(rad_mean >= E->viscosity.lv_min_radius) break;
        }
        nz_min[m] = e;

        /* find index of radius corresponding to lv_max_radius */
        for(e=E->lmesh.elz; e>=1; e--) {
            rad_mean = 0.5 * (E->sx[m][3][E->ien[m][e].node[1]] +
                              E->sx[m][3][E->ien[m][e].node[8]]);
            if(rad_mean <= E->viscosity.lv_max_radius) break;
        }
        nz_max[m] = e;
    }



    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        for(k=1; k<=E->lmesh.elx*E->lmesh.ely; k++) {
            for(i=nz_min[m]; i<=nz_max[m]; i++) {
                e = (k-1)*E->lmesh.elz + i;

                rad_mean = 0.5 * (E->sx[m][3][E->ien[m][e].node[1]] +
                                  E->sx[m][3][E->ien[m][e].node[8]]);

                /* loop over elements below e */
                for(ii=i; ii>=nz_min[m]; ii--) {
                    ee = (k-1)*E->lmesh.elz + ii;

                    /* if ee has tracers in it */
                    if(E->trace.ntracer_flavor[m][flavor][ee] > 0) {
                        F[e] = E->viscosity.lv_reduction;
                        break;
                    }
                }
            }
        }
    }


    /** debug **
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(e=1; e<=E->lmesh.nel; e++)
            fprintf(stderr, "lv_reduction: %d %e\n", e, F[e]);
    */

    return;
}


void stress_2_inv(E,eedot,SQRT,iterate)
    struct All_variables *E;
    float **eedot;
    int SQRT,iterate;
{

    void exchange_node_d();
    void get_rtf_at_ppts();
    void velo_from_element();
    void construct_c3x3matrix_el();
    void velo_from_element_d();
    int i,j,k,e,node,snode,m,nel2;
    
    double VV[4][9],Szz,Sxx,Syy,Sxy,Sxz,Szy;
    double Vxyz1,Vxyz2,Vxyz3,Vxyz4,Vxyz5,Vxyz6;
    double Sxyz1,Sxyz2,Sxyz3,Sxyz4,Sxyz5,Sxyz6;
    double sinaa,cosaa,ct,el_volume,alpha,onep_alpha,visc0,visc1,lamdaO,maxwell0,esmuO,a,b,rtf[4][9];
    double edot[4][4],dudx[4][4];
    double old_Skk;
    double lambda, elambda0, S2kk, mu, esmu, trace_edot;
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

   alpha = E->viscosity.alpha;
   onep_alpha = 1.0 - E->viscosity.alpha;

//fprintf(E->fp_out,"in_stress_2_inv %d %d\n",E->monitor.solution_cycles,E->viscosity.iterate);

  for(m=1;m<=E->sphere.caps_per_proc;m++) {
  for(e=1;e<=E->lmesh.nel;e++)  {

    if (iterate==0) {
        edot[1][1] = E->S2xx[m][e]; 
        edot[2][2] = E->S2yy[m][e]; 
        edot[3][3] = E->S2zz[m][e];
        edot[1][2] = E->S2xy[m][e];
        edot[1][3] = E->S2xz[m][e];
        edot[2][3] = E->S2zy[m][e];
        /*        E->S2xx[m][e] = E->Maxwelltime[m][e]*E->S2xx[m][e]; 
                E->S2yy[m][e] = E->Maxwelltime[m][e]*E->S2yy[m][e]; 
                E->S2zz[m][e] = E->Maxwelltime[m][e]*E->S2zz[m][e]; 
                E->S2xy[m][e] = E->Maxwelltime[m][e]*E->S2xy[m][e]; 
                E->S2xz[m][e] = E->Maxwelltime[m][e]*E->S2xz[m][e]; 
                E->S2zy[m][e] = E->Maxwelltime[m][e]*E->S2zy[m][e]; 
        */      
    }
    else {
      gnxx = E->gNX[m][e].ppt;
      get_rtf_at_ppts(E,m,E->mesh.levmax,e,rtf);
      if ((e-1)%E->lmesh.elz==0)
        construct_c3x3matrix_el(E,e,&Cc,&Ccx,E->mesh.levmax,m,1);
      velo_from_element_d(E,VV,m,e,sphere_key);

      // only for one vpt for each element

      sinaa = sin(rtf[1][1]);
      cosaa = cos(rtf[1][1]);
      ct = cosaa/sinaa;

      for(i=1;i<=ends;i++)   {
        for(k=1;k<=dims;k++)   {
          Sxyz1 = VV[k][i]*rtf[3][1]*
                 (gnxx[GNPXINDEX(0,i,1)]*Cc.ppt[BPINDEX(1,k,i,1)]
                 +E->N.ppt[GNPINDEX(i,1)]*Ccx.ppt[BPXINDEX(1,k,1,i,1)]
                 +E->N.ppt[GNPINDEX(i,1)]*Cc.ppt[BPINDEX(3,k,i,1)]);
          Sxyz2 = VV[k][i]*rtf[3][1]*
                (E->N.ppt[GNPINDEX(i,1)]*Cc.ppt[BPINDEX(1,k,i,1)]*ct
                +E->N.ppt[GNPINDEX(i,1)]*Cc.ppt[BPINDEX(3,k,i,1)]
                +(gnxx[GNPXINDEX(1,i,1)]*Cc.ppt[BPINDEX(2,k,i,1)]
                 +E->N.ppt[GNPINDEX(i,1)]*Ccx.ppt[BPXINDEX(2,k,2,i,1)])/sinaa);
          Sxyz3 = VV[k][i]*gnxx[GNPXINDEX(2,i,1)]*Cc.ppt[BPINDEX(3,k,i,1)];
          Sxyz4 = VV[k][i]*rtf[3][1]*
                (gnxx[GNPXINDEX(0,i,1)]*Cc.ppt[BPINDEX(2,k,i,1)]
                +E->N.ppt[GNPINDEX(i,1)]*Ccx.ppt[BPXINDEX(2,k,1,i,1)]
                -ct*Cc.ppt[BPINDEX(2,k,i,1)]*E->N.ppt[GNPINDEX(i,1)]
                +(E->N.ppt[GNPINDEX(i,1)]*Ccx.ppt[BPXINDEX(1,k,2,i,1)]
                 +gnxx[GNPXINDEX(1,i,1)]*Cc.ppt[BPINDEX(1,k,i,1)])/sinaa);
          Sxyz5 = VV[k][i]*rtf[3][1]*
                (gnxx[GNPXINDEX(2,i,1)]*Cc.ppt[BPINDEX(1,k,i,1)]/rtf[3][1]
                +(E->N.ppt[GNPINDEX(i,1)]*Ccx.ppt[BPXINDEX(3,k,1,i,1)]
                 +gnxx[GNPXINDEX(0,i,1)]*Cc.ppt[BPINDEX(3,k,i,1)]
                 -E->N.ppt[GNPINDEX(i,1)]*Cc.ppt[BPINDEX(1,k,i,1)]));
          Sxyz6 = VV[k][i]*
                (gnxx[GNPXINDEX(2,i,1)]*Cc.ppt[BPINDEX(2,k,i,1)]
                -rtf[3][1]*E->N.ppt[GNPINDEX(i,1)]*Cc.ppt[BPINDEX(2,k,i,1)]
                +rtf[3][1]/sinaa*(E->N.ppt[GNPINDEX(i,1)]*Ccx.ppt[BPXINDEX(3,k,2,i,1)]+gnxx[GNPXINDEX(1,i,1)]*Cc.ppt[BPINDEX(3,k,i,1)]));
        }
      }

      /*
      Concern: For iteration == 1, i.e., the second iteration, the EVolder is 1 for all elements.
                Is this good? Should we use Eold only for iteration == 1?
                e.g. if (iteration == 1) visc1 = E->EVold[m][e];
      */
      visc1=onep_alpha*E->EVolder[m][e]+alpha*E->EVold[m][e];
      a = E->esmu_o[E->mesh.levmax][m][e]/(2.0*visc1)*E->advection.timestep; 
      visc0 = E->esmu_o[E->mesh.levmax][m][e]/(1.0+a);
      maxwell0=(1.0-a)/(1.0+a);

      // Note at this point, the add_viscoelastic is not applied, so be careful with the variables like EVI, Maxwelltime
      if(!E->ve_data_cont.compressible){  // for incompressible
        edot[1][1] = 2.0*visc0*Sxyz1 + maxwell0*E->S2xx[m][e]; 
        edot[2][2] = 2.0*visc0*Sxyz2 + maxwell0*E->S2yy[m][e]; 
        edot[3][3] = 2.0*visc0*Sxyz3 + maxwell0*E->S2zz[m][e];
        edot[1][2] =    visc0*Sxyz4  + maxwell0*E->S2xy[m][e];
        edot[1][3] =    visc0*Sxyz5  + maxwell0*E->S2xz[m][e];
        edot[2][3] =    visc0*Sxyz6  + maxwell0*E->S2zy[m][e];
      } else {  // for compressible
        lambda = E->elambda_o[E->mesh.levmax][m][e]; // at this point, elambda is still lambda (non-dimensionlized)
        mu = E->esmu_o[E->mesh.levmax][m][e]; // at this point, esmu is still mu (non-dimensionlized)
        elambda0 = (lambda + (lambda + (2.0/3.0)* mu)*a)/(1.0+a);
        esmu = 2.0/3.0*a/(1.0+a);

        S2kk = E->S2xx[m][e] + E->S2yy[m][e] + E->S2zz[m][e]; // trace of previous stress

        edot[1][1] = 2.0*visc0*Sxyz1 + elambda0 * (Sxyz1 + Sxyz2 + Sxyz3) 
                        + maxwell0*E->S2xx[m][e] + esmu * S2kk;
        edot[2][2] = 2.0*visc0*Sxyz2 + elambda0 * (Sxyz1 + Sxyz2 + Sxyz3) 
                        + maxwell0*E->S2yy[m][e] + esmu * S2kk;
        edot[3][3] = 2.0*visc0*Sxyz3 + elambda0 * (Sxyz1 + Sxyz2 + Sxyz3) 
                        + maxwell0*E->S2zz[m][e] + esmu * S2kk;
        edot[1][2] =    visc0*Sxyz4  + maxwell0*E->S2xy[m][e];
        edot[1][3] =    visc0*Sxyz5  + maxwell0*E->S2xz[m][e];
        edot[2][3] =    visc0*Sxyz6  + maxwell0*E->S2zy[m][e];

      }
    }

    if(E->ve_data_cont.compressible){
      // find deviatoric part of edot for calculating the second invariant
      trace_edot = edot[1][1] + edot[2][2] + edot[3][3];
      edot[1][1] -= trace_edot/3.0;
      edot[2][2] -= trace_edot/3.0;
      edot[3][3] -= trace_edot/3.0;
    }

    eedot[m][e] = edot[1][1]*edot[1][1] + edot[1][2]*edot[1][2]*2.0
                + edot[2][2]*edot[2][2] + edot[2][3]*edot[2][3]*2.0
                + edot[3][3]*edot[3][3] + edot[1][3]*edot[1][3]*2.0;
  }    /* end for el */
  }     /* end for m */

if (SQRT)
  for(m=1;m<=E->sphere.caps_per_proc;m++)   
     for(e=1;e<=E->lmesh.nel;e++)
       eedot[m][e] = sqrt(0.5*eedot[m][e]);
else
  for(m=1;m<=E->sphere.caps_per_proc;m++)   
     for(e=1;e<=E->lmesh.nel;e++)
       eedot[m][e] = 0.5*eedot[m][e];

    return; 
  }
