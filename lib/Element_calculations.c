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
/* Functions to assemble the element k matrices and the element f vector.
   Note that for the regular grid case the calculation of k becomes repetitive
   to the point of redundancy. */

#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "material_properties.h"

/* else, PGI would complain */
void construct_side_c3x3matrix_el(struct All_variables *,int ,
				  struct CC *,struct CCX *,
				  int ,int ,int ,int );
void construct_c3x3matrix(struct All_variables *);
void construct_c3x3matrix_el (struct All_variables *,int ,struct CC *,
			      struct CCX *,int ,int ,int );
void assemble_div_u(struct All_variables *,
                    double **, double **, int );
void get_elt_tr(struct All_variables *, int , int , double [24], int );
void get_elt_tr_pseudo_surf(struct All_variables *, int , int , double [24], int );


void add_force(struct All_variables *E, int e, double elt_f[24], int m)
{
  const int dims=E->mesh.nsd;
  const int ends=enodes[E->mesh.nsd];
  int a, a1, a2, a3, p, node;

  for(a=1;a<=ends;a++)          {
    node = E->ien[m][e].node[a];
    p=(a-1)*dims;
    a1=E->id[m][node].doff[1];
    E->F[m][a1] += elt_f[p];
    a2=E->id[m][node].doff[2];
    E->F[m][a2] += elt_f[p+1];
    a3=E->id[m][node].doff[3];
    E->F[m][a3] += elt_f[p+2];
  }
}



/* ================================================================
   Function to assemble the global  F vector.
                     +
   Function to get the global H vector (mixed method driving terms)
   ================================================================ */

void assemble_forces(struct All_variables *E, int count)
 //    struct All_variables *E;
 //    int count;
{
  double elt_f[24],elt_ps[24],temp1;
  int m,a,e,i;

  void get_buoyancy();
  void get_elt_f();
  void get_elt_ps();
  void get_elt_tr();
  void strip_bcs_from_residual();
  double global_vdot();

  void load_boundary_forces();
  void get_potential();
  void append_boundary_forces();
  void append_multi_boundary_forces();
  void get_oceanload();
  void get_dynamic_oceanload();
 
  void parallel_process_termination();
  void get_boundary_forces();

  double zero=0;
  int p,node,a1,a2,a3;

  static double *F0[NCS];
  static int been_here = 0;

  const int neq=E->lmesh.neq;
  const int nel=E->lmesh.nel;
  const int lev=E->mesh.levmax;
  const int ends = enodes[E->mesh.nsd];
  const int dims = E->mesh.nsd;
  const int n=loc_mat_size[dims];

  double CPU_time0(),time;


//  get_buoyancy(E,E->buoyancy);
//
   if (been_here==0)   {  // initialize array F0
        for(m=1;m<=E->sphere.caps_per_proc;m++)   {
            F0[m]= (double *)malloc((neq+1)*sizeof(double));
            for(a=0;a<neq;a++)  F0[m][a] = zero;
        }
        been_here ++;
    }

  for (i=0;i<24;i++)
     elt_f[i] = 0.0;

// Set slice.surf/botm[2] for use in get_potential:
   get_boundary_forces(E,count);

      //Sets E->potential[1][m][j] and E->potential[0][m][j]
   if (E->ve_data_cont.SELFG)
      get_potential(E,count);

       // calculate E->slice.oceanload according to sea level equation
       //  (based on E->potential[0] and E->slice.surf[3])
	if (E->ve_data_cont.SLE)
      		get_dynamic_oceanload(E,count);

    // Get height of load above geoid, and put it in
    //   E->slice.load[1][m][j]  and   E->slice.load[3][m][j]
    //>>>> modify by tao:
    if(E->ve_data_cont.compressible){
			append_multi_boundary_forces(E);
    	}
	else{
   			append_boundary_forces(E); 
		}

	//<<<< end tao

  if (count==0)
   for(m=1;m<=E->sphere.caps_per_proc;m++)    {
        for(a=0;a<neq;a++)  F0[m][a] = zero;
        for (e=1;e<=nel;e++)  {

//            get_elt_f(E,e,elt_f,1,m);
              get_elt_ps(E,e,elt_ps,1,m);

            for(a=1;a<=ends;a++)  {
                node = E->ien[m][e].node[a];
                p=(a-1)*dims;
                a1=E->id[m][node].doff[1];
                F0[m][a1] += /*elt_f[p] +*/ elt_ps[p];
                a2=E->id[m][node].doff[2];
                F0[m][a2] += /*elt_f[p+1] +*/ elt_ps[p+1];
                a3=E->id[m][node].doff[3];
                F0[m][a3] += /*elt_f[p+2] +*/ elt_ps[p+2];
                }
             }
          }     /* end for m, if (count==0) */


     // for elements that are at the surface or CMB or density interfaces:
  for(m=1;m<=E->sphere.caps_per_proc;m++)    {
    for(a=0;a<neq;a++)
      E->F[m][a] = 0.0;
    for (e=1;e<=nel;e++)  {
           for (p=0;p<n;p++) elt_f[p] = zero;
		   if(E->ve_data_cont.compressible){
				load_boundary_forces(E,e,elt_f,m); 
				for(a=1;a<=ends;a++)		  {
					 node = E->ien[m][e].node[a];
					 p=(a-1)*dims;
					 a1=E->id[m][node].doff[1];
					 E->F[m][a1] += elt_f[p];
					 a2=E->id[m][node].doff[2];
					 E->F[m][a2] += elt_f[p+1];
					 a3=E->id[m][node].doff[3];
					 E->F[m][a3] += elt_f[p+2];
				 }
		   	}
		   else{
		   		if (   e%E->lmesh.elz==0 		// element at surface
				   || (e-1)%E->lmesh.elz==0  )	{	// element at CMB
		           load_boundary_forces(E,e,elt_f,m); // uses slice.load[1] and [3]
		           for(a=1;a<=ends;a++)          {
		                node = E->ien[m][e].node[a];
		                p=(a-1)*dims;
		                a1=E->id[m][node].doff[1];
		                E->F[m][a1] += elt_f[p];
		                a2=E->id[m][node].doff[2];
		                E->F[m][a2] += elt_f[p+1];
		                a3=E->id[m][node].doff[3];
		                E->F[m][a3] += elt_f[p+2];
		            } 
			  }
		   	}
      } 
  }       /* end for m */

  E->monitor.fdotf = sqrt(global_vdot(E, E->F, E->F, lev));
  for(m=1;m<=E->sphere.caps_per_proc;m++)
    for(a=0;a<neq;a++) {
        E->F[m][a] = E->F[m][a] + F0[m][a];
     }

  (E->solver.exchange_id_d)(E, E->F, lev);
  strip_bcs_from_residual(E,E->F,lev);

  /* compute the norm of E->F */
  temp1 = sqrt(global_vdot(E, F0, F0, lev));

  if(E->parallel.me==0) {
      fprintf(stderr, "Momentum equation force %.9e %.9e\n",
              E->monitor.fdotf,temp1);
      fprintf(E->fp, "Momentum equation force %.9e %.9e\n",
              E->monitor.fdotf,temp1);
  }


  return;
}


void assemble_forces_pseudo_surf(E,penalty)
     struct All_variables *E;
     int penalty;
{
  double elt_f[24];
  int m,a,e,i;

  void get_buoyancy();
  void get_elt_f();
  void get_elt_tr_pseudo_surf();
  void strip_bcs_from_residual();

  const int neq=E->lmesh.neq;
  const int nel=E->lmesh.nel;
  const int lev=E->mesh.levmax;

  get_buoyancy(E,E->buoyancy);

  for(m=1;m<=E->sphere.caps_per_proc;m++)    {

    for(a=0;a<neq;a++)
      E->F[m][a] = 0.0;

    for (e=1;e<=nel;e++)  {
      get_elt_f(E,e,elt_f,1,m);
      add_force(E, e, elt_f, m);
    }

    /* for traction bc */
    for(i=1; i<=E->boundary.nel; i++) {
      e = E->boundary.element[m][i];

      for(a=0;a<24;a++) elt_f[a] = 0.0;
      for(a=SIDE_BEGIN; a<=SIDE_END; a++)
	get_elt_tr_pseudo_surf(E, i, a, elt_f, m);

      add_force(E, e, elt_f, m);
    }
  }       /* end for m */

  (E->solver.exchange_id_d)(E, E->F, lev);
  strip_bcs_from_residual(E,E->F,lev);
  return;
}


/*==============================================================
  Function to supply the element strain-displacement matrix Ba at velocity
  quadrature points, which is used to compute element stiffness matrix
  ==============================================================  */

void get_ba(struct Shape_function *N, struct Shape_function_dx *GNx,
       struct CC *cc, struct CCX *ccx, double rtf[4][9],
       int dims, double ba[9][9][4][7])
{
    int k, a, n;
    const int vpts = VPOINTS3D;
    const int ends = ENODES3D;

    double ra[9], isi[9], ct[9];
    double gnx0, gnx1, gnx2, shp, cc1, cc2, cc3;

    for(k=1;k<=vpts;k++) {
    ra[k] = rtf[3][k];
    isi[k] = 1.0 / sin(rtf[1][k]);
    ct[k] = cos(rtf[1][k]) * isi[k];
    }

    for(a=1;a<=ends;a++)
        for(k=1;k<=vpts;k++) {
            gnx0 = GNx->vpt[GNVXINDEX(0,a,k)];
            gnx1 = GNx->vpt[GNVXINDEX(1,a,k)];
            gnx2 = GNx->vpt[GNVXINDEX(2,a,k)];
            shp  = N->vpt[GNVINDEX(a,k)];
            for(n=1;n<=dims;n++) {
                cc1 = cc->vpt[BVINDEX(1,n,a,k)];
                cc2 = cc->vpt[BVINDEX(2,n,a,k)];
                cc3 = cc->vpt[BVINDEX(3,n,a,k)];

        ba[a][k][n][1] = ( gnx0 * cc1
                                   + shp * ccx->vpt[BVXINDEX(1,n,1,a,k)]
                                   + shp * cc3 ) * ra[k];

        ba[a][k][n][2] = ( shp * cc1 * ct[k]
                                   + shp * cc3
                                   + ( gnx1 * cc2
                                       + shp * ccx->vpt[BVXINDEX(2,n,2,a,k)] )
                                   * isi[k] ) * ra[k];

        ba[a][k][n][3] = gnx2 * cc3;

        ba[a][k][n][4] = ( gnx0 * cc2
                                   + shp * ccx->vpt[BVXINDEX(2,n,1,a,k)]
                                   - shp * cc2 * ct[k]
                                   + ( gnx1 * cc1
                                       + shp * ccx->vpt[BVXINDEX(1,n,2,a,k)] )
                                   * isi[k] ) * ra[k];

        ba[a][k][n][5] = gnx2 * cc1
                               + ( gnx0 * cc3
                                   + shp * ( ccx->vpt[BVXINDEX(3,n,1,a,k)]
                                             - cc1 ) ) * ra[k];

        ba[a][k][n][6] = gnx2 * cc2
                               - ra[k] * shp * cc2
                               + ( gnx1 * cc3
                                   + shp * ccx->vpt[BVXINDEX(3,n,2,a,k)] )
                               * isi[k] * ra[k];
        }
        }

    return;
}


/*==============================================================
  Function to supply the element strain-displacement matrix Ba at pressure
  quadrature points, which is used to compute strain rate
  ==============================================================  */

void get_ba_p(struct Shape_function *N, struct Shape_function_dx *GNx,
              struct CC *cc, struct CCX *ccx, double rtf[4][9],
              int dims, double ba[9][9][4][7])
{
    int k, a, n;
    const int ppts = PPOINTS3D;
    const int ends = ENODES3D;

    double ra[9], isi[9], ct[9];
    double gnx0, gnx1, gnx2, shp, cc1, cc2, cc3;

    for(k=1;k<=ppts;k++) {
        ra[k] = rtf[3][k];
        isi[k] = 1.0 / sin(rtf[1][k]);
        ct[k] = cos(rtf[1][k]) * isi[k];
    }

    for(k=1;k<=ppts;k++)
        for(a=1;a<=ends;a++) {
            gnx0 = GNx->ppt[GNPXINDEX(0,a,k)];
            gnx1 = GNx->ppt[GNPXINDEX(1,a,k)];
            gnx2 = GNx->ppt[GNPXINDEX(2,a,k)];
            shp  = N->ppt[GNPINDEX(a,k)];
            for(n=1;n<=dims;n++) {
                cc1 = cc->ppt[BPINDEX(1,n,a,k)];
                cc2 = cc->ppt[BPINDEX(2,n,a,k)];
                cc3 = cc->ppt[BPINDEX(3,n,a,k)];

        ba[a][k][n][1] = ( gnx0 * cc1
                           + shp * ccx->ppt[BPXINDEX(1,n,1,a,k)]
                           + shp * cc3 ) * ra[k];

        ba[a][k][n][2] = ( shp * cc1 * ct[k]
                           + shp * cc3
                           + ( gnx1 * cc2
                               + shp * ccx->ppt[BPXINDEX(2,n,2,a,k)] )
                           * isi[k] ) * ra[k];

        ba[a][k][n][3] = gnx2 * cc3;

        ba[a][k][n][4] = ( gnx0 * cc2
                           + shp * ccx->ppt[BPXINDEX(2,n,1,a,k)]
                           - shp * cc2 * ct[k]
                           + ( gnx1 * cc1
                               + shp * ccx->ppt[BPXINDEX(1,n,2,a,k)] )
                           * isi[k] ) * ra[k];

        ba[a][k][n][5] = gnx2 * cc1
                       + ( gnx0 * cc3
                           + shp * ( ccx->ppt[BPXINDEX(3,n,1,a,k)]
                                     - cc1 ) ) * ra[k];

        ba[a][k][n][6] = gnx2 * cc2
                       - ra[k] * shp * cc2
                       + ( gnx1 * cc3
                           + shp * ccx->ppt[BPXINDEX(3,n,2,a,k)] )
                       * isi[k] * ra[k];
            }
        }
    return;
}



/*==============================================================
  Function to supply the element k matrix for a given element e.
  ==============================================================  */

void get_elt_k(E,el,elt_k,lev,m,iconv)
     struct All_variables *E;
     int el,m;
     double elt_k[24*24];
     int lev, iconv;
{
    double bdbmu[4][4];
    int pn,qn,ad,bd;

    int a,b,i,j,i1,j1,k;
    double rtf[4][9],W[9],W1[9],W2[9];
    double stress_scaling;

    const double two = 2.0;
    const double two_thirds = 2.0/3.0;

    void get_rtf_at_vpts();

    double ba[9][9][4][7]; /* integration points,node,3x6 matrix */

    const int nn=loc_mat_size[E->mesh.nsd];
    const int vpts = VPOINTS3D;
    const int ends = ENODES3D;
    const int dims=E->mesh.nsd;

    stress_scaling = E->data.density*E->data.grav_acc*E->sphere.dradius/E->ve_data_cont.shear_mod;
    get_rtf_at_vpts(E, m, lev, el, rtf);

    if (iconv || (el-1)%E->lmesh.ELZ[lev]==0)
      construct_c3x3matrix_el(E,el,&E->element_Cc,&E->element_Ccx,lev,m,0);

    /* Note N[a].gauss_pt[n] is the value of shape fn a at the nth gaussian
       quadrature point. Nx[d] is the derivative wrt x[d]. */

    for(k=1;k<=vpts;k++) {
      W[k]=g_point[k].weight[dims-1]*E->GDA[lev][m][el].vpt[k]*E->EVI[lev][m][(el-1)*vpts+k];
    }

   if(E->ve_data_cont.compressible)  {
     for(k=1;k<=vpts;k++) {
       W1[k]=g_point[k].weight[dims-1]*E->GDA[lev][m][el].vpt[k]*E->elambda[lev][m][el];
       W2[k]=g_point[k].weight[dims-1]*E->GDA[lev][m][el].vpt[k]*E->erho[lev][m][el]*E->egrav[lev][m][el]*stress_scaling;
       }
    }

    get_ba(&(E->N), &(E->GNX[lev][m][el]), &E->element_Cc, &E->element_Ccx,
           rtf, E->mesh.nsd, ba);


  for(a=1;a<=ends;a++)
    for(b=a;b<=ends;b++)   {
      bdbmu[1][1]=bdbmu[1][2]=bdbmu[1][3]=
      bdbmu[2][1]=bdbmu[2][2]=bdbmu[2][3]=
      bdbmu[3][1]=bdbmu[3][2]=bdbmu[3][3]=0.0;

      for(i=1;i<=dims;i++)
        for(j=1;j<=dims;j++)
          for(k=1;k<=VPOINTS3D;k++)
              bdbmu[i][j] += W[k] * ( two * ( ba[a][k][i][1]*ba[b][k][j][1] +
                                              ba[a][k][i][2]*ba[b][k][j][2] +
                                              ba[a][k][i][3]*ba[b][k][j][3] ) +
                                      ba[a][k][i][4]*ba[b][k][j][4] +
                                      ba[a][k][i][5]*ba[b][k][j][5] +
                                      ba[a][k][i][6]*ba[b][k][j][6] );

/*      if(E->control.inv_gruneisen != 0)
        for(i=1;i<=dims;i++)
          for(j=1;j<=dims;j++)
            for(k=1;k<=VPOINTS3D;k++)
                bdbmu[i][j] -= W[k] * two_thirds *
                    ( ba[a][k][i][1] + ba[a][k][i][2] + ba[a][k][i][3] ) *
                    ( ba[b][k][j][1] + ba[b][k][j][2] + ba[b][k][j][3] );
*/

      if(E->ve_data_cont.compressible)  {
        for(i=1;i<=dims;i++)
          for(j=1;j<=dims;j++)
            for(k=1;k<=VPOINTS3D;k++)  {
// add in the lambda related terms for compressible rheology 
                        bdbmu[i][j] += W1[k]*(ba[a][k][i][1]+ba[a][k][i][2]+ba[a][k][i][3])*
                                             (ba[b][k][j][1]+ba[b][k][j][2]+ba[b][k][j][3]);
// add in the -( div_w*rho*g*v_r + w_r*rho*g*div_v ) terms  
                        bdbmu[i][j] -= W2[k]*(ba[a][k][i][1]+ba[a][k][i][2]+ba[a][k][i][3])*
                                 E->N.vpt[GNVINDEX(b,k)]*E->element_Cc.vpt[BVINDEX(3,j,b,k)];
                        bdbmu[i][j] -= W2[k]*(ba[b][k][j][1]+ba[b][k][j][2]+ba[b][k][j][3])*
                                 E->N.vpt[GNVINDEX(a,k)]*E->element_Cc.vpt[BVINDEX(3,i,a,k)];
   
         	}
         }


                /**/
      ad=dims*(a-1);
      bd=dims*(b-1);

      pn=ad*nn+bd;
      qn=bd*nn+ad;

      elt_k[pn       ] = bdbmu[1][1] ; /* above */
      elt_k[pn+1     ] = bdbmu[1][2] ;
      elt_k[pn+2     ] = bdbmu[1][3] ;
      elt_k[pn+nn    ] = bdbmu[2][1] ;
      elt_k[pn+nn+1  ] = bdbmu[2][2] ;
      elt_k[pn+nn+2  ] = bdbmu[2][3] ;
      elt_k[pn+2*nn  ] = bdbmu[3][1] ;
      elt_k[pn+2*nn+1] = bdbmu[3][2] ;
      elt_k[pn+2*nn+2] = bdbmu[3][3] ;

      elt_k[qn       ] = bdbmu[1][1] ; /* below diag */
      elt_k[qn+1     ] = bdbmu[2][1] ;
      elt_k[qn+2     ] = bdbmu[3][1] ;
      elt_k[qn+nn    ] = bdbmu[1][2] ;
      elt_k[qn+nn+1  ] = bdbmu[2][2] ;
      elt_k[qn+nn+2  ] = bdbmu[3][2] ;
      elt_k[qn+2*nn  ] = bdbmu[1][3] ;
      elt_k[qn+2*nn+1] = bdbmu[2][3] ;
      elt_k[qn+2*nn+2] = bdbmu[3][3] ;
                /**/

      } /*  Sum over all the a,b's to obtain full  elt_k matrix */


    return;
}


/* =============================================
   General calling function for del_squared:
   according to whether it should be element by
   element or node by node.
   ============================================= */

void assemble_del2_u(E,u,Au,level,strip_bcs)
     struct All_variables *E;
     double **u,**Au;
     int level;
     int strip_bcs;
{
  void e_assemble_del2_u();
  void n_assemble_del2_u();

  if(E->control.NMULTIGRID||E->control.NASSEMBLE)
    n_assemble_del2_u(E,u,Au,level,strip_bcs);
  else
    e_assemble_del2_u(E,u,Au,level,strip_bcs);

  return;
}

/* ======================================
   Assemble del_squared_u vector el by el
   ======================================   */

void e_assemble_del2_u(E,u,Au,level,strip_bcs)
  struct All_variables *E;
  double **u,**Au;
  int level;
  int strip_bcs;

{
  int  e,i,a,b,a1,a2,a3,ii,m,nodeb;
  void strip_bcs_from_residual();

  const int n=loc_mat_size[E->mesh.nsd];
  const int ends=enodes[E->mesh.nsd];
  const int dims=E->mesh.nsd;
  const int nel=E->lmesh.NEL[level];
  const int neq=E->lmesh.NEQ[level];

  for (m=1;m<=E->sphere.caps_per_proc;m++)   {
    for(i=0;i<neq;i++)
      Au[m][i] = 0.0;

    for(e=1;e<=nel;e++)   {
      for(a=1;a<=ends;a++) {
	ii = E->IEN[level][m][e].node[a];
	a1 = E->ID[level][m][ii].doff[1];
	a2 = E->ID[level][m][ii].doff[2];
	a3 = E->ID[level][m][ii].doff[3];
	for(b=1;b<=ends;b++) {
	        nodeb = E->IEN[level][m][e].node[b];
	        ii = (a*n+b)*dims-(dims*n+dims);
			/* i=1, j=1,2 */
	          /* i=1, j=1,2,3 */
		Au[m][a1] +=
		        E->elt_k[level][m][e].k[ii] *
			u[m][E->ID[level][m][nodeb].doff[1]]
		      + E->elt_k[level][m][e].k[ii+1] *
			u[m][E->ID[level][m][nodeb].doff[2]]
		      + E->elt_k[level][m][e].k[ii+2] *
			u[m][E->ID[level][m][nodeb].doff[3]];
		/* i=2, j=1,2,3 */
		Au[m][a2] +=
		        E->elt_k[level][m][e].k[ii+n] *
			u[m][E->ID[level][m][nodeb].doff[1]]
		      + E->elt_k[level][m][e].k[ii+n+1] *
			u[m][E->ID[level][m][nodeb].doff[2]]
		      + E->elt_k[level][m][e].k[ii+n+2] *
			u[m][E->ID[level][m][nodeb].doff[3]];
		/* i=3, j=1,2,3 */
		Au[m][a3] +=
		        E->elt_k[level][m][e].k[ii+n+n] *
			u[m][E->ID[level][m][nodeb].doff[1]]
		      + E->elt_k[level][m][e].k[ii+n+n+1] *
			u[m][E->ID[level][m][nodeb].doff[2]]
		      + E->elt_k[level][m][e].k[ii+n+n+2] *
			u[m][E->ID[level][m][nodeb].doff[3]];

 	    }         /* end for loop b */
        }             /* end for loop a */

       }          /* end for e */
     }         /* end for m  */

    (E->solver.exchange_id_d)(E, Au, level);

  if(strip_bcs)
     strip_bcs_from_residual(E,Au,level);

  return; }


/* ======================================================
   Assemble Au using stored, nodal coefficients.
   ====================================================== */

void n_assemble_del2_u(E,u,Au,level,strip_bcs)
     struct All_variables *E;
     double **u,**Au;
     int level;
     int strip_bcs;
{
    int m, e,i;
    int eqn1,eqn2,eqn3;

    double UU,U1,U2,U3;
    void strip_bcs_from_residual();
    void add_restoring();
	void add_restoring_comp();

    int *C;
    higher_precision *B1,*B2,*B3;

    const int neq=E->lmesh.NEQ[level];
    const int nno=E->lmesh.NNO[level];
    const int dims=E->mesh.nsd;
    const int max_eqn = dims*14;


  for (m=1;m<=E->sphere.caps_per_proc;m++)  {

     for(e=0;e<=neq;e++)
	Au[m][e]=0.0;

     if (strip_bcs){
	 	if(E->ve_data_cont.compressible)
			add_restoring_comp(E,level,m,Au,u);
	 	else
        	add_restoring(E,level,m,Au,u);
     	}
	 
     u[m][neq] = 0.0;

     for(e=1;e<=nno;e++)     {

       eqn1=E->ID[level][m][e].doff[1];
       eqn2=E->ID[level][m][e].doff[2];
       eqn3=E->ID[level][m][e].doff[3];

       U1 = u[m][eqn1];
       U2 = u[m][eqn2];
       U3 = u[m][eqn3];

       C=E->Node_map[level][m] + (e-1)*max_eqn;
       B1=E->Eqn_k1[level][m]+(e-1)*max_eqn;
       B2=E->Eqn_k2[level][m]+(e-1)*max_eqn;
       B3=E->Eqn_k3[level][m]+(e-1)*max_eqn;

       for(i=3;i<max_eqn;i++)  {
	  UU = u[m][C[i]];
  	  Au[m][eqn1] += B1[i]*UU;
  	  Au[m][eqn2] += B2[i]*UU;
  	  Au[m][eqn3] += B3[i]*UU;
       }
       for(i=0;i<max_eqn;i++)
          Au[m][C[i]] += B1[i]*U1+B2[i]*U2+B3[i]*U3;

       }     /* end for e */
     }     /* end for m */

     (E->solver.exchange_id_d)(E, Au, level);

//    if (strip_bcs)
	strip_bcs_from_residual(E,Au,level);

    return;
}


void build_diagonal_of_K(E,el,elt_k,level,m)
     struct All_variables *E;
     int level,el,m;
     double elt_k[24*24];

{
    int a,a1,a2,p,node;

    const int n=loc_mat_size[E->mesh.nsd];
    const int dims=E->mesh.nsd;
    const int ends=enodes[E->mesh.nsd];

    for(a=1;a<=ends;a++) {
	    node=E->IEN[level][m][el].node[a];
	    /* dirn 1 */
	    a1 = E->ID[level][m][node].doff[1];
	    p=(a-1)*dims;
	    E->BI[level][m][a1] += elt_k[p*n+p];

	    /* dirn 2 */
	    a2 = E->ID[level][m][node].doff[2];
	    p=(a-1)*dims+1;
	    E->BI[level][m][a2] += elt_k[p*n+p];

	    /* dirn 3 */
	    a1 = E->ID[level][m][node].doff[3];
	    p=(a-1)*dims+2;
	    E->BI[level][m][a1] += elt_k[p*n+p];
            }

  return;
}

void build_diagonal_of_Ahat(E)
    struct All_variables *E;
{
    double assemble_dAhatp_entry();

    double BU;
    int m,e,npno,neq,level;

 for (level=E->mesh.gridmin;level<=E->mesh.gridmax;level++)

   for (m=1;m<=E->sphere.caps_per_proc;m++)    {

     npno = E->lmesh.NPNO[level];
     neq=E->lmesh.NEQ[level];

     for(e=1;e<=npno;e++)
	E->BPI[level][m][e]=1.0;

     if(!E->control.precondition)
	return;

     for(e=1;e<=npno;e++)  {
	BU=assemble_dAhatp_entry(E,e,level,m);
	if(BU != 0.0)
	    E->BPI[level][m][e] = 1.0/BU;
	else
	    E->BPI[level][m][e] = 1.0;
        }
     }

    return;
}


/* =====================================================
   Assemble grad(rho_ref*ez)*V element by element.
   Note that the storage is not zero'd before assembling.
   =====================================================  */

void assemble_c_u(struct All_variables *E,
                  double **U, double **result, int level)
{
    int e,j1,j2,j3,p,a,b,m;

    const int nel = E->lmesh.NEL[level];
    const int ends = enodes[E->mesh.nsd];
    const int dims = E->mesh.nsd;
    const int npno = E->lmesh.NPNO[level];

    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(a=1;a<=ends;a++) {
            p = (a-1)*dims;
            for(e=1;e<=nel;e++) {
                b = E->IEN[level][m][e].node[a];
                j1= E->ID[level][m][b].doff[1];
                j2= E->ID[level][m][b].doff[2];
                j3= E->ID[level][m][b].doff[3];
                result[m][e] += E->elt_c[level][m][e].c[p  ][0] * U[m][j1]
                              + E->elt_c[level][m][e].c[p+1][0] * U[m][j2]
                              + E->elt_c[level][m][e].c[p+2][0] * U[m][j3];
            }
        }

    return;
}



/* =====================================================
   Assemble div(rho_ref*V) = div(V) + grad(rho_ref*ez)*V
   element by element
   =====================================================  */

void assemble_div_rho_u(struct All_variables *E,
                        double **U, double **result, int level)
{
    void assemble_div_u();
    assemble_div_u(E, U, result, level);
    assemble_c_u(E, U, result, level);

    return;
}


/* ==========================================
   Assemble a div_u vector element by element
   ==========================================  */

void assemble_div_u(struct All_variables *E,
                    double **U, double **divU, int level)
{
    int e,j1,j2,j3,p,a,b,m;

    const int nel=E->lmesh.NEL[level];
    const int ends=enodes[E->mesh.nsd];
    const int dims=E->mesh.nsd;
    const int npno=E->lmesh.NPNO[level];

  for(m=1;m<=E->sphere.caps_per_proc;m++)
    for(e=1;e<=npno;e++)
	divU[m][e] = 0.0;

  for(m=1;m<=E->sphere.caps_per_proc;m++)
       for(a=1;a<=ends;a++)   {
	  p = (a-1)*dims;
          for(e=1;e<=nel;e++) {
	    b = E->IEN[level][m][e].node[a];
	    j1= E->ID[level][m][b].doff[1];
            j2= E->ID[level][m][b].doff[2];
	    j3= E->ID[level][m][b].doff[3];
	    divU[m][e] += E->elt_del[level][m][e].g[p  ][0] * U[m][j1]
	                + E->elt_del[level][m][e].g[p+1][0] * U[m][j2]
	                + E->elt_del[level][m][e].g[p+2][0] * U[m][j3];
	    }
	 }

    return;
}


/* ==========================================
   Assemble a grad_P vector element by element
   ==========================================  */

void assemble_grad_p(E,P,gradP,lev)
     struct All_variables *E;
     double **P,**gradP;
     int lev;

{
  int m,e,i,j1,j2,j3,p,a,b,nel,neq;
  void strip_bcs_from_residual();

  const int ends=enodes[E->mesh.nsd];
  const int dims=E->mesh.nsd;

  for(m=1;m<=E->sphere.caps_per_proc;m++)  {

    nel=E->lmesh.NEL[lev];
    neq=E->lmesh.NEQ[lev];

    for(i=0;i<neq;i++)
      gradP[m][i] = 0.0;

    for(e=1;e<=nel;e++) {

	if(0.0==P[m][e])
	    continue;

	for(a=1;a<=ends;a++)       {
	     p = (a-1)*dims;
	     b = E->IEN[lev][m][e].node[a];
	     j1= E->ID[lev][m][b].doff[1];
	     j2= E->ID[lev][m][b].doff[2];
	     j3= E->ID[lev][m][b].doff[3];
		        /*for(b=0;b<ploc_mat_size[E->mesh.nsd];b++)  */
             gradP[m][j1] += E->elt_del[lev][m][e].g[p  ][0] * P[m][e];
             gradP[m][j2] += E->elt_del[lev][m][e].g[p+1][0] * P[m][e];
             gradP[m][j3] += E->elt_del[lev][m][e].g[p+2][0] * P[m][e];
	     }
        }       /* end for el */
     }       /* end for m */

  (E->solver.exchange_id_d)(E, gradP,  lev); /*  correct gradP   */


  strip_bcs_from_residual(E,gradP,lev);

return;
}


double assemble_dAhatp_entry(E,e,level,m)
     struct All_variables *E;
     int e,level,m;

{
    int i,j,p,a,b,node,npno;
    void strip_bcs_from_residual();

    double gradP[81],divU;

    const int ends=enodes[E->mesh.nsd];
    const int dims=E->mesh.nsd;

    npno=E->lmesh.NPNO[level];

    for(i=0;i<81;i++)
	gradP[i] = 0.0;

    divU=0.0;

    for(a=1;a<=ends;a++) {
      p = (a-1)*dims;
      node = E->IEN[level][m][e].node[a];
      j=E->ID[level][m][node].doff[1];
      gradP[p] += E->BI[level][m][j]*E->elt_del[level][m][e].g[p][0];

      j=E->ID[level][m][node].doff[2];
      gradP[p+1] += E->BI[level][m][j]*E->elt_del[level][m][e].g[p+1][0];

      j=E->ID[level][m][node].doff[3];
      gradP[p+2] += E->BI[level][m][j]*E->elt_del[level][m][e].g[p+2][0];
      }


    /* calculate div U from the same thing .... */

    /* only need to run over nodes with non-zero grad P, i.e. the ones in
       the element accessed above, BUT it is only necessary to update the
       value in the original element, because the diagonal is all we use at
       the end ... */

    for(b=1;b<=ends;b++) {
      p = (b-1)*dims;
      divU +=E->elt_del[level][m][e].g[p][0] * gradP[p];
      divU +=E->elt_del[level][m][e].g[p+1][0] * gradP[p+1];
      divU +=E->elt_del[level][m][e].g[p+2][0] * gradP[p+2];
      }

return(divU);  }


/*==============================================================
  Function to supply the element c matrix for a given element e.
  ==============================================================  */

void get_elt_c(struct All_variables *E, int el,
               higher_precision elt_c[24][1], int lev, int m)
{

    int p, a, i, j, nz;
    double temp, beta, rho_avg, x[4];

    double rho[9];

    const int dims = E->mesh.nsd;
    const int ends = enodes[dims];

    if ((el-1)%E->lmesh.ELZ[lev]==0)
        construct_c3x3matrix_el(E,el,&E->element_Cc,&E->element_Ccx,lev,m,1);

    temp = p_point[1].weight[dims-1] * E->GDA[lev][m][el].ppt[1];

    switch (E->refstate.choice) {
    case 1:
        /* the reference state is computed by rho=exp((1-r)Di/gamma) */
        /* so d(rho)/dr/rho == -Di/gamma */

        beta = - E->control.disptn_number * E->control.inv_gruneisen;

        for(a=1;a<=ends;a++) {
            for (i=1;i<=dims;i++) {
                x[i] = E->N.ppt[GNPINDEX(a,1)]
                    * E->element_Cc.ppt[BPINDEX(3,i,a,1)];
            }
            p=dims*(a-1);
            elt_c[p  ][0] = -x[1] * temp * beta;
            elt_c[p+1][0] = -x[2] * temp * beta;
            elt_c[p+2][0] = -x[3] * temp * beta;
        }
        break;
    default:
        /* compute d(rho)/dr/rho from rho(r) */

        for(a=1;a<=ends;a++) {
            j = E->IEN[lev][m][el].node[a];
            nz = (j - 1) % E->lmesh.noz + 1;
            rho[a] = E->refstate.rho[nz];
        }

        rho_avg = 0;
        for(a=1;a<=ends;a++) {
            rho_avg += rho[a];
        }
        rho_avg /= ends;

        for(a=1;a<=ends;a++) {
            for (i=1;i<=dims;i++) {
                x[i] = rho[a] * E->GNX[lev][m][el].ppt[GNPXINDEX(2,a,1)]
                    * E->N.ppt[GNPINDEX(a,1)]
                    * E->element_Cc.ppt[BPINDEX(3,i,a,1)];
            }
            p=dims*(a-1);
            elt_c[p  ][0] = -x[1] * temp / rho_avg;
            elt_c[p+1][0] = -x[2] * temp / rho_avg;
            elt_c[p+2][0] = -x[3] * temp / rho_avg;
        }

    }

    return;
}


/*==============================================================
  Function to supply the element g matrix for a given element e.
  ==============================================================  */

void get_elt_g(E,el,elt_del,lev,m)
     struct All_variables *E;
     int el,m;
     higher_precision elt_del[24][1];
     int lev;

{
   void get_rtf_at_ppts();
   int p,a,i;
   double ra,ct,si,x[4],rtf[4][9];
   double temp;

   const int dims=E->mesh.nsd;
   const int ends=enodes[dims];

   /* Special case, 4/8 node bilinear cartesian square/cube element -> 1 pressure point */

   if ((el-1)%E->lmesh.ELZ[lev]==0)
      construct_c3x3matrix_el(E,el,&E->element_Cc,&E->element_Ccx,lev,m,1);

   get_rtf_at_ppts(E, m, lev, el, rtf);

   temp=p_point[1].weight[dims-1] * E->GDA[lev][m][el].ppt[1];

   ra = rtf[3][1];
   si = 1.0/sin(rtf[1][1]);
   ct = cos(rtf[1][1])*si;

   for(a=1;a<=ends;a++)      {
     for (i=1;i<=dims;i++)
       x[i]=E->GNX[lev][m][el].ppt[GNPXINDEX(2,a,1)]*E->element_Cc.ppt[BPINDEX(3,i,a,1)]
        + 2.0*ra*E->N.ppt[GNPINDEX(a,1)]*E->element_Cc.ppt[BPINDEX(3,i,a,1)]
        + ra*(E->GNX[lev][m][el].ppt[GNPXINDEX(0,a,1)]*E->element_Cc.ppt[BPINDEX(1,i,a,1)]
        +E->N.ppt[GNPINDEX(a,1)]*E->element_Ccx.ppt[BPXINDEX(1,i,1,a,1)]
        +ct*E->N.ppt[GNPINDEX(a,1)]*E->element_Cc.ppt[BPINDEX(1,i,a,1)]
        +si*(E->GNX[lev][m][el].ppt[GNPXINDEX(1,a,1)]*E->element_Cc.ppt[BPINDEX(2,i,a,1)]
        +E->N.ppt[GNPINDEX(a,1)]*E->element_Ccx.ppt[BPXINDEX(2,i,2,a,1)]));

     p=dims*(a-1);
     elt_del[p  ][0] = -x[1] * temp;
     elt_del[p+1][0] = -x[2] * temp;
     elt_del[p+2][0] = -x[3] * temp;

      /* fprintf (E->fp,"B= %d %d %g %g %g %g %g\n",el,a,E->GDA[lev][m][el].ppt[1],E->GNX[lev][m][el].ppt[GNPXINDEX(0,a,1)],E->GNX[lev][m][el].ppt[GNPXINDEX(1,a,1)],elt_del[p][0],elt_del[p+1][0]);
      */
     }

   return;
 }


// for both compressible and incompressible, contributions to the forcing vector from pre-stress 
// term and previous step displacement U 
void get_elt_ps(struct All_variables *E, int el,double elt_f[24], int bcs,int m)
//     struct All_variables *E;
//     int el,m;
//     double elt_f[24];
//     int bcs;

{

   int aid,i,p,a,b,d,j,k,q,es;
    int node[5],back_front,got_elt_k,nodea,nodeb;
    unsigned int type;

    double force[9],force_at_gs[9],stress[9],elt_k[24*24];
    double vector[4],magnitude;
    double sinaa[9],ct[9],tmp,rtf[4][9],ss[9][9];
    double VV[4][9],Ur[8],div_w[4][9][9],skk[9],stress_scaling,maxwelltime,esmu,erho,egrav,con;

    void construct_c3x3matrix_el();
    void get_rtf_at_vpts();
    void parallel_process_termination();

    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
    const int n=loc_mat_size[dims];
    const int ends=enodes[dims];
    const int vpts=vpoints[dims];
    const int lev=E->mesh.levmax;
    const int sphere_key=1;
    double *gnxx,*gnda;    

   stress_scaling = E->data.density*E->data.grav_acc*E->sphere.dradius/E->ve_data_cont.shear_mod;

   if (E->ve_data_cont.SELFG)
      con = E->ve_data_cont.Rsg;  // comment by tao: seems not used
   else 
      con = 0;

    gnxx = E->gNX[m][el].vpt;
    gnda = E->gDA[m][el].vpt;
    es = (el-1)/E->lmesh.elz + 1;

    get_rtf_at_vpts(E, m, lev, el, rtf);

    for (j=1;j<=vpts;j++) {
      sinaa[j] = sin(rtf[1][j]);
      ct[j] = cos(rtf[1][j])/sinaa[j];
      }

	maxwelltime = E->Maxwelltime[m][el]; // comment by tao: before I accidentally defined it in compressible case, making it =0 for incompressible case.

   if (E->ve_data_cont.compressible)   {
      if ((el-1)%E->lmesh.elz==0)
         construct_c3x3matrix_el(E,el,&E->element_Cc,&E->element_Ccx,E->mesh.levmax,m,0);

    /* compute the velocity field at element nodes, and then U_r at the guassian points, also compute B_1k+B_2k+B_3k at node i and gp j */
      for (i=1;i<=ends;i++) {
        VV[1][i] = E->sphere.cap[m].total_V[1][E->ien[m][el].node[i]];
        VV[2][i] = E->sphere.cap[m].total_V[2][E->ien[m][el].node[i]];
        VV[3][i] = E->sphere.cap[m].total_V[3][E->ien[m][el].node[i]];
      }

      esmu = E->esmu[E->mesh.levmax][m][el];
      erho = E->erho[E->mesh.levmax][m][el];
      egrav= E->egrav[E->mesh.levmax][m][el];

      for (j=1;j<=vpts;j++) {
       //sinaa[j] = sin(rtf[1][j]);
       //ct[j] = cos(rtf[1][j])/sinaa[j];
        Ur[j] = 0.0;
        skk[j] = E->Sxx[m][(el-1)*vpts+j] + E->Syy[m][(el-1)*vpts+j] + E->Szz[m][(el-1)*vpts+j];
        for (i=1;i<=ends;i++) {
          for (k=1;k<=dims;k++){
            Ur[j] += VV[k][i]*E->N.vpt[GNVINDEX(i,j)]*E->element_Cc.vpt[BVINDEX(3,k,i,j)];

            div_w[k][i][j] = rtf[3][j]*
                 (gnxx[GNVXINDEX(0,i,j)]*E->element_Cc.vpt[BVINDEX(1,k,i,j)]
                 +E->N.vpt[GNVINDEX(i,j)]*E->element_Ccx.vpt[BVXINDEX(1,k,1,i,j)]
                 +E->N.vpt[GNVINDEX(i,j)]*E->element_Cc.vpt[BVINDEX(3,k,i,j)])
                          +rtf[3][j]*
                (E->N.vpt[GNVINDEX(i,j)]*E->element_Cc.vpt[BVINDEX(1,k,i,j)]*ct[j]
                +E->N.vpt[GNVINDEX(i,j)]*E->element_Cc.vpt[BVINDEX(3,k,i,j)]
                +(gnxx[GNVXINDEX(1,i,j)]*E->element_Cc.vpt[BVINDEX(2,k,i,j)]
                +E->N.vpt[GNVINDEX(i,j)]*E->element_Ccx.vpt[BVXINDEX(2,k,2,i,j)])/sinaa[j])
                          +gnxx[GNVXINDEX(2,i,j)]*E->element_Cc.vpt[BVINDEX(3,k,i,j)];
          }
        }
      }
    }

    for(p=0;p<n;p++) elt_f[p] = 0.0;

    for(j=1;j<=vpts;j++)   {
       ss[1][j] = E->Sxx[m][(el-1)*vpts+j]*maxwelltime;
       ss[2][j] = E->Syy[m][(el-1)*vpts+j]*maxwelltime;
       ss[3][j] = E->Szz[m][(el-1)*vpts+j]*maxwelltime;
       ss[4][j] = E->Sxy[m][(el-1)*vpts+j]*maxwelltime;
       ss[5][j] = E->Sxz[m][(el-1)*vpts+j]*maxwelltime;
       ss[6][j] = E->Szy[m][(el-1)*vpts+j]*maxwelltime;
       if (E->ve_data_cont.compressible)   {
         ss[1][j] += skk[j]*esmu;
         ss[2][j] += skk[j]*esmu;
         ss[3][j] += skk[j]*esmu;
         }
     }

   for(a=1;a<=ends;a++)  {
      p= dims*(a-1);
      for(j=1;j<=vpts;j++)     {
           elt_f[p  ] -= (  ss[1][j]*gnxx[GNVXINDEX(0,a,j)]
                          + ss[2][j]*ct[j]*E->N.vpt[GNVINDEX(a,j)]
                          + ss[4][j]/sinaa[j]*gnxx[GNVXINDEX(1,a,j)]
                          + ss[5][j]*(gnxx[GNVXINDEX(2,a,j)]/rtf[3][j]
                                            - E->N.vpt[GNVINDEX(a,j)]) )
                         * rtf[3][j]*gnda[j]*g_point[j].weight[dims-1] ;
           elt_f[p+1] -= ( ss[2][j]*gnxx[GNVXINDEX(1,a,j)]/sinaa[j]+
                           ss[4][j]*(gnxx[GNVXINDEX(0,a,j)]-ct[j]
                                     *E->N.vpt[GNVINDEX(a,j)])+
                           ss[6][j]*(gnxx[GNVXINDEX(2,a,j)]
                                     /rtf[3][j]-E->N.vpt[GNVINDEX(a,j)])) 
                         * rtf[3][j]*gnda[j]*g_point[j].weight[dims-1];
           elt_f[p+2] -= ((ss[1][j]+ss[2][j])*E->N.vpt[GNVINDEX(a,j)]+
                           ss[3][j]*gnxx[GNVXINDEX(2,a,j)]/rtf[3][j]+
                           ss[5][j]*gnxx[GNVXINDEX(0,a,j)]+
                           ss[6][j]/sinaa[j]*gnxx[GNVXINDEX(1,a,j)])
                         * rtf[3][j]*gnda[j]*g_point[j].weight[dims-1];
                         
        if (E->ve_data_cont.compressible)   {
       /* The 'pre' terms from the additional terms due to compressibility in stiffness matrix */
       /*  w_r*rho*g*div_v + div_w*rho*g*v_r  */                        
           elt_f[p  ] += stress_scaling*erho*egrav*gnda[j]*g_point[j].weight[dims-1]*
             (E->Ekk[m][(el-1)*vpts+j]*E->N.vpt[GNVINDEX(a,j)]*E->element_Cc.vpt[BVINDEX(3,1,a,j)]
             +div_w[1][a][j]*Ur[j]);
           elt_f[p+1] += stress_scaling*erho*egrav*gnda[j]*g_point[j].weight[dims-1]*
             (E->Ekk[m][(el-1)*vpts+j]*E->N.vpt[GNVINDEX(a,j)]*E->element_Cc.vpt[BVINDEX(3,2,a,j)]
             +div_w[2][a][j]*Ur[j]);
           elt_f[p+2] += stress_scaling*erho*egrav*gnda[j]*g_point[j].weight[dims-1]*
             (E->Ekk[m][(el-1)*vpts+j]*E->N.vpt[GNVINDEX(a,j)]*E->element_Cc.vpt[BVINDEX(3,3,a,j)]
             +div_w[3][a][j]*Ur[j]);
           }
         }
 
       }
   return;
 }


// for incompressible VE's pre-stress contribution to the forcing vector
//void get_elt_ps_old(E,el,elt_f,bcs,m)
//     struct All_variables *E;
//     int el,m;
//     double elt_f[24];
//     int bcs;
//
//{
//
//   int aid,i,p,a,b,d,j,k,q,es;
//    int node[5],back_front,got_elt_k,nodea,nodeb;
//    unsigned int type;
//
//    double force[9],force_at_gs[9],stress[9],elt_k[24*24];
//    double vector[4],magnitude;
//    double sinaa[9],ct[9],tmp,rtf[4][9],ss[9][9];
//
//    void construct_c3x3matrix_el();
//    void get_rtf_at_vpts();
//
//    static struct CC Cc;
//    static struct CCX Ccx;
//
//    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
//    const int n=loc_mat_size[dims];
//    const int ends=enodes[dims];
//    const int vpts=vpoints[dims];
//    const int lev=E->mesh.levmax;
//    const int sphere_key=1;
//    double *gnxx,*gnda;    
//
//    gnxx = E->gNX[m][el].vpt;
//    gnda = E->gDA[m][el].vpt;
//    es = (el-1)/E->lmesh.elz + 1;
//
//    get_rtf_at_vpts(E, m, lev, el, rtf);
//
//    for(p=0;p<n;p++) elt_f[p] = 0.0;
//
//        for(j=1;j<=vpts;j++)   {
//            ss[1][j] = E->Sxx[m][(el-1)*vpts+j]*E->Maxwelltime[m][el];
//            ss[2][j] = E->Syy[m][(el-1)*vpts+j]*E->Maxwelltime[m][el];
//            ss[3][j] = E->Szz[m][(el-1)*vpts+j]*E->Maxwelltime[m][el];
//            ss[4][j] = E->Sxy[m][(el-1)*vpts+j]*E->Maxwelltime[m][el];
//            ss[5][j] = E->Sxz[m][(el-1)*vpts+j]*E->Maxwelltime[m][el];
//            ss[6][j] = E->Szy[m][(el-1)*vpts+j]*E->Maxwelltime[m][el];
//            sinaa[j] = sin(rtf[1][j]);
//            ct[j] = cos(rtf[1][j])/sinaa[j];
//        }
//
//
//        for(a=1;a<=ends;a++)  {
//            p= dims*(a-1);
//            for(j=1;j<=vpts;j++)     {
//                elt_f[p  ] -= (  ss[1][j]*gnxx[GNVXINDEX(0,a,j)]
//                               + ss[2][j]*ct[j]*E->N.vpt[GNVINDEX(a,j)]
//                               + ss[4][j]/sinaa[j]*gnxx[GNVXINDEX(1,a,j)]
//                               + ss[5][j]*(gnxx[GNVXINDEX(2,a,j)]/rtf[3][j]
//                                            - E->N.vpt[GNVINDEX(a,j)]) )
//                              * rtf[3][j]*gnda[j]*g_point[j].weight[dims-1] ;
//                elt_f[p+1] -= ( ss[2][j]*gnxx[GNVXINDEX(1,a,j)]/sinaa[j]+
//                                ss[4][j]*(gnxx[GNVXINDEX(0,a,j)]-ct[j]
//                                          *E->N.vpt[GNVINDEX(a,j)])+
//                                ss[6][j]*(gnxx[GNVXINDEX(2,a,j)]
//                                          /rtf[3][j]-E->N.vpt[GNVINDEX(a,j)])) 
//                              * rtf[3][j]*gnda[j]*g_point[j].weight[dims-1];
//               elt_f[p+2] -= ((ss[1][j]+ss[2][j])*E->N.vpt[GNVINDEX(a,j)]+
//                                ss[3][j]*gnxx[GNVXINDEX(2,a,j)]/rtf[3][j]+
//                                ss[5][j]*gnxx[GNVXINDEX(0,a,j)]+
//                                ss[6][j]/sinaa[j]*gnxx[GNVXINDEX(1,a,j)])
//                              * rtf[3][j]*gnda[j]*g_point[j].weight[dims-1];
//               }
// 
//             }
//
//   return;
// }


/*=================================================================
  Function to create the element force vector (allowing for velocity b.c.'s)
  ================================================================= */

void get_elt_f(E,el,elt_f,bcs,m)
     struct All_variables *E;
     int el,m;
     double elt_f[24];
     int bcs;

{

  int i,p,a,b,j,k,q,es;
  int got_elt_k,nodea,nodeb;
  unsigned int type;
  const unsigned int vbc_flag[] = {0, VBX, VBY, VBZ};

  double force[9],force_at_gs[9],elt_k[24*24];


  const int dims=E->mesh.nsd;
  const int n=loc_mat_size[dims];
  const int ends=enodes[dims];
  const int vpts=vpoints[dims];

  es = (el-1)/E->lmesh.elz + 1;

  if ((el-1)%E->lmesh.elz==0)
      construct_c3x3matrix_el(E,el,&E->element_Cc,&E->element_Ccx,E->mesh.levmax,m,0);

  for(p=0;p<n;p++) elt_f[p] = 0.0;

  for(p=1;p<=ends;p++)
    force[p] = E->buoyancy[m][E->ien[m][el].node[p]];

  for(j=1;j<=vpts;j++)       {   /*compute force at each int point */
    force_at_gs[j] = 0.0;
    for(k=1;k<=ends;k++)
      force_at_gs[j] += force[k] * E->N.vpt[GNVINDEX(k,j)] ;
    }

  for(i=1;i<=dims;i++)  {
    for(a=1;a<=ends;a++)  {
      nodea=E->ien[m][el].node[a];
      p= dims*(a-1)+i-1;

      for(j=1;j<=vpts;j++)     /*compute sum(Na(j)*F(j)*det(j)) */
        elt_f[p] += force_at_gs[j] * E->N.vpt[GNVINDEX(a,j)]
           *E->gDA[m][el].vpt[j]*g_point[j].weight[dims-1]
           *E->element_Cc.vpt[BVINDEX(3,i,a,j)];

	  /* imposed velocity terms */

      if(bcs)  {
        got_elt_k = 0;
        for(j=1;j<=dims;j++) {
	  type=vbc_flag[j];
          for(b=1;b<=ends;b++) {
            nodeb=E->ien[m][el].node[b];
            if ((E->node[m][nodeb]&type)&&(E->sphere.cap[m].VB[j][nodeb]!=0.0)){
              if(!got_elt_k) {
                get_elt_k(E,el,elt_k,E->mesh.levmax,m,1);
                got_elt_k = 1;
                }
              q = dims*(b-1)+j-1;
              if(p!=q) {
                elt_f[p] -= elt_k[p*n+q] * E->sphere.cap[m].VB[j][nodeb];
                }
              }
            }  /* end for b */
          }      /* end for j */
        }      /* end if for if bcs */

      }
    } /*  Complete the loops for a,i  	*/



  return;
}


/*=================================================================
  Function to create the element force vector due to stress b.c.
  ================================================================= */

void get_elt_tr(struct All_variables *E, int bel, int side, double elt_tr[24], int m)
{

	const int dims=E->mesh.nsd;
	const int ends1=enodes[dims-1];
	const int oned = onedvpoints[dims];

	struct CC Cc;
	struct CCX Ccx;

	const unsigned sbc_flag[4] = {0,SBX,SBY,SBZ};

	double traction[4][5],traction_at_gs[4][5], value, tmp;
	int j, b, p, k, a, nodea, d;
	int el = E->boundary.element[m][bel];
	int flagged;
	int found = 0;

	const float rho = E->data.density;
	const float g = E->data.grav_acc;
	const float R = 6371000.0;
	const float eta = E->data.ref_viscosity;
	const float kappa = E->data.therm_diff;
	const float factor = 1.0e+00;
	int nodeas;

	if(E->control.side_sbcs)
		for(a=1;a<=ends1;a++)  {
			nodea = E->ien[m][el].node[ sidenodes[side][a] ];
			for(d=1;d<=dims;d++) {
				value = E->sbc.SB[m][side][d][ E->sbc.node[m][nodea] ];
				flagged = (E->node[m][nodea] & sbc_flag[d]) && (value);
				found |= flagged;
				traction[d][a] = ( flagged ? value : 0.0 );
			}
		}
	else {
		/* if side_sbcs is false, only apply sbc on top and bottom surfaces */
		if(side == SIDE_BOTTOM || side == SIDE_TOP) {
			for(a=1;a<=ends1;a++)  {
				nodea = E->ien[m][el].node[ sidenodes[side][a] ];
				for(d=1;d<=dims;d++) {
					value = E->sphere.cap[m].VB[d][nodea];
					flagged = (E->node[m][nodea] & sbc_flag[d]) && (value);
					found |= flagged;
					traction[d][a] = ( flagged ? value : 0.0 );
				}
			}
		}
	}

	/* skip the following computation if no sbc_flag is set
	   or value of sbcs are zero */
	if(!found) return;

	/* compute traction at each int point */
	construct_side_c3x3matrix_el(E,el,&Cc,&Ccx,
				     E->mesh.levmax,m,0,side);

	for(k=1;k<=oned;k++)
		for(d=1;d<=dims;d++) {
			traction_at_gs[d][k] = 0.0;
			for(j=1;j<=ends1;j++)
				traction_at_gs[d][k] += traction[d][j] * E->M.vpt[GMVINDEX(j,k)] ;
		}

	for(j=1;j<=ends1;j++) {
		a = sidenodes[side][j];
		for(d=1;d<=dims;d++) {
			p = dims*(a-1)+d-1;
			for(k=1;k<=oned;k++) {
				tmp = 0.0;
				for(b=1;b<=dims;b++)
					tmp += traction_at_gs[b][k] * Cc.vpt[BVINDEX(b,d,a,k)];

				elt_tr[p] += tmp * E->M.vpt[GMVINDEX(j,k)]
					* E->boundary.det[m][side][k][bel] * g_1d[k].weight[dims-1];

			}
		}
	}
}

void get_elt_tr_pseudo_surf(struct All_variables *E, int bel, int side, double elt_tr[24], int m)
{

	const int dims=E->mesh.nsd;
	const int ends1=enodes[dims-1];
	const int oned = onedvpoints[dims];

	struct CC Cc;
	struct CCX Ccx;

	const unsigned sbc_flag[4] = {0,SBX,SBY,SBZ};

	double traction[4][5],traction_at_gs[4][5], value, tmp;
	int j, b, p, k, a, nodea, d;
	int el = E->boundary.element[m][bel];
	int flagged;
	int found = 0;

	const float rho = E->data.density;
	const float g = E->data.grav_acc;
	const float R = 6371000.0;
	const float eta = E->data.ref_viscosity;
	const float kappa = E->data.therm_diff;
	const float factor = 1.0e+00;
	int nodeas;

	if(E->control.side_sbcs)
		for(a=1;a<=ends1;a++)  {
			nodea = E->ien[m][el].node[ sidenodes[side][a] ];
			for(d=1;d<=dims;d++) {
				value = E->sbc.SB[m][side][d][ E->sbc.node[m][nodea] ];
				flagged = (E->node[m][nodea] & sbc_flag[d]) && (value);
				found |= flagged;
				traction[d][a] = ( flagged ? value : 0.0 );
			}
		}
	else {
		if( side == SIDE_TOP && E->parallel.me_loc[3]==E->parallel.nprocz-1 && (el%E->lmesh.elz==0)) {
			for(a=1;a<=ends1;a++)  {
				nodea = E->ien[m][el].node[ sidenodes[side][a] ];
				nodeas = E->ien[m][el].node[ sidenodes[side][a] ]/E->lmesh.noz;
				traction[1][a] = 0.0;
				traction[2][a] = 0.0;
				traction[3][a] = -1.0*factor*rho*g*(R*R*R)/(eta*kappa)
					*(E->slice.freesurf[m][nodeas]+E->sphere.cap[m].V[3][nodea]*E->advection.timestep);
				if(E->parallel.me==11 && nodea==3328)
					fprintf(stderr,"traction=%e vnew=%e timestep=%e coeff=%e\n",traction[3][a],E->sphere.cap[m].V[3][nodea],E->advection.timestep,-1.0*factor*rho*g*(R*R*R)/(eta*kappa));
				found = 1;
#if 0
				if(found && E->parallel.me==1)
					fprintf(stderr,"me=%d bel=%d el=%d side=%d TOP=%d a=%d sidenodes=%d ien=%d noz=%d nodea=%d traction=%e %e %e\n",
						E->parallel.me,bel,el,side,SIDE_TOP,a,sidenodes[side][a],
						E->ien[m][el].node[ sidenodes[side][a] ],E->lmesh.noz,
						nodea,traction[1][a],traction[2][a],traction[3][a]);

#endif
			}
		}
		else {
			for(a=1;a<=ends1;a++)  {
				nodea = E->ien[m][el].node[ sidenodes[side][a] ];
				for(d=1;d<=dims;d++) {
					value = E->sphere.cap[m].VB[d][nodea];
					flagged = (E->node[m][nodea] & sbc_flag[d]) && (value);
					found |= flagged;
					traction[d][a] = ( flagged ? value : 0.0 );
				}
			}
		}
	}

	/* skip the following computation if no sbc_flag is set
	   or value of sbcs are zero */
	if(!found) return;

	/* compute traction at each int point */
	construct_side_c3x3matrix_el(E,el,&Cc,&Ccx,
				     E->mesh.levmax,m,0,side);

	for(k=1;k<=oned;k++)
		for(d=1;d<=dims;d++) {
			traction_at_gs[d][k] = 0.0;
			for(j=1;j<=ends1;j++)
				traction_at_gs[d][k] += traction[d][j] * E->M.vpt[GMVINDEX(j,k)] ;
		}

	for(j=1;j<=ends1;j++) {
		a = sidenodes[side][j];
		for(d=1;d<=dims;d++) {
			p = dims*(a-1)+d-1;
			for(k=1;k<=oned;k++) {
				tmp = 0.0;
				for(b=1;b<=dims;b++)
					tmp += traction_at_gs[b][k] * Cc.vpt[BVINDEX(b,d,a,k)];

				elt_tr[p] += tmp * E->M.vpt[GMVINDEX(j,k)]
					* E->boundary.det[m][side][k][bel] * g_1d[k].weight[dims-1];

			}
		}
	}
}


/* =================================================================
 subroutine to get augmented lagrange part of stiffness matrix
================================================================== */

void get_aug_k(E,el,elt_k,level,m)
     struct All_variables *E;
     int el,m;
     double elt_k[24*24];
     int level;
{
     int i,p[9],a,b,nodea,nodeb;
     double Visc;

     const int n=loc_mat_size[E->mesh.nsd];
     const int ends=enodes[E->mesh.nsd];
     const int vpts=vpoints[E->mesh.nsd];
     const int dims=E->mesh.nsd;

     Visc = 0.0;
     for(a=1;a<=vpts;a++) {
	  p[a] = (a-1)*dims;
	  Visc += E->EVI[level][m][(el-1)*vpts+a];
       }
     Visc = Visc/vpts;

     for(a=1;a<=ends;a++) {
        nodea=E->IEN[level][m][el].node[a];
        for(b=1;b<=ends;b++) {
           nodeb=E->IEN[level][m][el].node[b];      /* for Kab dims*dims  */
	   i = (a-1)*n*dims+(b-1)*dims;
	   elt_k[i  ] += Visc*E->control.augmented*
	              E->elt_del[level][m][el].g[p[a]][0]*
		      E->elt_del[level][m][el].g[p[b]][0];   /*for 11 */
	   elt_k[i+1] += Visc*E->control.augmented*
	              E->elt_del[level][m][el].g[p[a]][0]*
		      E->elt_del[level][m][el].g[p[b]+1][0];  /* for 12 */
	   elt_k[i+n] += Visc*E->control.augmented*
	              E->elt_del[level][m][el].g[p[a]+1][0]*
		      E->elt_del[level][m][el].g[p[b]][0];    /* for 21 */
	   elt_k[i+n+1] += Visc*E->control.augmented*
	              E->elt_del[level][m][el].g[p[a]+1][0]*
		      E->elt_del[level][m][el].g[p[b]+1][0];  /* for 22 */

           if(3==dims) {
	       elt_k[i+2] += Visc*E->control.augmented*
	              E->elt_del[level][m][el].g[p[a]][0]*
		      E->elt_del[level][m][el].g[p[b]+2][0];  /* for 13 */
	       elt_k[i+n+2] += Visc*E->control.augmented*
	              E->elt_del[level][m][el].g[p[a]+1][0]*
		      E->elt_del[level][m][el].g[p[b]+2][0];  /* for 23 */
	       elt_k[i+n+n] += Visc*E->control.augmented*
	              E->elt_del[level][m][el].g[p[a]+2][0]*
		      E->elt_del[level][m][el].g[p[b]][0];    /* for 31 */
	       elt_k[i+n+n+1] += Visc*E->control.augmented*
	              E->elt_del[level][m][el].g[p[a]+2][0]*
		      E->elt_del[level][m][el].g[p[b]+1][0];  /* for 32 */
	       elt_k[i+n+n+2] += Visc*E->control.augmented*
	              E->elt_del[level][m][el].g[p[a]+2][0]*
		      E->elt_del[level][m][el].g[p[b]+2][0];  /* for 33 */
               }
           }
       }

   return;
   }


/* version */
/* $Id: Element_calculations.c 13196 2008-10-29 23:17:11Z tan2 $ */

/* End of file  */
