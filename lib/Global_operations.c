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
#include <mpi.h>

#include <math.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"
#ifdef ALLOW_ELLIPTICAL
double theta_g(double , struct All_variables *);
#endif

void calc_cbase_at_tp(float , float , float *);

/* ===============================================
   strips horizontal average from nodal field X.
   Assumes orthogonal mesh, otherwise, horizontals
   aren't & another method is required.
   =============================================== */

void remove_horiz_ave(E,X,H,store_or_not)
     struct All_variables *E;
     double **X, *H;
     int store_or_not;

{
    int m,i,j,k,n,nox,noz,noy;
    void return_horiz_ave();

    const int dims = E->mesh.nsd;

    noy = E->lmesh.noy;
    noz = E->lmesh.noz;
    nox = E->lmesh.nox;

    return_horiz_ave(E,X,H);

  for(m=1;m<=E->sphere.caps_per_proc;m++)
    for(k=1;k<=noy;k++)
      for(j=1;j<=nox;j++)
	for(i=1;i<=noz;i++) {
            n = i+(j-1)*noz+(k-1)*noz*nox;
            X[m][n] -= H[i];
	}

   return;
}


void remove_horiz_ave2(struct All_variables *E, double **X)
{
    double *H;

    H = (double *)malloc( (E->lmesh.noz+1)*sizeof(double));
    remove_horiz_ave(E, X, H, 0);
    free ((void *) H);
}


void return_horiz_ave(E,X,H)
     struct All_variables *E;
     double **X, *H;
{
  const int dims = E->mesh.nsd;
  int m,i,j,k,d,nint,noz,nox,noy,el,elz,elx,ely,j1,j2,i1,i2,k1,k2,nproc;
  int top,lnode[5], sizeofH, noz2,iroot;
  double *Have,*temp,aa[5];
  struct Shape_function1 M;
  struct Shape_function1_dA dGamma;
  void get_global_1d_shape_fn();

  sizeofH = (2*E->lmesh.noz+2)*sizeof(double);

  Have = (double *)malloc(sizeofH);
  temp = (double *)malloc(sizeofH);

  noz = E->lmesh.noz;
  noy = E->lmesh.noy;
  elz = E->lmesh.elz;
  elx = E->lmesh.elx;
  ely = E->lmesh.ely;
  noz2 = 2*noz;

  for (i=1;i<=elz;i++)  {
    temp[i] = temp[i+noz] = 0.0;
    temp[i+1] = temp[i+1+noz] = 0.0;
    top = 0;
    if (i==elz) top = 1;
    for (m=1;m<=E->sphere.caps_per_proc;m++)
      for (k=1;k<=ely;k++)
        for (j=1;j<=elx;j++)     {
          el = i + (j-1)*elz + (k-1)*elx*elz;
          get_global_1d_shape_fn(E,el,&M,&dGamma,top,m);

          lnode[1] = E->ien[m][el].node[1];
          lnode[2] = E->ien[m][el].node[2];
          lnode[3] = E->ien[m][el].node[3];
          lnode[4] = E->ien[m][el].node[4];

          for(nint=1;nint<=onedvpoints[E->mesh.nsd];nint++)   {
            for(d=1;d<=onedvpoints[E->mesh.nsd];d++)
              temp[i] += X[m][lnode[d]] * E->M.vpt[GMVINDEX(d,nint)]
                          * dGamma.vpt[GMVGAMMA(0,nint)];
            temp[i+noz] += dGamma.vpt[GMVGAMMA(0,nint)];
            }

          if (i==elz)  {
            lnode[1] = E->ien[m][el].node[5];
            lnode[2] = E->ien[m][el].node[6];
            lnode[3] = E->ien[m][el].node[7];
            lnode[4] = E->ien[m][el].node[8];

            for(nint=1;nint<=onedvpoints[E->mesh.nsd];nint++)   {
              for(d=1;d<=onedvpoints[E->mesh.nsd];d++)
                temp[i+1] += X[m][lnode[d]] * E->M.vpt[GMVINDEX(d,nint)]
                          * dGamma.vpt[GMVGAMMA(1,nint)];
              temp[i+1+noz] += dGamma.vpt[GMVGAMMA(1,nint)];
              }

            }   /* end of if i==elz    */
          }   /* end of j  and k, and m  */
     }        /* Done for i */

  MPI_Allreduce(temp,Have,noz2+1,MPI_DOUBLE,MPI_SUM,E->parallel.horizontal_comm);

  for (i=1;i<=noz;i++) {
    if(Have[i+noz] != 0.0)
       H[i] = Have[i]/Have[i+noz];
    }
 /* if (E->parallel.me==0)
    for(i=1;i<=noz;i++)
      fprintf(stderr,"area %d %d %g\n",E->parallel.me,i,Have[i+noz]);
*/
  free ((void *) Have);
  free ((void *) temp);

  return;
  }

void return_horiz_ave_f(E,X,H)
     struct All_variables *E;
     float **X, *H;
{
  const int dims = E->mesh.nsd;
  int m,i,j,k,d,nint,noz,nox,noy,el,elz,elx,ely,j1,j2,i1,i2,k1,k2,nproc;
  int top,lnode[5], sizeofH, noz2,iroot;
  float *Have,*temp,aa[5];
  struct Shape_function1 M;
  struct Shape_function1_dA dGamma;
  void get_global_1d_shape_fn();

  sizeofH = (2*E->lmesh.noz+2)*sizeof(float);

  Have = (float *)malloc(sizeofH);
  temp = (float *)malloc(sizeofH);

  noz = E->lmesh.noz;
  noy = E->lmesh.noy;
  elz = E->lmesh.elz;
  elx = E->lmesh.elx;
  ely = E->lmesh.ely;
  noz2 = 2*noz;

  for (i=1;i<=elz;i++)  {
    temp[i] = temp[i+noz] = 0.0;
    temp[i+1] = temp[i+1+noz] = 0.0;
    top = 0;
    if (i==elz) top = 1;
    for (m=1;m<=E->sphere.caps_per_proc;m++)
      for (k=1;k<=ely;k++)
        for (j=1;j<=elx;j++)     {
          el = i + (j-1)*elz + (k-1)*elx*elz;
          get_global_1d_shape_fn(E,el,&M,&dGamma,top,m);

          lnode[1] = E->ien[m][el].node[1];
          lnode[2] = E->ien[m][el].node[2];
          lnode[3] = E->ien[m][el].node[3];
          lnode[4] = E->ien[m][el].node[4];

          for(nint=1;nint<=onedvpoints[E->mesh.nsd];nint++)   {
            for(d=1;d<=onedvpoints[E->mesh.nsd];d++)
              temp[i] += X[m][lnode[d]] * E->M.vpt[GMVINDEX(d,nint)]
                          * dGamma.vpt[GMVGAMMA(0,nint)];
            temp[i+noz] += dGamma.vpt[GMVGAMMA(0,nint)];
            }

          if (i==elz)  {
            lnode[1] = E->ien[m][el].node[5];
            lnode[2] = E->ien[m][el].node[6];
            lnode[3] = E->ien[m][el].node[7];
            lnode[4] = E->ien[m][el].node[8];

            for(nint=1;nint<=onedvpoints[E->mesh.nsd];nint++)   {
              for(d=1;d<=onedvpoints[E->mesh.nsd];d++)
                temp[i+1] += X[m][lnode[d]] * E->M.vpt[GMVINDEX(d,nint)]
                          * dGamma.vpt[GMVGAMMA(1,nint)];
              temp[i+1+noz] += dGamma.vpt[GMVGAMMA(1,nint)];
              }

            }   /* end of if i==elz    */
          }   /* end of j  and k, and m  */
     }        /* Done for i */

  MPI_Allreduce(temp,Have,noz2+1,MPI_FLOAT,MPI_SUM,E->parallel.horizontal_comm);

  for (i=1;i<=noz;i++) {
    if(Have[i+noz] != 0.0)
       H[i] = Have[i]/Have[i+noz];
    }
 /* if (E->parallel.me==0)
    for(i=1;i<=noz;i++)
      fprintf(stderr,"area %d %d %g\n",E->parallel.me,i,Have[i+noz]);
*/
  free ((void *) Have);
  free ((void *) temp);

  return;
  }


/******* RETURN ELEMENTWISE HORIZ AVE ********************************/
/*                                                                   */
/* This function is similar to return_horiz_ave in the citcom code   */
/* however here, elemental horizontal averages are given rather than */
/* nodal averages. Also note, here is average per element            */

void return_elementwise_horiz_ave(E,X,H)
     struct All_variables *E;
     double **X, *H;
{

  int m,i,j,k,d,noz,noy,el,elz,elx,ely,nproc;
  int sizeofH;
  int elz2;
  double *Have,*temp;

  sizeofH = (2*E->lmesh.elz+2)*sizeof(double);

  Have = (double *)malloc(sizeofH);
  temp = (double *)malloc(sizeofH);

  noz = E->lmesh.noz;
  noy = E->lmesh.noy;
  elz = E->lmesh.elz;
  elx = E->lmesh.elx;
  ely = E->lmesh.ely;
  elz2 = 2*elz;

  for (i=0;i<=(elz*2+1);i++)
  {
    temp[i]=0.0;
  }

  for (i=1;i<=elz;i++)
  {
    for (m=1;m<=E->sphere.caps_per_proc;m++)
    {
      for (k=1;k<=ely;k++)
      {
        for (j=1;j<=elx;j++)
        {
          el = i + (j-1)*elz + (k-1)*elx*elz;
          temp[i] += X[m][el]*E->ECO[E->mesh.levmax][m][el].area;
          temp[i+elz] += E->ECO[E->mesh.levmax][m][el].area;
        }
      }
    }
  }



/* determine which processors should get the message from me for
               computing the layer averages */

  MPI_Allreduce(temp,Have,elz2+1,MPI_DOUBLE,MPI_SUM,E->parallel.horizontal_comm);

  for (i=1;i<=elz;i++) {
    if(Have[i+elz] != 0.0)
       H[i] = Have[i]/Have[i+elz];
    }


  free ((void *) Have);
  free ((void *) temp);

  return;
}

float return_bulk_value(E,Z,average)
     struct All_variables *E;
     float **Z;
     int average;

{
    int n,i,j,k,el,m;
    float volume,integral,volume1,integral1;

    const int vpts = vpoints[E->mesh.nsd];
    const int ends = enodes[E->mesh.nsd];

    volume1=0.0;
    integral1=0.0;

    for (m=1;m<=E->sphere.caps_per_proc;m++)
       for (el=1;el<=E->lmesh.nel;el++)  {

	  for(j=1;j<=vpts;j++)
	    for(i=1;i<=ends;i++) {
		n = E->ien[m][el].node[i];
		volume1 += E->N.vpt[GNVINDEX(i,j)] * E->gDA[m][el].vpt[j];
		integral1 += Z[m][n] * E->N.vpt[GNVINDEX(i,j)] * E->gDA[m][el].vpt[j];
                }

          }


    MPI_Allreduce(&volume1  ,&volume  ,1,MPI_FLOAT,MPI_SUM,E->parallel.world);
    MPI_Allreduce(&integral1,&integral,1,MPI_FLOAT,MPI_SUM,E->parallel.world);

    if(average && volume != 0.0)
 	   integral /= volume;

    return((float)integral);
}

/************ RETURN BULK VALUE_D *****************************************/
/*                                                                        */
/* Same as return_bulk_value but allowing double instead of float.        */
/* I think when integer average =1, volume average is returned.           */
/*         when integer average =0, integral is returned.           */


double return_bulk_value_d(E,Z,average)
     struct All_variables *E;
     double **Z;
     int average;

{
    int n,i,j,el,m;
    double volume,integral,volume1,integral1;

    const int vpts = vpoints[E->mesh.nsd];
    const int ends = enodes[E->mesh.nsd];

    volume1=0.0;
    integral1=0.0;

    for (m=1;m<=E->sphere.caps_per_proc;m++)
       for (el=1;el<=E->lmesh.nel;el++)  {

          for(j=1;j<=vpts;j++)
            for(i=1;i<=ends;i++) {
                n = E->ien[m][el].node[i];
                volume1 += E->N.vpt[GNVINDEX(i,j)] * E->gDA[m][el].vpt[j];
                integral1 += Z[m][n] * E->N.vpt[GNVINDEX(i,j)] * E->gDA[m][el].vpt[j];
            }

       }


    MPI_Allreduce(&volume1  ,&volume  ,1,MPI_DOUBLE,MPI_SUM,E->parallel.world);
    MPI_Allreduce(&integral1,&integral,1,MPI_DOUBLE,MPI_SUM,E->parallel.world);

    if(average && volume != 0.0)
           integral /= volume;

    return((double)integral);
}

/* ================================================== */
float find_max_horizontal(E,Tmax)
struct All_variables *E;
float Tmax;
{
 float ttmax;

 MPI_Allreduce(&Tmax,&ttmax,1,MPI_FLOAT,MPI_MAX,E->parallel.horizontal_comm);

 return(ttmax);
 }

/* ================================================== */
void sum_across_surf_sph2(E,sphc,sphs)
struct All_variables *E;
double *sphc,*sphs;
{
 int jumpp,total,j,d;
 double *sphcs,*temp;

 temp = (double *) malloc((E->sphere.hindice*2+3)*sizeof(double));
 sphcs = (double *) malloc((E->sphere.hindice*2+3)*sizeof(double));

 /* pack */
 jumpp = E->sphere.hindice;
 total = E->sphere.hindice*2+3;
 for (j=0;j<E->sphere.hindice;j++)   {
   sphcs[j] = sphc[j];
   sphcs[j+jumpp] = sphs[j];
 }
 sphcs[total-1] = sphs[E->sphere.hindice];

 /* sum across processors in horizontal direction */
 MPI_Allreduce(sphcs,temp,total,MPI_DOUBLE,MPI_SUM,E->parallel.horizontal_comm);

 /* unpack */
 for (j=0;j<E->sphere.hindice;j++)   {
   sphc[j] = temp[j];
   sphs[j] = temp[j+jumpp];
 }
 sphs[E->sphere.hindice]=temp[total-1];

 free((void *)temp);
 free((void *)sphcs);

 return;
}

/* ================================================== */

/* ================================================== */
void sum_across_surf_sph1(E,sphc,sphs)
struct All_variables *E;
float *sphc,*sphs;
{
 int jumpp,total,j,d;
 float *sphcs,*temp;

 temp = (float *) malloc((E->sphere.hindice*2+3)*sizeof(float));
 sphcs = (float *) malloc((E->sphere.hindice*2+3)*sizeof(float));

 /* pack */
 jumpp = E->sphere.hindice;
 total = E->sphere.hindice*2+3;
 for (j=0;j<E->sphere.hindice;j++)   {
   sphcs[j] = sphc[j];
   sphcs[j+jumpp] = sphs[j];
 }
 sphcs[total-1] = sphs[E->sphere.hindice];

 /* sum across processors in horizontal direction */
 MPI_Allreduce(sphcs,temp,total,MPI_FLOAT,MPI_SUM,E->parallel.horizontal_comm);

 /* unpack */
 for (j=0;j<E->sphere.hindice;j++)   {
   sphc[j] = temp[j];
   sphs[j] = temp[j+jumpp];
 }
 sphs[E->sphere.hindice]=temp[total-1];

 free((void *)temp);
 free((void *)sphcs);

 return;
}

/* ================================================== */


float global_fvdot(E,A,B,lev)
   struct All_variables *E;
   float **A,**B;
   int lev;

{
  int m,i,neq;
  float prod, temp,temp1;

  neq=E->lmesh.NEQ[lev];

  temp = 0.0;
  temp1 = 0.0;
  prod = 0.0;
  for (m=1;m<=E->sphere.caps_per_proc;m++)  {
    neq=E->lmesh.NEQ[lev];
    temp1 = 0.0;
    for (i=0;i<neq;i++)
      temp += A[m][i]*B[m][i];

    for (i=1;i<=E->parallel.Skip_neq[lev][m];i++)
       temp1 += A[m][E->parallel.Skip_id[lev][m][i]]*B[m][E->parallel.Skip_id[lev][m][i]];

    temp -= temp1;

    }

  MPI_Allreduce(&temp, &prod,1,MPI_FLOAT,MPI_SUM,E->parallel.world);

  return (prod);
}


double kineticE_radial(E,A,lev)
   struct All_variables *E;
   double **A;
   int lev;

{
  int m,i,neq;
  double prod, temp,temp1;

    temp = 0.0;
    prod = 0.0;

  for (m=1;m<=E->sphere.caps_per_proc;m++)  {
    neq=E->lmesh.NEQ[lev];
    temp1 = 0.0;
    for (i=0;i<neq;i++)
      if ((i+1)%3==0)
        temp += A[m][i]*A[m][i];

    for (i=1;i<=E->parallel.Skip_neq[lev][m];i++)
      if ((E->parallel.Skip_id[lev][m][i]+1)%3==0)
        temp1 += A[m][E->parallel.Skip_id[lev][m][i]]*A[m][E->parallel.Skip_id[lev][m][i]];

    temp -= temp1;

    }

  MPI_Allreduce(&temp, &prod,1,MPI_DOUBLE,MPI_SUM,E->parallel.world);

  return (prod);
}

double global_vdot_e(E,A,B,lev)
   struct All_variables *E;
   double **A,**B;
   int lev;

{
  int m,a,i,node,e,j,neq,nel;
  double prod, temp,temp1;

  const int ends=enodes[E->mesh.nsd];

    temp = 0.0;
    prod = 0.0;

  for (m=1;m<=E->sphere.caps_per_proc;m++)  {
    nel=E->lmesh.nel;
    temp1 = 0.0;
    for (e=0;e<=nel;e++) {
      if (E->viscosity.sdepv_expt[E->mat[m][e]-1]>1.01)  {   // for non-Newt 
        for (a=1;a<=ends;a++) {
          node = E->ien[m][e].node[a];
          for (j=1;j<=E->mesh.nsd;j++) {
            i = E->id[m][node].doff[j];
            temp += A[m][i]*B[m][i];
            }
          }
        }

       }
    }

  MPI_Allreduce(&temp, &prod,1,MPI_DOUBLE,MPI_SUM,E->parallel.world);

  return (prod);
   }

double global_vdot(E,A,B,lev)
   struct All_variables *E;
   double **A,**B;
   int lev;

{
  int m,i,neq;
  double prod, temp,temp1;

    temp = 0.0;
    prod = 0.0;

  for (m=1;m<=E->sphere.caps_per_proc;m++)  {
    neq=E->lmesh.NEQ[lev];
    temp1 = 0.0;
    for (i=0;i<neq;i++)
      temp += A[m][i]*B[m][i];

    for (i=1;i<=E->parallel.Skip_neq[lev][m];i++)
       temp1 += A[m][E->parallel.Skip_id[lev][m][i]]*B[m][E->parallel.Skip_id[lev][m][i]];

    temp -= temp1;

    }

  MPI_Allreduce(&temp, &prod,1,MPI_DOUBLE,MPI_SUM,E->parallel.world);

  return (prod);
}


double global_pdot(E,A,B,lev)
   struct All_variables *E;
   double **A,**B;
   int lev;

{
  int i,m,npno;
  double prod, temp;

  npno=E->lmesh.NPNO[lev];

  temp = 0.0;
  prod = 0.0;
  for (m=1;m<=E->sphere.caps_per_proc;m++)  {
    npno=E->lmesh.NPNO[lev];
    for (i=1;i<=npno;i++)
      temp += A[m][i]*B[m][i];
    }

  MPI_Allreduce(&temp, &prod,1,MPI_DOUBLE,MPI_SUM,E->parallel.world);

  return (prod);
}


/* return ||V||^2 */
double global_v_norm2(struct All_variables *E,  double **V)
{
    int i, m, d;
    int eqn1, eqn2, eqn3;
    double prod, temp;

    temp = 0.0;
    prod = 0.0;
    for (m=1; m<=E->sphere.caps_per_proc; m++)
        for (i=1; i<=E->lmesh.nno; i++) {
            eqn1 = E->id[m][i].doff[1];
            eqn2 = E->id[m][i].doff[2];
            eqn3 = E->id[m][i].doff[3];
            /* L2 norm  */
            temp += (V[m][eqn1] * V[m][eqn1] +
                     V[m][eqn2] * V[m][eqn2] +
                     V[m][eqn3] * V[m][eqn3]) * E->NMass[m][i];
        }

    MPI_Allreduce(&temp, &prod, 1, MPI_DOUBLE, MPI_SUM, E->parallel.world);

    return (prod/E->mesh.volume);
}


/* return ||P||^2 */
double global_p_norm2(struct All_variables *E,  double **P)
{
    int i, m;
    double prod, temp;

    temp = 0.0;
    prod = 0.0;
    for (m=1; m<=E->sphere.caps_per_proc; m++)
        for (i=1; i<=E->lmesh.npno; i++) {
            /* L2 norm */
            temp += P[m][i] * P[m][i] * E->eco[m][i].area;
        }

    MPI_Allreduce(&temp, &prod, 1, MPI_DOUBLE, MPI_SUM, E->parallel.world);

    return (prod/E->mesh.volume);
}


/* return ||A||^2, where A_i is \int{div(u) d\Omega_i} */
double global_div_norm2(struct All_variables *E,  double **A)
{
    int i, m;
    double prod, temp;

    temp = 0.0;
    prod = 0.0;
    for (m=1; m<=E->sphere.caps_per_proc; m++)
        for (i=1; i<=E->lmesh.npno; i++) {
            /* L2 norm of div(u) */
            temp += A[m][i] * A[m][i] / E->eco[m][i].area;

            /* L1 norm */
            /*temp += fabs(A[m][i]);*/
        }

    MPI_Allreduce(&temp, &prod, 1, MPI_DOUBLE, MPI_SUM, E->parallel.world);

    return (prod/E->mesh.volume);
}


double global_tdot_d(E,A,B,lev)
   struct All_variables *E;
   double **A,**B;
   int lev;

{
  int i,nno,m;
  double prod, temp;

  nno=E->lmesh.NNO[lev];
  temp = 0.0;
  prod = 0.0;
  for (m=1;m<=E->sphere.caps_per_proc;m++)  {
    nno=E->lmesh.NNO[lev];
    for (i=1;i<=nno;i++)
    if (/*1 ||*/ !(E->NODE[lev][m][i] & SKIP))
      temp += i/E->lmesh.noz * A[m][i] ; // *B[m][i];
    }
  MPI_Allreduce(&temp, &prod,1,MPI_DOUBLE,MPI_SUM,E->parallel.world);

  return (prod);
  }

double global_tdot_d_noskip(E,A,B,lev)
   struct All_variables *E;
   double **A,**B;
   int lev;

{
  int i,nno,m;
  double prod, temp;

  nno=E->lmesh.NNO[lev];
  temp = 0.0;
  prod = 0.0;
  for (m=1;m<=E->sphere.caps_per_proc;m++)  {
    nno=E->lmesh.NNO[lev];
    for (i=1;i<=nno;i++)
    if (1 || !(E->NODE[lev][m][i] & SKIP))
      temp += i/E->lmesh.noz * A[m][i] ; // *B[m][i];
    }
  MPI_Allreduce(&temp, &prod,1,MPI_DOUBLE,MPI_SUM,E->parallel.world);

  return (prod);
  }

float global_tdot(E,A,B,lev)
   struct All_variables *E;
   float **A,**B;
   int lev;

{
  int i,nno,m;
  float prod, temp;


  temp = 0.0;
  prod = 0.0;
  for (m=1;m<=E->sphere.caps_per_proc;m++)  {
    nno=E->lmesh.NNO[lev];
    for (i=1;i<=nno;i++)
      if (!(E->NODE[lev][m][i] & SKIP))
        temp += A[m][i]*B[m][i];
    }

  MPI_Allreduce(&temp, &prod,1,MPI_FLOAT,MPI_SUM,E->parallel.world);

  return (prod);
  }


float global_fmin(E,a)
   struct All_variables *E;
   float a;
{
  float temp;
  MPI_Allreduce(&a, &temp,1,MPI_FLOAT,MPI_MIN,E->parallel.world);
  return (temp);
  }

double global_dmax(E,a)
   struct All_variables *E;
   double a;
{
  double temp;
  MPI_Allreduce(&a, &temp,1,MPI_DOUBLE,MPI_MAX,E->parallel.world);
  return (temp);
  }


float global_fmax(E,a)
   struct All_variables *E;
   float a;
{
  float temp;
  MPI_Allreduce(&a, &temp,1,MPI_FLOAT,MPI_MAX,E->parallel.world);
  return (temp);
  }

double Tmaxd(E,T)
  struct All_variables *E;
  double **T;
{
  double global_dmax(),temp,temp1;
  int i,m;

  temp = -10.0;
  for (m=1;m<=E->sphere.caps_per_proc;m++)
    for(i=1;i<=E->lmesh.nno;i++)
      temp = max(T[m][i],temp);

  temp1 = global_dmax(E,temp);
  return (temp1);
  }


float Tmaxf(E,T)
  struct All_variables *E;
  float **T;
{
  float global_fmax(),temp,temp1;
  int i,m;

  temp = -10.0;
  for (m=1;m<=E->sphere.caps_per_proc;m++)
    for(i=1;i<=E->lmesh.nno;i++)
      temp = max(T[m][i],temp);

  temp1 = global_fmax(E,temp);
  return (temp1);
  }


double  vnorm_nonnewt(E,dU,U,lev)
  struct All_variables *E;
  double **dU,**U;
  int lev;
{
 double temp1,temp2,dtemp,temp;
 int a,e,i,m,node;
 const int dims = E->mesh.nsd;
 const int ends = enodes[dims];
 const int nel=E->lmesh.nel;

 dtemp=0.0;
 temp=0.0;
for (m=1;m<=E->sphere.caps_per_proc;m++)
  for (e=1;e<=nel;e++)
   /*if (E->mat[m][e]==1)*/
     for (i=1;i<=dims;i++)
       for (a=1;a<=ends;a++) {
	 node = E->IEN[lev][m][e].node[a];
         dtemp += dU[m][ E->ID[lev][m][node].doff[i] ]*
                  dU[m][ E->ID[lev][m][node].doff[i] ];
         temp += U[m][ E->ID[lev][m][node].doff[i] ]*
                 U[m][ E->ID[lev][m][node].doff[i] ];
         }


  MPI_Allreduce(&dtemp, &temp2,1,MPI_DOUBLE,MPI_SUM,E->parallel.world);
  MPI_Allreduce(&temp, &temp1,1,MPI_DOUBLE,MPI_SUM,E->parallel.world);

  temp1 = sqrt(temp2/temp1);

  return (temp1);
}


void sum_across_depth_sph1(E,sphc,sphs)
     struct All_variables *E;
     double  *sphc,*sphs;
{
    int jumpp,total,j;

    double *sphcs,*temp;

    if (E->parallel.nprocz > 1)  {
	total = E->sphere.hindice*2;
	temp = (double *) malloc(total*sizeof(double));
	sphcs = (double *) malloc(total*sizeof(double));

	/* pack sphc[] and sphs[] into sphcs[] */
	jumpp = E->sphere.hindice;
	for (j=0;j<E->sphere.hindice;j++)   {
	    sphcs[j] = sphc[j];
	    sphcs[j+jumpp] = sphs[j];
	}

	/* sum across processors in z direction */
	MPI_Allreduce(sphcs, temp, total, MPI_DOUBLE, MPI_SUM,
		      E->parallel.vertical_comm);

	/* unpack */
	for (j=0;j<E->sphere.hindice;j++)   {
	    sphc[j] = temp[j];
	    sphs[j] = temp[j+jumpp];
	}

	free(temp);
	free(sphcs);
    }


    return;
}


/* ================================================== */
/* ================================================== */
void broadcast_vertical(struct All_variables *E,
                        double *sphc, double *sphs,
                        int root)
{
    int jumpp, total, j;
    double *temp;

    if(E->parallel.nprocz == 1) return;

    jumpp = E->sphere.hindice;
    total = E->sphere.hindice*2;
    temp = (double *) malloc(total*sizeof(double));

    if (E->parallel.me_loc[3] == root) {
        /* pack */
        for (j=0; j<E->sphere.hindice; j++)   {
            temp[j] = sphc[j];
            temp[j+jumpp] = sphs[j];
        }
    }

    MPI_Bcast(temp, total, MPI_DOUBLE, root, E->parallel.vertical_comm);

    if (E->parallel.me_loc[3] != root) {
        /* unpack */
        for (j=0; j<E->sphere.hindice; j++)   {
            sphc[j] = temp[j];
            sphs[j] = temp[j+jumpp];
        }
    }

    free((void *)temp);

    return;
}


/*
 * remove rigid body rotation or angular momentum from the velocity
 */

void remove_rigid_rot(struct All_variables *E)
{
    void velo_from_element_d();
    double myatan();
    double wx, wy, wz, v_theta, v_phi, cos_t,sin_t,sin_f, cos_f,frd;
    double vx[9], vy[9], vz[9];
    double r, t, f, efac,tg;
    float cart_base[9];
    double exyz[4], fxyz[4];

    int m, e, i, k, j, node;
    const int lev = E->mesh.levmax;
    const int nno = E->lmesh.nno;
    const int ends = ENODES3D;
    const int ppts = PPOINTS3D;
    const int vpts = VPOINTS3D;
    const int sphere_key = 1;
    double VV[4][9];
    double rot, fr, tr;
    double tmp, moment_of_inertia, rho;
    double ref_density;


    if(E->control.remove_angular_momentum) {
        moment_of_inertia = tmp = 0;
        for (i=1;i<=E->lmesh.elz;i++)
        {
          if(!E->ve_data_cont.compressible)
              ref_density = 0.5*(E->refstate.rho[i] + E->refstate.rho[i+1]);
          else if(E->ve_data_cont.compressible)
              ref_density = E->erho[lev][1][i];

          tmp += (8.0*M_PI/15.0)*
//                0.5*(E->refstate.rho[i] + E->refstate.rho[i+1])*
                  // E->erho[lev][1][i] *
                  ref_density *
                (pow(E->sx[1][3][i+1],5.0) - pow(E->sx[1][3][i],5.0));

        }

        MPI_Allreduce(&tmp, &moment_of_inertia, 1, MPI_DOUBLE,
                      MPI_SUM, E->parallel.vertical_comm);
    } else {
         /* no need to weight in rho(r) here. */
        moment_of_inertia = (8.0*M_PI/15.0)*
            (pow(E->sphere.ro,5.0) - pow(E->sphere.ri,5.0));
    }

    /* compute and add angular momentum components */
    
    exyz[1] = exyz[2] = exyz[3] = 0.0;
    fxyz[1] = fxyz[2] = fxyz[3] = 0.0;
    
    
    for (m=1;m<=E->sphere.caps_per_proc;m++) {
      for (e=1;e<=E->lmesh.nel;e++) {
#ifdef ALLOW_ELLIPTICAL
	t = theta_g(E->eco[m][e].centre[1],E);
#else
	t = E->eco[m][e].centre[1];
#endif
	f = E->eco[m][e].centre[2];
	r = E->eco[m][e].centre[3];
	
	cos_t = cos(t);sin_t = sin(t);
	sin_f = sin(f);cos_f = cos(f);
	
	/* get Cartesian, element local velocities */
	velo_from_element_d(E,VV,m,e,sphere_key);
	for (j=1;j<=ppts;j++)   {
	  vx[j] = 0.0;vy[j] = 0.0;
	}
	for (j=1;j<=ppts;j++)   {
	  for (i=1;i<=ends;i++)   {
	    vx[j] += VV[1][i]*E->N.ppt[GNPINDEX(i,j)]; 
	    vy[j] += VV[2][i]*E->N.ppt[GNPINDEX(i,j)]; 
	  }
	}

        wx = -r*vy[1];
        wy =  r*vx[1];

        if(E->control.remove_angular_momentum) {
            int nz = (e-1) % E->lmesh.elz + 1;
            if(!E->ve_data_cont.compressible)
              rho = 0.5 * (E->refstate.rho[nz] + E->refstate.rho[nz+1]);
            else if(E->ve_data_cont.compressible)
              rho = E->erho[lev][m][e];
        } else {
            rho = 1;
        }
	exyz[1] += (wx*cos_t*cos_f - wy*sin_f) * E->eco[m][e].area * rho;
	exyz[2] += (wx*cos_t*sin_f + wy*cos_f) * E->eco[m][e].area * rho;
	exyz[3] -= (wx*sin_t                 ) * E->eco[m][e].area * rho;
      }
    } /* end cap */
    
    MPI_Allreduce(exyz,fxyz,4,MPI_DOUBLE,MPI_SUM,E->parallel.world);
    
    fxyz[1] = fxyz[1] / moment_of_inertia;
    fxyz[2] = fxyz[2] / moment_of_inertia;
    fxyz[3] = fxyz[3] / moment_of_inertia;
    
    rot = sqrt(fxyz[1]*fxyz[1] + fxyz[2]*fxyz[2] + fxyz[3]*fxyz[3]);
    fr = myatan(fxyz[2], fxyz[1]);
    tr = acos(fxyz[3] / rot);
    
    if (E->parallel.me==0) {
        if(E->control.remove_angular_momentum) {
            fprintf(E->fp,"Angular momentum: rot=%e tr=%e fr=%e\n",rot,tr*180/M_PI,fr*180/M_PI);
            fprintf(stderr,"Angular momentum: rot=%e tr=%e fr=%e\n",rot,tr*180/M_PI,fr*180/M_PI);
        } else {
            fprintf(E->fp,"Rigid rotation: rot=%e tr=%e fr=%e\n",rot,tr*180/M_PI,fr*180/M_PI);
            fprintf(stderr,"Rigid rotation: rot=%e tr=%e fr=%e\n",rot,tr*180/M_PI,fr*180/M_PI);
        }
    }
    /*
      remove rigid rotation 
    */
#ifdef ALLOW_ELLIPTICAL
    for (m=1;m<=E->sphere.caps_per_proc;m++)  {
      for (node=1;node<=nno;node++)   {
	/* cartesian velocity = omega \cross r  */
	vx[0] = fxyz[2]* E->x[m][3][node] - fxyz[3]*E->x[m][2][node];
	vx[1] = fxyz[3]* E->x[m][1][node] - fxyz[1]*E->x[m][3][node];
	vx[2] = fxyz[1]* E->x[m][2][node] - fxyz[2]*E->x[m][1][node];
	/* project into theta, phi */
	calc_cbase_at_node(m,node,cart_base,E);
	v_theta = vx[0]*cart_base[3] + vx[1]*cart_base[4] + vx[2]*cart_base[5] ;
	v_phi   = vx[0]*cart_base[6] + vx[1]*cart_base[7];
	E->sphere.cap[m].V[1][node] -= v_theta;
	E->sphere.cap[m].V[2][node] -= v_phi;
      }
    }
#else
    sin_t = sin(tr) * rot;
    cos_t = cos(tr) * rot;
    for (m=1;m<=E->sphere.caps_per_proc;m++)  {
      for (node=1;node<=nno;node++)   {
	frd = fr - E->sx[m][2][node];
	v_theta = E->sx[m][3][node] * sin_t * sin(frd);
	v_phi =   E->sx[m][3][node] * 
	  (  E->SinCos[lev][m][0][node] * cos_t - E->SinCos[lev][m][2][node]  * sin_t * cos(frd) );
	
	E->sphere.cap[m].V[1][node] -= v_theta;
	E->sphere.cap[m].V[2][node] -= v_phi;
//        E->U[m][E->id[m][node].doff[1]] = E->sphere.cap[m].V[1][node];
//        E->U[m][E->id[m][node].doff[2]] = E->sphere.cap[m].V[2][node];
             // now U and V are identical
      }
    }
#endif

    if (E->ve_data_cont.compressible){
      // compute Net rotation rate of surface
      double moment_of_inertia_surface = 8 * M_PI * pow(E->sphere.ro,4.0) / 3.0; 
      double phi, theta;
      double total_surface_integral();
      // assume radius of 1.0
      // assert(E->sphere.ro == 1.0);
      // double omega_surface[4];  // rotation rate of surface, see E->ve_data_cont.Omega_surface
      for (m=1;m<=E->sphere.caps_per_proc;m++)  
        for (j=1;j<=E->lmesh.nsf;j++)  { 
            i=j*E->lmesh.noz;
            phi = E->sx[m][2][i];
            theta = E->sx[m][1][i];
            cos_t = cos(theta); sin_t = sin(theta);
            sin_f = sin(phi); cos_f = cos(phi);
            wx = -E->sphere.cap[m].V[2][i] * E->sphere.ro;
            wy =  E->sphere.cap[m].V[1][i] * E->sphere.ro;
            E->Xsurf[1][m][j] = wx*cos_t*cos_f - wy*sin_f; // x-component of (r x V)
            E->Xsurf[2][m][j] = wx*cos_t*sin_f + wy*cos_f; // y-component of (r x V)
            E->Xsurf[3][m][j] = -1 * wx*sin_t;                  // z-component of (r x V)
        }
      
      // compute surface integral
      E->ve_data_cont.Omega_surface[0] = pow(E->sphere.ro,2.0) * total_surface_integral(E,E->Xsurf[1], 1); // 1 for surface
      E->ve_data_cont.Omega_surface[1] = pow(E->sphere.ro,2.0) * total_surface_integral(E,E->Xsurf[2], 1); // 
      E->ve_data_cont.Omega_surface[2] = pow(E->sphere.ro,2.0) * total_surface_integral(E,E->Xsurf[3], 1); // 
      E->ve_data_cont.Omega_surface[0] /= moment_of_inertia_surface;
      E->ve_data_cont.Omega_surface[1] /= moment_of_inertia_surface;
      E->ve_data_cont.Omega_surface[2] /= moment_of_inertia_surface;
      E->ve_data_cont.Omega_surface[0] *= 180.0/M_PI / E->ve_data_cont.tau_in_years;  // convert to degrees per year, seems need (E->ve_data_cont.tau*E->advection.timestep)
      E->ve_data_cont.Omega_surface[1] *= 180.0/M_PI / E->ve_data_cont.tau_in_years;  
      E->ve_data_cont.Omega_surface[2] *= 180.0/M_PI / E->ve_data_cont.tau_in_years;  
      // written to topo_s

    }

    return;

}

/****************************************************************************
 * truncate_Ylm_expansion(E,X3,l_max)
 * Computes the spherical harmonic expansion of surface field X3[m][j] 
 * (where m is in sphere.caps_per_proc and j is a surface node index),
 * then re-sums that expansion up to l_max.
 ***************************************************************************/
void truncate_Ylm_expansion(E,X3,l_max)
    struct All_variables *E;
    double **X3;
    int l_max;
{
    int m,j,i,p;
    int ll,mm;
    int es,el,k,k1;
    double area, temp0[5] ;
    double *sphc,*sphs;
    void sphere_expansion();
    void parallel_process_sync();
    void sum_across_surf_sph2();
    double CPU_time0(),time;
    FILE *fp;
    static int been_here=0;
    int verbose=1;   // write out data from SLE calculations
    //

    if (E->parallel.me_loc[3] != E->parallel.nprocz-1)  return;

    if (E->parallel.me == 0)  time = CPU_time0();
    if (!been_here){
        been_here = 1;
    }

    sphc = E->sphere.sphc[0] ;
    sphs = E->sphere.sphs[0] ;

    // perform the sph harmonic decomposition:
    for (i=0;i<E->sphere.hindice;i++) sphc[i] = sphs[i] = 0.0 ;
    area = 0.0;
    for (j=1;j<=E->sphere.caps_per_proc;j++)
    for (es=1;es<=E->lmesh.snel;es++) {
        el = (E->lmesh.elz*es) ;
        for(k=1;k<=onedvpoints[E->mesh.nsd];k++)   {
            temp0[k] = 0.0;
            for(k1=1;k1<=onedvpoints[E->mesh.nsd];k1++)   
                temp0[k] +=  X3[j][E->sien[j][es].node[k1]] 
                           * E->M.vpt[GMVINDEX(k1,k)];
        }

        area += E->gDA1[j][es].vpt[5];
        for(k=1;k<=onedvpoints[E->mesh.nsd];k++)
            temp0[k] = temp0[k]*E->gDA1[j][es].vpt[k];

        //for (ll=1; ll<=l_max; ll++)   
        for (ll=0; ll<=l_max; ll++)   
        for (mm=0; mm <= ll;  mm++) { 
            p = E->sphere.hindex[ll][mm];
            for(k=1;k<=onedvpoints[E->mesh.nsd];k++)   {
                sphc[p]+= temp0[k] *E->Sph_Harm_Tblcs[j][mm][es].cs[k]
                                   *E->Sph_Harm_Tbllm[j][p][es].lm[k];
                sphs[p]+= temp0[k] *E->Sph_Harm_Tblcs[j][mm][es].sn[k]
                                   *E->Sph_Harm_Tbllm[j][p][es].lm[k];
            }
        }   // end l,m loop
    }       // end surface-element loop
    sphs[E->sphere.hindice] = area;
    sum_across_surf_sph2(E,sphc,sphs);
    area = 4.0*M_PI/sphs[E->sphere.hindice];
    for (i=0;i<E->sphere.hindice;i++)    {
        sphc[i] = sphc[i]*area;
        sphs[i] = sphs[i]*area;
    }

    // the above spherical harm. expansion can also be done like this:
    //sphere_expansion (E,1,X3,sphc,sphs);
    // (but it takes ll from 1 to E->sphere.output_llmax)
    parallel_process_sync(E);

    // now re-sum back into X3:
    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (j=1;j<=E->lmesh.nsf;j++) {
        X3[m][j] = 0.0; 
        for (ll=0; ll<=l_max; ll++)   
        //for (ll=1; ll<=l_max; ll++)   
        for (mm=0; mm <= ll;  mm++) { 
            p = E->sphere.hindex[ll][mm]; 
            X3[m][j] += E->Tbl_lm[m][p][j] *
                        ( E->Tbl_cs[m][mm][j]* sphc[p] 
                        + E->Tbl_sn[m][mm][j]* sphs[p] ) ;
        } // end l,m loop
     }    // end all-node loop

    // print some messages
    if (verbose && E->parallel.me==0)
    for (i=0;i<2;i++) {
        fp = (i==0)? stderr : E->fp ;
        fprintf(fp,"Performed spherical harmonic truncation in %g seconds.\n",
                   CPU_time0()-time); 
        fflush(fp);
    }

    return;
}

/****************************************************************************
 * total_surface_integral(E,X3,ic)
 * Computes the surface integral of X3[m][j] 
 *   (where m is in sphere.caps_per_proc and j is a surface node index)
 * Specifically it returns  integral( X3.dO )  where dO is an element of
 *   solid angle, and the integration is over all 4pi.
 * ic indicates surface (ic=1) or CMB (ic=0).
 ***************************************************************************/
double total_surface_integral(E,X3,ic)
    struct All_variables *E;
    double **X3;
    int ic;  // ic=1 for surface, ic=0 for CMB
{
    const int dims = E->mesh.nsd;
    const int lev = E->mesh.levmax;
    int e,n,m,i,j,k,d,nint,nox,noy,el,elx,ely,elz,j1,j2,i1,i2,k1,k2,nproc;
    int top,lnode[5],iroot;
    double Have[3],temp[3],integral;
    struct Shape_function1 M;
    struct Shape_function1_dA dGamma;

    nox = E->lmesh.nox;
    noy = E->lmesh.noy;
    elx = E->lmesh.elx;
    ely = E->lmesh.ely;
    elz = E->lmesh.elz;

    temp[0] = temp[1] = temp[2] = 0.0;
    Have[0] = Have[1] = Have[2] = 0.0;
    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (k=1;k<=ely;k++)
    for (j=1;j<=elx;j++)     {
        el = j + (k-1)*elx;

        for(d=1;d<=onedvpoints[E->mesh.nsd];d++)
            lnode[d] = E->sien[m][el].node[d];

        if (ic==0)   {
            for(nint=1;nint<=onedvpoints[E->mesh.nsd];nint++)
            for(d=1;d<=onedvpoints[E->mesh.nsd];d++)
                temp[1] += X3[m][lnode[d]] * E->M.vpt[GMVINDEX(d,nint)] 
                          * E->gDA0[m][el].vpt[nint];
            temp[2] += E->gDA0[m][el].vpt[5];
        }
        else if (ic==1)  {
            for(nint=1;nint<=onedvpoints[E->mesh.nsd];nint++)
            for(d=1;d<=onedvpoints[E->mesh.nsd];d++)
                temp[1] += X3[m][lnode[d]] * E->M.vpt[GMVINDEX(d,nint)] 
                          * E->gDA1[m][el].vpt[nint];
            temp[2] += E->gDA1[m][el].vpt[5];
        }
    }   /* end of j  and k, and m  */

    /* determine which processors should get the message from me for 
       computing the layer averages */

  MPI_Allreduce(temp,Have,3,MPI_DOUBLE,MPI_SUM,E->parallel.horizontal_comm);

    if(Have[2] != 0.0) 
        integral = Have[1]/Have[2];

    integral *= 4.0 * M_PI;  // convert surface average to total integral

    return (integral); 
}

/* =============================================================================
 * load_to_CM(E,X_surf,X_cmb,add_incr_CM)
 * =============================================================================
 * Calculate the incremental center-of-mass due to the incremental loads at
 * the surface and CMB (X_surf, and X_cmb), then add a load to these fields
 * which keeps the CM of the earth-load system at the origin of the coordinate
 * system.  
 *
 * E->data.CM_incr will be affected depending on the value of add_incr_CM:
 *    0 => E->data.CM_incr is not changed
 *    1 => the computed CM is put into E->data.CM_incr 
 *    2 => the computed CM is added to E->data.CM_incr 
 *
 * This routine and the next two were added by Archie (August 2003) so that we
 * may compute the degree-1 response.
 * Usage notes:
 *    The X arrays must be non-dimensional stresses at surface and CMB.
 *    X_surf should include the load.
 *    The X arrays should have zero average.
 *    This routine requires that both the surface and cmb arrays are on
 *      the same processor (ie, only compatable with nprocz==1).
 * =============================================================================
 */
void load_to_CM(E,X_surf,X_cmb, add_incr_CM)
    struct All_variables *E;
    double **X_surf;
    double **X_cmb;
    int add_incr_CM ; // flag to control if incr CM is added to total CM
{
    const int dims = E->mesh.nsd;
    const int lev = E->mesh.levmax;
    int e,n,m,i,j,k,d,d1,el,nproc;
    int lnode[5];
    double temp[3],tempx[5],tempy[5],tempz[5],temps[5],tempb[5],temp0[3];
    //
    static int been_here=0;
    double tmp0, cm[3];      // x,y,z coordinates of CM
    double t0, f0, r0; // theta,phi,radius coordinates of CM
    static double c_surf, c_cmb, earth_mass ;
    double xy_mag, load_height_surf, load_height_cmb;
    double t1,f1,rel_cos;
    double x_surf, x_cmb ;
    FILE *fp;
    
    double mn,mx;
    double **X;

    if (!been_here) {
        been_here = 1;
        earth_mass =  4.0*M_PI / 3.0 * pow(E->sphere.dradius,3) 
                     *(   pow(E->sphere.ri,3) * E->data.density_below
                        + pow(E->sphere.ro,3) * E->data.density      
                        - pow(E->sphere.ri,3) * E->data.density   );
        // scale factors:
        c_surf  = 1.0/ E->ve_data_cont.surf_scaling ; // nondim stress -> nondim height
        c_surf *= E->sphere.dradius ;    // -> dim'ful height 
        c_surf *= E->data.density ;      // -> dim'ful mass/area
        c_surf *= ( E->sphere.dradius*E->sphere.dradius ) / earth_mass; 
        c_surf *= E->sphere.ro ;
        if (1) {  // these two methods are identical
            c_cmb   = 1.0/ E->ve_data_cont.botm_scaling ;
            c_cmb  *= E->sphere.dradius ;    // -> dim'ful height 
            c_cmb  *= (E->data.density_below - E->data.density) ;
        } else {
            c_cmb   = E->ve_data_cont.shear_mod ;   // rho.g.h
            c_cmb  *= 1.0/E->data.grav_acc ;
        }
        c_cmb  *= ( E->sphere.dradius*E->sphere.dradius ) / earth_mass; 
        c_cmb  *= E->sphere.ri ;
        // initialize total CM:
        E->ve_data_cont.CM[0] = 0.0 ;
        E->ve_data_cont.CM[1] = 0.0 ;
        E->ve_data_cont.CM[2] = 0.0 ;
    }

    cm[0] = cm[1] = cm[2] = 0.0;
    temp[0] = temp[1] = temp[2] = 0.0;
    temp0[0] = temp0[1] = temp0[2] = 0.0;

    /* Calculate the integrals
     *      x_cm = Integral( (X_s*r_s + X_b*r_b) * sin(t) * cos(f) * dO )
     *      y_cm = Integral( (X_s*r_s + X_b*r_b) * sin(t) * cos(f) * dO )
     *      z_cm = Integral( (X_s*r_s + X_b*r_b) * cos(t) * dO )
     *  where r_* is the radius at the surface (s) and cmb (b),
     *        X_* is the mass per unit-solid-angle at the surface and cmb,
     *        dO is d-omega (element of solid angle),
     *        t and f are theta and phi.
     *  note: this technique actually calculates Integral( ... r^2 dO),
     *        but since we're using mass/area, for X_*, we've already divided
     *        by r^2 in the constants c_surf and c_cmb.
     */
    for (m=1; m<=E->sphere.caps_per_proc; m++)
    for (k=1; k<=E->lmesh.ely; k++)
    for (j=1; j<=E->lmesh.elx; j++) {
        el = j + (k-1)*E->lmesh.elx;
        for(d=1;d<=onedvpoints[E->mesh.nsd];d++) {
            lnode[d] = E->sien[m][el].node[d];
            tempx[d] = 0.0;
            tempy[d] = 0.0;
            tempz[d] = 0.0;
            temps[d] = 0.0;
            tempb[d] = 0.0;
            }

        for(d =1; d <=onedvpoints[E->mesh.nsd]; d++) 
        for(d1=1; d1<=onedvpoints[E->mesh.nsd]; d1++) {
            i = lnode[d1]*E->lmesh.noz;
            tempx[d] += E->X[E->mesh.levmax][m][1][i]*E->M.vpt[GMVINDEX(d1,d)];
            tempy[d] += E->X[E->mesh.levmax][m][2][i]*E->M.vpt[GMVINDEX(d1,d)];
            tempz[d] += E->X[E->mesh.levmax][m][3][i]*E->M.vpt[GMVINDEX(d1,d)];
            temps[d] += X_surf[m][lnode[d1]]*E->M.vpt[GMVINDEX(d1,d)];
            tempb[d] += X_cmb[m][lnode[d1]]*E->M.vpt[GMVINDEX(d1,d)];
            }

        if (E->parallel.me_loc[3] == E->parallel.nprocz-1) {
          for(d =1; d <=onedvpoints[E->mesh.nsd]; d++)   { 
            x_surf = c_surf * temps[d] 
                      * E->gDA1[m][el].vpt[d];

            temp[0] +=  x_surf * tempx[d];       // x
            temp[1] +=  x_surf * tempy[d];       // y
            temp[2] +=  x_surf * tempz[d];       // z
          }
        }
        if (E->parallel.me_loc[3] == 0) {
          for(d =1; d <=onedvpoints[E->mesh.nsd]; d++)   { 
            x_cmb  = c_cmb  * tempb[d] 
                      * E->gDA0[m][el].vpt[d];

            temp0[0] += x_cmb * tempx[d];       // x
            temp0[1] += x_cmb * tempy[d];       // y
            temp0[2] += x_cmb * tempz[d];       // z
          }
        }

    }  // end loop over surface elements

 if (add_incr_CM>0)  {
       fprintf(E->fp_out,"CMtemp %g %g %g %d\n",temp[0],temp[1],temp[2],E->parallel.me);
       fprintf(E->fp_out,"CMtemp %g %g %g %d\n",temp0[0],temp0[1],temp0[2],E->parallel.me);
  }

    temp[0] = temp[0] + temp0[0];
    temp[1] = temp[1] + temp0[1];
    temp[2] = temp[2] + temp0[2];


        // put totals of temp into cm:
        MPI_Allreduce(temp,cm,3,MPI_DOUBLE,MPI_SUM,E->parallel.world);

 if (add_incr_CM>0)  {
       fprintf(E->fp_out,"CMtemp %g %g %g %d\n",temp[0],temp[1],temp[2],E->parallel.me);
       fprintf(E->fp_out,"CM %g %g %g %d %d %d\n",cm[0],cm[1],cm[2],E->parallel.me,E->monitor.solution_cycles,add_incr_CM);
       fflush(E->fp_out);
  }
    // now 'cm' has the x,y,z distances of the CM from the origin

    if (add_incr_CM==1) 
        for (i=0;i<3;i++) E->ve_data_cont.CM_incr[i] = cm[i] ;
    else if (add_incr_CM==2) 
        for (i=0;i<3;i++) E->ve_data_cont.CM_incr[i] += cm[i] ;
    // If add_incr_CM is 0, do nothing.

    // convert CM coordinates to spherical (r0,t0,f0):
    xy_mag = sqrt( cm[0]*cm[0] + cm[1]*cm[1] );
    r0 = sqrt( cm[2]*cm[2] + xy_mag*xy_mag );
    if (r0) {
        t0 = acos( cm[2]/r0 );
        f0 = acos( cm[0]/xy_mag );
        if ( cm[1]<0 ) f0 = 2.0*M_PI - f0 ;
    } 
    else t0 = f0 = 0.0 ;

    // Finally, apply additional surface loads to X_surf/cmb so that the CM is
    // at the origin:
    load_height_surf = r0 * E->ve_data_cont.surf_scaling  ;
    load_height_cmb  = r0 * E->ve_data_cont.botm_scaling  ;
    if (E->parallel.me_loc[3]==E->parallel.nprocz-1) {
      for (m=1;m<=E->sphere.caps_per_proc;m++)
      for (j=1;j<=E->lmesh.nsf;j++)    {
        i = j*E->lmesh.noz;
        t1 = E->SX[E->mesh.levmax][m][1][i];
        f1 = E->SX[E->mesh.levmax][m][2][i];
        rel_cos = -1.0 * (cos(t1)*cos(t0) + sin(t1)*sin(t0)*cos(f0-f1)) ;
        X_surf[m][j] += load_height_surf * rel_cos ;
        }
    }
    if (E->parallel.me_loc[3]==0) {
      for (m=1;m<=E->sphere.caps_per_proc;m++)
      for (j=1;j<=E->lmesh.nsf;j++)    {
        i = j*E->lmesh.noz;
        t1 = E->SX[E->mesh.levmax][m][1][i];
        f1 = E->SX[E->mesh.levmax][m][2][i];
        rel_cos = -1.0 * (cos(t1)*cos(t0) + sin(t1)*sin(t0)*cos(f0-f1)) ;
        X_cmb[m][j]  += load_height_cmb  * rel_cos ;
        }
    }

    return; 
}



// add by tao: follow geruo, calculate deg 1 deformation in calculate_potential_comp.
// icon always 1 when called.
void load_to_CM_grav_comp(struct All_variables *E, int load_only, int add_incr_CM)
{
	void calculate_potential_comp_OR_correct_CM(struct All_variables *E, double ** all_potential, int load_only, int icon, int Only_Correct_CM, int add_incr_CM);

  if(E->parallel.me==0)
	  printf("debug deg-1, now with the load_to_CM_grav_comp implementation\n");
	
	calculate_potential_comp_OR_correct_CM(E, E->slice_ve.all_potential, load_only, 1, 1, add_incr_CM);

}

// use Geruo's method to compute degree-1 gravitational potential at the surface
// that is then used to determine the CM

void load_to_CM_grav(E,X_surf,X_cmb, add_incr_CM)
    struct All_variables *E;
    double **X_surf;
    double **X_cmb;
    int add_incr_CM ; // flag to control if incr CM is added to total CM
{
    void calculate_potential_deg1_2();

    const int dims = E->mesh.nsd;
    const int lev = E->mesh.levmax;
    int e,n,m,i,j,k,d,d1,el,nproc;
    int lnode[5];
    double deg1phi[4], temp[3],tempx[5],tempy[5],tempz[5],temps[5],tempb[5],temp0[3];
    //
    static int been_here=0;
    double tmp0, cm[3];      // x,y,z coordinates of CM
    double t0, f0, r0; // theta,phi,radius coordinates of CM
    static double c_surf, c_cmb, earth_mass ;
    double xy_mag, load_height_surf, load_height_cmb;
    double t1,f1,rel_cos;
    double x_surf, x_cmb ;
    FILE *fp;
    
    double mn,mx;
    double **X;

    if (!been_here) {
        been_here = 1;
        earth_mass =  4.0*M_PI / 3.0 * pow(E->sphere.dradius,3) 
                     *(   pow(E->sphere.ri,3) * E->data.density_below
                        + pow(E->sphere.ro,3) * E->data.density      
                        - pow(E->sphere.ri,3) * E->data.density   );
        earth_mass = earth_mass/(E->data.density*pow(E->sphere.dradius,3));
        // initialize total CM:
        E->ve_data_cont.CM[0] = 0.0 ;
        E->ve_data_cont.CM[1] = 0.0 ;
        E->ve_data_cont.CM[2] = 0.0 ;
    }

    calculate_potential_deg1_2(E, X_surf, X_cmb);

    cm[0] = -E->ve_data_cont.CM_pot[1]*sqrt(12.0*M_PI)/earth_mass;
    cm[1] = -E->ve_data_cont.CM_pot[2]*sqrt(12.0*M_PI)/earth_mass;
    cm[2] =  E->ve_data_cont.CM_pot[0]*sqrt(12.0*M_PI)/earth_mass;

 if (add_incr_CM>0)  {
       fprintf(E->fp_out,"CM %g %g %g %d %d %d\n",cm[0],cm[1],cm[2],E->parallel.me,E->monitor.solution_cycles,add_incr_CM);
       fflush(E->fp_out);
  }
    // now 'cm' has the x,y,z distances of the CM from the origin

    if (add_incr_CM==1) 
        for (i=0;i<3;i++) E->ve_data_cont.CM_incr[i] = cm[i] ;
    else if (add_incr_CM==2) 
        for (i=0;i<3;i++) E->ve_data_cont.CM_incr[i] += cm[i] ;
    // If add_incr_CM is 0, do nothing.

    // convert CM coordinates to spherical (r0,t0,f0):
    xy_mag = sqrt( cm[0]*cm[0] + cm[1]*cm[1] );
    r0 = sqrt( cm[2]*cm[2] + xy_mag*xy_mag );
    if (r0) {
        t0 = acos( cm[2]/r0 );
        f0 = acos( cm[0]/xy_mag );
        if ( cm[1]<0 ) f0 = 2.0*M_PI - f0 ;
    } 
    else t0 = f0 = 0.0 ;

    // Finally, apply additional surface loads to X_surf/cmb so that the CM is
    // at the origin:
    load_height_surf = r0 * E->ve_data_cont.surf_scaling  ;
    load_height_cmb  = r0 * E->ve_data_cont.botm_scaling  ;
    if (E->parallel.me_loc[3]==E->parallel.nprocz-1) {
      for (m=1;m<=E->sphere.caps_per_proc;m++)
      for (j=1;j<=E->lmesh.nsf;j++)    {
        i = j*E->lmesh.noz;
        t1 = E->SX[E->mesh.levmax][m][1][i];
        f1 = E->SX[E->mesh.levmax][m][2][i];
        rel_cos = -1.0 * (cos(t1)*cos(t0) + sin(t1)*sin(t0)*cos(f0-f1)) ;
        X_surf[m][j] += load_height_surf * rel_cos ;
        }
    }
    if (E->parallel.me_loc[3]==0) {
      for (m=1;m<=E->sphere.caps_per_proc;m++)
      for (j=1;j<=E->lmesh.nsf;j++)    {
        i = j*E->lmesh.noz;
        t1 = E->SX[E->mesh.levmax][m][1][i];
        f1 = E->SX[E->mesh.levmax][m][2][i];
        rel_cos = -1.0 * (cos(t1)*cos(t0) + sin(t1)*sin(t0)*cos(f0-f1)) ;
        X_cmb[m][j]  += load_height_cmb  * rel_cos ;
        }
    }

    return; 
 }


/* =============================================================================
 * calculate_CM(E,X_surf,X_cmb)
 * =============================================================================
 * Just calculate and report CM displacement (used for debugging and 
 * developement).  This is identical to load_to_CM(), except that it 
 * doesn't change the surface arrays or E->data.CM*.
 * =============================================================================
 */
void calculate_CM(E,X_surf,X_cmb)
    struct All_variables *E;
    double **X_surf;
    double **X_cmb;
{
    const int dims = E->mesh.nsd;
    const int lev = E->mesh.levmax;
    int e,n,m,i,j,k,d,d1,el,nproc;
    int lnode[5];
    double temp[3];
    //
    static int been_here=0;
    double tmp0, cm[3];      // x,y,z coordinates of CM
    double t0, f0, r0; // theta,phi,radius coordinates of CM
    static double c_surf, c_cmb, earth_mass ;
    double xy_mag, load_height_surf, load_height_cmb;
    double t1,f1,rel_cos;
    double x_surf, x_cmb ;
    FILE *fp;

    if (!been_here) {
        been_here = 1;
        //c = sqrt((4.0*M_PI)/3.0) ;
        earth_mass =  4.0*M_PI / 3.0 * pow(E->sphere.dradius,3) 
                     *(   pow(E->sphere.ri,3) * E->data.density_below
                        + pow(E->sphere.ro,3) * E->data.density      
                        - pow(E->sphere.ri,3) * E->data.density   );
        // scale factors:
        c_surf  = 1.0/ E->ve_data_cont.surf_scaling ; // nondim stress -> nondim height
        c_surf *= E->sphere.dradius ;    // -> dim'ful height 
        c_surf *= E->data.density ;      // -> dim'ful mass/area
        c_surf *= ( E->sphere.dradius*E->sphere.dradius ) / earth_mass; 
        c_surf *= E->sphere.ro ;
        c_cmb   = E->ve_data_cont.shear_mod ;   // rho.g.h
        c_cmb  *= 1.0/E->data.grav_acc ;
        c_cmb  *= ( E->sphere.dradius*E->sphere.dradius ) / earth_mass; 
        c_cmb  *= E->sphere.ri ;
        // initialize total CM:
        E->ve_data_cont.CM[0] = 0.0 ;
        E->ve_data_cont.CM[1] = 0.0 ;
        E->ve_data_cont.CM[2] = 0.0 ;
    }

    cm[0] = cm[1] = cm[2] = 0.0;
    temp[0] = temp[1] = temp[2] = 0.0;

    /* Calculate the integrals
     *      x_cm = Integral( (X_s*r_s + X_b*r_b) * sin(t) * cos(f) * dO )
     *      y_cm = Integral( (X_s*r_s + X_b*r_b) * sin(t) * cos(f) * dO )
     *      z_cm = Integral( (X_s*r_s + X_b*r_b) * cos(t) * dO )
     *  where r_* is the radius at the surface (s) and cmb (b),
     *        X_* is the mass per unit-solid-angle at the surface and cmb,
     *        dO is d-omega (element of solid angle),
     *        t and f are theta and phi.
     *  note: this technique actually calculates Integral( ... r^2 dO),
     *        but we've divided by r^2 in the constants c_surf and c_cmb.
     */
    for (m=1; m<=E->sphere.caps_per_proc; m++)
    for (k=1; k<=E->lmesh.ely; k++)
    for (j=1; j<=E->lmesh.elx; j++) {
        el = j + (k-1)*E->lmesh.elx;
        for(d=1;d<=onedvpoints[E->mesh.nsd];d++)
            lnode[d] = E->sien[m][el].node[d];
        for(d =1; d <=onedvpoints[E->mesh.nsd]; d++) 
        for(d1=1; d1<=onedvpoints[E->mesh.nsd]; d1++) {
            i = lnode[d1]*E->lmesh.noz;
            t1 = E->SX[E->mesh.levmax][m][1][i];
            f1 = E->SX[E->mesh.levmax][m][2][i];
            x_surf = c_surf * X_surf[m][lnode[d1]] 
                      * E->gDA1[m][el].vpt[d] * E->M.vpt[GMVINDEX(d1,d)] ;
            x_cmb  = c_cmb  *  X_cmb[m][lnode[d1]] 
                      * E->gDA0[m][el].vpt[d] * E->M.vpt[GMVINDEX(d1,d)] ;

            temp[0] += ( x_surf + x_cmb ) * sin(t1) * cos(f1) ;  // x
            temp[1] += ( x_surf + x_cmb ) * sin(t1) * sin(f1) ;  // y
            temp[2] += ( x_surf + x_cmb ) * cos(t1) ;            // z
        }

    }  // end loop over surface elements

        // put totals of temp into cm:
   MPI_Allreduce(temp,cm,3,MPI_DOUBLE,MPI_SUM,E->parallel.horizontal_comm);

    // now 'cm' has the x,y,z distances of the CM from the origin

    for (i=0;i<3;i++) E->ve_data_cont.CM_incr[i] = cm[i] ;
    
    xy_mag = sqrt( cm[0]*cm[0] + cm[1]*cm[1] );
    r0 = sqrt( cm[2]*cm[2] + xy_mag*xy_mag );
    if (r0) {
        t0 = acos( cm[2]/r0 );
        f0 = acos( cm[0]/xy_mag );
        if ( cm[1]<0 ) f0 = 2.0*M_PI - f0 ;
    } 
    else t0 = f0 = 0.0 ;

    if (E->parallel.me==0)
    for (i=0;i<2;i++) {
        fp = (i==0)? stderr : E->fp ;
        fprintf(fp,"CM_incr report (step %d )  th: %g  fi: %g  r: %g\n",
                    E->monitor.solution_cycles,
                    t0*180./M_PI, f0*180./M_PI, r0 );
        fprintf(fp,"CM report      (step %d )  x: %g  y: %g  z: %g\n",
                    E->monitor.solution_cycles,
                    E->ve_data_cont.CM[0], E->ve_data_cont.CM[1], E->ve_data_cont.CM[2] );
        fflush(fp);
    }

    return; 
}

/* =============================================================================
 * shift_to_CM(E,X)
 * =============================================================================
 * Takes the surface field X and subtracts the current incremental center of
 * mass location (stored in E->data.CM_incr).  This keeps the field X in a
 * reference frame where the origin is at the center of mass of the earth-load
 * system.  At the end of a timestep, there has been some CM motion (which is
 * removed by placing extra degree-1 loads in load_to_CM). So this routine is
 * invoked to remove that little motion from the field X.
 * Usage notes:
 *    The field X should be nondimensional distances.
 * =============================================================================
 */
void shift_to_CM(E,X)
    struct All_variables *E;
    double **X;
{
    int m,j,i;
    double t1,f1, rel_cos;
    double xy_mag;
    double t0,f0,r0 ; // coordinates of CM

    xy_mag = sqrt(  E->ve_data_cont.CM_incr[0]*E->ve_data_cont.CM_incr[0] 
                  + E->ve_data_cont.CM_incr[1]*E->ve_data_cont.CM_incr[1] );
    r0 = sqrt( E->ve_data_cont.CM_incr[2]*E->ve_data_cont.CM_incr[2] + xy_mag*xy_mag );
    if (r0) {
        t0 = acos( E->ve_data_cont.CM_incr[2]/r0 );
        f0 = acos( E->ve_data_cont.CM_incr[0]/xy_mag );
        if ( E->ve_data_cont.CM_incr[1]<0 ) f0 = 2.0*M_PI - f0 ;
    } 
    else t0 = f0 = 0.0 ;

    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (j=1;j<=E->lmesh.nsf;j++)    {
       i = j*E->lmesh.noz;
       t1 = E->SX[E->mesh.levmax][m][1][i];
       f1 = E->SX[E->mesh.levmax][m][2][i];
       rel_cos = cos(t1)*cos(t0) + sin(t1)*sin(t0)*cos(f0-f1) ;
       X[m][j] -= r0 * rel_cos ;
    }
    return;
}

// similar to Geruo's to modify displacement field U as well, in addition to 
// topography at surface XS and CMB XB

void shift_U_to_CM(struct All_variables * E,double ** XS,double ** XB)
 //   struct All_variables *E;
 //   double **XS, **XB;
{
    int m,j,i,k;
    double t1,f1, rel_cos;
    double xy_mag, u1,u2,u3;
    double v1,v2,v3;
    double t0,f0,r0 ; // coordinates of CM

    xy_mag = sqrt(  E->ve_data_cont.CM_incr[0]*E->ve_data_cont.CM_incr[0] 
                  + E->ve_data_cont.CM_incr[1]*E->ve_data_cont.CM_incr[1] );
    r0 = sqrt( E->ve_data_cont.CM_incr[2]*E->ve_data_cont.CM_incr[2] + xy_mag*xy_mag );
    if (r0) {
        t0 = acos( E->ve_data_cont.CM_incr[2]/r0 );
        f0 = acos( E->ve_data_cont.CM_incr[0]/xy_mag );
        if ( E->ve_data_cont.CM_incr[1]<0 ) f0 = 2.0*M_PI - f0 ;
    } 
    else t0 = f0 = 0.0 ;

    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (j=1;j<=E->lmesh.nsf;j++)    
    for (k=1;k<=E->lmesh.noz;k++)    {
       i = k+(j-1)*E->lmesh.noz;
       t1 = E->SX[E->mesh.levmax][m][1][i];
       f1 = E->SX[E->mesh.levmax][m][2][i];
       rel_cos = -sin(t1)*cos(t0) + cos(t1)*sin(t0)*cos(f0-f1);
       u1 = E->sphere.cap[m].V[1][i] - r0 * rel_cos ;

       rel_cos = sin(t0)*sin(f0-f1);
       u2 = E->sphere.cap[m].V[2][i] - r0 * rel_cos ;

       rel_cos = cos(t1)*cos(t0) + sin(t1)*sin(t0)*cos(f0-f1) ;
       u3 = E->sphere.cap[m].V[3][i] - r0 * rel_cos ;


/*
    if ( k==E->lmesh.noz && E->parallel.me_loc[3]==E->parallel.nprocz-1) {
       E->sphere.cap[m].total_VS[1][j] += E->sphere.cap[m].V[1][i];
       E->sphere.cap[m].total_VS[2][j] += E->sphere.cap[m].V[2][i];
       E->sphere.cap[m].total_VS[3][j] += E->sphere.cap[m].V[3][i];  
     }
*/
	    if ( k==E->lmesh.noz && E->parallel.me_loc[3]==E->parallel.nprocz-1) {
	       E->sphere.cap[m].total_VS[1][j] += u1;
	       E->sphere.cap[m].total_VS[2][j] += u2;
	       E->sphere.cap[m].total_VS[3][j] += u3;
	     }

/*       if ( k==E->lmesh.noz && E->parallel.me_loc[3]==E->parallel.nprocz-1) {
          v1 = E->sphere.cap[m].V[1][i];
          v2 = E->sphere.cap[m].V[2][i];
          v3 = E->sphere.cap[m].V[3][i];
          v1 = u1;
          v2 = u2;
          v3 = u3;
// total_V[1]-[3] are in x, y and z Cartesian
          E->sphere.cap[m].total_VC[1][j] += sin(t1)*cos(f1)*v3 + cos(t1)*cos(f1)*v1 - sin(f1)*v2;
          E->sphere.cap[m].total_VC[2][j] += sin(t1)*sin(f1)*v3 + cos(t1)*sin(f1)*v1 + cos(f1)*v2;
          E->sphere.cap[m].total_VC[3][j] += cos(t1)*v3         - sin(t1)*v1;
          v1 = E->sphere.cap[m].total_VC[1][j];
          v2 = E->sphere.cap[m].total_VC[2][j];
          v3 = E->sphere.cap[m].total_VC[3][j];
          E->sphere.cap[m].total_VS[3][j] = sin(t1)*cos(f1)*v1 + sin(t1)*sin(f1)*v2 + cos(t1)*v3;
          E->sphere.cap[m].total_VS[1][j] = cos(t1)*cos(f1)*v1 + cos(t1)*sin(f1)*v2 - sin(t1)*v3;
          E->sphere.cap[m].total_VS[2][j] = -sin(f1)*v1         + cos(f1)*v2;  
          }
*/

       E->sphere.cap[m].V[1][i] = u1;
       E->sphere.cap[m].V[2][i] = u2; 
       E->sphere.cap[m].V[3][i] = u3;


       if ( k==E->lmesh.noz && E->parallel.me_loc[3]==E->parallel.nprocz-1) {
          XS[m][j] -= r0 * rel_cos ;
          }
       if ( k==1 && E->parallel.me_loc[3]==0)  {
          XB[m][j] -= r0 * rel_cos ;
          }
    }
    return;
}

void shift_V_to_CM(E)
    struct All_variables *E;
{
    int m,j,i;
    double t1,f1, rel_cos;
    double xy_mag;
    double t0,f0,r0 ; // coordinates of CM

    xy_mag = sqrt(  E->ve_data_cont.CM_incr[0]*E->ve_data_cont.CM_incr[0]
                  + E->ve_data_cont.CM_incr[1]*E->ve_data_cont.CM_incr[1] );
    r0 = sqrt( E->ve_data_cont.CM_incr[2]*E->ve_data_cont.CM_incr[2] + xy_mag*xy_mag );
    if (r0) {
        t0 = acos( E->ve_data_cont.CM_incr[2]/r0 );
        f0 = acos( E->ve_data_cont.CM_incr[0]/xy_mag );
        if ( E->ve_data_cont.CM_incr[1]<0 ) f0 = 2.0*M_PI - f0 ;
    }
    else t0 = f0 = 0.0 ;

    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (j=1;j<=E->lmesh.nno;j++)    {
       t1 = E->SX[E->mesh.levmax][m][1][j];
       f1 = E->SX[E->mesh.levmax][m][2][j];
       // r - direction:
       rel_cos = cos(t1)*cos(t0) + sin(t1)*sin(t0)*cos(f0-f1) ;
       E->sphere.cap[m].total_V[3][j] -= r0 * rel_cos ;
       // t - direction:
       rel_cos = -sin(t1)*cos(t0) + cos(t1)*sin(t0)*cos(f0-f1);
       E->sphere.cap[m].total_V[1][j] -= r0 * rel_cos ;
       // f - direction:
       rel_cos = sin(t0)*sin(f0-f1);
       E->sphere.cap[m].total_V[2][j] -= r0 * rel_cos ;
    }
    return;
}


void remove_average(E,X3,ic)
    struct All_variables *E;
    double **X3;
    int ic;
{
    const int dims = E->mesh.nsd;
    const int lev = E->mesh.levmax;
    int e,n,m,i,j,k,d,nint,nox,noy,el,elx,ely,elz,j1,j2,i1,i2,k1,k2,nproc;
    int top,lnode[5],iroot;
    double Have[3],temp[3],average;
    struct Shape_function1 M;
    struct Shape_function1_dA dGamma;

    nox = E->lmesh.nox;
    noy = E->lmesh.noy;
    elx = E->lmesh.elx;
    ely = E->lmesh.ely;
    elz = E->lmesh.elz;

    temp[0] = temp[1] = temp[2] = 0.0;
    Have[0] = Have[1] = Have[2] = 0.0;
    for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (k=1;k<=ely;k++)
            for (j=1;j<=elx;j++)     {
                el = j + (k-1)*elx;

                for(d=1;d<=onedvpoints[E->mesh.nsd];d++)
                    lnode[d] = E->sien[m][el].node[d];

                if (ic==0)   {
                    for(nint=1;nint<=onedvpoints[E->mesh.nsd];nint++)
                        for(d=1;d<=onedvpoints[E->mesh.nsd];d++)
                            temp[1] += X3[m][lnode[d]] * E->M.vpt[GMVINDEX(d,nint)] 
                                * E->gDA0[m][el].vpt[nint];
                    temp[2] += E->gDA0[m][el].vpt[5];
                }
                else if (ic==1)  {
                    for(nint=1;nint<=onedvpoints[E->mesh.nsd];nint++)
                        for(d=1;d<=onedvpoints[E->mesh.nsd];d++)
                            temp[1] += X3[m][lnode[d]] * E->M.vpt[GMVINDEX(d,nint)] 
                                * E->gDA1[m][el].vpt[nint];
                    temp[2] += E->gDA1[m][el].vpt[5];
                }

            }   /* end of j  and k, and m  */

    /* determine which processors should get the message from me for 
       computing the layer averages */

   MPI_Allreduce(temp,Have,3,MPI_DOUBLE,MPI_SUM,E->parallel.horizontal_comm);

    if(Have[2] != 0.0) 
        average = Have[1]/Have[2];

    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(k=1;k<=noy;k++)      
            for(j=1;j<=nox;j++)     {
                n = j+(k-1)*nox;
                X3[m][n] -= average;  
            }


    return; 
}


/* =============================================================================
 * polar_wander_effects (E,X_surf,X_cmb,potential_surf,potential_cmb)
 * =============================================================================
 * (This routine is only called by calculate_potential.)
 * Adds an incremental centrifugal potential (due to the surface masses passed
 * in) to the potential fields passed in. The centrifugal potential is the
 * potential due to motion of the axis of rotation (polar wander) which is due
 * to mass variations from loading and deformation. It is given by
 *      -  R^2 * omega^2 * sin(t)*cos(t)*[ m1 cos(f) + m2 sin(f)  ]
 * where omega is 2pi/(24*3600)
 *       m1, m2  give the incremental pole motion:  w = omega*(m1,m2,1+m3)
 *       t,f are theta and phi
 *       R is the radius of the surface or cmb, for potential there
 * The fields X_surf and X_cmb (passed in from calculate_potential) are in terms
 * of nondimensional height of boundary, where X_surf includes the mass of the
 * load. The fields potential_surf and potential_cmb are in terms of dimless
 * potential.
 * Written by Archie, March 2004.
 * =============================================================================
 */
void polar_wander_effects(E,X_surf,X_cmb,potential_surf,potential_cmb)
    struct All_variables *E;
    double **X_surf;
    double **X_cmb;
    double **potential_surf;
    double **potential_cmb;
{
    void find_rotation_axis();
    void find_rotation_axis_grav();
    static int been_here=0;
    double M[3]; // space for M1 and M2 (the rotation axis)
    FILE *fp;
    int i,m,n;
    double t,f;  // theta, phi
    double C,A;   // unperturbed polar & equatorial moments of inertia
    double scale_pot, omega, centri_pot, potential[3],prod[3],temp;
    static double *oldpotl_surf[NCS], *oldpotl_cmb[NCS];
    static double c_surf, c_cmb;   // scale factors for potential 
    static int tstep=0;       // tstep,M1,M2 are for computing the total
    static double M0=0.0;     //      polar wander for testing purposes
    static double M1=0.0;     // 

    if (!been_here) {
        been_here = 1;
        // Concerning units:  First, a nondimensional M1 and M2 are calculated
        // in find_rotation_axis, then this becomes an MKS potential when
        // multiplied by R^2.omega^2 (where R is the radius of the boundary
        // you're working on).
        // To dimensionalize a nondim'l potential (such as E->potential or
        // E->incr_potential), multiply by 4.pi.G.rho.R^2 (where rho is the
        // density change over the boundary, and R is its radius).  Since the
        // cent. pot. we calculate below is dim'ful (MKS), we will divide it by
        // that factor (the R^2  part is cancelled out).
        omega = 2.0*M_PI/(24.0*3600.0);  // rotation_rate
        scale_pot = 1.0/( 4.0*M_PI * E->data.grav_const * E->data.density );
        c_surf = -1.0 * scale_pot * omega*omega ;
        c_cmb = c_surf * E->sphere.ri*E->sphere.ri/(E->sphere.ro*E->sphere.ro);
        E->ve_data_cont.PW_incr[0] = 0.0;
        E->ve_data_cont.PW_incr[1] = 0.0;
        
       for (m=1;m<=E->sphere.caps_per_proc;m++)   {
          oldpotl_surf[m] =(double*)malloc((E->lmesh.nsf+2)*sizeof(double));
          oldpotl_cmb[m] = (double*)malloc((E->lmesh.nsf+2)*sizeof(double));
          for (n=1;n<=E->lmesh.nsf;n++)  { 
             oldpotl_surf[m][n] = 0.0;
             oldpotl_cmb[m][n] = 0.0;
        }

        }

    }

    // get (nondim) M1 and M2 (the incr polar position), in M[0] and M[1]:

    find_rotation_axis(E,X_surf,X_cmb,M);     // Archie's
//    find_rotation_axis_grav(E,X_surf,X_cmb,M);  // Geruo's

//     beta (>0 and <1) is relaxation or damping parameter to help iteration 
    M0 = E->viscosity.beta*M[0] + (1.0-E->viscosity.beta)*E->ve_data_cont.PW_incr[0]; 
    M1 = E->viscosity.beta*M[1] + (1.0-E->viscosity.beta)*E->ve_data_cont.PW_incr[1]; 

    // add new centrifugal potential into incr_potential:

    potential[0]=potential[1]=0.0;
    prod[0]=prod[1]=0.0;
    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (n=1;n<=E->lmesh.nsf;n++)   {
        t = E->SX[E->mesh.levmax][m][1][n*E->lmesh.NOZ[E->mesh.levmax]];
        f = E->SX[E->mesh.levmax][m][2][n*E->lmesh.NOZ[E->mesh.levmax]];
        centri_pot = sin(t) * cos(t) * ( M0*cos(f) + M1*sin(f) ) ;
        if (E->parallel.me_loc[3]==E->parallel.nprocz-1)   {
                potential_surf[m][n] += c_surf * centri_pot ;
                potential[0] += potential_surf[m][n]*potential_surf[m][n];
                temp = potential_surf[m][n]-oldpotl_surf[m][n];
                potential[1] += temp*temp;
                }
        if (E->parallel.me_loc[3]==0)   {
                potential_cmb [m][n] += c_cmb  * centri_pot ;
                potential[0] += potential_cmb[m][n]*potential_cmb[m][n];
                temp = potential_cmb[m][n]-oldpotl_cmb[m][n];
                potential[1] += temp*temp;
                }
      }


   MPI_Allreduce(potential,prod,3,MPI_DOUBLE,MPI_SUM,E->parallel.world);

   E->ve_data_cont.potential_vary_PW = sqrt(prod[1]/prod[0]);

 if (E->parallel.me==0) fprintf(E->fp_out,"pw step=%d %g %g %g %g %g\n",E->monitor.solution_cycles,E->viscosity.beta,M[0],M[1],E->ve_data_cont.PW_incr[0],E->ve_data_cont.PW_incr[1]);
   E->ve_data_cont.PW_incr[0] = M[0];
   E->ve_data_cont.PW_incr[1] = M[1];


    for (m=1;m<=E->sphere.caps_per_proc;m++)
    for (n=1;n<=E->lmesh.nsf;n++)   {
        if (E->parallel.me_loc[3]==E->parallel.nprocz-1)
                oldpotl_surf[m][n] = potential_surf[m][n];
        if (E->parallel.me_loc[3]==0)   
                oldpotl_cmb[m][n] = potential_cmb[m][n];
      }

    return; 
}


// using Geruo's method: from (2,1) gravitational potential 
// that was computed in calculating CM_

void find_rotation_axis_grav(E,X_surf,X_cmb, M)
    struct All_variables *E;
    double **X_surf;
    double **X_cmb;
    double M[3];      // space for M1 and M2
{
    void calculate_potential_deg1_2();
    static int been_here=0;
    static double c_surf, c_cmb, earth_mass ;
    static double M_conv ;
    double omega; // earth's rotation rate

    if (!been_here)  {
        been_here = 1;
        omega = 2.0*M_PI/(24.0*3600.0);  // rotation_rate
        M_conv = (3.0*E->data.grav_const*E->data.density) /
          (omega*omega*E->ve_data_cont.kf)*5.0*sqrt(4.0*M_PI/15.0);
    }

    calculate_potential_deg1_2(E, X_surf, X_cmb);

    M[0] = M_conv*E->ve_data_cont.PW_pot[0];
    M[1] = M_conv*E->ve_data_cont.PW_pot[1];

    return;
  }


void find_rotation_axis(E,X_surf,X_cmb, M)
    struct All_variables *E;
    double **X_surf;
    double **X_cmb;
    double M[3];      // space for M1 and M2
{
    static int been_here=0;
    static double c_surf, c_cmb, earth_mass ;
    static double I_to_M_conv ;  // 1/(C-A)   converts inertia to rotation axis
    double omega; // earth's rotation rate
    double temp[3];
    int e,n,m,i,j,k,d,d1,el,nproc;
    int lnode[5];
    double t1,f1,rel_cos;
    double x_surf, x_cmb ;
    double density_surf,density_cmb;
    double I[3];      // space for moments of inertia I13 and I23
    // variables for parallel integration:

    if (!been_here)  {
        been_here = 1;
        // C and A here are the polar and equatorial moments of intertia.
        // These quantities should give a C-A slightly larger than the
        // hydrostatic limit from rotation effects in order to have a stable
        // rotation axis (or, equivalently, to remove the s=0 mode of the
        // body-tide k love number).  To do this, we take C and A from Lambeck
        // 1980 _Earth's Variable Rotation_ (p27), and then mulitply by CApert
        // (given in the input file). 
        // I_to_M_conv is 1/(C-A), with MKS units.
        //
        //I_to_M_conv = 1.0/ ( ( 0.3306 - 0.3295 )   // C-A
        //                    *( 5.974e24 * 6.37814e6*6.37814e6 ) // M.R^2 
        //                    * E->control.CApert ); // non-hydrostatic effect
   	density_surf = 1.0;
    	density_cmb = (E->data.density_below-E->data.density)/E->data.density;

        omega = 2.0*M_PI/(24.0*3600.0);  // rotation_rate
        I_to_M_conv = ( 3.0 * E->data.grav_const * E->data.density) /
                      (omega*omega*E->ve_data_cont.kf);

        // The fields X_surf and X_cmb are the surface masses, passed in
        c_surf  = -1.0; 
        c_surf  = c_surf/E->ve_data_cont.surf_scaling*density_surf;    //turn load to surface mass density, sigma
        c_cmb   = -E->sphere.ri*E->sphere.ri/(E->sphere.ro*E->sphere.ro);  
        c_cmb  = c_cmb/E->ve_data_cont.botm_scaling*density_cmb;      //turn load to surface mass density, sigma
    }


    /* Calculate the integrals given in the comment above.
     * Note: this technique actually calculates Integral( f r^2 dO), where f is
     * c_surf*X_surf (for example), and r is nondimensinal.
     */
    I[0] = I[1] = I[2] = 0.0;
    temp[0] = temp[1] = temp[2] = 0.0;
    for (m=1; m<=E->sphere.caps_per_proc; m++)
    for (k=1; k<=E->lmesh.ely; k++)
    for (j=1; j<=E->lmesh.elx; j++) {
        el = j + (k-1)*E->lmesh.elx;
        for(d=1;d<=onedvpoints[E->mesh.nsd];d++)
            lnode[d] = E->sien[m][el].node[d];
        for(d =1; d <=onedvpoints[E->mesh.nsd]; d++) 
        for(d1=1; d1<=onedvpoints[E->mesh.nsd]; d1++) {
            i = lnode[d1]*E->lmesh.noz;
            t1 = E->SX[E->mesh.levmax][m][1][i];
            f1 = E->SX[E->mesh.levmax][m][2][i];
            x_surf = x_cmb = 0.0;
            if (E->parallel.me_loc[3]==E->parallel.nprocz-1) 
                x_surf = c_surf * X_surf[m][lnode[d1]] 
                      * E->gDA1[m][el].vpt[d] * E->M.vpt[GMVINDEX(d1,d)] ;
            if (E->parallel.me_loc[3]==0) 
                x_cmb  = c_cmb  *  X_cmb[m][lnode[d1]] 
                      * E->gDA0[m][el].vpt[d] * E->M.vpt[GMVINDEX(d1,d)] ;

            temp[0] += ( x_surf + x_cmb ) * sin(t1)*cos(t1)*cos(f1)  ;
            temp[1] += ( x_surf + x_cmb ) * sin(t1)*cos(t1)*sin(f1)  ;
            //temp[2] += ( x_surf + x_cmb ) * cos(t1) ;         // unused
        }
    }  // end loop over surface elements

    // gather all processor results into final integral, I:
    MPI_Allreduce(temp,I,3,MPI_DOUBLE,MPI_SUM,E->parallel.world);

    // Finally, convert to M:
    for (i=0;i<=2;i++)  M[i] = I[i] * I_to_M_conv ;

    return;
}


/* =============================================================================
 * find_rotation_axis (E,X_surf,X_cmb, M)
 * =============================================================================
 * (This routine is only called by polar_wander_effects.)
 * Calculates two components of the inertia tensor of the earth and converts
 * them to M, the incremental pole position, for purposes of polar wander
 * calculations. This routine is called only by polar_wander_effects(), above.
 *    I_13 = - Integral( (X_s*r_s^4 + X_b*r_b^4) *sin(t)*cos(t)*cos(f)* dO )
 *    I_23 = - Integral( (X_s*r_s^4 + X_b*r_b^4) *sin(t)*cos(t)*sin(f)* dO )
 *  where r_* is the radius at the surface (s) and cmb (b),
 *        X_* is the mass per unit-solid-angle at the surface and cmb,
 *        dO is d-omega (element of solid angle),
 *        t and f are theta and phi.
 * The returned values are unitless. Dimensionalize by mulitplying by M.R^2.
 * This code is largely copied from load_to_CM().
 * Written by Archie, March 2004.
 * =============================================================================
 */

void find_rotation_axis1(E,X_surf,X_cmb, M)
    struct All_variables *E;
    double **X_surf;
    double **X_cmb;
    double M[3];      // space for M1 and M2
{
    static int been_here=0;
    static double c_surf, c_cmb, earth_mass ;
    static double I_to_M_conv ;  // 1/(C-A)   converts inertia to rotation axis
    double omega; // earth's rotation rate
    double temp[3];
    int e,n,m,i,j,k,d,d1,el,nproc;
    int lnode[5];
    double t1,f1,rel_cos;
    double x_surf, x_cmb ;
    double I[3];      // space for moments of inertia I13 and I23
    // variables for parallel integration:

    if (!been_here)  {
        been_here = 1;
        // C and A here are the polar and equatorial moments of intertia.
        // These quantities should give a C-A slightly larger than the
        // hydrostatic limit from rotation effects in order to have a stable
        // rotation axis (or, equivalently, to remove the s=0 mode of the
        // body-tide k love number).  To do this, we take C and A from Lambeck
        // 1980 _Earth's Variable Rotation_ (p27), and then mulitply by CApert
        // (given in the input file). 
        // I_to_M_conv is 1/(C-A), with MKS units.
        //
        //I_to_M_conv = 1.0/ ( ( 0.3306 - 0.3295 )   // C-A
        //                    *( 5.974e24 * 6.37814e6*6.37814e6 ) // M.R^2 
        //                    * E->control.CApert ); // non-hydrostatic effect

        omega = 2.0*M_PI/(24.0*3600.0);  // rotation_rate
        I_to_M_conv = ( 3.0 * E->data.grav_const ) /
                      (omega*omega * pow(E->sphere.dradius,5) * E->ve_data_cont.kf);

        // The fields X_surf and X_cmb are the surface masses, passed in
        // (from calculate_potential) as nondimensional height of boundary.
        // Scale factors convert this nondim height to
        //      mass-per-nondim-area * radius^2
        c_surf  = -E->sphere.dradius ;    // -> dim'ful height 
        c_surf *= E->data.density ;       // -> dim'ful mass/area
        c_surf *= (E->sphere.dradius*E->sphere.dradius) ; // nondim area
        c_surf *= (E->sphere.dradius*E->sphere.dradius) ; // dim'ful radius^2
        // repeat for cmb:
        c_cmb   = -E->sphere.dradius ;  
        c_cmb  *= (E->data.density_below - E->data.density) ;
        c_cmb  *= ( E->sphere.dradius*E->sphere.dradius ) ;
        c_cmb  *= ( E->sphere.dradius*E->sphere.dradius ) * pow(E->sphere.ri,2);
        // Now we can integrate c_surf*X_surf+c_cmb*X_cmb to get dim'ful moment
        // of inertia, and then get nondimensional M with I_to_M_conv.
    }

    /* Calculate the integrals given in the comment above.
     * Note: this technique actually calculates Integral( f r^2 dO), where f is
     * c_surf*X_surf (for example), and r is nondimensinal.
     */
    I[0] = I[1] = I[2] = 0.0;
    temp[0] = temp[1] = temp[2] = 0.0;
    for (m=1; m<=E->sphere.caps_per_proc; m++)
    for (k=1; k<=E->lmesh.ely; k++)
    for (j=1; j<=E->lmesh.elx; j++) {
        el = j + (k-1)*E->lmesh.elx;
        for(d=1;d<=onedvpoints[E->mesh.nsd];d++)
            lnode[d] = E->sien[m][el].node[d];
        for(d =1; d <=onedvpoints[E->mesh.nsd]; d++) 
        for(d1=1; d1<=onedvpoints[E->mesh.nsd]; d1++) {
            i = lnode[d1]*E->lmesh.noz;
            t1 = E->SX[E->mesh.levmax][m][1][i];
            f1 = E->SX[E->mesh.levmax][m][2][i];
            x_surf = x_cmb = 0.0;
            if (E->parallel.me_loc[3]==E->parallel.nprocz-1) 
                x_surf = c_surf * X_surf[m][lnode[d1]] 
                      * E->gDA1[m][el].vpt[d] * E->M.vpt[GMVINDEX(d1,d)] ;
            if (E->parallel.me_loc[3]==0) 
                x_cmb  = c_cmb  *  X_cmb[m][lnode[d1]] 
                      * E->gDA0[m][el].vpt[d] * E->M.vpt[GMVINDEX(d1,d)] ;

            temp[0] += ( x_surf + x_cmb ) * sin(t1)*cos(t1)*cos(f1)  ;
            temp[1] += ( x_surf + x_cmb ) * sin(t1)*cos(t1)*sin(f1)  ;
            //temp[2] += ( x_surf + x_cmb ) * cos(t1) ;         // unused
        }
    }  // end loop over surface elements

    // gather all processor results into final integral, I:
    MPI_Allreduce(temp,I,3,MPI_DOUBLE,MPI_SUM,E->parallel.world);

    // Finally, convert to M:
    for (i=0;i<=2;i++)  M[i] = I[i] * I_to_M_conv ;

    return;
}

/* ================================================== */
void sum_across_surface(E,data,total)
struct All_variables *E;
float *data;
int total;
{
 int j,d;
 float *temp;

 temp = (float *)malloc((total+1)*sizeof(float));
 MPI_Allreduce(data,temp,total,MPI_FLOAT,MPI_SUM,E->parallel.horizontal_comm);

 for (j=0;j<total;j++) {
   data[j] = temp[j];
 }

 free((void *)temp);

 return;
}


double global_dmin(E,a)
    struct All_variables *E;
    double a;
{
    double temp,temp1,a1;
    a1 = a;
    MPI_Allreduce(&a1, &temp1,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    temp = temp1;
    return (temp);
}
