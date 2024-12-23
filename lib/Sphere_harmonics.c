/* Functions relating to the building and use of mesh locations ... */


#include <math.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"
#include <stdlib.h>

static void compute_sphereh_table(struct All_variables *);

/*   ======================================================================
     ======================================================================  */

void set_sphere_harmonics(E)
     struct All_variables *E;

{
    int m,node,ll,mm,i,j;

    i=0;
    for (ll=0;ll<=E->output.llmax;ll++)
        for (mm=0;mm<=ll;mm++)   {
            E->sphere.hindex[ll][mm] = i;
            i++;
        }

    E->sphere.hindice = i;

    /* spherical harmonic coeff (0=cos, 1=sin)
       for surface topo, cmb topo and geoid */
    for (i=0;i<=1;i++)   {
        E->sphere.harm_geoid[i]=(float*)malloc(E->sphere.hindice*sizeof(float));
        E->sphere.harm_geoid_from_bncy[i]=(float*)malloc(E->sphere.hindice*sizeof(float));
        E->sphere.harm_geoid_from_bncy_botm[i]=(float*)malloc(E->sphere.hindice*sizeof(float));
        E->sphere.harm_geoid_from_tpgt[i]=(float*)malloc(E->sphere.hindice*sizeof(float));
        E->sphere.harm_geoid_from_tpgb[i]=(float*)malloc(E->sphere.hindice*sizeof(float));
        E->sphere.harm_tpgt[i]=(float*)malloc(E->sphere.hindice*sizeof(float));
        E->sphere.harm_tpgb[i]=(float*)malloc(E->sphere.hindice*sizeof(float));
    }

    compute_sphereh_table(E);

    return;
}

/* =====================================
   Generalized Legendre polynomials
   =====================================*/
double modified_plgndr_a(int l, int m, double t)
{
    int i,ll;
    double x,fact1,fact2,fact,pll,pmm,pmmp1,somx2,plgndr;
    const double three=3.0;
    const double two=2.0;
    const double one=1.0;

    x = cos(t);
    pmm=one;
    if(m>0) {
        somx2=sqrt((one-x)*(one+x));
        fact1= three;
        fact2= two;
        for (i=1;i<=m;i++)   {
            fact=sqrt(fact1/fact2);
            pmm = -pmm*fact*somx2;
            fact1+=  two;
            fact2+=  two;
        }
    }

    if (l==m)
        plgndr = pmm;
    else  {
        pmmp1 = x*sqrt(two*m+three)*pmm;
        if(l==m+1)
            plgndr = pmmp1;
        else   {
            for (ll=m+2;ll<=l;ll++)  {
                fact1= sqrt((4.0*ll*ll-one)*(double)(ll-m)/(double)(ll+m));
                fact2= sqrt((2.0*ll+one)*(ll-m)*(ll+m-one)*(ll-m-one)
                            /(double)((two*ll-three)*(ll+m)));
                pll = ( x*fact1*pmmp1-fact2*pmm)/(ll-m);
                pmm = pmmp1;
                pmmp1 = pll;
            }
            plgndr = pll;
        }
    }

    plgndr /= sqrt(4.0*M_PI);

    if (m!=0) plgndr *= sqrt(two);

    return plgndr;
}

void sphere_expansion_output(E,ic,TG,sphc,sphs,ii,filen)
     struct All_variables *E;
     double **TG,*sphc,*sphs;
     int ii,ic;     //ic=1 for surface and ic=0 for CMB
     char * filen;
{
    void sphere_expansion_VE();
    void print_spectral();
    void parallel_process_sync();
    void parallel_process_termination();
    int proc_loc;

    proc_loc = E->parallel.me_loc[3];

    sphere_expansion_VE(E,ic,TG,sphc,sphs,E->output.llmax);

    print_spectral(E,sphc,sphs,proc_loc,ii,filen);

//    parallel_process_sync(E);

   return;
  }


/* =========================================================
   expand the field TG into spherical harmonics
   ========================================================= */
void sphere_expansion(E,TG,sphc,sphs)
     struct All_variables *E;
     float **TG,*sphc,*sphs;
{
    int el,nint,d,p,i,m,j,es,mm,ll,rand();
    void sum_across_surf_sph1();

    for (i=0;i<E->sphere.hindice;i++)    {
        sphc[i] = 0.0;
        sphs[i] = 0.0;
    }

    for (m=1;m<=E->sphere.caps_per_proc;m++)
        for (es=1;es<=E->lmesh.snel;es++)   {

            for (ll=0;ll<=E->output.llmax;ll++)
                for (mm=0; mm<=ll; mm++)   {

                    p = E->sphere.hindex[ll][mm];

                    for(nint=1;nint<=onedvpoints[E->mesh.nsd];nint++)   {
                        for(d=1;d<=onedvpoints[E->mesh.nsd];d++)   {
                            j = E->sien[m][es].node[d];
                            sphc[p] += TG[m][E->sien[m][es].node[d]]
                                * E->sphere.tablesplm[m][j][p]
                                * E->sphere.tablescosf[m][j][mm]
                                * E->M.vpt[GMVINDEX(d,nint)]
                                * E->surf_det[m][nint][es];
                            sphs[p] += TG[m][E->sien[m][es].node[d]]
                                * E->sphere.tablesplm[m][j][p]
                                * E->sphere.tablessinf[m][j][mm]
                                * E->M.vpt[GMVINDEX(d,nint)]
                                * E->surf_det[m][nint][es];
                        }
                    }

                }       /* end for ll and mm  */

        }

    sum_across_surf_sph1(E,sphc,sphs);

    return;
}


void debug_sphere_expansion(struct All_variables *E)
{
    /* expand temperature field (which should be a sph. harm. load)
     * and output the expansion coeff. to stderr
     */
    int m, i, j, k, p, node;
    int ll, mm;
    float *TT[NCS], *sph_harm[2];

    for(m=1;m<=E->sphere.caps_per_proc;m++)
        TT[m] = (float *) malloc ((E->lmesh.nsf+1)*sizeof(float));

    /* sin coeff */
    sph_harm[0] = (float*)malloc(E->sphere.hindice*sizeof(float));
    /* cos coeff */
    sph_harm[1] = (float*)malloc(E->sphere.hindice*sizeof(float));

    for(k=1;k<=E->lmesh.noz;k++)  {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(i=1;i<=E->lmesh.noy;i++)
                for(j=1;j<=E->lmesh.nox;j++)  {
                    node= k + (j-1)*E->lmesh.noz + (i-1)*E->lmesh.nox*E->lmesh.noz;
                    p = j + (i-1)*E->lmesh.nox;
                    TT[m][p] = E->T[m][node];
                }

        /* expand TT into spherical harmonics */
        sphere_expansion(E, TT, sph_harm[0], sph_harm[1]);

        /* only the first nprocz CPU needs output */
        if(E->parallel.me < E->parallel.nprocz) {
            for (ll=0;ll<=E->output.llmax;ll++)
                for (mm=0; mm<=ll; mm++)   {
                    p = E->sphere.hindex[ll][mm];
                    fprintf(stderr, "T expanded layer=%d ll=%d mm=%d -- %12g %12g\n",
                            k+E->lmesh.nzs-1, ll, mm,
                            sph_harm[0][p], sph_harm[1][p]);
                }
        }
    }

    return;
}


/* ==================================================*/
/* ==================================================*/
static void  compute_sphereh_table(struct All_variables *E)
  //   struct All_variables *E;
{
    double modified_plgndr_a();
   void get_sphere_expansion_coff();
   void get_potential_coff();

    int m,node,ll,mm,i,j,p,e,lev;
    double t,f,mmf;
	
	if(1){  // comment by tao: here I allocate it even for compressible case. since it is also used in output(postprocess) subroutines.
		  E->sphere.sphc[0] = (double *)malloc((E->sphere.hindice+3)*sizeof(double));
		  E->sphere.sphc[1] = (double *)malloc((E->sphere.hindice+3)*sizeof(double));
		  E->sphere.sphs[0] = (double *)malloc((E->sphere.hindice+3)*sizeof(double));
		  E->sphere.sphs[1] = (double *)malloc((E->sphere.hindice+3)*sizeof(double));

		  E->sphere.sphc[2] = (double *)malloc((E->sphere.hindice+3)*sizeof(double));
		  E->sphere.sphc[3] = (double *)malloc((E->sphere.hindice+3)*sizeof(double));
		  E->sphere.sphc[4] = (double *)malloc((E->sphere.hindice+3)*sizeof(double));
		  E->sphere.sphc[5] = (double *)malloc((E->sphere.hindice+3)*sizeof(double));
		  E->sphere.sphs[2] = (double *)malloc((E->sphere.hindice+3)*sizeof(double));
		  E->sphere.sphs[3] = (double *)malloc((E->sphere.hindice+3)*sizeof(double));
		  E->sphere.sphs[4] = (double *)malloc((E->sphere.hindice+3)*sizeof(double));
		  E->sphere.sphs[5] = (double *)malloc((E->sphere.hindice+3)*sizeof(double));

          E->sphere.sphc_init_load = (double *)malloc((E->sphere.hindice+3)*sizeof(double));
          E->sphere.sphs_init_load = (double *)malloc((E->sphere.hindice+3)*sizeof(double));
		}
	if(E->ve_data_cont.compressible){
		  E->sphere.sphc_c = (double **)malloc((E->mesh.noz+3)*sizeof(double));
		  E->sphere.sphs_c = (double **)malloc((E->mesh.noz+3)*sizeof(double));
		  E->sphere.sphc_v = (double **)malloc((E->mesh.noz+3)*sizeof(double));
		  E->sphere.sphs_v = (double **)malloc((E->mesh.noz+3)*sizeof(double));
		  for (i=0;i<=E->mesh.noz;i++) {  // comment by tao: add by tao
		    E->sphere.sphc_c[i] = (double *)malloc((E->sphere.hindice+3)*sizeof(double));
		    E->sphere.sphs_c[i] = (double *)malloc((E->sphere.hindice+3)*sizeof(double));
		    E->sphere.sphc_v[i] = (double *)malloc((E->sphere.hindice+3)*sizeof(double));
		    E->sphere.sphs_v[i] = (double *)malloc((E->sphere.hindice+3)*sizeof(double));
 		   }
		}
	
  for (m=1;m<=E->sphere.caps_per_proc;m++)   {
    for (mm=0;mm<=E->output.llmax;mm++)  {
      E->Sph_Harm_Tblcs[m][mm] = (struct SPH1 *) malloc((E->lmesh.snel+2)*sizeof(struct SPH1));
      E->Tbl_cs[m][mm] = (float *) malloc((E->lmesh.nsf+2)*sizeof(float));
      E->Tbl_sn[m][mm] = (float *) malloc((E->lmesh.nsf+2)*sizeof(float));
      }
    for (p=0;p<E->sphere.hindice;p++)    {
      E->Sph_Harm_Tbllm[m][p] = (struct SPH2 *) malloc((E->lmesh.snel+2)*sizeof(struct SPH2));
      E->Tbl_lm[m][p] = (float *) malloc((E->lmesh.nsf+2)*sizeof(float));
      }
    }

 lev = E->mesh.levmax;
 for(m=1;m<=E->sphere.caps_per_proc;m++)  {
    for (e=1;e<=E->lmesh.snel;e++)   {
       get_sphere_expansion_coff(E,e,lev,m);
       }
    get_potential_coff(E,lev,m);
    }
 

 for(m=1;m<=E->sphere.caps_per_proc;m++)  {
        E->sphere.tablesplm[m]   = (double **) malloc((E->lmesh.nsf+1)*sizeof(double*));
        E->sphere.tablescosf[m] = (double **) malloc((E->lmesh.nsf+1)*sizeof(double*));
        E->sphere.tablessinf[m] = (double **) malloc((E->lmesh.nsf+1)*sizeof(double*));

        for (i=1;i<=E->lmesh.nsf;i++)   {
            E->sphere.tablesplm[m][i]= (double *)malloc((E->sphere.hindice)*sizeof(double));
            E->sphere.tablescosf[m][i]= (double *)malloc((E->output.llmax+1)*sizeof(double));
            E->sphere.tablessinf[m][i]= (double *)malloc((E->output.llmax+1)*sizeof(double));
        }
    }

  for(m=1;m<=E->sphere.caps_per_proc;m++)  {
        for (j=1;j<=E->lmesh.nsf;j++)  {
            node = j*E->lmesh.noz;
            f=E->sx[m][2][node];
            t=E->sx[m][1][node];
            for (mm=0;mm<=E->output.llmax;mm++)   {
	      mmf = (double)(mm)*f;
                E->sphere.tablescosf[m][j][mm] = cos( mmf );
                E->sphere.tablessinf[m][j][mm] = sin( mmf );
            }

            for (ll=0;ll<=E->output.llmax;ll++)
                for (mm=0;mm<=ll;mm++)  {
                    p = E->sphere.hindex[ll][mm];
                    E->sphere.tablesplm[m][j][p] = modified_plgndr_a(ll,mm,t) ;
                }
        }
    }

    return;
}

/* =========================================================
  ========================================================= */
 void sphere_expansion_VE(E,ic,TG,sphc,sphs,llmax)
 struct All_variables *E;
 double **TG,*sphc,*sphs;
 int ic,llmax;
 {
  void parallel_process_termination();
 void sum_across_surf_sph2();
 int mm,ll;
 double modified_plgndr_a(),t[5],f[5],temp,sphere_h();

 int hindice,j1,el,a,b,i,j,k1,p,q,nint,d,k,n,es;
 double area,temp0[5],temp1[5],temp2[5],temp3[5];

  void get_global_1d_shape_fn_2();

  struct Shape_function1 GM;
  struct Shape_function1_dA dGammax;

    const int lev=E->mesh.levmax;
    const int dims=E->mesh.nsd,dofs=E->mesh.dof;

   hindice = (llmax+2)*(llmax+1)/2;

   for (i=0;i<E->sphere.hindice;i++)    {
      sphc[i] = 0.0;
      sphs[i] = 0.0;
      }

   j1 = (ic==1)?4:0;
   area = 0.0;

   for (j=1;j<=E->sphere.caps_per_proc;j++)
     for (es=1;es<=E->lmesh.snel;es++) {

       el = (ic==1)?(E->lmesh.elz*es):(E->lmesh.elz*(es-1)+1);

/*       get_global_1d_shape_fn_2(E,el,&GM,&dGammax,ic,lev,j);

       for(k1=1;k1<=onedvpoints[E->mesh.nsd];k1++)   {
         t[k1] = E->sx[1][j][E->ien[j][el].node[k1+j1]];
         f[k1] = E->sx[2][j][E->ien[j][el].node[k1+j1]];
         }
*/

       for(k=1;k<=onedvpoints[E->mesh.nsd];k++)   {
         temp0[k] = 0.0;
         for(k1=1;k1<=onedvpoints[E->mesh.nsd];k1++)   
           temp0[k] += TG[j][E->sien[j][es].node[k1]] * E->M.vpt[GMVINDEX(k1,k)];
         }

       if (ic==1)  {
         area += E->gDA1[j][es].vpt[5];
         for(k=1;k<=onedvpoints[E->mesh.nsd];k++)
           temp0[k] = temp0[k]*E->gDA1[j][es].vpt[k];
         }
       else if (ic==0)  {
         area += E->gDA0[j][es].vpt[5];
         for(k=1;k<=onedvpoints[E->mesh.nsd];k++)
           temp0[k] = temp0[k]*E->gDA0[j][es].vpt[k];
         }

       for (ll=0;ll<=llmax;ll++)
         for (mm=0; mm<=ll; mm++)   {
           p = E->sphere.hindex[ll][mm];

           for(k=1;k<=onedvpoints[E->mesh.nsd];k++)   {
             temp1[k] = temp2[k] = temp3[k] = 0.0;
/*             for(k1=1;k1<=onedvpoints[E->mesh.nsd];k1++)   {
                temp1[k] += cos(mm*f[k1]) * E->M.vpt[GMVINDEX(k1,k)];
                temp2[k] += sin(mm*f[k1]) * E->M.vpt[GMVINDEX(k1,k)];
                temp3[k] += modified_plgndr_a(ll,mm,t[k1]) * E->M.vpt[GMVINDEX(k1,k)];
                }
*/

             sphc[p]+=temp0[k]*E->Sph_Harm_Tblcs[j][mm][es].cs[k]*E->Sph_Harm_Tbllm[j][p][es].lm[k];
             sphs[p]+=temp0[k]*E->Sph_Harm_Tblcs[j][mm][es].sn[k]*E->Sph_Harm_Tbllm[j][p][es].lm[k];

             }
           }   /* end for ll and mm*/

       }  


   sphs[E->sphere.hindice] = area;

   sum_across_surf_sph2(E,sphc,sphs);

   area = sphs[E->sphere.hindice];

   area = 4.0*M_PI/area;

   for (i=0;i<E->sphere.hindice;i++)    {
      sphc[i] = sphc[i]*area;
      sphs[i] = sphs[i]*area;
      }

 return;
 }

