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
#include "citcom_init.h"
#include "interuption.h"
#include "output.h"
#include "parallel_related.h"
#include "checkpoints.h"

extern int Emergency_stop;

void solver_init(struct All_variables *E);

int main(argc,argv)
     int argc;
     char **argv;

{	/* Functions called by main*/
  void general_equations_of_motions_solver();
  void global_default_values();
  void read_instructions();
  void initial_setup();
  void initial_conditions();
  void post_processing();
  void parallel_process_termination();

  void process_new_velocity();
  void pickup_dt();
  void apply_new_loads();
  void get_iceModel();
  void get_static_oceanload();

  void read_velocity_boundary_from_file();
  void read_rayleigh_from_file();
  void read_mat_from_file();
  void read_temperature_boundary_from_file();

  void output_finalize();

  float cpu_time_on_vp_it;

  int cpu_total_seconds,k,need_init_sol;
  double CPU_time0(),time,initial_time,start_time;

  struct All_variables *E;
  MPI_Comm world;

  MPI_Init(&argc,&argv); /* added here to allow command-line input */

  if (argc < 2)   {
    fprintf(stderr,"Usage: %s PARAMETERFILE\n", argv[0]);
    parallel_process_termination();
  }



  /* this section reads input, allocates memory, and set some initial values;
   * replaced by CitcomS.Controller.initialize() and
   * CitcomS.Solver.initialize() in Pyre. */
  world = MPI_COMM_WORLD;
  E = citcom_init(&world); /* allocate global E and do initializaion here */

  /* define common aliases for full/regional functions */
  solver_init(E);

  start_time = time = CPU_time0();

  /* Global interuption handling routine defined once here */
  set_signal();

  /* default values for various parameters */
  global_default_values(E);

  /* read input parameters from file */
  read_instructions(E, argv[1]);

  /* create mesh, setup solvers etc. */
  initial_setup(E);


  cpu_time_on_vp_it = CPU_time0();
  initial_time = cpu_time_on_vp_it - time;
  if (E->parallel.me == 0)  {
    fprintf(stderr,"Input parameters taken from file '%s'\n",argv[1]);
    fprintf(stderr,"Initialization complete after %g seconds\n\n",initial_time);
    fprintf(E->fp,"Initialization complete after %g seconds\n\n",initial_time);
    fflush(E->fp);
  }



  /* this section sets the initial condition;
   * replaced by CitcomS.Controller.launch() ->
   * CitcomS.Solver.launch() in Pyre. */
  if (E->control.post_p) {
      /* the initial condition is from previous checkpoint */
      read_checkpoint(E);

      /* the program will finish after post_processing */
      post_processing(E);
      (E->problem_output)(E, E->monitor.solution_cycles);

      citcom_finalize(E, 0);
  }

  if (E->control.restart) {
      /* the initial condition is from previous checkpoint */
      read_checkpoint(E);
      need_init_sol = 0;
  }
  else {
      /* regular init, or read T from file only */
      initial_conditions(E);
      need_init_sol = 1;	/*  */
  }

    if (E->parallel.me == 0) {
      fprintf(E->fp,"step %d %g %g need_init_sol %d\n",E->monitor.solution_cycles,E->monitor.elapsed_time,E->advection.timestep,need_init_sol);
      fprintf(stderr,"step %d %g %g need_init_sol %d\n",E->monitor.solution_cycles,E->monitor.elapsed_time,E->advection.timestep,need_init_sol);
      fprintf(E->fp,"step %d %g %g %g stage %d %d\n",E->monitor.solution_cycles,E->monitor.elapsed_time,E->advection.timestep,E->advection.next_timestep,E->ve_data_cont.stage,E->ve_data_cont.DIRECT);
      fprintf(stderr,"step %d %g %g %g stage %d %d\n",E->monitor.solution_cycles,E->monitor.elapsed_time,E->advection.timestep,E->advection.next_timestep,E->ve_data_cont.stage,E->ve_data_cont.DIRECT);
      }

  if(need_init_sol){
    /* find first solution */
       general_equations_of_motions_solver(E);
       process_new_velocity(E, E->monitor.solution_cycles);
//printf("after first timestep: after process_new_velocity\n");
  }

  /* stop the computation if only computes elastic problem */
  if (E->control.stokes)  {
    citcom_finalize(E, 0);
  }
  /* information about simulation time and wall clock time */
  output_time(E, E->monitor.solution_cycles);
  if(!E->control.restart)	/* if we have not restarted, print new
				   checkpoint, else leave as is to
				   allow reusing directories */
    output_checkpoint(E);


  while ( E->control.keep_going   &&  (Emergency_stop == 0) ) {

    E->monitor.solution_cycles++;
    E->advection.timesteps = E->monitor.solution_cycles;

    if(E->monitor.solution_cycles>E->control.print_convergence)
      E->control.print_convergence=1;

    if(E->monitor.solution_cycles>E->advection.max_timesteps)
      E->control.keep_going = 0;
    else
      E->control.keep_going = 1;
    pickup_dt(E);
    E->monitor.elapsed_time += E->advection.timestep;
    if (E->parallel.me == 0) {
      fprintf(E->fp,"step %d %g %g %g stage %d %d\n",E->monitor.solution_cycles,E->monitor.elapsed_time,E->advection.timestep,E->advection.next_timestep,E->ve_data_cont.stage,E->ve_data_cont.DIRECT); fflush(E->fp);
      fprintf(stderr,"step %d %g %g %g stage %d %d\n",E->monitor.solution_cycles,E->monitor.elapsed_time,E->advection.timestep,E->advection.next_timestep,E->ve_data_cont.stage,E->ve_data_cont.DIRECT);
      }

    if (E->ve_data_cont.Heaviside==2)   {
            get_iceModel(E,E->slice_ve.iceload[0]);
            if (E->ve_data_cont.SLE) get_static_oceanload(E);
        }
    apply_new_loads(E);
    general_equations_of_motions_solver(E);
    process_new_velocity(E,E->monitor.solution_cycles);

    cpu_total_seconds = CPU_time0()-start_time;
    if (cpu_total_seconds > E->control.record_all_until)  {
      E->control.keep_going = 0;
    }

    /* information about simulation time and wall clock time */
    output_time(E, E->monitor.solution_cycles);

    /* print checkpoint every checkpoint_frequency, unless we have restarted,
       then, we would like to avoid overwriting 
    */
    if ( ((E->monitor.solution_cycles % E->control.checkpoint_frequency)==0) &&
	 ((!E->control.restart) || (E->monitor.solution_cycles != E->monitor.solution_cycles_init))){
	output_checkpoint(E);
    }

    if (E->parallel.me == 0)  {
      fprintf(E->fp,"CPU total = %g & CPU = %g for step %d g_nonl_iteration=%d time = %.4e dt = %.4e  stage %d DIRECT %d\n",CPU_time0()-start_time,CPU_time0()-time,E->monitor.solution_cycles,E->viscosity.iterate,E->monitor.elapsed_time,E->advection.timestep,E->ve_data_cont.stage,E->ve_data_cont.DIRECT);

      time = CPU_time0();
    }

  }


  /* this section prints time accounting;
   * no counterpart in pyre */
  if (E->parallel.me == 0)  {
    fprintf(stderr,"cycles=%d\n",E->monitor.solution_cycles);
    cpu_time_on_vp_it=CPU_time0()-cpu_time_on_vp_it;
    fprintf(stderr,"Average cpu time taken for velocity step = %f\n",
	    cpu_time_on_vp_it/((float)(E->monitor.solution_cycles-E->control.restart)));
    fprintf(E->fp,"Initialization overhead = %f\n",initial_time);
    fprintf(E->fp,"Average cpu time taken for velocity step = %f\n",
	    cpu_time_on_vp_it/((float)(E->monitor.solution_cycles-E->control.restart)));
  }
  citcom_finalize(E, 0);
  return(0);

}
