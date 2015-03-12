/*
 * multi-steady-state-boundaries.c - fast simulation to identify slip time
 *   series and identify the steady state use MSER. Simultaneously runs
 *   multiple simulations to determine the true transient time.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "stretchy_buffer.h"
#include "load_data.h"
#include "KISSrand.h"
#define CODE_DIRECTORY "/home/dcmertens/projects/2015/03-06-Kuramoto-Slips/"

typedef struct {
	int time_step;
	int boundary;
	int cumulative_slips;
} slip_event;

/* Open the filename and read one double at a time. Returns a double stretchybuffer */
double * load_data (char * filename_stem);
/* Run the time stepper */
slip_event * run_sim(double * om, double * th, double * th_diff, int N_time_steps,
	int starting_time_step, double K, slip_event * slips
);
/* identifies the slip index (not time step) of the first non-transient slip */
int mser_first_steady_slip (slip_event * slips);
/* Identifies the most recent slip event that actually changes the status of a boundary */
int time_step_of_last_boundary_change(int start, slip_event * slips, int * is_boundary, int t_change);

int main(int argc, char **argv) {
	if (argc != 4) {
		printf("You must give three arguments\n");
		exit(1);
	}
	
	/* Make sure the working directory is in a clean state */
	char command[160];
	sprintf(command, "perl " CODE_DIRECTORY "add_boundaries_provenance.pl %s %s", argv[1], argv[2]);
//	if (system(command) != 0) exit(1);
	
	/* Args are value of K, number of time steps, number of simultaneous simulations. */
	double K = atof(argv[1]);
	int N_time_steps = atoi(argv[2]);
	int N_sims = atoi(argv[3]);
	int N_desired_steady_state_steps = N_time_steps;
	printf("Simultaneously running %d simulations at K=%f with initial run of %d time steps\n",
		N_sims, K, N_time_steps);
	
	/* Load omegas from file */
	double * omegas = load_data("population");
	int N_osc = sb_count(omegas);
	
	/* Create the random initial positions, allocate theta differences,
	 * allocate slip and boundary tracking */
	init_KISS();
	double ** thetas = malloc(sizeof(double*) * N_sims);
	double ** theta_diffs = malloc(sizeof(double*) * N_sims);
	int ** is_boundary = malloc(sizeof(int*) * N_sims);
	slip_event ** slips = malloc(sizeof(slip_event*) * N_sims);
	int sim, i;
	for (sim = 0; sim < N_sims; sim++) {
		/* Allocate memory */
		thetas[sim] = malloc(sizeof(double) * N_osc);
		theta_diffs[sim] = malloc(sizeof(double) * (N_osc-1));
		is_boundary[sim] = malloc(sizeof(int) * (N_osc-1));
		slips[sim] = 0;
		
		/* Initialize random thetas, compute dthetas, set zero boundary counts */
		thetas[sim][0] = 2 * M_PI * (dKISS() - 0.5);
		for (i = 1; i < N_osc; i++) {
			thetas[sim][i] = 2 * M_PI * (dKISS() - 0.5);
			
			theta_diffs[sim][i-1] = thetas[sim][i] - thetas[sim][i-1] - M_PI_2;
			if (theta_diffs[sim][i-1] < -M_PI) theta_diffs[sim][i-1] += 2 * M_PI;
			if (theta_diffs[sim][i-1] > M_PI) theta_diffs[sim][i-1] -= 2 * M_PI;
			
			is_boundary[sim][i-1] = 0;
		}
	}
	/* these need to be allocated, but do not need to be initialized */
	int * steady_state_slip_offsets = malloc(sizeof(int) * N_sims);
	
	/* Run until we have the requested number of steady-state time steps */
	int N_steady_state_steps = 0;
	int total_steps_so_far = 0;
	int max_steady_state_time_step;
	while(N_steady_state_steps < N_desired_steady_state_steps) {
		int N_steps_this_round = N_desired_steady_state_steps - N_steady_state_steps;
		
		/* run each simulation for the given number of steps; use MSER to get steady-state
		 * time step */
		max_steady_state_time_step = 0;
		for (sim = 0; sim < N_sims; sim++) {
			/* time steps */
			slips[sim] = run_sim(omegas, thetas[sim], theta_diffs[sim],
				N_steps_this_round, total_steps_so_far, K, slips[sim]);
			/* mser */
			int slip_ss_offset = mser_first_steady_slip(slips[sim]);
			steady_state_slip_offsets[sim] = slip_ss_offset;
			if (slips[sim][slip_ss_offset].time_step > max_steady_state_time_step) {
				max_steady_state_time_step = slips[sim][slip_ss_offset].time_step;
			}
		}
		
		/* Always double the steady state steps, just to be absolutely sure */
		max_steady_state_time_step *= 2;
		
		/* Track how long we've simulated */
		total_steps_so_far += N_steps_this_round;
		N_steady_state_steps = total_steps_so_far - max_steady_state_time_step;
	}
	
	/* Keep running until we have the requested number of steps without any
	 * new boundaries added. Start by noting the boundaries identified since we
	 * started steady state, and the time at which they were added. */
	int * curr_ss_slip_offset = malloc(sizeof(int) * N_sims);
	int * init_ss_slip_offset = malloc(sizeof(int) * N_sims);
	int * t_change = malloc(sizeof(int) * N_sims);
	int max_t_change = 0;
	for (sim = 0; sim < N_sims; sim++) {
		
		/* Find first slip event after steady state time step */
		init_ss_slip_offset[sim] = curr_ss_slip_offset[sim]
			= get_slip_offset_after_t_step(slips[sim], max_steady_state_time_step);
		
		/* Examine steady state flips; calculate time of last change */
		t_change[sim] = time_step_of_last_boundary_change(curr_ss_slip_offset[sim],
			slips[sim], is_boundary[sim], 0);
		if (t_change[sim] > max_t_change) max_t_change = t_change[sim];
	}
	
	int round;
	for (round = 2; round < N_steady_state_steps / 2; round++) {
		/* As it gets longer, we become more tolerant. */
		int N_steps_this_round = N_desired_steady_state_steps / round;
		if (total_steps_so_far - max_t_change > N_steps_this_round) break;
		
		/* Run some more */
		printf("Boundary round %d\n", round);
		for (sim = 0; sim < N_sims; sim++) {
			/* Update our current steady state slip offset so that we
			 * only examine t_changes for new slip data */
			curr_ss_slip_offset[sim] = sb_count(slips[sim]);
			
			/* run it */
			slips[sim] = run_sim(omegas, thetas[sim], theta_diffs[sim],
				N_steps_this_round, total_steps_so_far, K, slips[sim]);
			
			/* Look for time of latest change */
			t_change[sim] = time_step_of_last_boundary_change(curr_ss_slip_offset[sim],
				slips[sim], is_boundary[sim], t_change[sim]);
			
			/* Track max time */
			if (t_change[sim] > max_t_change) max_t_change = t_change[sim];
		}
		
		/* Note the total time and obtain the time of the last change */
		total_steps_so_far += N_steps_this_round;
	}
	
	/* All done; save stuff */
	char out_filename[160];
	FILE * out_fh;
	for (sim = 0; sim < N_sims; sim++) {
		/* when done, save the transients and steady-state slips separately */
		sprintf(out_filename, "transient-slips,sim%d,K=%f.txt", sim, K);
		out_fh = fopen(out_filename, "w");
		int prev_cumu = 0;
		for (i = 0; i < init_ss_slip_offset[sim]; i++) {
			fprintf(out_fh, "%d %d %d\n", slips[sim][i].time_step, slips[sim][i].boundary,
				slips[sim][i].cumulative_slips - prev_cumu);
			prev_cumu = slips[sim][i].cumulative_slips;
		}
		fclose(out_fh);
	
		/* save the steady-state values, too */
		sprintf(out_filename, "steady-state-slips,sim%d,K=%f.txt", sim, K);
		out_fh = fopen(out_filename, "w");
		for (; i < sb_count(slips[sim]); i++) {
			fprintf(out_fh, "%d %d %d\n", slips[sim][i].time_step, slips[sim][i].boundary,
				slips[sim][i].cumulative_slips - prev_cumu);
			prev_cumu = slips[sim][i].cumulative_slips;
		}
		fclose(out_fh);
		
		/* save the per-simulation coherent regions */
		sprintf(out_filename, "boundaries,sim%d,K=%f.txt", sim, K);
		out_fh = fopen(out_filename, "w");
		int left, right;
		for (left = 0; left < N_osc-1; left++) {
			if (is_boundary[sim][left] == 1) continue;
			for (right = left + 1; right < N_osc; right++) {
				if (is_boundary[sim][right-1]) break;
			}
			fprintf(out_fh, "%d-%d\n", left, right-1);
			left = right - 1; /* less one in prep for the for's increment */
		}
		fclose(out_fh);
	}
	
	/* save the cross-simulation coherent regions */
	sprintf(out_filename, "boundaries,K=%f.txt", K);
	out_fh = fopen(out_filename, "w");
	int left, right, left_is_good;
	for (left = 0; left < N_osc-1; left++) {
		/* See if the left is good across simulations */
		left_is_good = 1;
		for (sim = 0; sim < N_sims; sim++) {
			if (is_boundary[sim][left] == 1) left_is_good = 0;
		}
		if (!left_is_good) continue;
		
		/* find the furthest right that is good across simulations */
		for (right = left + 1; right < N_osc; right++) {
			for (sim = 0; sim < N_sims; sim++) {
				if (is_boundary[sim][right-1]) goto PRINT_RESULTS;
			}
		}
		PRINT_RESULTS: fprintf(out_fh, "%d-%d\n", left, right-1);
		left = right - 1; /* less one in prep for the for's increment */
	}
	fclose(out_fh);
	
	
	/* Save some information about this simulation */
	sprintf(out_filename, "slip-details,K=%f.txt", K);
	out_fh = fopen(out_filename, "w");
	
	fprintf(out_fh, "Total time steps: %d\n", total_steps_so_far);
	
	fprintf(out_fh, "Per-simulation transient time steps:");
	for (sim = 0; sim < N_sims; sim++) {
		fprintf(out_fh, " %d", slips[sim][steady_state_slip_offsets[sim]].time_step);
	}
	fprintf(out_fh, "\n");
	
	fprintf(out_fh, "Per-simulation transient slips:");
	for (sim = 0; sim < N_sims; sim++) {
		fprintf(out_fh, " %d", init_ss_slip_offset[sim] + 1);
	}
	fprintf(out_fh, "\n");
	
	fprintf(out_fh, "Per-simulation steady-state slips:");
	for (sim = 0; sim < N_sims; sim++) {
		fprintf(out_fh, " %d", sb_count(slips[sim]) - init_ss_slip_offset[sim]);
	}
	fprintf(out_fh, "\n");
	
	
	fclose(out_fh);
	
	/* Clean up */
	sb_free(omegas);
	for (sim = 0; sim < N_sims; sim++) {
		free(thetas[sim]);
		free(theta_diffs[sim]);
		free(is_boundary[sim]);
		sb_free(slips[sim]);
	}
	free(thetas);
	free(theta_diffs);
	free(is_boundary);
	free(slips);
	free(steady_state_slip_offsets);
	free(init_ss_slip_offset);
	free(curr_ss_slip_offset);
	free(t_change);
	return (0);
}

/* Compute the velocity for the given positions, or an Euler step from them
 * using the given velocity. */
void compute_vel(double * om, double * th, double K, double * prev_vel, double dt, double * out_vel);

slip_event * run_sim(double * om, double * th, double * th_diff, int N_time_steps,
	int starting_time_step, double K, slip_event * slips
) {
	int N_osc = sb_count(om);
	int update_tick = 1000000 / N_osc;
	int slip_dir;
	
	/* Will track the cumulative number of slips. Start at zero, or where
	 * we left off. */
	int cumulative_slips = 0;
	if (slips) cumulative_slips = slips[sb_count(slips)-1].cumulative_slips;
	
	/* A simulation constant */
	double dt = 0.01;
	
	/* Allocate a couple of working vectors */
	double * v1 = malloc(sizeof(double) * N_osc);
	double * v2 = malloc(sizeof(double) * N_osc);
	double * v3 = malloc(sizeof(double) * N_osc);
	double * v4 = malloc(sizeof(double) * N_osc);
	double * th_vel = malloc(sizeof(double) * N_osc);
	
	int i, j;
	for (i = 0; i < N_time_steps; i++) {
		if (i % update_tick == 0) {
			printf("\rStep %d/%d",
				i + starting_time_step, N_time_steps + starting_time_step);
			fflush(stdout);
		}
		
		/* Compute the velocity using rk4 */
		compute_vel(om, th, K, NULL, dt, v1);
		compute_vel(om, th, K, v1, dt/2.0, v2);
		compute_vel(om, th, K, v2, dt/2.0, v3);
		compute_vel(om, th, K, v3, dt, v4);
		
		/* Compute the velocities, update positions, and keep them within -pi
		 * and pi */
		for (j = 0; j < N_osc; j++) {
			th_vel[j] = (v1[j] + 2.0 * v2[j] + 2.0 * v3[j] + v4[j]) / 6.0;
			th[j] += dt * th_vel[j];
			if (th[j] > M_PI) th[j] -= 2.0 * M_PI;
			else if (th[j] < -M_PI) th[j] += 2.0 * M_PI;
		}
		
		/* Update relative phase positions, identifying 2pi slips in the process */
		for (j = 0; j < N_osc - 1; j++) {
			th_diff[j] += dt * (th_vel[j + 1] - th_vel[j]);
			slip_dir = 0;
			if (th_diff[j] > M_PI) slip_dir = 1;
			else if (th_diff[j] < -M_PI) slip_dir = -1;
			if (slip_dir != 0) {
				/* track the slip */
				cumulative_slips += slip_dir;
				slip_event new_event;
				new_event.time_step = i + starting_time_step;
				new_event.boundary = j;
				new_event.cumulative_slips = cumulative_slips;
				sb_push(slips, new_event);
				
				/* keep the difference within 2pi */
				th_diff[j] -= slip_dir * 2.0 * M_PI;
			}
		}
		
		/* Fix relative phase positions every 1000 steps so they don't get
		 * out of sync (ha!) with the real phase positions. */
		if (i % 1000 == 0) {
			double curr_diff;
			for (j = 0; j < N_osc - 1; j++) {
				curr_diff = th[j+1] - th[j] - M_PI_2;
				if (curr_diff > M_PI) curr_diff -= 2.0 * M_PI;
				else if (curr_diff < -M_PI) curr_diff += 2.0 * M_PI;
				/* If there is a big difference, it's due to 2pi wrapping.
				 * Figure out which way we're going to unwind and print a slip
				 * for the process. */
				slip_dir = 0;
				if (curr_diff - th_diff[j] < -6) slip_dir = 1;
				else if (curr_diff - th_diff[j] > 6) slip_dir = -1;
				if (slip_dir != 0) {
					cumulative_slips += slip_dir;
					slip_event new_event;
					new_event.time_step = i + starting_time_step;
					new_event.boundary = j;
					new_event.cumulative_slips = cumulative_slips;
					sb_push(slips, new_event);
				}
				
				th_diff[j] = curr_diff;
			}
		}
	}
	printf("\n");
	
	free(th_vel);
	free(v1);
	free(v2);
	free(v3);
	free(v4);
	
	return slips;
}

/* Compute the velocity for the given positions, or an Euler step from them
 * using the given velocity. */
void compute_vel(double * om, double * init_th, double K, double * prev_vel, double dt, double * out_vel) {
	int N_osc = sb_count(om);
	int i;
	double r, psi;
	double Rx = 0;
	double Ry = 0;
	double * curr_th;
	
	if (prev_vel) {
		/* Temporarily (ab)use out_vel's allocation to hold the positions of the quasi step */
		curr_th = out_vel;
		for (i = 0; i < N_osc; i++) curr_th[i] = init_th[i] + dt * prev_vel[i];
	}
	else {
		/* use the initial positions */
		curr_th = init_th;
	}
	
	/* Calculate r and psi for the quasi-step */
	for (i = 0; i < N_osc; i++) {
		Rx += cos(curr_th[i]);
		Ry += sin(curr_th[i]);
	}
	r = sqrt(Rx*Rx + Ry*Ry) / (double)N_osc;
	psi = atan2(Ry, Rx);
	
	/* Update the contents of out_vel to be the velocity at the quasi-step */
	for (i = 0; i < N_osc; i++) out_vel[i] = om[i] + r * K * sin(psi - curr_th[i]);
}

int mser_first_steady_slip (slip_event * slips) {
	/* Perform the MSER analysis, assuming a steady state of a constant
	 * slope (as opposed to a constant value) */
	int i, j, N_data;
	double N = 0;
	double S_x = 0;
	double S_y = 0;
	double S_xx = 0;
	double S_xy = 0;
	double S_yy = 0;
	double x, y;
	
	/* Accumulate information for the first (last) ten points; do not
	 * calculate MSER scores for them, however */
	int N_slips = sb_count(slips);
	for (i = N_slips - 1; i > N_slips - 11; i--) {
		x = slips[i].time_step * 0.01;
		y = slips[i].cumulative_slips;
		N++;
		S_x += x;
		S_y += y;
		S_xx += x*x;
		S_xy += x*y;
		S_yy += y*y;
	}
	double slope = (S_xy * N - S_x * S_y) / (S_xx * N - S_x * S_x);
	double y0 = (S_xy - slope * S_xx) / S_x;
	
	/* keep track of the best score and its index */
	double best_MSER_score 
		= (N * y0*y0 + 2 * y0 * slope * S_x
		- 2 * y0 * S_y + slope*slope * S_xx
		 - 2 * slope * S_xy + S_yy) / N/N;
	int best_slip = i;
	
	/* Now compute the MSER scores for the remaining points */
	for (; i >= 0; i--) {
		x = slips[i].time_step * 0.01;
		y = slips[i].cumulative_slips;
		N++;
		S_x += x;
		S_y += y;
		S_xx += x*x;
		S_xy += x*y;
		S_yy += y*y;
		
		slope = (S_xy * N - S_x * S_y) / (S_xx * N - S_x * S_x);
		y0 = (S_xy - slope * S_xx) / S_x;
		
		double curr_MSER_score
			= (N * y0*y0 + 2 * y0 * slope * S_x
			- 2 * y0 * S_y + slope*slope * S_xx
			 - 2 * slope * S_xy + S_yy) / N/N;
		if (curr_MSER_score < best_MSER_score) {
			best_MSER_score = curr_MSER_score;
			best_slip = i;
		}
	}
	return best_slip;
}

/* Perform a binary search to find first slip event on or after the given
 * time step. */
int get_slip_offset_after_t_step(slip_event * slips, int min_t_step) {
	int left = 0;
	int right = sb_count(slips) - 1;
	while(right - left > 1) {
		int mid = (left + right) / 2;
		int mid_t_step = slips[mid].time_step;
		if (mid_t_step == min_t_step) return mid;
		if (mid_t_step < min_t_step) left = mid;
		else right = mid;
	}
	return right;
}

int time_step_of_last_boundary_change(int start, slip_event * slips, int * is_boundary, int t_change) {
	int curr_boundary;
	int i;
	for (i = start; i < sb_count(slips); i++) {
		curr_boundary = slips[i].boundary;
		/* track the time at which any boundary changes from zero... */
		if (is_boundary[curr_boundary] == 0) t_change = slips[i].time_step;
		is_boundary[curr_boundary] = 1;
	}
	return t_change;
}
