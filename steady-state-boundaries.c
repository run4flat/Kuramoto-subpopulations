/*
 * steady-state-slips.c - fast simulation to identify slip time series
 *   and identify the steady state use MSER.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "stretchy_buffer.h"
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
	if (argc != 3) {
		printf("You must give two arguments\n");
		exit(1);
	}
	
	/* Make sure the working directory is in a clean state */
//	char command[160];
//	sprintf(command, "perl " CODE_DIRECTORY "add_slips_provenance.pl %s %s", argv[1], argv[2]);
//	if (system(command) != 0) exit(1);
	
	/* Args are input population file, value of K, and number of time steps. */
	double K = atof(argv[1]);
	int N_time_steps = atoi(argv[2]);
	int N_desired_steady_state_steps = N_time_steps;
	printf("Simulating at K=%f with initial run of %d time steps\n", K, N_time_steps);
	
	/* Load omegas and thetas from file */
	double * omegas = load_data("population");
	double * thetas = load_data("positions");
	int N = sb_count(omegas);
	if (sb_count(thetas) != N) {
		printf("Omegas and thetas have different numbers of oscillators!\n");
		exit(1);
	}
	
	/* Build theta differences */
	double * theta_diffs = malloc(sizeof(double) * (N-1));
	int * is_boundary = malloc(sizeof(int) * (N-1));
	int i;
	for (i = 1; i < N; i++) {
		theta_diffs[i-1] = thetas[i] - thetas[i-1] - M_PI_2;
		if (theta_diffs[i-1] < -M_PI) theta_diffs[i-1] += 2 * M_PI;
		if (theta_diffs[i-1] > M_PI) theta_diffs[i-1] -= 2 * M_PI;
		
		is_boundary[i-1] = 0;
	}
	
	/* Run until we have the requested number of steady-state time steps */
	slip_event * slips = 0;
	int N_steady_state_steps = 0;
	int total_steps_so_far = 0;
	int steady_state_slip_offset;
	while(N_steady_state_steps < N_desired_steady_state_steps) {
		int N_steps_this_round = N_desired_steady_state_steps - N_steady_state_steps;
		
		/* Run the expected number of steps to steady state */
		slips = run_sim(omegas, thetas, theta_diffs, N_steps_this_round,
			total_steps_so_far, K, slips);
		total_steps_so_far += N_steps_this_round;
		
		/* mser */
		steady_state_slip_offset = mser_first_steady_slip(slips);
		N_steady_state_steps
			= total_steps_so_far - slips[steady_state_slip_offset].time_step;
	}
	
	/* Keep running until we have the requested number of steps without any
	 * new boundaries added. Start by noting the boundaries identified since we
	 * started steady state, and the time at which they were added. */
	int t_change = time_step_of_last_boundary_change(steady_state_slip_offset,
		slips, is_boundary, 0);
	
	int round;
	for (round = 2; round < N_steady_state_steps / 2; round++) {
		int boundary_check_start = sb_count(slips);
		
		/* As it gets longer, we become more tolerant. */
		int N_steps_this_round = N_desired_steady_state_steps / round;
		if (total_steps_so_far - t_change > N_steps_this_round) break;
		
		/* Run some more */
		printf("Boundary round %d\n", round);
		slips = run_sim(omegas, thetas, theta_diffs, N_steps_this_round,
			total_steps_so_far, K, slips);
		
		/* Note the total time and obtain the time of the last change */
		total_steps_so_far += N_steps_this_round;
		t_change = time_step_of_last_boundary_change(boundary_check_start, slips, is_boundary, t_change);
	}
	
	/* when done, save the transients and steady-state slips separately */
	char out_filename[160];
	sprintf(out_filename, "transient-slips,K=%f.txt", K);
	FILE * out_fh = fopen(out_filename, "w");
	int prev_cumu = 0;
	for (i = 0; i < steady_state_slip_offset; i++) {
		fprintf(out_fh, "%d %d %d\n", slips[i].time_step, slips[i].boundary,
			slips[i].cumulative_slips - prev_cumu);
		prev_cumu = slips[i].cumulative_slips;
	}
	fclose(out_fh);
	
	/* save the steady-state values, too */
	sprintf(out_filename, "steady-state-slips,K=%f.txt", K);
	out_fh = fopen(out_filename, "w");
	for (; i < sb_count(slips); i++) {
		fprintf(out_fh, "%d %d %d\n", slips[i].time_step, slips[i].boundary,
			slips[i].cumulative_slips - prev_cumu);
		prev_cumu = slips[i].cumulative_slips;
	}
	fclose(out_fh);
	
	/* save the coherent regions */
	sprintf(out_filename, "boundaries,K=%f.txt", K);
	out_fh = fopen(out_filename, "w");
	int left, right;
	for (left = 0; left < N-1; left++) {
		if (is_boundary[left] == 1) continue;
		for (right = left + 1; right < N; right++) {
			if (is_boundary[right-1]) break;
		}
		fprintf(out_fh, "%d-%d\n", left, right-1);
		left = right - 1; /* less one in prep for the for's increment */
	}
	fclose(out_fh);
	
	/* Save some information about this simulation */
	sprintf(out_filename, "slip-details,K=%f.txt", K);
	out_fh = fopen(out_filename, "w");
	fprintf(out_fh, "Total time steps: %d\n", total_steps_so_far);
	fprintf(out_fh, "Transient time steps: %d\n",
		slips[steady_state_slip_offset-1].time_step);
	fprintf(out_fh, "total slips: %d\n", sb_count(slips));
	fprintf(out_fh, "Transient slips: %d\n", steady_state_slip_offset-1);
	fclose(out_fh);
	
	/* Clean up */
	sb_free(omegas);
	sb_free(thetas);
	free(theta_diffs);
	sb_free(slips);
	free(is_boundary);
	return (0);
}

/* Given a stem (like "population") loads raw data from "population.bin"
 * if it exists, or text data from "population.txt" if it exists, or
 * prints a warning and returns null. */
double * load_data (char * filename_stem) {
	char full_filename[160];
	double * data = 0;
	FILE * in_fh;
	
	/* Try the binary file */
	sprintf(full_filename, "%s.bin", filename_stem);
	if (in_fh = fopen(full_filename, "r")) {
		double datum;
		int bytes_read;
		while(1) {
			/* Read on number at a time */
			bytes_read = fread(&datum, sizeof(double), 1, in_fh);
			/* Exit when there's nothing else to read */
			if (bytes_read < sizeof(double)) break;
			/* Add the data we read and advance to the next */
			sb_push(data, datum);
		}
		fclose(in_fh);
		return data;
	}
	
	/* Try the text file */
	sprintf(full_filename, "%s.txt", filename_stem);
	if (in_fh = fopen(full_filename, "r")) {
		double datum;
		while(!feof(in_fh)) {
			if(fscanf(in_fh, "%lf", &datum) == 1) sb_push(data, datum);
		}
		fclose(in_fh);
		return data;
	}
	
	/* Print a warning about trouble */
	printf("Unable to open file [%s]\n", filename_stem);
	return 0;
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
void compute_vel(double * om, double * th, double K, double * prev_vel, double dt, double * out_vel) {
	int N_osc = sb_count(om);
	int i;
	double r, psi;
	double Rx = 0;
	double Ry = 0;
	
	/* Temporarily abuse out_vel to hold the positions of the quasi step */
	if (prev_vel) {
		for (i = 0; i < N_osc; i++) out_vel[i] = th[i] + dt * prev_vel[i];
	}
	else {
		for (i = 0; i < N_osc; i++) out_vel[i] = th[i];
	}
	
	/* Calculate r and psi for the quasi-step */
	for (i = 0; i < N_osc; i++) {
		Rx += cos(out_vel[i]);
		Ry += sin(out_vel[i]);
	}
	r = sqrt(Rx*Rx + Ry*Ry) / (double)N_osc;
	psi = atan2(Ry, Rx);
	
	/* Update the contents of out_vel to be the velocity at the quasi-step */
	for (i = 0; i < N_osc; i++) out_vel[i] = om[i] + r * K * sin(psi - out_vel[i]);
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
