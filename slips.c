/*
 * slips.c - fast simulation to identify slip time series.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "stretchy_buffer.h"
#include "load_data.h"
#define CODE_DIRECTORY "/home/dcmertens/projects/2015/03-06-Kuramoto-Slips/"

/* Open the filename and read one double at a time. Returns a double stretchybuffer */
double * load_data (char * filename_stem);
/* Run the time stepper */
void run_sim(double * om, double * th, double * th_diff, int N_time_steps, double K);

int main(int argc, char **argv) {
	if (argc != 3) {
		printf("You must give two arguments\n");
		exit(1);
	}
	
	/* Make sure the working directory is in a clean state */
	char command[160];
	sprintf(command, "perl " CODE_DIRECTORY "add_slips_provenance.pl %s %s", argv[1], argv[2]);
	if (system(command) != 0) exit(1);
	
	/* Args are input population file, value of K, and number of time steps. */
	double K = atof(argv[1]);
	int N_time_steps = atoi(argv[2]);
	printf("Simulating at K=%f for %d time steps\n", K, N_time_steps);
	
	/* Load omegas from file */
	double * omegas = load_data("population");
	double * thetas = load_data("positions");
	int N = sb_count(omegas);
	if (sb_count(thetas) != N) {
		printf("Omegas and thetas have different numbers of oscillators!\n");
		exit(1);
	}
	
	/* Build theta differences */
	double * theta_diffs = malloc(sizeof(double) * (N-1));
	int i;
	for (i = 1; i < N; i++) {
		theta_diffs[i-1] = thetas[i] - thetas[i-1] - M_PI_2;
		if (theta_diffs[i-1] < -M_PI) theta_diffs[i-1] += 2 * M_PI;
		if (theta_diffs[i-1] > M_PI) theta_diffs[i-1] -= 2 * M_PI;
	}
	
	/* Run a bunch of time steps */
	run_sim(omegas, thetas, theta_diffs, N_time_steps, K);
	
	/* Save the final theta positions in binary form for potential later use */
	char out_filename[160];
	sprintf(out_filename, "final-thetas,K=%s.bin", argv[1]);
	FILE * out_fh = fopen(out_filename, "w");
	fwrite(thetas, sizeof(double), sb_count(thetas), out_fh);
	fclose(out_fh);
	
	/* Clean up */
	sb_free(omegas);
	sb_free(thetas);
	free(theta_diffs);
	return (0);
}

/* Compute the velocity for the given positions, or an Euler step from them
 * using the given velocity. */
void compute_vel(double * om, double * th, double K, double * prev_vel, double dt, double * out_vel);

void run_sim(double * om, double * th, double * th_diff, int N_time_steps, double K) {
	int N_osc = sb_count(om);
	int update_tick = 1000000 / N_osc;
	
	/* A simulation constant */
	double dt = 0.01;
	
	/* Allocate a couple of working vectors */
	double * v1 = malloc(sizeof(double) * N_osc);
	double * v2 = malloc(sizeof(double) * N_osc);
	double * v3 = malloc(sizeof(double) * N_osc);
	double * v4 = malloc(sizeof(double) * N_osc);
	double * th_vel = malloc(sizeof(double) * N_osc);
	
	char filename[80];
	sprintf(filename, "slips,K=%f.txt", K);
	FILE * out_fh = fopen(filename, "w");
	
	int i, j;
	for (i = 0; i < N_time_steps; i++) {
		if (i % update_tick == 0) printf("\rStep %d/%d", i, N_time_steps);
		fflush(stdout);
		
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
			if (th_diff[j] > M_PI) {
				fprintf(out_fh, "%d %d 1\n", i, j);
				th_diff[j] -= 2.0 * M_PI;
			}
			else if (th_diff[j] < -M_PI) {
				fprintf(out_fh, "%d %d -1\n", i, j);
				th_diff[j] += 2.0 * M_PI;
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
				if (curr_diff - th_diff[j] > 6) {
					fprintf(out_fh, "%d %d -1\n", i, j);
				}
				else if (curr_diff - th_diff[j] < -6) {
					fprintf(out_fh, "%d %d 1\n", i, j);
				}
				else if (curr_diff - th_diff[j] > M_PI
					|| curr_diff - th_diff[j] < -M_PI
				) {
					printf("Excessively large dtheta difference of %f at i=%d, j=%d\n",
						curr_diff - th_diff[j], i, j);
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
	fclose(out_fh);
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

