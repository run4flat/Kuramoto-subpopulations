/*
 * mser.c - computes the time at which the phase slip vs time reaches equilibrium
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "stretchy_buffer.h"

void save_results(char * filename, int offset, double trans_time, double slope, double y0);

int main(int argc, char **argv) {
	if (argc != 3) {
		printf("You must give a filename and a boolean to indicate signed or unsigned analysis\n");
		exit(1);
	}
	
	/* Assume the first argument is a filename; get the signedness */
	int is_signed = atoi(argv[2]);
	int * slip_time_indices = 0;
	int * cumulative_N_slips = 0;
	int cumulative = 0;
	
	/* open the file and get the time of each slip */
	FILE * in_fh = fopen(argv[1], "r");
	int new_time_index, new_boundary_index, new_sign;
	while (!feof(in_fh)) {
		fscanf(in_fh, "%d %d %d", &new_time_index, &new_boundary_index, &new_sign);
		sb_push(slip_time_indices, new_time_index);
		cumulative += is_signed ? new_sign : 1;
		sb_push(cumulative_N_slips, cumulative);
	}
	fclose(in_fh);
	
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
	
	/* Accumulate information for the first ten points; do not calculate
	 * MSER scores for them, however */
	N_data = sb_count(slip_time_indices);
	for (i = N_data - 1; i > N_data - 11; i--) {
		x = slip_time_indices[i] * 0.01;
		y = cumulative_N_slips[i];
		N++;
		S_x += x;
		S_y += y;
		S_xx += x*x;
		S_xy += x*y;
		S_yy += y*y;
	}
	double slope = (S_xy * N - S_x * S_y) / (S_xx * N - S_x * S_x);
	double y0 = (S_xy - slope * S_xx) / S_x;
	
	/* keep track of the best score, its index, slope, and intercept */
	double best_MSER_score 
		= (N * y0*y0 + 2 * y0 * slope * S_x
		- 2 * y0 * S_y + slope*slope * S_xx
		 - 2 * slope * S_xy + S_yy) / N/N;
	int best_MSER_index = i;
	double best_slope = slope;
	double best_y0 = y0;
	
	/* Now compute the MSER score for the remaining points */
	for (; i >= 0; i--) {
		x = slip_time_indices[i] * 0.01;
		y = cumulative_N_slips[i];
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
			best_MSER_index = i;
			best_slope = slope;
			best_y0 = y0;
		}
	}
	
	/* Add this information to the list of transient data */
	save_results(argv[1], best_MSER_index,
		slip_time_indices[best_MSER_index] * 0.01, best_slope, best_y0);
	/* print the result, too */
	printf("Best MSER score occurs at index %d, time %f, with slope %g and intercept %f\n",
		best_MSER_index, slip_time_indices[best_MSER_index] * 0.01, best_slope, best_y0);
}

void save_results(char * filename, int offset, double trans_time, double slope, double y0) {
	FILE * in_fh;
	FILE * out_fh;
	out_fh = fopen("transients.tmp", "w");
	if (!out_fh) {
		printf("!!! Unable to save results!!!\n");
		return;
	}
	
	char curr_filename[160];
	int curr_offset;
	double curr_trans_time, curr_slope, curr_y0;
	
	/* Copy previous results, excluding any results for this file */
	if (in_fh = fopen("transients.txt", "r")) {
		while(!feof(in_fh)) {
			int N_read = fscanf(in_fh, "%s %d %lf %lf %lf", curr_filename, &curr_offset,
				&curr_trans_time, &curr_slope, &curr_y0);
			if (N_read != 5) break;
			if (strcmp(curr_filename, filename) != 0) {
				fprintf(out_fh, "%s %d %lf %lf %lf\n", curr_filename, curr_offset,
					curr_trans_time, curr_slope, curr_y0);
			}
			else {
				printf("Replacing previous result for %s\n", filename);
			}
		}
		fclose(in_fh);
	}
	
	/* Add this file's results */
	fprintf(out_fh, "%s %d %lf %lf %lf\n", filename, offset, trans_time,
		slope, y0);
	fclose(out_fh);
	
	/* Shuffle files around */
	remove("transients.txt");
	rename("transients.tmp", "transients.txt");
}
