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
		int N_read;
		while(1) {
			/* Read on number at a time */
			N_read = fread(&datum, sizeof(double), 1, in_fh);
			/* Exit when there's nothing else to read */
			if (N_read != 1) break;
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

