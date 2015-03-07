use strict;
use warnings;

# Creates a new population and set of initial positions, overwriting if
# necessary. Also overwrites the provenance file with new material.

use PDL;
my $N_osc = shift || die "You must give a number of oscillators\n";

# Make sure the working directory is clean
die "Please commit your changes before creating a new population\n"
	if !!`git diff --shortstat`;

# Create the population
zeros($N_osc)->grandom->qsort->wcols("population.txt");

# Create the initial positions
my $pi = atan2(-1, 0);
my $thetas = zeros($N_osc)->random * 2 * $pi;
$thetas -= $pi;
$thetas->wcols("positions.txt");

# Start a new provenance file
my $commit_hash = `git log -1 --pretty=%h`;
open my $out_fh, '>', 'provenance.txt';
print $out_fh "gen_pop.pl created population and positions at " . localtime . " with code from commit $commit_hash";
close $out_fh;
