use strict;
use warnings;

# Creates a new population and set of initial positions, overwriting if
# necessary. Also overwrites the provenance file with new material.

use PDL;
my $N_osc = shift || die "You must give a number of oscillators\n";

# Make sure the simulation code is in a clean state
require '/home/dcmertens/projects/2015/03-06-Kuramoto-Slips/RepoCheck.pm';

# Create the population
zeros($N_osc)->grandom->qsort->wcols("population.txt");

# Create the initial positions
my $pi = atan2(-1, 0);
my $thetas = zeros($N_osc)->random * 2 * $pi;
$thetas -= $pi;
$thetas->wcols("positions.txt");

$Repo::hash = $Repo::hash; # avoid stupid warning

# Start a new provenance file
open my $out_fh, '>', 'provenance.txt';
print $out_fh "gen_pop.pl created population and positions at " . localtime . " with code from commit $Repo::hash";
close $out_fh;
