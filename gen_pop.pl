use strict;
use warnings;

use PDL;
my $N_osc = shift || die "You must give a number of oscillators\n";
zeros($N_osc)->grandom->qsort->wcols("population.txt");

my $pi = atan2(-1, 0);
my $thetas = zeros($N_osc)->random * 2 * $pi;
$thetas -= $pi;
$thetas->wcols("positions.txt");
