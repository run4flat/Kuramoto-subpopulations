use strict;
use warnings;

# Called by slips.c to add information to the provenance file based on
# the current simulation.

# Make sure the simulation code is in a clean state
require '/home/dcmertens/projects/2015/03-06-Kuramoto-Slips/RepoCheck.pm';

# Get the simulation parameters
my ($K, $N_time_steps) = @ARGV;

$Repo::hash = $Repo::hash; # avoid stupid warning

# Add the provenance method
open my $out_fh, '>>', 'provenance.txt';
print $out_fh "slips.c ran a simulation at " . localtime . " with K=$K and $N_time_steps time steps, using code from commit $Repo::hash";
close $out_fh;
