# Uses results from steady-state-slips.c to determine the subpopulation
# structure as a function of time. Plots the number of identified
# subpopulations as a function of time.
use strict;
use warnings;

die "You must give at least one input file\n" unless @ARGV;

use PDL;
use PDL::Graphics::Prima::Simple;

# Get the number of oscillatos
my $N_osc = 0;
open my $in_fh, '<', 'population.txt';
$N_osc++ while(<$in_fh>);
close $in_fh;

# Create one data set per data file
my %datasets;

for my $file (@ARGV) {
	my (@times, @N_subpops);
	my $N_slips = zeros($N_osc - 1);
	open my $in_fh, '<', $file;
	while(my $line = <$in_fh>) {
		if ($line =~ /(\d+)\s+(\d+)\s+(-?\d+)/) {
			my ($time_offset, $boundary, $dir) = ($1, $2, $3);
			$N_slips->slice([$boundary]) += $dir;
			push @times, $time_offset * 0.01;
			push @N_subpops, sum($N_slips != 0);
		}
		else {
			print "Weird data for file $file on line $.\n";
		}
	}
	close $in_fh;
	
	my $xs = pdl(\@times);
	my $ys = pdl(\@N_subpops);
	$datasets{'-' . $file} = ds::Pair($xs, $ys, plotTypes => ppair::Lines);
}
plot(%datasets);
