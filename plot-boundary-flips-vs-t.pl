use strict;
use warnings;

die "You must give an input file\n" unless @ARGV >= 1;
my $file = shift;
my $range = shift || '0:-1';

# Get the number of oscillatos
my $N_osc = 0;
open my $in_fh, '<', 'population.txt';
$N_osc++ while(<$in_fh>);
close $in_fh;

use PDL;
use PDL::Graphics::Prima::Simple;
my $is_tracking = zeros($N_osc);
$is_tracking->slice($range)++;

# Track the time and flips for each boundary
open $in_fh, '<', $file or die "Could not open $file\n";
my (@flip_counts, @times);
while (my $line = <$in_fh>) {
	if ($line =~ /(\d+)\s+(\d+)\s+(-?\d+)/) {
		my ($time_offset, $boundary, $dir) = ($1, $2, $3);
		next unless $is_tracking->at($boundary);
		
		# Get the current list of times and counts
		my $curr_counts = $flip_counts[$boundary];
		my $curr_times = $times[$boundary];
		if ($curr_counts) {
			# Add another entry
			push @$curr_counts, $curr_counts->[-1] + $dir;
			push @$curr_times, $time_offset * 0.01;
		}
		else {
			# Create a new list if it's really new
			$flip_counts[$boundary] = [$dir];
			$times[$boundary] = [$time_offset * 0.01];
		}
	}
	else {
		print "Weird data for file $file on line $.\n";
	}
}

# Convert boundary arrays into vectors
my %datasets;
for my $i (0 .. $#flip_counts) {
	next unless $flip_counts[$i];
	$datasets{"-$i"} = ds::Pair(pdl($times[$i]), pdl($flip_counts[$i]),
		plotType => ppair::Lines);
}
plot(%datasets);
