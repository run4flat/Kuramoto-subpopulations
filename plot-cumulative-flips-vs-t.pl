use strict;
use warnings;

my $is_signed = 0;
if ($ARGV[0] and $ARGV[0] eq 'signed') {
	$is_signed = 1;
	shift @ARGV;
}

die "You must give at least one input file\n" unless @ARGV;

use PDL;
use PDL::Graphics::Prima::Simple;

# Create one data set per data file
my %datasets;
for my $file (@ARGV) {
	my ($xs, $dir) = rcols($file, 0, 2);
	$xs *= 0.01;
	my $ys = $is_signed ? $dir->cumusumover : $xs->sequence;
	$datasets{'-' . $file} = ds::Pair($xs, $ys, plotTypes => ppair::Lines);
}
plot(%datasets);

__END__
my $S = $xs->nelem;
my $S_x = sumover($xs);
my $S_y = sumover($ys);
my $S_xx = sumover($xs*$xs);
my $S_xy = sumover($xs*$ys);
my $slope = ($S_xy * $S - $S_x * $S_y) / ($S_xx * $S - $S_x * $S_x);
my $y0 = ($S_xy - $slope * $S_xx) / $S_x;

my $dy = $ys - $y0 - $slope * $xs;
line_plot($xs, $dy);
