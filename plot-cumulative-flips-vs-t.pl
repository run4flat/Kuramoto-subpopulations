use strict;
use warnings;
use PDL;
my $infile = shift @ARGV or die "You must give an input file\n";
my $xs = rcols($infile) * 0.01;
my $ys = $xs->sequence;
use PDL::Graphics::Prima::Simple;
plot(
	-data => ds::Pair($xs, $ys, plotTypes => [
		ppair::Lines, #ppair::TrendLines(color => cl::LightRed)
	]),
);

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
