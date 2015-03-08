use strict;
use warnings;

# A very simple module that checks that the repository with the simulation
# code is checked-in.

my $CODE_DIRECTORY = "/home/dcmertens/projects/2015/03-06-Kuramoto-Slips";

# This can be run from a directory that is not the code directory. Temporarily
# switch to the code directory to perform the git check.
use Cwd;
my $cwd = getcwd;
chdir $CODE_DIRECTORY or die "Unable to chdir into code directory\n";

# make sure the repository is good
die "Please commit your changes to the codebase before continuing.\n"
	if !!`git diff --shortstat`;

# Get the latest commit hash
$Repo::hash = `git log -1 --pretty=%h`;

# All done, switch back
chdir $cwd;

# Finish with a true statement.
1;
