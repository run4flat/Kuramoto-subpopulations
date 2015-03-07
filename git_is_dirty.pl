use strict;
use warnings;

# Checks if the current git repo is in a clean state.
die "Please commit your changes before continuing.\n"
	if !!`git diff --shortstat`;
exit(0);
