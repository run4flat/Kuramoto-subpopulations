/* From GoodPracticeRNG.pdf */
unsigned int devrand(void)
{
	int fn;
	unsigned int r;
	fn = fopen("/dev/urandom", "r");
	if (fn == -1)
		exit(-1); /* Failed! */
	if (read(fn, &r, 4) != 4)
		exit(-1); /* Failed! */
	close(fn);
	return r;
}
static unsigned int x = 123456789,y = 362436000,z = 521288629,c = 7654321; /* Seed variables */
/* Initialise KISS generator using /dev/urandom */
void init_KISS() {
	x = devrand();
	while (!(y = devrand())); /* y must not be zero! */
	z = devrand();
	/* We don’t really need to set c as well but let's anyway... */
	/* NOTE: offset c by 1 to avoid z=c=0 */
	c = devrand() % 698769068 + 1; /* Should be less than 698769069 */
	printf("Seeding with x=%d, y=%d, z=%d, c=%d\n", x, y, z, c);
}
unsigned int KISS() {
	unsigned long long t, a = 698769069ULL;
	x = 69069*x+12345;
	y ^= (y<<13); y ^= (y>>17); y ^= (y<<5); /* y must never be set to zero! */
	t = a*z+c; c = (t>>32); /* Also avoid setting z=c=0! */
	return x+y+(z=t);
}
#define unsigned_int_max 4294967295.0
double dKISS() {
	return ((double)KISS() / unsigned_int_max);
}
