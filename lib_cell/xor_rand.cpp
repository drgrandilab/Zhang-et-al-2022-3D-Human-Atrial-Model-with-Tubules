#include "xor_rand.hpp"

xor_rand::xor_rand(unsigned int seed, unsigned int id) {

	reset(seed, id);
}

double xor_rand::gen_rand() {
	unsigned int t = (xx ^ (xx << 11));
	xx = yy;
	yy = zz;
	zz = ww;
	return ( ww = (ww ^ (ww >> 19)) ^ (t ^ (t >> 8)) ) / (double)(UINT_MAX);
}

unsigned int xor_rand::gen_rand_uint() {
	unsigned int t = (xx ^ (xx << 11));
	xx = yy;
	yy = zz;
	zz = ww;
	return ( ww = (ww ^ (ww >> 19)) ^ (t ^ (t >> 8)) );
}

void xor_rand::reset(unsigned int seed, unsigned int id) {

	xx = 123456789 + id + seed;
	yy = 362436069 + id * 100 + seed * 10;
	zz = 521288629 + id * 1000 + seed * 100;
	ww = 88675123 + id * 10000 + seed * 1000;

	for (int i = 0; i < 1000; i++)
		gen_rand();
}