#ifndef XOR_RAND_HPP
#define XOR_RAND_HPP

#define _USE_MATH_DEFINES
#include <cmath>
#include <climits>
#include <stdint.h>

class xor_rand
{
public:
	xor_rand(unsigned int seed = 0, unsigned int id = 0);
	~xor_rand() {};

	double gen_rand();
	void reset(unsigned int seed, unsigned int id);

	unsigned int gen_rand_uint();

	unsigned long int xx, yy, zz, ww;
};
#endif