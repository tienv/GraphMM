#include "rand.h"
using namespace std;

/* Period parameters */
#define Nseed 624
#define MM 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

/* Tempering parameters */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)


/* Initializing the array with a seed */
void Rand::MT_sgenrand(unsigned int seed, unsigned int *mt, int & mti)
{
    int i;
    
    for (i = 0; i < Nseed; i++) {
        mt[i] = seed & 0xffff0000;
        seed = 69069 * seed + 1;
        mt[i] |= (seed & 0xffff0000) >> 16;
        seed = 69069 * seed + 1;
    }
    mti = Nseed;
}

unsigned int Rand::MT_genrand()
{
    unsigned int *mt = seeds + 1; /* the array for the state vector  */
    int mti = Nseed + 1; /* mti==N+1 means mt[N] is not initialized */

    unsigned int y;
    static unsigned int mag01[2]={0x0, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    mti = seeds[0];

    if (mti >= Nseed) { /* generate N words at one time */
        int kk;
    
        if (mti == Nseed+1)   /* if sgenrand() has not been called, */
            MT_sgenrand(4357, mt, mti); /* a default initial seed is used   */

        for (kk = 0; kk < Nseed - MM; kk++) {
            y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
            mt[kk] = mt[kk+MM] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (; kk < Nseed - 1; kk++) {
            y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
            mt[kk] = mt[kk+(MM-Nseed)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[Nseed-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
        mt[Nseed-1] = mt[MM-1] ^ (y >> 1) ^ mag01[y & 0x1];

        mti = 0;
    }

    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);
    seeds[0] = mti;

    return (y); 
}

double Rand::MT_genrand_real()
{
    return (double(MT_genrand()) * 2.3283064365386963e-10);
}
unsigned int Rand::MT_genrand_int(unsigned int n)
{
    return (MT_genrand() % n);
}

/* Functions to check rand class (compare with R) */
// void PrintVec2(const vector<unsigned int> & V, const char *s)
// {
//     int n;
//     n = V.size();
//     cout << s << "\n";
//     for (int i=0; i<n; i++)
//     {
//         cout << (int)V[i] << " ";
//     }
//     cout << "\n";
// }
// 
// int main()
// {
//     vector<unsigned int> seeds;
//     ifstream fin("LastSeed.txt");
//     if (fin.is_open())
//     {
//         ReadVec(seeds, fin);
//         fin.close();
//     }
//     cout.precision(15);
//     cout << "seed[0]: " << (int)seeds[0] << " seed[624]: " << (int)seeds[624] << endl;
//     Rand * rand_item = new Rand(seeds);
//     cout << "Normal: " << "\n";
//     for (int i=0; i<50; i++)
//     {
//         cout << rand_item->runif() << endl;
//         cout << my_rnorm(2.5, 1.5, rand_item) << endl;
//         cout << my_rbeta(50.0, 20.0, rand_item) << endl;
//     }
//     vector<unsigned int> s(625);
//     rand_item->getSeed(s);
//     cout << "seed[0]: " << (int)s[0] << " seed[50]: "<< (int)s[50] << " seed[624]: " << (int)s[624] << endl;
// }
// 
