#ifndef RAND_H_
#define RAND_H_
#include <vector>

class Rand {
    /* A version of Mersenne Twister, code is extracted from R/src/main/RNG.c */
    private:
        unsigned int seeds[625];
        
    public:
        Rand(const std::vector<unsigned int> & s) {
            for (int i=0; i<625; i++)
                seeds[i] = s[i];
            }
        ~Rand() {}
    
        void setSeed(const std::vector<unsigned int> & s) { 
            for (int i=0; i<625; i++)
                seeds[i] = s[i];
            }
    
        void getSeed(std::vector<unsigned int> & s) { 
            for (int i=0; i<625; i++)
                s[i] = seeds[i];
            }
        
        double runif() {
            // Returns a pseudo-random number between 0 and 1.
            double ans = MT_genrand_real();
            /* ensure 0 and 1 are never returned */
            if(ans <= 0.0) return 0.5*2.328306437080797e-10;
            if((1.0 - ans) <= 0.0) return 1.0 - 0.5*2.328306437080797e-10;
            return ans;
        }
        unsigned int rint(unsigned int n) {
            return (MT_genrand_int(n));
        }
    private:
        void MT_sgenrand(unsigned int seed, unsigned int *mt, int & mti);
        unsigned int MT_genrand();
        double MT_genrand_real();
        unsigned int MT_genrand_int(unsigned int n);
};

#endif
