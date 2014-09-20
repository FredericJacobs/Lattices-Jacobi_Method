//
//  main.cpp
//  Jacobi-LLL-Benchmarking
//
//  Created by Frederic Jacobs on 20/09/14.
//  Copyright (c) 2014 EPFL. All rights reserved.
//

#include <iostream>
#include <fplll.h>
#include "lib/lattice_gen/generate_random.cpp"

int main(int argc, const char * argv[]) {
    
    /**
     *  Anja recommended that I use the SVP Challenge Generator
     */
    
    mat_ZZ randomLattice = genSVPChallengeLattice();
    
    std::cout << "Yo, lattices: " << randomLattice;
    
    /**
     *  Try Jacobian Lattices
     */
    
    
    /**
     *  Try LLL Reduction
     */
    NTL::ZZ det;
    std::cout << "Reducted Lattice: " << LLL(det, randomLattice);

    
    /**
     *  Compare results
     */
    
    
    std::cout << "Hello, World!\n";
    return 0;
}
