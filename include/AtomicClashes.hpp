/* Purpose: Fix atomic clashes */

#ifndef AtomicClashes_HPP
#define AtomicClashes_HPP 1

#include <math.h>
#include <stdbool.h>
#include <vector>
#include <cstdlib>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <map>

#include "GeometryTools.hpp"
#include "PDBParser.hpp"

using namespace std;

int fix_clashes (ModelUnit &pep);

// A
// Atomic order:     P   OP1  OP2   O5'  C5' C4'   O4'  C3'  O3'  C2'  O2'  C1'   N9   C8   N7     C5     C6   N6    N1    C2    N3    C4
static double A_vdw[22] = {1.8, 1.4, 1.4, 1.4, 1.7, 1.7, 1.4, 1.7, 1.4, 1.7, 1.4, 1.7, 1.55, 1.75, 1.55, 1.75, 1.75, 1.55, 1.55, 1.75, 1.55, 1.75};

// C
// Atomic order:     P   OP1  OP2   O5'  C5'  C4'  O4'  C3'  O3'  C2'  O2'  C1'   N1   C2   O2    N3   C4    N4    C5   C6
static double C_vdw[20] = {1.8, 1.4, 1.4, 1.4, 1.7, 1.7, 1.4, 1.7, 1.4, 1.7, 1.4, 1.7, 1.55, 1.75, 1.4, 1.55, 1.75, 1.55, 1.75, 1.75};

// G
// Atomic order:     P   OP1  OP2   O5'  C5'  C4'  O4'  C3'  O3'  C2'  O2'  C1'   N9   C8    N7    C5   C6     O6    N1   C2    N2    N3    C4
static double G_vdw[23] = {1.8, 1.4, 1.4, 1.4, 1.7, 1.7, 1.4, 1.7, 1.4, 1.7, 1.4, 1.7, 1.55, 1.75, 1.55, 1.75, 1.75, 1.4, 1.55, 1.75, 1.55, 1.55, 1.75};

// U
// Atomic order:     P   OP1  OP2   O5'  C5'  C4'  O4'  C3'  O3'  C2'  O2'  C1'   N1   C2   O2    N3    C4   O4    C5   C6
static double U_vdw[20] = {1.8, 1.4, 1.4, 1.4, 1.7, 1.7, 1.4, 1.7, 1.4, 1.7, 1.4, 1.7, 1.55, 1.75, 1.4, 1.55, 1.75, 1.4, 1.75, 1.75};


// The order of each array is the standard PDB order, with hydrogen atoms omitted. VDW radii are in units of Angstroms.

#endif