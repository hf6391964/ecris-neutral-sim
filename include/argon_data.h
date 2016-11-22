#pragma once

#include "element_data.h"

const ElementData ARGON_DATA = {
    18,  // Z
    39.962384,  // mass in eV
    { 
      15.8, 2.532, { -2.672, 2.543, -0.769, 0.008, 0.006 },  // low energy parameters
      4.131,  // impact energy limit
      15.8, 4.337, { 3.092, -21.253, 14.626, 0.018, 0.031 }  // high energy parameters
    }
};

