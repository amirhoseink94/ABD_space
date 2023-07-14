#ifndef UTIL_H
#define UTIL_H

#include <set>
#include <iostream>
#include <stdio.h>
#include <armadillo>
#include <bits/stdc++.h>

#include "Particle.h"
#include "Segment.h"
#include "Face.h"
#include "Body.h"

using namespace std;
using namespace arma;


namespace ABD
{
	class Util
	{
	public:
		static double distance_VF(Particle& P, Face& F);
		static bool CCD_VF(vec A, vec B, Face& F, double& toi);
	};
}
#endif
