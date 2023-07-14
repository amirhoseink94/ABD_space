// openGL libraries
#include <GL/gl.h>
#include <GLFW/glfw3.h>
#include <GL/glut.h>

// my libraries
#include "../include/Util.h"


// c++ libraries
#include <vector>
#include <cmath>
#include <algorithm>




// namespaces
using namespace std;
using namespace arma;

namespace ABD
{
  double Util::distance_VF(Particle& P, Face& F)
  {
    vec z(12, fill::zeros);
    //dQ_P = z;
    //dQ_F = z;
    vec AP = P.pos - F.x->pos;
    vec BP = P.pos - F.y->pos;
    vec CP = P.pos - F.z->pos;


    vec AB = F.y->pos - F.x->pos;
    vec BC = F.z->pos - F.y->pos;
    vec CA = F.x->pos - F.z->pos;

    // check for on plan
    double v_distance = abs(dot(F.n, P.pos - F.x->pos));
    //dQ_P = derivitives_VF_dV(F.n, P.J, v_distance, ddQ_P );
    double dist = 0;

    //mat dH = F.y->J - F.x->J;
    //mat dJ = F.z->J - F.x->J;
    //mat dK = F.z->J - F.y->J;

    //dQ_F = derivitives_VF_dF(dH, dJ, F.x->J, Q_F, P.pos, F.x->pos, v_distance, ddQ_F);

    int sgn_AB = (dot(F.n_AB, P.pos)-dot(F.n_AB, F.x->pos) >= 0) ? 1 : 0;
    int sgn_BC = (dot(F.n_BC, P.pos)-dot(F.n_BC, F.y->pos) >= 0) ? 1 : 0;
    int sgn_CA = (dot(F.n_CA, P.pos)-dot(F.n_CA, F.z->pos) >= 0) ? 1 : 0;

    if(sgn_AB == 1 && sgn_BC == 1 && sgn_CA == 1)
    {
      // in this case, we do not need to do anything extra, the point is in the
      // easy position to calculate everything. S, yay!

      dist = v_distance;

    }
    else if(sgn_AB == 0 && sgn_BC == 1 && sgn_CA == 1)
    {
      // in this case, the point P is outside of the triangel, i.e., out of AB,
      // so we have to look into an extra length, and then add it fo dQ_F, and ....
      int n_A = (norm_dot(AP, AB)>0)?1:0;
      int n_B = (norm_dot(BP, -AB)>0)?1:0;


      if(n_A == 1 && n_B == 1)
      {

        double h_distance = -(dot(F.n_AB, P.pos) - dot(F.n_AB, F.x->pos));
        dist =  sqrt(v_distance*v_distance + h_distance*h_distance);
      }
      else if(n_A == 0 && n_B == 1)
      {
        dist = norm(AP);
      }
      else if(n_A == 1 && n_B == 0)
      {
        dist = norm(BP);
      }
    }
    else if(sgn_AB == 1 && sgn_BC == 0 && sgn_CA == 1)
    {
      int n_B = (norm_dot(BP, BC)>0)?1:0;
      int n_C = (norm_dot(CP, -BC)>0)?1:0;
      if(n_B == 1 && n_C == 1)
      {
        double h_distance = -(dot(F.n_BC, P.pos) - dot(F.n_BC, F.y->pos));
        dist = sqrt(v_distance*v_distance + h_distance*h_distance);
      }
      else if(n_B == 0 && n_C == 1)
      {
        dist = norm(BP);
      }
      else if(n_B == 1 && n_C == 0)
      {
        dist = norm(CP);
      }
    }
    else if(sgn_AB == 1 && sgn_BC == 1 && sgn_CA == 0)
    {
      int n_C = (norm_dot(CP, CA)>0)?1:0;
      int n_A = (norm_dot(AP, -CA)>0)?1:0;
      if(n_C == 1 && n_A == 1)
      {
        double h_distance = -(dot(F.n_CA, P.pos) - dot(F.n_CA, F.z->pos));
        dist = sqrt(v_distance*v_distance + h_distance*h_distance);
      }
      else if(n_C == 0 && n_A == 1)
      {
        dist = norm(CP);
      }
      else if(n_C == 1 && n_A == 0)
      {
        dist = norm(AP);
      }
    }
    else if(sgn_AB == 1 && sgn_BC == 0 && sgn_CA == 0)
    {
      dist = norm(CP);
    }
    else if(sgn_AB == 0 && sgn_BC == 1 && sgn_CA == 0)
    {
      dist = norm(AP);
    }
    else if(sgn_AB == 0 && sgn_BC == 0 && sgn_CA == 1)
    {
      dist = norm(BP);
    }
    return dist;
  }

  bool Util::CCD_VF(vec A, vec B, Face& F, double& toi)
  {
    vec AB = B - A;

    double c = dot(F.n, F.x->pos);

    double sign_A = dot(F.n, A) - c;
    double sign_B = dot(F.n, B) - c;

    if(abs(sign_A) < 1e-10)
    {
      vec z(3, fill::zeros);
      Particle p(A, z, 1, z);
      double dist = distance_VF(p, F);
      if(dist < 1e-10)
      {
        toi = 0;
        return true;
      }
      else
      {
        toi = -1;
        return false;
      }
    }
    if(abs(sign_B) < 1e-10)
    {
      vec z(3, fill::zeros);
      Particle p(B, z, 1, z);
      double dist = distance_VF(p, F);
      if(dist < 1e-10)
      {
        toi = 1;
        return true;
      }
      else
      {
        toi = -1;
        return false;
      }
    }

    if(sign_A * sign_B < 0)
    {
      double t = (c - dot(F.n, A)) / dot(F.n, AB);
      vec C = A + t * AB;
      vec z(3, fill::zeros);
      Particle p(C, z, 1, z);
      double dist = distance_VF(p, F);
      if(dist < 1e-10)
      {
        toi = t;
        return true;
      }
      else
      {
        toi = -1;
        return false;
      }
    }
    else
    {
      toi = -1;
      return false;
    }
  }
}
