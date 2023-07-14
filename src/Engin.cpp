// openGL libraries
#include <GL/gl.h>
#include <GLFW/glfw3.h>
#include <GL/glut.h>

// my libraries
#include "../include/Engin.h"
#include "../include/Util.h"

// c++ libraries
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>


// namespaces
using namespace std;
using namespace arma;

namespace ABD
{
  void Engin::add_object(Body* b)
  {
    if(b->movable)
      dynamic_objects.push_back(b);
    else
      static_objects.push_back(b);
  }

// =============================================================================
  void Engin::init()
  {
    N = dynamic_objects.size();

    Q.zeros(12 * N);
    dQ.zeros(12 * N);
    H.zeros(12 * N, 12 * N);
    int i = 0;
    for(auto itr = dynamic_objects.begin(); itr!=dynamic_objects.end(); itr++)
    {
        Q.subvec(i, i+11) = (*itr)->q;
        i = i+12;
    }
    Dl_t = 1./32.;
    G = 1.1;
    big_mass = 100;
    sp_value = 100;
    distance_tresh_hold = 3;
    distance_tresh_hold_contanct = 0.1;

    energy_tresh_hold = 0.01;
  }

// =============================================================================
  vec Engin::calculate_next_Q()
  {
    // update q_mad for all of the obkects
    for(auto itr = dynamic_objects.begin(); itr!=dynamic_objects.end(); itr++)
    {
      (*itr)->update_q_mad();
    }

    vec Q_p = Q;


    vec Q_current_vec = Q;

    vec dQ_col(12*N, fill::zeros);

    // previous energy
    double energy_p = 0; // calculate_energy(Q_current_vec);
    double energy_sp = 0;
    double energy_collision = 0; //calculate_energy_collision(Q_current_vec, dQ_col); //calculating the energy of collision

    //cout<<dQ_col<<"--\n"<<energy_collision<<endl;
    energy_p += energy_collision;
    double energy = 0;
    // -----

    vec step;

    double gamma = 1;
    int counter = 0;
    int s_counter = 0;
    vec Q_temp;
    bool check_first = true;
    double old_step;
    do
    {
      // previous energy
      energy_p = calculate_energy(Q_current_vec);

      //calculating the energy of collision

      energy_p += energy_collision;

      if(check_first)
      {
        energy_sp = energy_p;
      }
      if(abs(energy_p-energy_sp)<1e-9)
        s_counter++;
      else
      {
        energy_sp = energy_p;
      }
      if(s_counter>200)
        break;
      //cout<<"start the loop "<<s_counter<<endl;
      int loop = 0;
      //calculate te global derivation
      int i = 0;
      int j = 0;
      for(auto itr = dynamic_objects.begin(); itr!=dynamic_objects.end(); itr++)
      {
        dQ.subvec(i, i+11) = (*itr)->E_der( Q_current_vec.subvec(i, i+11) );
        i = i+12;
        j++;
      }

      dQ += dQ_col; // derivatio of collision of the current Q



      // calculate the global hessian
      i = 0;
      j = 0;
      for(auto itr = dynamic_objects.begin(); itr!=dynamic_objects.end(); itr++)
      {
        H.submat(i, i, i+11, i+11) = (*itr)->E_der_der( Q_current_vec.subvec(i, i+11) );
        i = i+12;
        j++;
      }
      // Hessian at this stage will remain the same, because we chose
      // the collision potential vert carefully
      H = PSPD(H);


      step = -1 * inv(H) * dQ;
      if(check_first)
      {
        check_first = false;
        old_step = (1.0 / Dl_t) * norm(step, "inf");
      }



      //cout<<"step is: "<<endl<<step<<endl;
      gamma = 1;
      Q_temp = Q_p + gamma * step;

      double toi = collision(Q_p, Q_temp);

      if(toi != 2)
        gamma = toi;


      int repeat_inner_loop=0;
      //finding the best local place that we can go with this derivite
      do
      {
        loop++;
        // we should go for the collision, chck for wich gamma ther is no collision
        Q_temp = Q_p + gamma * step;

        gamma = gamma * 0.5;
        energy = calculate_energy(Q_temp);


        // with respect to the new position Q_temp, we have to calculate the new energy_p
        // the overal energy(movements and orthoganility) is calculated,
        // we just need to calculate the energy of collisio based on Q_temp
        if(abs(energy - energy_p)<1e-9)
        {
          Q_current_vec = Q_temp;
          //cout<<"==="<<endl<<energy<<"|"<<energy_p<<endl;
          break;
        }
        //cout<<"==="<<endl<<energy<<"|"<<energy_p<<endl;
        if(abs(energy - energy_p)<1e-9)
          repeat_inner_loop++;
        if(repeat_inner_loop>=100)
          break;
      } while(energy > energy_p && loop<1024 );

      energy_p = energy;

      Q_p = Q_current_vec = Q_temp;

      double w_step = (1.0 / Dl_t) * norm(step, "inf");
      if(abs(w_step-old_step)<1e-9)
        counter++;
      else
      {
        old_step = w_step;
      }
      if(counter>200)
        break;
    } while( (1.0 / Dl_t) * norm(step, "inf") > energy_tresh_hold );
    //cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<the search is over: "<<
    //      counter<<" >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
    return Q_current_vec;
  }

// =============================================================================

  void Engin::init_particles(int n)
  {
    //particles.resize(n);
    vec z(3, fill::zeros);
    int up_range = 50;
    for(int i=0; i<n;i++)
    {
      vec pos(3, fill::zeros);
      vec v = arma::randu<arma::vec>(3) * 2;
      random_device rd;   // Obtain a random seed from the hardware
    	mt19937 gen(rd());  // Seed the random number generator
    	uniform_int_distribution<int> dis(0, up_range);  // Define the range


      pos(0) = dis(gen) - up_range/2;
      pos(1) = dis(gen) - up_range/2;
      pos(2) = dis(gen) - up_range/2;
      Particle* temp = new Particle(pos, z, 1, z);
      //temp->history.push_back(pos);
      particles.push_back(temp);

    }
    return;
  }

  void Engin::render_field()
  {
    vec pos_new;
    vec z(3, fill::zeros);

    for(auto itr = particles.begin(); itr != particles.end(); itr++)
    {
      (*itr)->F = (*itr)->F_list[0] + (*itr)->F_list[1] + (*itr)->F_list[2];
      (*itr)->F_list[0] = z;
      (*itr)->F_list[1] = z;
      (*itr)->F_list[2] = z;
      pos_new = (*itr)->pos +
                    Dl_t * (*itr)->v +
                    Dl_t * Dl_t * (*itr)->F;

      //double t = collision_particles((*itr)->pos, pos_new);
      double t = collision_particles_dynamic((*itr)->pos, pos_new);
      if(t!=2)
      {
        pos_new = (*itr)->pos + t * (pos_new - (*itr)->pos);
        delete (*itr);
        vec pos(3, fill::zeros);
        vec v = arma::randu<arma::vec>(3) * 2;
        vec z(3, fill::zeros);
        int up_range = 50;
        random_device rd;   // Obtain a random seed from the hardware
        mt19937 gen(rd());  // Seed the random number generator
        uniform_int_distribution<int> dis(0, up_range);  // Define the range

        pos(0) = dis(gen) - up_range/2;
        pos(1) = dis(gen) - up_range/2;
        pos(2) = dis(gen) - up_range/2;
        (*itr) = new Particle(pos, z, 1, z);
      }
      else
      {
        (*itr)->v = pos_new - (*itr)->pos;
        (*itr)->pos = pos_new;
        (*itr)->history.insert((*itr)->history.begin(), pos_new);
        if((*itr)->history.size() > 100)
        {
          (*itr)->history.pop_back();
        }
      }
    }
  }

  void Engin::apply_force_to_particles_dynamic()
  {
    int j = 0;
    // on each particle

    for(unsigned int i_particles = 0;
                     i_particles < particles.size();
                     i_particles++)
    { // for begin on i_points
      // on each object, which here is only one

      for(unsigned int j_faces = 0;
                       j_faces < dynamic_objects[j]->faces.size();
                       j_faces++)
      { // for begins on j_faces
        // apply gravity forces here

        // first check if the disance to the face is not zeros

        double dist = Util::distance_VF(*particles[i_particles]
                                       ,dynamic_objects[j]->faces[j_faces]);

        // cif the dist !=0 , then we can apply some force to it
        if(dist > 1e-20)
        {
          // x point
          vec dx_ji = particles[i_particles]->pos -
                       dynamic_objects[j]->faces[j_faces].x->pos;
          double length_x = norm(dx_ji);
          double gf = 0;
          if(length_x > distance_tresh_hold)
            gf = G * particles[i_particles]->m *
                   dynamic_objects[j]->faces[j_faces].x->m / (length_x * length_x);
          vec nx_ji = dx_ji / length_x;
          vec gFx_ji = gf * nx_ji / dynamic_objects[j]->faces[j_faces].x->part;
          particles[i_particles]->F_list[0] += -gFx_ji;

          // y point
          vec dy_ji = particles[i_particles]->pos -
                       dynamic_objects[j]->faces[j_faces].y->pos;
          double length_y = norm(dy_ji);
          gf = 0;
          if(length_y > distance_tresh_hold)
            gf = G * particles[i_particles]->m *
                   dynamic_objects[j]->faces[j_faces].y->m / (length_y * length_y);
          vec ny_ji = dy_ji / length_y;
          vec gFy_ji = gf * ny_ji / dynamic_objects[j]->faces[j_faces].y->part;
          particles[i_particles]->F_list[0] += -gFy_ji;

          // z point
          vec dz_ji = particles[i_particles]->pos -
                       dynamic_objects[j]->faces[j_faces].z->pos;
          double length_z = norm(dz_ji);
          gf = 0;
          if(length_z > distance_tresh_hold)
            gf = G * particles[i_particles]->m *
                   dynamic_objects[j]->faces[j_faces].z->m / (length_z * length_z);
          vec nz_ji = dz_ji / length_z;
          vec gFz_ji = gf * nz_ji / dynamic_objects[j]->faces[j_faces].z->part;
          particles[i_particles]->F_list[0] += -gFz_ji;
        }
      } // end for on j_points

    } //end for on i_point

    return;
  }
  void Engin::apply_force_to_particles()
  {
    int j = 0;
    // on each particle
    for(unsigned int i_particles = 0;
                     i_particles < particles.size();
                     i_particles++)
    { // for begin on i_points
      // on each object, which here is only one
      for(unsigned int j_faces = 0;
                       j_faces < static_objects[j]->faces.size();
                       j_faces++)
      { // for begins on j_faces
        // apply gravity forces here

        // first check if the disance to the face is not zeros

        double dist = Util::distance_VF(*particles[i_particles]
                                       ,static_objects[j]->faces[j_faces]);

        // cif the dist !=0 , then we can apply some force to it
        if(dist > 1e-20)
        {
          // x point
          vec dx_ji = particles[i_particles]->pos -
                       static_objects[j]->faces[j_faces].x->pos;
          double length_x = norm(dx_ji);
          double gf = 0;
          if(length_x > distance_tresh_hold)
            gf = G * particles[i_particles]->m *
                   static_objects[j]->faces[j_faces].x->m / (length_x * length_x);
          vec nx_ji = dx_ji / length_x;
          vec gFx_ji = gf * nx_ji / static_objects[j]->faces[j_faces].x->part;
          particles[i_particles]->F_list[0] += -gFx_ji;

          // y point
          vec dy_ji = particles[i_particles]->pos -
                       static_objects[j]->faces[j_faces].y->pos;
          double length_y = norm(dy_ji);
          gf = 0;
          if(length_y > distance_tresh_hold)
            gf = G * particles[i_particles]->m *
                   static_objects[j]->faces[j_faces].y->m / (length_y * length_y);
          vec ny_ji = dy_ji / length_y;
          vec gFy_ji = gf * ny_ji / static_objects[j]->faces[j_faces].y->part;
          particles[i_particles]->F_list[0] += -gFy_ji;

          // z point
          vec dz_ji = particles[i_particles]->pos -
                       static_objects[j]->faces[j_faces].z->pos;
          double length_z = norm(dz_ji);
          gf = 0;
          if(length_z > distance_tresh_hold)
            gf = G * particles[i_particles]->m *
                   static_objects[j]->faces[j_faces].z->m / (length_z * length_z);
          vec nz_ji = dz_ji / length_z;
          vec gFz_ji = gf * nz_ji / static_objects[j]->faces[j_faces].z->part;
          particles[i_particles]->F_list[0] += -gFz_ji;
        }
      } // end for on j_points
    } //end for on i_point
    return;
  }
  void Engin::draw_field()
  {
    for(auto itr = particles.begin(); itr != particles.end(); itr++)
    {
      glColor3ub(170, 170, 120);
      for(auto itr2 = ((*itr)->history.begin()); itr2 != ((*itr)->history.end()-1); itr2++)
      {
        glBegin(GL_LINES);
        glVertex3f((*itr2)[0], (*itr2)[1], (*itr2)[2]);
        glVertex3f((*(itr2+1))[0], (*(itr2+1))[1], (*(itr2+1))[2]);
        glEnd();
      }
    }
  }

// =============================================================================
  mat Engin::PSPD(const mat& X)
  {
    vec eigval;
  	mat eigvec;

  	eig_sym(eigval, eigvec, X);
    //cout<<"the init is:"<<endl<<eigval<<endl;
    for(unsigned long int i=0;i<eigval.size();i++)
  	{
  		if(eigval(i)<0)
  			eigval(i) = 0.000001;
  	}

  	mat C = eigvec * diagmat(eigval) * eigvec.t();
    return C;
  }

// =============================================================================
  double Engin::calculate_energy(const vec& X)
  {
    double global_e = 0;
    int i = 0;
    for(auto itr = dynamic_objects.begin(); itr!=dynamic_objects.end(); itr++)
    {
      global_e += (*itr)->calculate_energy(X.subvec(i, i+11));
      i = i+12;
    }
    return global_e;
  }

// =============================================================================
  void Engin::apply_tranformation()
  {
    // we update the each q here
    int i=0;
    for(auto itr = dynamic_objects.begin(); itr!=dynamic_objects.end(); itr++)
    {
      (*itr)->q = Q.subvec(i, i+11);
      (*itr)->apply_tranformation();
      i+=12;
    }
  }

// =============================================================================
  void Engin::apply_tranformation_temp(const vec& Q_temp)
  {
    // we update the each q here
    int i=0;
    for(auto itr = dynamic_objects.begin(); itr!=dynamic_objects.end(); itr++)
    {
      (*itr)->q = Q_temp.subvec(i, i+11);
      (*itr)->apply_tranformation();
      i+=12;
    }
  }
  // =============================================================================
    void Engin::apply_force_dynamic()
    {
      for(unsigned int i = 0; i<dynamic_objects.size(); i++)
      {

        for(unsigned int i_points = 0;
                         i_points < dynamic_objects[i]->points.size();
                         i_points++)
        { // for begin on i_points

          // apply gravity forces here
          vec d_ji = dynamic_objects[i]->points[i_points]->pos;
          double length = norm(d_ji);
          double gf = 0;
          if(length > 0.1)
            gf = G * dynamic_objects[i]->points[i_points]->m * big_mass;
          vec n_ji = d_ji / length;
          vec gF_ji = gf * n_ji;
          dynamic_objects[i]->points[i_points]->F_list[0] += -gF_ji;

          vec h = d_ji/length * sp_value;
          vec sp_force(3, fill::zeros);
          sp_force(0) = - h(0);
          sp_force(2) = h(2);
          //dynamic_objects[i]->points[i_points]->F_list[2] += sp_force;
        } //end for on i_point

      } //end for on i for other dynamic objects

      for(auto itr = dynamic_objects.begin(); itr!= dynamic_objects.end(); itr++)
      {
        (*itr)->apply_force();
      }
      return;
    }
// =============================================================================
  void Engin::apply_force()
  {
    for(unsigned int i = 0; i<dynamic_objects.size(); i++)
    {
      // looking at other dynamic objects
      for(unsigned int j=i+1; j<dynamic_objects.size(); j++)
      {
        // each Particle of i with each particle of j
        for(unsigned int i_points = 0;
                         i_points < dynamic_objects[i]->points.size();
                         i_points++)
        { // for begin on i_points
          for(unsigned int j_points = 0;
                           j_points < dynamic_objects[j]->points.size();
                           j_points++)
          { // for begins on j_points
            // apply gravity forces here
            vec d_ji = dynamic_objects[i]->points[i_points]->pos - dynamic_objects[j]->points[j_points]->pos;
            double length = norm(d_ji);
            double gf = 0;
            if(length > 0.1)
              gf = G * dynamic_objects[i]->points[i_points]->m * dynamic_objects[j]->points[j_points]->m / (length * length);
            vec n_ji = d_ji / length;
            vec gF_ji = gf * n_ji;
            dynamic_objects[i]->points[i_points]->F_list[0] += -gF_ji;
            dynamic_objects[j]->points[j_points]->F_list[0] += gF_ji;
            int charging = dynamic_objects[i]->points[i_points]->charge * dynamic_objects[j]->points[j_points]->charge;
            /*if( charging != 0)
            {
              if( charging == 1) // pushing away
              {
                dynamic_objects[i]->points[i_points]->F_list[1] += n_ji * 20;
                dynamic_objects[j]->points[j_points]->F_list[1] -= n_ji * 20;
              }
              else // pulling
              {
                dynamic_objects[i]->points[i_points]->F_list[1] -= n_ji * 20;
                dynamic_objects[j]->points[j_points]->F_list[1] += n_ji * 20;
              }
            }*/
          } // end for on j_points
        } //end for on i_point
      } // end for on j

      // ==

      for(unsigned int j=0; j<static_objects.size(); j++)
      {
        // each Particle of i with each particle of j
        for(unsigned int i_points = 0;
                         i_points < dynamic_objects[i]->points.size();
                         i_points++)
        { // for begin on i_points
          for(unsigned int j_points = 0;
                           j_points < static_objects[j]->points.size();
                           j_points++)
          { // for begins on j_points
            // apply gravity forces here
            vec d_ji = dynamic_objects[i]->points[i_points]->pos - static_objects[j]->points[j_points]->pos;
            double length = norm(d_ji);
            double gf = 0;
            if(length > 0.01)
              gf = G * dynamic_objects[i]->points[i_points]->m * static_objects[j]->points[j_points]->m / (length * length);
            vec n_ji = d_ji / length;
            vec gF_ji = gf * n_ji;
            dynamic_objects[i]->points[i_points]->F_list[0] += -gF_ji;
            //dynamic_objects[j]->points[j_points]->F_list[0] += gF_ji;
            //int charging = dynamic_objects[i]->points[i_points]->charge * dynamic_objects[j]->points[j_points]->charge;
            /*if( charging != 0)
            {
              if( charging == 1) // pushing away
              {
                dynamic_objects[i]->points[i_points]->F_list[1] += n_ji * 20;
                dynamic_objects[j]->points[j_points]->F_list[1] -= n_ji * 20;
              }
              else // pulling
              {
                dynamic_objects[i]->points[i_points]->F_list[1] -= n_ji * 20;
                dynamic_objects[j]->points[j_points]->F_list[1] += n_ji * 20;
              }
            }*/
          } // end for on j_points
        } //end for on i_point
      } // end for on j

    } //end for on i for other dynamic objects

    for(auto itr = dynamic_objects.begin(); itr!= dynamic_objects.end(); itr++)
    {
      (*itr)->apply_force();
    }
    return;
  }
// =============================================================================
  void Engin::rotate_dynamic()
  {
    double angle = 0.1;  // Rotation angle in degrees
    double radians = datum::pi * angle / 180.0;  // Convert angle to radians

    mat r = {{cos(radians), 0.0, sin(radians)},
               {0.0, 1.0, 0.0},
               {-sin(radians), 0.0, cos(radians)}};
    for(auto itr = dynamic_objects.begin(); itr!= dynamic_objects.end(); itr++)
    {
      for(auto it = (*itr)->points.begin(); it != (*itr)->points.end(); it++)
      {
        (*it)->pos = r * (*it)->pos;
      }
    }
    return;
  }
// =============================================================================
  void Engin::apply_force_null()
  {
    for(auto itr = dynamic_objects.begin(); itr!= dynamic_objects.end(); itr++)
    {
      (*itr)->apply_force();
    }
    return;
  }
// =============================================================================
  void Engin::apply_force_face()
  {
    for(unsigned int i = 0; i<dynamic_objects.size(); i++)
    {
      // looking at other dynamic objects
      for(unsigned int j=i+1; j<dynamic_objects.size(); j++)
      {
        // each Particle of i with each particle of j
        for(unsigned int i_points = 0;
                         i_points < dynamic_objects[i]->points.size();
                         i_points++)
        { // for begin on i_points
          for(unsigned int j_points = 0;
                           j_points < dynamic_objects[j]->points.size();
                           j_points++)
          { // for begins on j_points
            // apply gravity forces here
            vec d_ji = dynamic_objects[i]->points[i_points]->pos - dynamic_objects[j]->points[j_points]->pos;
            double length = norm(d_ji);
            double gf = 0;
            if(length > 0.1)
              gf = G * dynamic_objects[i]->points[i_points]->m * dynamic_objects[j]->points[j_points]->m / (length * length);
            vec n_ji = d_ji / length;
            vec gF_ji = gf * n_ji;
            dynamic_objects[i]->points[i_points]->F_list[0] += -gF_ji;
            dynamic_objects[j]->points[j_points]->F_list[0] += gF_ji;
            //int charging = dynamic_objects[i]->points[i_points]->charge * dynamic_objects[j]->points[j_points]->charge;
            /*if( charging != 0)
            {
              if( charging == 1) // pushing away
              {
                dynamic_objects[i]->points[i_points]->F_list[1] += n_ji * 20;
                dynamic_objects[j]->points[j_points]->F_list[1] -= n_ji * 20;
              }
              else // pulling
              {
                dynamic_objects[i]->points[i_points]->F_list[1] -= n_ji * 20;
                dynamic_objects[j]->points[j_points]->F_list[1] += n_ji * 20;
              }
            }*/
          } // end for on j_points
        } //end for on i_point
      } // end for on j

      // ==

      for(unsigned int j=0; j<static_objects.size(); j++)
      {
        // each Particle of i with each particle of j
        for(unsigned int i_points = 0;
                         i_points < dynamic_objects[i]->points.size();
                         i_points++)
        { // for begin on i_points
          for(unsigned int j_faces = 0;
                           j_faces < static_objects[j]->faces.size();
                           j_faces++)
          { // for begins on j_faces
            // apply gravity forces here

            // first check if the disance to the face is not zeros

            double dist = Util::distance_VF(*dynamic_objects[i]->points[i_points]
                                           ,static_objects[j]->faces[j_faces]);

            // cif the dist !=0 , then we can apply some force to it
            if(dist > 1e-6)
            {
              // x point
              vec dx_ji = dynamic_objects[i]->points[i_points]->pos -
                           static_objects[j]->faces[j_faces].x->pos;
              double length_x = norm(dx_ji);
              double gf = 0;
              if(length_x > distance_tresh_hold)
                gf = G * dynamic_objects[i]->points[i_points]->m *
                       static_objects[j]->faces[j_faces].x->m / (length_x * length_x);
              vec nx_ji = dx_ji / length_x;
              vec gFx_ji = gf * nx_ji / static_objects[j]->faces[j_faces].x->part;
              dynamic_objects[i]->points[i_points]->F_list[0] += -gFx_ji;

              // y point
              vec dy_ji = dynamic_objects[i]->points[i_points]->pos -
                           static_objects[j]->faces[j_faces].y->pos;
              double length_y = norm(dy_ji);
              gf = 0;
              if(length_y > distance_tresh_hold)
                gf = G * dynamic_objects[i]->points[i_points]->m *
                       static_objects[j]->faces[j_faces].y->m / (length_y * length_y);
              vec ny_ji = dy_ji / length_y;
              vec gFy_ji = gf * ny_ji / static_objects[j]->faces[j_faces].y->part;
              dynamic_objects[i]->points[i_points]->F_list[0] += -gFy_ji;

              // z point
              vec dz_ji = dynamic_objects[i]->points[i_points]->pos -
                           static_objects[j]->faces[j_faces].z->pos;
              double length_z = norm(dz_ji);
              gf = 0;
              if(length_z > distance_tresh_hold)
                gf = G * dynamic_objects[i]->points[i_points]->m *
                       static_objects[j]->faces[j_faces].z->m / (length_z * length_z);
              vec nz_ji = dz_ji / length_z;
              vec gFz_ji = gf * nz_ji / static_objects[j]->faces[j_faces].z->part;
              dynamic_objects[i]->points[i_points]->F_list[0] += -gFz_ji;
            }

            /*if(dist < distance_tresh_hold_contanct)
            {
              vec P = dynamic_objects[i]->points[i_points]->pos;
              vec A = static_objects[j]->faces[j_faces].x->pos;
              vec B = static_objects[j]->faces[j_faces].y->pos;
              vec C = static_objects[j]->faces[j_faces].z->pos;

              vec AP = P-A;
              vec BP = P-B;
              vec CP = P-C;
              cout<<"we had a bloody contact!"<<endl;
              double distance_potential = -pow(dist - distance_tresh_hold, 2) * log(dist / distance_tresh_hold);
              double distance_potential_2 = -log(dist / distance_tresh_hold);
              vec der = -log(norm(AP)) * AP + -log(norm(BP)) * BP + -log(norm(CP)) * CP;
              der = 100 * der / norm(der);
              der = distance_potential_2 * der;
              dynamic_objects[i]->points[i_points]->F_list[1] += der;
              cout<<der<<"--"<<endl;
              cout<<dynamic_objects[i]->points[i_points]->F_list[0]<<endl;
            }*/

          } // end for on j_points
        } //end for on i_point
      } // end for on j

    } //end for on i for other dynamic objects

    for(auto itr = dynamic_objects.begin(); itr!= dynamic_objects.end(); itr++)
    {
      (*itr)->apply_force();
    }
    return;
  }
// =============================================================================
  double Engin::collision(vec Q1, vec Q2)
  {
    double toi = 2;

    for(unsigned int i = 0; i<dynamic_objects.size(); i++)
    {

      for(unsigned int j=0; j<static_objects.size(); j++)
      {
        // each Particle of i with each particle of j

        for(unsigned int i_points = 0;
                         i_points < dynamic_objects[i]->points.size();
                         i_points++)
        { // for begin on i_points
          for(unsigned int j_faces = 0;
                           j_faces < static_objects[j]->faces.size();
                           j_faces++)
          { // for begins on j_faces

              vec pos1 = dynamic_objects[i]->points[i_points]->J * Q1.subvec(i, i+11);
              vec pos2 = dynamic_objects[i]->points[i_points]->J * Q2.subvec(i, i+11);

              double t;
              if(Util::CCD_VF(pos1, pos2, static_objects[j]->faces[j_faces], t))
              {
                if(t<toi)
                {
                 toi = t;
                 cout<<"we have an impact"<<endl;
                }
              }
          } // end for on j_points
        } //end for on i_point
      } // end for on j
    } //end for on i for other dynamic objects
    return toi;
  }
// =============================================================================
  double Engin::collision_particles(vec pos1, vec pos2)
  {
    double toi = 2;
    int j = 0;
    for(unsigned int j_faces = 0;
                     j_faces < static_objects[j]->faces.size();
                     j_faces++)
    { // for begins on j_faces
      double t;
      if(Util::CCD_VF(pos1, pos2, static_objects[j]->faces[j_faces], t))
      {
        if(t<toi)
         toi = t;
      }
    } // end for on j_points
    return toi;
  }
  double Engin::collision_particles_dynamic(vec pos1, vec pos2)
  {
    double toi = 2;
    int j = 0;
    for(unsigned int j_faces = 0;
                     j_faces < dynamic_objects[j]->faces.size();
                     j_faces++)
    { // for begins on j_faces
      double t;
      if(Util::CCD_VF(pos1, pos2, dynamic_objects[j]->faces[j_faces], t))
      {
        if(t<toi)
         toi = t;
      }
    } // end for on j_points
    return toi;
  }
}
