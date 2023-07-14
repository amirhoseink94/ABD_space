#ifndef ENGIN_H
#define ENGIN_H

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


namespace ABD{
	class Engin
	{
		public:
		// All objects, as the body type, that would be nice :D
		vector<Body*> dynamic_objects;
		vector<Body*> static_objects;
		unsigned int N;
		vec Q;												// stacked Q =(q_1, q_2,..., q_\ell)
		vec dQ;												// derivitive of Q
		mat H;												// hessian mattix

		double Dl_t;									// time step, TODO: moving it to a better plave
		double G;

		double big_mass;
		double sp_value;

		double energy_tresh_hold;
		double distance_tresh_hold;
		double distance_tresh_hold_contanct;



		// initialization phase needs these two:
		void init();										// the initilization,

		void add_object(Body*);

		// First step:

		void apply_force();							// apply external forces to each object, for now, it is just the particles applyforce

		void apply_force_face();
		void apply_force_dynamic();
		void apply_force_null();

		vec calculate_next_Q();

		// Thirs spte: moving
		void update_A_p();							// update the affine trnasformation A and translation p based on q, it seems we do not need this
		void apply_tranformation();

		void apply_tranformation_temp(const vec&);

		double calculate_energy(const vec&);
		double calculate_energy_collision(const vec&, vec&);

		double collision(vec, vec);
		double collision_top(vec, vec);
		double collision_particles(vec, vec);
		double collision_particles_dynamic(vec, vec);

		// tools
		mat PSPD(const mat&);					// normalizing the hessian matrix to be semipositive
		void rotate_dynamic();


		// extracing field
		void init_particles(int);
		void render_field();
		void apply_force_to_particles();
		void apply_force_to_particles_dynamic();
		void render_particles();
		void draw_field();
		//void collision_particle(vec, vec);
		private:
		vector<Particle*> particles;

	};


	static double distance_VF(Particle& P, Face& F);
}
#endif
