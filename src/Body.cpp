// openGL libraries
#include <GL/gl.h>
#include <GLFW/glfw3.h>
#include <GL/glut.h>

// my libraries
#include "../include/Body.h"

// c++ libraries
#include <vector>

// namespaces
using namespace std;
using namespace arma;

namespace ABD{

  Body::Body(bool move, SHAPE sh, vec cntr, double r, double weight)
  {
    //load_from_file(SHAPE sh);
    radius = r;
    W = weight;
    this->sh = sh;
    switch (sh)
    {
      case SPHERE:
      {
        load_from_file_sphere(cntr);
        init(move);
        break;
      }
      case PLANE:
      {
        cout<<"we are making a plane"<<endl;
        load_from_file_plane(cntr);
        init(move);
        break;
      }
      case PLANE_STIP:
      {
        cout<<"we are making a stip plane"<<endl;
        load_from_file_plane_stip(cntr, 30);
        init(move);
        break;
      }
      case PLANE_CURVED:
      {
        cout<<"we are making a plane curved"<<endl;
        load_from_file_plane_curved(cntr);
        init(move);
        break;
      }
      case MOB:
      {
        cout<<"we are making mobius"<<endl;
        load_from_file_MOB(cntr);
        init(move);
        break;
      }
    }

  }
  void Body::init(bool move)
  {
    cout<<"starting the init with: "<<move<<endl;
    movable = move;
    if(movable)
      set_up_dynamic();
    k = 1000;
    v = 1;
    Dl_t = 1./32.;
    Dl_t_r = 24;
    energy_tresh_hold = 0.05;
    cout<<"BODY is ready"<<endl;
  }


  void Body::set_up_dynamic()
  {
    M = Mat<double>(12,12, fill::zeros);
    cout<<"we are here"<<endl;
    for(long unsigned int i=0; i<points.size(); i++)
    {
      M = M + (points[i]->m) * points[i]->J_T * points[i]->J;
    }
    cout<<M<<endl;
    if(movable)
      M_inv = inv(M,inv_opts::no_ugly	);

    A = Mat<double>(3,3, fill::eye);
    p = vec(3, fill::zeros);

    //A(1,1) = sqrt(3)/2;
    //A(2,1) = 0.5;

    //A(1,2) = -0.5;;
    //A(2,2) = A(1,1);

    q = vec(12, fill::zeros);
    q_dot = vec(12, fill::zeros);
    q_mad = vec(12, fill::zeros);

    q(0) = p(0);
    q(1) = p(1);
    q(2) = p(2);

    // first column
    q(3) = A(0,0);
    q(4) = A(1,0);
    q(5) = A(2,0);

    // second column
    q(6) = A(0,1);
    q(7) = A(1,1);
    q(8) = A(2,1);

    // third column
    q(9) = A(0,2);
    q(10) = A(1,2);
    q(11) = A(2,2);

  }

  void Body::update_A_p()
  {
    p(0) = q(0);
    p(1) = q(1);
    p(2) = q(2);

    A(0,0) = q(3);
    A(1,0) = q(4);
    A(2,0) = q(5);

    A(0,1) = q(6);
    A(1,1) = q(7);
    A(2,1) = q(8);

    A(0,2) = q(9);
    A(1,2) = q(10);
    A(2,2) = q(11);
  }


  void Body::update_q_mad()
  {
    vec F_tot(12, fill::zeros);
    for(auto itr=points.begin(); itr!=points.end(); itr++)
    {
      F_tot = F_tot + (*itr)->J_T * (*itr)->F;
    }

    q_mad = q + (Dl_t * q_dot) + (Dl_t * Dl_t * M_inv * F_tot);

  }

  vec Body::E_der(const Mat<double>& X)
  {
    // Note, we have to update q_mad whenever we want to use the derivation and hessian
    vec e_der = 1 * (X - q_mad);

    vec a[3]={vec(3, fill::zeros), vec(3, fill::zeros), vec(3, fill::zeros)};
    a[0](0) = X[3];
    a[0](1) = X[4];
    a[0](2) = X[5];

    a[1](0) = X[6];
    a[1](1) = X[7];
    a[1](2) = X[8];

    a[2](0) = X[9];
    a[2](1) = X[10];
    a[2](2) = X[11];

    vec b[3]={vec(3, fill::zeros), vec(3, fill::zeros), vec(3, fill::zeros)};
    b[0] = 4 * (dot(a[0], a[0]) - 1) * a[0] + 4 * dot(a[0], a[1]) * a[1] + 4 * dot(a[0], a[2]) * a[2];

    b[1] = 4 * (dot(a[1], a[1]) - 1) * a[1] + 4 * dot(a[1], a[2]) * a[2] + 4 * dot(a[1], a[0]) * a[0];

    b[2] = 4 * (dot(a[2], a[2]) - 1) * a[2] + 4 * dot(a[2], a[0]) * a[0] + 4 * dot(a[2], a[1]) * a[1];

    e_der[3] = e_der[3] + k * v * b[0](0);
    e_der[4] = e_der[4] + k * v * b[0](1);
    e_der[5] = e_der[5] + k * v * b[0](2);

    e_der[6] = e_der[6] + k * v * b[1](0);
    e_der[7] = e_der[7] + k * v * b[1](1);
    e_der[8] = e_der[8] + k * v * b[1](2);

    e_der[9] = e_der[9] + k * v * b[2](0);
    e_der[10] = e_der[10] + k * v * b[2](1);
    e_der[11] = e_der[11] + k * v * b[2](2);

    //cout<<e_der<<endl;


    return e_der;
  }
  Mat<double> Body::E_der_der(const Mat<double>& X)
  {
    Mat<double> e_der_der(12, 12, fill::eye);

    Mat<double> id(3,3, fill::eye);

    // pay attention to this detil: update the q_mad before using this function

    e_der_der = 2 * e_der_der;

    Mat<double> a[3]={Mat<double>(3,1, fill::zeros), Mat<double>(3,1, fill::zeros), Mat<double>(3,1, fill::zeros)};
    a[0](0) = X[3];
    a[0](1) = X[4];
    a[0](2) = X[5];

    a[1](0) = X[6];
    a[1](1) = X[7];
    a[1](2) = X[8];

    a[2](0) = X[9];
    a[2](1) = X[10];
    a[2](2) = X[11];

    Mat<double> b[3][3];
    //b[0][0] = Mat<double>(3,3, fill::zeros);

    b[0][0] = 4 * (2 * a[0] * a[0].t() + 1 * id * (dot(a[0], a[0]) - 1)) + 4 * (a[1] * a[1].t() + a[2] * a[2].t());
    b[1][1] = 4 * (2 * a[1] * a[1].t() + 1 * id * (dot(a[1], a[1]) - 1)) + 4 * (a[0] * a[0].t() + a[2] * a[2].t());
    b[2][2] = 4 * (2 * a[2] * a[2].t() + 1 * id * (dot(a[2], a[2]) - 1)) + 4 * (a[0] * a[0].t() + a[1] * a[1].t());

    b[0][1] = 4 * id * dot(a[0], a[1]) + 4 * a[1] * a[0].t();
    b[0][2] = 4 * id * dot(a[0], a[2]) + 4 * a[2] * a[0].t();

    b[1][0] = 4 * id * dot(a[1], a[0]) + 4 * a[0] * a[1].t();
    b[1][2] = 4 * id * dot(a[1], a[2]) + 4 * a[2] * a[1].t();

    b[2][0] = 4 * id * dot(a[2], a[0]) + 4 * a[0] * a[2].t();
    b[2][1] = 4 * id * dot(a[2], a[1]) + 4 * a[1] * a[2].t();

    for(int i=0; i<9; i+=3)
    {
      for(int j=0; j<9; j+=3)
      {
        for(int ii=0; ii<3; ii++)
        {
          for(int jj=0; jj<3; jj++)
          {
            e_der_der(i+ii+3, j+jj+3) = e_der_der(i+ii+3, j+jj+3) + k * v * b[i/3][j/3](ii, jj);
            //cout<<i+ii+3<<","<< j+jj+3<< "=" << "["<<i/3<<","<<j/3<<"]"<<"("<<ii<<","<<jj<<")"<<endl;
          }
        }
      }
    }

    return e_der_der;
  }

  Mat<double> Body::PSPD(const Mat<double>& X)
  {
    vec eigval;
  	mat eigvec;

  	eig_sym(eigval, eigvec, X);
    cout<<"the init is:"<<endl<<eigval<<endl;
    for(unsigned int i=0;i<eigval.size();i++)
  	{
  		if(eigval(i)<0)
  			eigval(i) = 0.000001;
  	}

  	mat C = eigvec * diagmat(eigval) * eigvec.t();
    return C;
  }
  double Body::calculate_energy(const Mat<double>& X)
  {
    Mat<double> s = X - q_mad;
    Mat<double> B(3, 3, fill::zeros);
    Mat<double> id(3,3, fill::eye);

    B(0,0) = X(3);
    B(1,0) = X(4);
    B(2,0) = X(5);

    B(0,1) = X(6);
    B(1,1) = X(7);
    B(2,1) = X(8);

    B(0,2) = X(9);
    B(1,2) = X(10);
    B(2,2) = X(11);

    double e = 0;

    Mat<double> stiff_poten = (B.t() * B - id);

    double stiff_pot_value = 1 * k * v * accu(stiff_poten % stiff_poten);
    e = 0.5 * accu(s%s) + stiff_pot_value;

    return e;
  }

  void Body::apply_force()
  {
    vec z(3, fill::zeros);
    for(unsigned long int i=0;i<points.size();i++)
    {
      //points[i]->F_list[1] = -points[i]->J * q_dot * points[i]->fr_value * points[i]->m;
      points[i]->F = 1 * points[i]->F_list[0] + 1 * points[i]->F_list[1] + points[i]->F_list[2];
      //cout<<points[i]->F<<endl;
      //int x;
      //cin>>x;
      points[i]->F_list[0] = z;
      points[i]->F_list[1] = z;
      points[i]->F_list[2] = z;
      //points[i]->F = z;
    }
    /*double x = points[36]->pos(0);
    double y = points[36]->pos(1);
    vec force(3, fill::zeros);
    double theta = atan(y/x);
    force(0) = -  radius * sin(theta);
    force(1) =   radius * cos(theta);
    cout<<force<<endl;

    points[36]->F = force;*/
  }

  void Body::apply_tranformation()
  {
    for(auto itr=points.begin(); itr!=points.end(); itr++)
    {
      (*itr)->pos = (*itr)->J * q ;
    }
  }




  //==============================================================================

  void Body::load_from_file_sphere(vec cntr)
  {
    ifstream reader("Generator/topologycal_points_sphere_small.txt");

    int N_p;
    reader>>N_p;
    vec Z(3, fill::zeros);

    points.resize(N_p);
    radius = 2;

    for(int i=0; i< N_p; i++)
    {
      int num;
      double x, y, z;
      reader>> num >> x >> y >> z;
      vec pos(3, fill::zeros);
      pos(0) = x;
      pos(1) = y;
      pos(2) = z;

      double l = norm(pos, 2);
      pos = (radius/l)*pos;
      pos = pos + cntr;
      points[i] = new Particle(pos, Z, W/N_p, Z);
      //points[i]->b = this;
    }
    center = new Particle(cntr, Z, 1, Z);
    int N_s;
    reader>>N_s;

    random_device rd;   // Obtain a random seed from the hardware
  	mt19937 gen(rd());  // Seed the random number generator
  	uniform_int_distribution<int> dis(0, 255);  // Define the range

    for(int i=0; i< N_s; i++)
    {
      int x, y;
      reader>> x >> y;
      Segment temp(points[x], points[y]);
      temp.R = dis(gen);
      temp.G = dis(gen);
      temp.B = dis(gen);
      cout<<x<<" "<<y<<endl;
      segments.push_back(temp);
    }

    int N_f;
    reader>>N_f;

    for(int i=0; i< N_f; i++)
    {
      int x, y, z;
      reader>> x >> y >> z;

      Face temp(points[x], points[y], points[z]);
      faces.push_back(temp);
    }
    cout<<"we are done"<<endl;
    reader.close();
  }

  void Body::load_from_file_MOB(vec cntr)
  {
    ifstream reader("/localhome/akazemin/Desktop/OptEnvironment/Generator/mob_saved_2.txt");

    int N_p;
    reader>>N_p;
    vec Z(3, fill::zeros);

    points.resize(N_p);
    //radius = 2;

    for(int i=0; i< N_p; i++)
    {
      int num;
      double x, y, z;
      reader>> num >> x >> y >> z;
      vec pos(3, fill::zeros);
      pos(0) = x;
      pos(1) = y;
      pos(2) = z;
      mat r_x(3,3, fill::eye);
      mat r_y(3,3, fill::eye);
      double angle = 50.0;  // Rotation angle in degrees
      double radians = datum::pi * angle / 180.0;  // Convert angle to radians

      mat r_z = {{cos(radians), -sin(radians), 0.0},
                 {sin(radians), cos(radians), 0.0},
                 {0.0, 0.0, 1.0}};





      pos = r_z * (pos)  + cntr;
      points[i] = new Particle(pos * 2, Z, W/N_p, Z);
      //points[i]->b = this;
    }
    center = new Particle(cntr, Z, 1, Z);
    int N_s;
    reader>>N_s;

    for(int i=0; i< N_s; i++)
    {
      int x, y;
      reader>> x >> y;
      Segment temp(points[x], points[y]);
      cout<<x<<" "<<y<<endl;
      segments.push_back(temp);
    }

    int N_f;
    reader>>N_f;

    for(int i=0; i< N_f; i++)
    {
      int x, y, z;
      reader>> x >> y >> z;
      points[x]->part++;
      points[y]->part++;
      points[z]->part++;
      Face temp(points[x], points[y], points[z]);
      faces.push_back(temp);
    }
    cout<<"we are done"<<endl;
    reader.close();
  }

  void Body::load_from_file_plane(vec cntr)
  {
    ifstream reader("Generator/topologycal_points_plane_small.txt");
    cout<<"start reading from file"<<endl;
    int N_p;
    reader>>N_p;
    vec Z(3, fill::zeros);

    points.resize(N_p);

    for(int i=0; i< N_p; i++)
    {
      int num;
      double x, y, z;
      reader>> num >> x >> y >> z;
      vec pos(3, fill::zeros);
      pos(0) = x;
      pos(1) = y;
      pos(2) = z;

      vec t(3, fill::zeros);
      t(0) = 0;
      t(1) = 0;
      t(2) = 0;

      //double l = norm(pos, 2);
      //pos = (10./l)*pos;
      pos = pos + cntr + t;
      //pos = pos * 10;
      points[i] = new Particle(pos, Z, 1, Z);
    }

    int N_s;
    reader>>N_s;
    for(int i=0; i< N_s; i++)
    {
      int x, y;
      reader>> x >> y;
      Segment temp(points[x], points[y]);

      segments.push_back(temp);
    }

    int N_f;
    reader>>N_f;
    for(int i=0; i< N_f; i++)
    {
      int x, y, z;
      reader>> x >> y >> z;
      Face temp(points[x], points[y], points[z]);

      faces.push_back(temp);
    }

    reader.close();

  }
  void Body::load_from_file_plane_curved(vec cntr)
  {
    ifstream reader("Generator/topologycal_points_plane_medium.txt");
    cout<<"start reading from file"<<endl;
    int N_p;
    reader>>N_p;
    vec Z(3, fill::zeros);

    points.resize(N_p);
    double r = 30;
    for(int i=0; i< N_p; i++)
    {
      int num;
      double x, y, z;
      reader>> num >> x >> y >> z;
      vec pos(3, fill::zeros);
      pos(0) = 1 * (x);
      pos(2) = 1 * (z);// - sqrt(r*r - x*x - y*y);
      pos(1) = - sqrt(r*r - pos(0)*pos(0) - pos(2)*pos(2)) + r;

      vec t(3, fill::zeros);


      //double l = norm(pos, 2);
      //pos = (10./l)*pos;
      pos = pos + cntr ;
      //pos = pos * 10;
      points[i] = new Particle(pos, Z, 1, Z);
    }

    int N_s;
    reader>>N_s;
    for(int i=0; i< N_s; i++)
    {
      int x, y;
      reader>> x >> y;
      Segment temp(points[x], points[y]);

      segments.push_back(temp);
    }

    int N_f;
    reader>>N_f;
    for(int i=0; i< N_f; i++)
    {
      int x, y, z;
      reader>> x >> y >> z;
      Face temp(points[x], points[y], points[z]);

      faces.push_back(temp);
    }

    reader.close();

  }
  void Body::load_from_file_plane_stip(vec cntr, double m)
  {
    ifstream reader("Generator/topologycal_points_plane_medium.txt");
    cout<<"start reading from file medium"<<endl;
    int N_p;
    reader>>N_p;
    vec Z(3, fill::zeros);

    points.resize(N_p);
    mat rot(3,3, fill::zeros);
    rot(0,0) = cos(m * M_PI / 180);
    rot(1,1) = rot(0,0);
    rot(2,2) = 1;
    rot(0,1) = -sin(m * M_PI / 180);
    rot(1,0) = -rot(0,1);

    for(int i=0; i< N_p; i++)
    {
      int num;
      double x, y, z;
      reader>> num >> x >> y >> z;
      vec pos(3, fill::zeros);
      pos(0) = 2 * x;
      pos(1) = y;
      pos(2) = 0.5 * z;
      pos = rot * pos;
      pos = pos + cntr;
      //pos = pos * 10;
      points[i] = new Particle(pos, Z, 1, Z);
      cout<<"particel number: "<<i<<" is done"<<endl;
    }

    int N_s;
    reader>>N_s;
    for(int i=0; i< N_s; i++)
    {
      int x, y;
      reader>> x >> y;
      Segment temp(points[x], points[y]);

      segments.push_back(temp);
    }

    int N_f;
    reader>>N_f;
    for(int i=0; i< N_f; i++)
    {
      int x, y, z;
      reader>> x >> y >> z;
      Face temp(points[x], points[y], points[z]);

      faces.push_back(temp);
    }

    reader.close();

  }

  void Body::reshape()
  {

    for(auto itr = points.begin(); itr!= points.end(); itr++)
    {
      double l = norm((*itr)->pos_0, 2);
      (*itr)->pos_0 = (10./l) * (*itr)->pos_0;
      (*itr)->pos = (*itr)->pos_0;
      //*itr).construct_J();
    }
  }

  void Body::draw_body()
  {
    for(long unsigned int i=0; i<segments.size(); i++)
    {
      if(sh == PLANE || sh == PLANE_STIP || sh == PLANE_CURVED)
        glColor3ub(150, 75, 0);
      if(sh == SPHERE)
        glColor3ub(segments[i].R, segments[i].G, segments[i].B);
      else
        glColor3ub(100, 100, 100);
      glBegin(GL_LINES);
      glVertex3f(this->segments[i].x->pos[0], this->segments[i].x->pos[1], this->segments[i].x->pos[2]);
      glVertex3f(this->segments[i].y->pos[0], this->segments[i].y->pos[1], this->segments[i].y->pos[2]);
      glEnd();
    }

    /*for(long unsigned int i=0; i<faces.size(); i++)
    {
      if(sh == PLANE || sh == PLANE_STIP || sh == PLANE_CURVED)
        glColor3ub(150, 75, 0);
      if(sh == SPHERE)
        glColor3ub(segments[i].R, segments[i].G, segments[i].B);
      else
        glColor3ub(100, 100, 100);
      glBegin(GL_LINES);
      glVertex3f(this->faces[i].x->pos[0], this->faces[i].x->pos[1], this->faces[i].x->pos[2]);
      glVertex3f(this->faces[i].y->pos[0], this->faces[i].y->pos[1], this->faces[i].y->pos[2]);
      glEnd();
      glBegin(GL_LINES);
      glVertex3f(this->faces[i].x->pos[0], this->faces[i].x->pos[1], this->faces[i].x->pos[2]);
      glVertex3f(this->faces[i].z->pos[0], this->faces[i].z->pos[1], this->faces[i].z->pos[2]);
      glEnd();
      glBegin(GL_LINES);
      glVertex3f(this->faces[i].z->pos[0], this->faces[i].z->pos[1], this->faces[i].z->pos[2]);
      glVertex3f(this->faces[i].y->pos[0], this->faces[i].y->pos[1], this->faces[i].y->pos[2]);
      glEnd();
    }*/

  }

  Body& Body::operator= (const Body& obj)
  {
    this->points = obj.points;
    this->segments = obj.segments;
    this->faces = obj.faces;

    return *this;
  }

  void Body::print()
  {
    cout<<"_________________________________________________________"<<endl;
    cout<<"________________BODY_____________________________________"<<endl;
    cout<<"000000000000000000"<<endl;
    for(long unsigned int i=0;i<points.size();i++)
    {
      cout<<*points[i]<<endl;
    }
    cout<<"%%%%%%%%%%%%%%%%%%"<<endl;
    for(long unsigned int i=0;i<segments.size();i++)
    {
      cout<<segments[i]<<endl;
    }
    cout<<"/_\\/_\\/_\\/_\\/_\\/_\\"<<endl;
    for(long unsigned int i=0;i<faces.size();i++)
    {
      cout<<faces[i]<<endl;
    }
    cout<<"________________BODY_____________________________________"<<endl;
    cout<<"---------------------------------------------------------"<<endl;

  }
}
