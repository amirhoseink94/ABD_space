// openGL libraries
#include <GL/gl.h>
#include <GLFW/glfw3.h>
#include <GL/glut.h>

// armadillo library
#include <armadillo>

// system and c++ libraries
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <time.h>
#include <vector>
#include <random>

// my libraries
#include "include/Particle.h"
#include "include/Segment.h"
#include "include/Face.h"
#include "include/Engin.h"
#include "include/Body.h"
#include "include/Util.h"

// opencv
#include <opencv2/opencv.hpp>

//namespaces
using namespace std;
using namespace arma;
using namespace ABD;
using namespace cv;



// --------------------------------- global variables --------------------------------

int widthp = 1024, height = 924;   // Size of the drawing area, to be set in reshape().
int frameNumber = 0;     // For use in animation.
int rotate = 0;

int zoom = 45.;
int scr_width = 1024;
int scr_height = 924;

int x_camera = 0.;
int y_camera = 0.;
int z_camera = 0.;

int x_rotate_angle = 0;
int y_rotate_angle = 0;

int azimuth = -300.;
int elevation = 0.;

bool right_button_flag = false;
bool left_button_flag = false;
int oldX = 0.;
int oldY = 0.;

// Animations settings
int frame_no = 0;
double dl_t = 1/2.;
double G = 0.0001;
double tresh_hold = 5;
// key handlers
bool UP_key_flag = false;
bool DOWN_key_flag = false;
bool RIGHT_key_flag = false;
bool LEFT_key_flag = false;

bool n_key_flag = false;
bool m_key_flag = false;
bool w_key_flag = false;
bool s_key_flag = false;
bool a_key_flag = false;
bool d_key_flag = false;
bool z_key_flag = false;
bool p_key_flag = false;
bool q_key_flag = false;

// functions for key handling
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
void key_handler();

// -----------------------------------------------------------------------------------

// rednders and drawing
void render();
void render_particles(vector<Particle>&);
void draw_cube();
void draw_space();
void draw_particles(vector<Particle>);

// tool functions
double distance(vec, vec);
double vec_length(const vec&);

void random_force(Body&);

// -----------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------
int main( int argc, char** argv)
{
	vec z(3, fill::zeros);
	vec pos(3, fill::zeros);
	pos[0] = 0;
	pos[1] = 10;
	pos[2] = -5;
	vec pos1(3, fill::zeros);
	pos1[0] = 0;
	pos1[1] = 0;
	pos1[2] = 0;
	vec pos2(3, fill::zeros);
	pos2[0] = 0;
	pos2[1] = 0;
	pos2[2] = 10;
	vec pos3(3, fill::zeros);
	pos3[0] = 10;
	pos3[1] = 10;
	pos3[2] = 10;

	vec pos4(3, fill::zeros);
	pos4[0] = 20;
	pos4[1] = 20;
	pos4[2] = 20;

	/*Particle p1(pos, z, 1, z);
	Particle p2(pos1, z, 1, z);
	Particle p3(pos2, z, 1, z);

	Face f(&p1, &p2, &p3);

	double t;

	cout<<Util::CCD_VF(pos4, pos3, f, t)<<endl;
	cout<<t<<endl;

	int xx;
	cin>>xx;
*/


	int up_range = 200;
	random_device rd;   // Obtain a random seed from the hardware
	mt19937 gen(rd());  // Seed the random number generator
	uniform_int_distribution<int> dis(0, up_range);  // Define the range


	Engin eng;
	/*vector<Body*> bodies;
	for(int i=0; i<200; i++)
	{
		vec pos(3, fill::zeros);
		int x = dis(gen) - up_range/2;
		int y = dis(gen) - up_range/2;
		int z = dis(gen) - up_range/2;
		pos[0] = x;
		pos[1] = y;
		pos[2] = z;
		Body *b = new Body(true, ABD::SPHERE, pos, 0.1, 0.1);
		bodies.push_back(b);
		eng.add_object(b);
	}*/
	Body mob(false, ABD::MOB, pos, 5, 100000);
	//random_force(b);
	Body c(true, ABD::SPHERE, pos1, 4, 1);
	//random_force(c);
	//Body d(true, ABD::SPHERE, pos2, 5, 1000);
	//cout<<"two objects are made"<<endl;

	eng.add_object(&mob);
	eng.add_object(&c);
	//eng.add_object(&d);*/

	eng.init();
	eng.init_particles(20);
	int xty;
	cin>>xty;

	if ( !glfwInit() )
		return -1;
	cout << "(the glfw just initilized!)" << endl;
	GLFWwindow* window = glfwCreateWindow(widthp, height, "test", NULL, NULL);
	if( !window )
	{
		glfwTerminate();
		return -1;
	}

	glfwMakeContextCurrent(window);
	glfwSwapInterval(1);

	// keyboard handles
	glfwSetKeyCallback(window, key_callback);

	//the main loop
	render();
	glfwSwapBuffers(window);

	vec Q_dot;

	VideoWriter outputVideo;
	outputVideo.open( "video.avi",  /*Video Name*/
			cv::VideoWriter::fourcc('M','J','P','G'), /* fourcc */
			32.0f,                      /*Yuchen: Frame Rate */
			cv::Size( widthp, height ),  /*Yuchen: Frame Size of the Video  */
			true);


	int frame_counter = 0;
	string frame_address = "image/image_frame_";
	bool action = false;
	bool action_field = false;
	while( !glfwWindowShouldClose( window) )
	{
		cout<<"we start the loop"<<endl;
		glfwPollEvents();
		key_handler();
		render();
		draw_space();
		for(auto itr = eng.dynamic_objects.begin(); itr != eng.dynamic_objects.end(); itr++)
			(*itr)->draw_body();

		for(auto itr = eng.static_objects.begin(); itr != eng.static_objects.end(); itr++)
			(*itr)->draw_body();
		//cout<<"we are here"<<endl;
		eng.draw_field();
		//b.draw_body();
		c.draw_body();
		//d.draw_body();

		if(z_key_flag)
		{
			action = true;
		}
		if(action)
		{
			cout<<"in the loop!"<<endl;

			eng.apply_force_face();
			//eng.apply_force_dynamic();
			//eng.apply_force_null();

			vec Q_next = eng.calculate_next_Q();

			Q_dot = (Q_next - eng.Q)/eng.Dl_t;
			//cout<<"Q_dot is: "<<endl<<Q_dot<<endl;

			int i=0;
			for(auto itr = eng.dynamic_objects.begin(); itr!=eng.dynamic_objects.end(); itr++)
	    {
	      (*itr)->q_dot = Q_dot.subvec(i, i+11);
				i+=12;
	    }

			eng.Q = Q_next;

			eng.apply_tranformation();




			cv::Mat pixels( /* num of rows */ height, /* num of cols */ widthp, CV_8UC3 );
			glReadPixels(0, 0, widthp, height, GL_RGB, GL_UNSIGNED_BYTE, pixels.data );
			cv::Mat cv_pixels( /* num of rows */ height, /* num of cols */ widthp, CV_8UC3 );

			for( int y=0; y<height; y++ ) for( int x=0; x<widthp; x++ )
			{
				cv_pixels.at<cv::Vec3b>(y,x)[2] = pixels.at<cv::Vec3b>(height-y-1,x)[0];
				cv_pixels.at<cv::Vec3b>(y,x)[1] = pixels.at<cv::Vec3b>(height-y-1,x)[1];
				cv_pixels.at<cv::Vec3b>(y,x)[0] = pixels.at<cv::Vec3b>(height-y-1,x)[2];
			}

			string frame_num_string = to_string(frame_counter);
			frame_num_string = string(5 - frame_num_string.length(), '0') + frame_num_string;
			cv::imwrite(frame_address+frame_num_string+".jpg", cv_pixels);
			frame_counter++;
			outputVideo << cv_pixels;
			cout<<"the loop ended"<<endl;
		}
		if(action_field)
		{
			cout<<"in the field loop!"<<endl;

			eng.apply_force_to_particles_dynamic();

			eng.render_field();

			eng.rotate_dynamic();

			cv::Mat pixels( /* num of rows */ height, /* num of cols */ widthp, CV_8UC3 );
			glReadPixels(0, 0, widthp, height, GL_RGB, GL_UNSIGNED_BYTE, pixels.data );
			cv::Mat cv_pixels( /* num of rows */ height, /* num of cols */ widthp, CV_8UC3 );

			for( int y=0; y<height; y++ ) for( int x=0; x<widthp; x++ )
			{
				cv_pixels.at<cv::Vec3b>(y,x)[2] = pixels.at<cv::Vec3b>(height-y-1,x)[0];
				cv_pixels.at<cv::Vec3b>(y,x)[1] = pixels.at<cv::Vec3b>(height-y-1,x)[1];
				cv_pixels.at<cv::Vec3b>(y,x)[0] = pixels.at<cv::Vec3b>(height-y-1,x)[2];
			}

			string frame_num_string = to_string(frame_counter);
			frame_num_string = string(5 - frame_num_string.length(), '0') + frame_num_string;
			cv::imwrite(frame_address+frame_num_string+".jpg", cv_pixels);
			frame_counter++;
			outputVideo << cv_pixels;
			cout<<"the loop ended"<<endl;
		}

		glfwSwapBuffers(window);
	}
	outputVideo.release();
	glfwTerminate();

	return 0;
}

// --------------------------------- rendering --------------------------------

// render camera and scene
void render()
{
    glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

    // reshape display window function, locating camera
    glMatrixMode(GL_PROJECTION); //these apply directy to my view, what I see in the clip!
    glLoadIdentity();
    gluPerspective(zoom, scr_width / scr_height, 1.0, 1000.0);


    // Use global parameters to rotate camera;
    glMatrixMode(GL_MODELVIEW); // this is all about the camera
    glLoadIdentity();
    gluLookAt(0, 0, 195.0, 0.0, 50.0, 0.0, 0.0, 1.0, 0.0);

    glTranslatef(-x_camera, -y_camera, -z_camera);
    glRotatef(x_rotate_angle, 1., 0., 0.);
    glRotatef(y_rotate_angle, 0., 1., 0.);

    // drawing the floor
    glColor3ub(255, 255, 255);

    // draw the bones
    glPushMatrix();
    glColor3ub(0, 255, 255);
    glPopMatrix();
}

// rendering objects at each time
void render_particles(vector<Particle>& particles)
{
	long unsigned int N = particles.size();
	for(long unsigned int i=0; i<N; i++)
	{
		vec F_tot(3, 1, fill::zeros);

		for(long unsigned int j=0; j<N; j++)
		{
			if( i!=j )
			{
				vec F = particles[j].pos - particles[i].pos;
				//cout<<"The force direction is: "<<i<<" "<<F<<endl;
				double len = vec_length(F);
				if(len < tresh_hold)
				{
					//cout<<"warning! we are in the dark age!"<<endl;
					continue;
				}
				F = F*(1./len);
				//cout<<"##"<<i<<len<<" "<<F<<endl;
				F = F*(G * particles[i].m * particles[j].m / pow(len, 2));
				F_tot = F_tot + F;
			}
		}
		//cout<<particles[i].F<<"|"<<F_tot<<endl;
		particles[i].F = F_tot;

	}
	/*particles[1].pos.x = 0;
	particles[1].pos.y = 50;
	particles[1].pos.y = 50;*/
	for(long unsigned int i=0; i<N; i++)
	{
		//Vec3D a =
	//	particles[i].pos = particles[i].pos + dl_t * particles[i].v;
		particles[i].v = particles[i].v + particles[i].F * ((1./particles[i].m) * dl_t);
		particles[i].pos = particles[i].pos + particles[i].v * (dl_t);
		//cout<<"pos: "<<i<<";"<<particles[i].pos<<endl;
		//cout<<"vel: "<<i<<";"<<particles[i].v<<endl;
	}

}

// --------------------------------- drawing objects --------------------------------

/*void draw_space()
{
    glBegin(GL_LINES);
    // x
    glColor3ub(255, 0, 0);
    glVertex3f(0, 0, 0);
    glVertex3f(500, 0, 0);
    // y
    glColor3ub(0, 255, 0);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 500, 0);
    // z
    glColor3ub(0, 0, 255);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0, 500);
    glEnd();
}*/
void draw_space()
{
    glBegin(GL_LINES);
    // x
    //glColor3ub(255, 0, 0);
    //glVertex3f(0, 0, 0);
    //glVertex3f(500, 0, 0);
    // y
    glColor3ub(0, 0, 150);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 500, 0);
    // z
    //glColor3ub(0, 0, 255);
    //glVertex3f(0, 0, 0);
    //glVertex3f(0, 0, 500);
    glEnd();
		for(int i=-500;i<500;i+=20)
		{
			glBegin(GL_LINES);
	    // x
	    glColor3ub(0, 0, 255);
	    glVertex3f(-500, 0, i);
	    glVertex3f(500, 0, i);
	    // z
	    glColor3ub(0, 0, 255);
	    glVertex3f(i, 0, -500);
	    glVertex3f(i, 0, 500);
	    glEnd();
		}
}

void draw_particles(vector<Particle> p)
{
	for(long unsigned int i=0; i<p.size(); i++)
	{
		glPointSize(2);
		glBegin(GL_POINTS);
			glColor3ub((i * 10)%255, 255, 255);
			glVertex3f(p[i].pos[0], p[i].pos[1], p[i].pos[2]);
		glEnd();
	}
}


void draw_cube()
{
	double r = 3;
	for(int i=0; i<=180; i++)
	{
		for(int j=0; j<360; j++)
		{
			glBegin(GL_POINTS);
				glColor3ub(255, i, j);
				glVertex3f(r * sin(i) * cos(j), r * sin(i) * sin(j), r * cos(i));
			glEnd();
		}
	}
}

// --------------------------------- Key handling --------------------------------
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if ( key == GLFW_KEY_Z )
    {
        if ( action == GLFW_PRESS )
            z_key_flag = true;
        else if ( action == GLFW_RELEASE )
            z_key_flag = false;
        else if ( action == GLFW_REPEAT)
            z_key_flag = true;
	}

    if ( key == GLFW_KEY_N )
    {
        if ( action == GLFW_PRESS )
            n_key_flag = true;
        else if ( action == GLFW_RELEASE )
            n_key_flag = false;
        else if ( action == GLFW_REPEAT)
            n_key_flag = true;
	}

	if ( key == GLFW_KEY_M )
    {
        if ( action == GLFW_PRESS )
            m_key_flag = true;
        else if ( action == GLFW_RELEASE )
            m_key_flag = false;
        else if ( action == GLFW_REPEAT)
            m_key_flag = true;
	}

	if ( key == GLFW_KEY_W )
    {
        if ( action == GLFW_PRESS )
            w_key_flag = true;
        else if ( action == GLFW_RELEASE )
            w_key_flag = false;
        else if ( action == GLFW_REPEAT)
            w_key_flag = true;
	}

	if ( key == GLFW_KEY_S )
    {
        if ( action == GLFW_PRESS )
            s_key_flag = true;
        else if ( action == GLFW_RELEASE )
            s_key_flag = false;
        else if ( action == GLFW_REPEAT)
            s_key_flag = true;
	}

	if ( key == GLFW_KEY_A )
    {
        if ( action == GLFW_PRESS )
            a_key_flag = true;
        else if ( action == GLFW_RELEASE )
            a_key_flag = false;
        else if ( action == GLFW_REPEAT)
            a_key_flag = true;
	}

	if ( key == GLFW_KEY_D )
    {
        if ( action == GLFW_PRESS )
            d_key_flag = true;
        else if ( action == GLFW_RELEASE )
            d_key_flag = false;
        else if ( action == GLFW_REPEAT)
            d_key_flag = true;
	}

	if ( key == GLFW_KEY_UP )
    {
        if ( action == GLFW_PRESS )
            UP_key_flag = true;
        else if ( action == GLFW_RELEASE )
            UP_key_flag = false;
        else if ( action == GLFW_REPEAT)
            UP_key_flag = true;
	}

	if ( key == GLFW_KEY_DOWN )
    {
        if ( action == GLFW_PRESS )
            DOWN_key_flag = true;
        else if ( action == GLFW_RELEASE )
            DOWN_key_flag = false;
        else if ( action == GLFW_REPEAT)
            DOWN_key_flag = true;
	}

	if ( key == GLFW_KEY_RIGHT )
    {
        if ( action == GLFW_PRESS )
            RIGHT_key_flag = true;
        else if ( action == GLFW_RELEASE )
            RIGHT_key_flag = false;
        else if ( action == GLFW_REPEAT)
            RIGHT_key_flag = true;
	}

	if ( key == GLFW_KEY_LEFT )
    {
        if ( action == GLFW_PRESS )
            LEFT_key_flag = true;
        else if ( action == GLFW_RELEASE )
            LEFT_key_flag = false;
        else if ( action == GLFW_REPEAT)
            LEFT_key_flag = true;
	}
}

void key_handler()
{
	if(n_key_flag)
		z_camera++;

	if(m_key_flag)
		z_camera--;


	if(w_key_flag)
		y_camera++;

	if(s_key_flag)
		y_camera--;

	if(a_key_flag)
		x_camera--;

	if(d_key_flag)
		x_camera++;

	if(UP_key_flag)
		x_rotate_angle++;

	if(DOWN_key_flag)
		x_rotate_angle--;

	if(LEFT_key_flag)
		y_rotate_angle++;

	if(RIGHT_key_flag)
		y_rotate_angle--;
}

double vec_length(const vec& a)
{
	double l = pow(a[0],2) + pow(a[1],2) + pow(a[2],2);
	return sqrt(l);
}

void random_force(Body& b)
{
	random_device rd;   // Obtain a random seed from the hardware
	mt19937 gen(rd());  // Seed the random number generator
	uniform_int_distribution<int> dis(0, b.points.size());  // Define the range
	int random_poins = dis(gen);  // Generate a random number between 0 and 255
	arma::arma_rng::set_seed_random();  // Set the random seed

  vec randomVec = arma::randu<arma::vec>(3);
	vec random_force = 10000 * randomVec;
	b.points[random_poins]->F_list[2] = random_force;
}
