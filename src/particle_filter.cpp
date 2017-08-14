/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <vector>
#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <sstream>
#include <string>
#include <iterator>
#include "map.h"
#include <map>

#include "particle_filter.h"

using namespace std;

#define NUM_PARTICLES (5)

// random number generator
static default_random_engine gen;



/**
 * init Initializes particle filter by initializing particles to Gaussian
 *   distribution around first position and all the weights to 1.
 *
 * @param x Initial x position [m] (simulated estimate from GPS)
 * @param y Initial y position [m]
 * @param theta Initial orientation [rad]
 * @param std[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
 *   standard deviation of yaw [rad]]
 */
void ParticleFilter::init(double x, double y, double theta, double std[]) {
  
  // number of particles
  num_particles = NUM_PARTICLES;

  // initial paramter noise
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  
  // initialization value for weights
  double weight = 1.0;
  
  // initialize particles
  for (int i = 0; i < num_particles; i++) {
    
    Particle p = {};
    p.id     = i;
    p.x      = dist_x(gen);
    p.y      = dist_y(gen);
    p.theta  = dist_theta(gen);
    p.weight = weight;
    particles.push_back(p);
    weights.push_back(weight);
  }
  
  // set the flag
  is_initialized = true;
  

}

/**
 * prediction Adds measurements to each particle and adds random Gaussian noise
 *
 * @param delta_t Time difference between predictions
 * @param std_pos
 * @param velocity Velocity of previous measurement
 * @param yaw_rate Yaw rate of previous measurement
 */
void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  
  double v_yaw       = velocity / yaw_rate;
  double dt_yaw      = delta_t * yaw_rate;
  double dt_velocity = delta_t * velocity;
  
  for (int i = 0; i < num_particles; ++i) {
    double new_x     = 0;
    double new_y     = 0;
    double new_theta = 0;
    
    if (yaw_rate == 0) {
      new_x     = particles[i].x + dt_velocity * cos(particles[i].theta);
      new_y     = particles[i].y + dt_velocity * sin(particles[i].theta);
      new_theta = particles[i].theta;
    }
    else {
      new_x     = particles[i].x + (v_yaw) * (sin(particles[i].theta + dt_yaw) - sin(particles[i].theta));
      new_y     = particles[i].y + (v_yaw) * (cos(particles[i].theta) - cos(particles[i].theta + dt_yaw));
      new_theta = particles[i].theta + dt_yaw;
    }
    
    // add noise
    normal_distribution<double> noise_x(new_x, std_pos[0]);
    normal_distribution<double> noise_y(new_y, std_pos[1]);
    normal_distribution<double> noise_theta(new_theta, std_pos[2]);
    
    particles[i].x     = noise_x(gen);
    particles[i].y     = noise_y(gen);
    particles[i].theta = noise_theta(gen);
    
  }
}

/**
 * updateWeights Updates the weights of each particle using a multi-variate Gaussian distribution
 *
 * @param sensor_range
 * @param std_landmark
 * @param observations
 * @param map_landmarks
 */
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
  
  //initialize
  int i;
  //num_particles = 5;
  Particle &particle = particles[i];
  
  double sumweight = 1.0;
  double prob      = 1.0;
  
  //initialize transfx, transfy
  for (int i = 0; i < observations.size(); i++) {
    transfx.push_back(0.0);
    transfy.push_back(0.0);
  }
  
  // must re-initialize weights
  for (int i = 0; i < num_particles; i++) {
    weights[i] = 1.0;
    particles[i].weight = 1.0;
  }
  
  // initialize pred_landmark
  for (int i = 0; i < 42; i++) {
    pred_landmark.push_back(70.0);
  }
  
  // initialize error
  for (int i = 0; i < num_particles; ++i) {
    double weight = 1.0;
    
    // Translating local coordinates to map coordinates
    for (int j = 0; j < observations.size(); j++) {
      
      transfx[j] = particles[i].x + ((observations[j].x) * (cos(particles[i].theta))) - ((observations[j].y) * (sin(particles[i].theta)));
      transfy[j] = particles[i].y + (observations[j].x * sin(particles[i].theta)) + (observations[j].y * cos(particles[i].theta));
      
      double min_error = 80.0;
      double error     = 90.0;
      
      for (int k = 0; k < 42; k++) {
        pred_landmark[k] = sqrt((map_landmarks.landmark_list[k].x_f - transfx[j])*(map_landmarks.landmark_list[k].x_f - transfx[j])+(map_landmarks.landmark_list[k].y_f - transfy[j])*(map_landmarks.landmark_list[k].y_f - transfy[j]));
        
        if (k > 0) {
          if (pred_landmark[k] < pred_landmark[k-1]) {
            error = pred_landmark[k];
          }
          else if (pred_landmark[k] > pred_landmark[k-1]) {
            error = pred_landmark[k-1];
          }
        }
        if (error < min_error) {
          min_error = error;
        }
      }
      
      prob = (1/(2*M_PI*std_landmark[0]*std_landmark[1]))* (exp((-min_error * min_error)/(2*std_landmark[0]*std_landmark[1])));
      
      particles[i].weight *= prob;
      weights[i] *= prob;
      
    }
      
    weight = weights[i];
    
    particle.weight = weight;
    particle.weight = particles[i].weight;
    sumweight = sumweight + weight;
    
    weights.push_back(weight);
  }
}

/**
 * resample Replaces particles with probability proportional to their weight
 */
void ParticleFilter::resample() {
  
  //implementing resampling wheel
  double sumweight = 0.0;
  
  for (int i = 0; i < num_particles; i++) {
    normprob.push_back(1);
  }
  
  for (int i = 0; i < num_particles; i++) {
    sumweight += particles[i].weight;
  }
  
  for (int i = 0; i < num_particles; i++) {
    normprob[i] = (particles[i].weight/sumweight) * num_particles + 0.5;
  }
  
  std::discrete_distribution<int>distribution{weights[0],weights[1],weights[2],weights[3],weights[4]};
  vector<Particle>resample_particles;
  
  for(int i = 0; i < num_particles; i++) {
    resample_particles.push_back(particles[distribution(gen)]);
  }
  particles = resample_particles;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
