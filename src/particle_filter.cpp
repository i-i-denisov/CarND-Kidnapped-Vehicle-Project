/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

std::vector<LandmarkObs> Transform_Obs (const std::vector<LandmarkObs> &observations,const Particle &particle)
{
  std::vector<LandmarkObs> transformed_obs=observations;
  for (int i=0;i<observations.size();i++)
    {
      transformed_obs[i].x=particle.x+cos(particle.theta)*observations[i].x-sin(particle.theta)*observations[i].y;
      transformed_obs[i].y=particle.y+sin(particle.theta)*observations[i].x+cos(particle.theta)*observations[i].y;
    }
return transformed_obs;
}

double Calculate_weights (const std::vector<LandmarkObs> &transformed_obs,const Map &map_landmarks,double std_landmark[])
{
  double weight=1.;
  //coords of associated map landmark
  double mu_x,mu_y;
  //loop through measurements
  //std::cout <<"calculating weights"<<std::endl;
  for (int i=0;i<transformed_obs.size();i++)
  {
    //getting coords of associated map landmark
    mu_x=map_landmarks.landmark_list[transformed_obs[i].id-1].x_f;
    mu_y=map_landmarks.landmark_list[transformed_obs[i].id-1].y_f;
    //calculating cumulative weight of given particle using transformed measuremtns
    //std::cout <<"observation "<<i<<std::endl;
    weight*=multiv_prob(std_landmark[0],std_landmark[1],transformed_obs[i].x,transformed_obs[i].y,mu_x,mu_y);
  }
  //std::cout <<"weight "<<weight<<std::endl;
  return weight;
}

double multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs,
                   double mu_x, double mu_y) {
  // calculate normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

  // calculate exponent
  double exponent;
  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
    
  // calculate weight using normalization terms and exponent
  double weight;
  weight = gauss_norm * exp(-exponent);
  //std::cout<<weight<<std::endl;  
  return weight;
  
}

void Associate_Obs(std::vector<LandmarkObs> &transformed_obs,const Map &map_landmarks)
{
  //id of a closest landmark on the map
  int id;
  //closest distance between currnet observation and any landmark
  double closest;
  //distance between current observation and current landmark
  double distance;
  for (int i=0;i<transformed_obs.size();i++)
  {
    closest=dist(transformed_obs[i].x,transformed_obs[i].y,map_landmarks.landmark_list[0].x_f,map_landmarks.landmark_list[0].y_f);
    id=0;
    for (int j=1;j<map_landmarks.landmark_list.size();j++)
    {
      distance=dist(transformed_obs[i].x,transformed_obs[i].y,map_landmarks.landmark_list[j].x_f,map_landmarks.landmark_list[j].y_f);
      if (distance<closest) 
      {
        id=map_landmarks.landmark_list[j].id_i;
        closest=distance;
      }
    }
  transformed_obs[i].id=id;
  //std::cout << "obs_x "<< transformed_obs[i].x<< " obs_y"<<transformed_obs[i].y<<std::endl;
  //std::cout <<"observation "<<i<< " associated landmark "<<id<<std::endl;
  }
}


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   *  Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   *  Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 10;  // TODO: Set the number of particles
  is_initialized=true;
  //creating randomizer device and creating noise distributions for x,y,theta
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<double> noise_x(0.0, std[0]);
  std::normal_distribution<double> noise_y(0.0, std[1]);
  std::normal_distribution<double> noise_theta(0.0, std[2]);

  //lopping through particles set and initializing each with gaussian noise
  for (int i=0;i<num_particles;i++)
  {
    Particle particle;
    particle.id=i;
    particle.x=x+noise_x(gen);
    particle.y=y+noise_y(gen);
    particle.theta=theta+noise_y(gen);
    particle.weight=1.;
    this->particles.push_back(particle);
  }
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
int number_of_particles=this->num_particles;
for(int i=0; i<number_of_particles;i++)
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<double> noise_x(0.0, std_pos[0]);
  std::normal_distribution<double> noise_y(0.0, std_pos[1]);
  std::normal_distribution<double> noise_theta(0.0, std_pos[2]);
  this->particles[i].x+=(delta_t*velocity*cos(yaw_rate)+noise_x(gen));
  this->particles[i].y+=(delta_t*velocity*sin(yaw_rate)+noise_y(gen));
  this->particles[i].theta=yaw_rate+noise_theta(gen);

}
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) 
{
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  int number_of_particles=this->num_particles;
  std::vector<LandmarkObs> transformed_obs=observations;
  double normalizer=0.;
  //loop through particles
  for (int i=0;i<number_of_particles;i++)
  {
      
      
      //Particle particle=this->particles[i];
      //debug
      std::cout <<"paricle "<<i<<std::endl;
      //std::cout <<"x "<<particle.x<<"y "<<particle.y<<std::endl;
      // for each particle transform observation from local to map coordinates
      transformed_obs=Transform_Obs(observations,this->particles[i]);
      //associate each observation with map entry
      Associate_Obs(transformed_obs,map_landmarks);
      //calculate weights
      this->particles[i].weight=Calculate_weights(transformed_obs,map_landmarks,std_landmark);
      std::cout<<"weigths updated "<< this->particles[i].weight<<std::endl;
      //TODO: particles outside the sensor range
      normalizer+=this->particles[i].weight;
  }
//TODO: normalize weights
if (normalizer==0){
  normalizer=1.;
  std::cout<<"failed to normalize weights"<<std::endl;
}
else {
  for (int i=0;i<number_of_particles;i++)
    {
    this->particles[i].weight=(this->particles[i].weight/normalizer);
    std::cout <<"Normalized weight "<<this->particles[i].weight<<std::endl;
    }
  }
}


void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}