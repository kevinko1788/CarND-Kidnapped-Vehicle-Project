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

using namespace std;

// declare a random engine to be used across multiple and various method calls
static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  // TODO: Set the number of particles
  num_particles = 100; 
  
  // normal distributions for sensor noise
  normal_distribution<double> dist_x (x, std[0]);
  normal_distribution<double> dist_y (y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);



  // Generate particles
  for (int i = 0; i < num_particles; i++){
    Particle particle;
    particle.id    = i;
    particle.x     = dist_x(gen);
    particle.y     = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight= 1.0;

    particles.push_back(particle);

  }
  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  //define normal distributions for sensor noise
  normal_distribution<float> dist_x (0, std_pos[0]);
  normal_distribution<float> dist_y (0, std_pos[1]);
  normal_distribution<float> dist_theta (0, std_pos[2]);

  for (int i = 0; i < num_particles; i++){
    // when yaw rate is very small (closed to 0)
    if(fabs(yaw_rate) < 0.00001){
      particles[i].x += velocity * delta_t * cos(particles[i].theta);
      particles[i].y += velocity * delta_t * sin(particles[i].theta);
    }
    // when yaw rate is != 0
    else{
      particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + (yaw_rate *  delta_t)) - sin(particles[i].theta));
      particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
      particles[i].theta += yaw_rate * delta_t;
    }
    // add noise
    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);
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
  unsigned int nObservations = observations.size();
  unsigned int nPredicted = predicted.size();

  for (unsigned int i = 0; i < nObservations; i++){

    // initializ big number for minimum distance
    double minDistance = numeric_limits<double>::max();

    // initialize map ID
    int mapID = -1;

    for (unsigned int j = 0; j < nPredicted; j++){
      // calculation or 
      double current_dist = dist(observations[i].x ,observations[i].y, predicted[i].x, predicted[j].y);
      // double xDistance    = observations[i].x - predicted[j].x;
      // double yDistance    = observations[i].y - predicted[j].y;
      // double current_dist = pow(xDistance,2) + pow(yDistance,2);


      if (current_dist < minDistance){
        minDistance = current_dist;
        mapID = predicted[j].id;
      }
    }
    observations[i].id = mapID;
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
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

  
  for (int i = 0; i < num_particles; i++){
    double particle_x     = particles[i].x;
    double particle_y     = particles[i].y;
    double particle_theta = particles[i].theta;
 
    vector<LandmarkObs> predictions;

    for (unsigned int j = 0; j<map_landmarks.landmark_list.size(); j++){
      //get id and x,y coordinates
      float lm_x = map_landmarks.landmark_list[j].x_f;
      float lm_y = map_landmarks.landmark_list[j].y_f;
      int   lm_id= map_landmarks.landmark_list[j].id_i;
    
      if(fabs(lm_x - particle_x) <= sensor_range && fabs(lm_y - particle_y) <= sensor_range){
        predictions.push_back(LandmarkObs{lm_id,lm_x,lm_y});
      }
    }
    // transform car observation coordinates to map coordinates
    vector<LandmarkObs> transformed_observations;

    for (unsigned int j = 0; j < observations.size(); j++){
      LandmarkObs transformed;
      transformed.id = j;
      transformed.x  = particle_x + (cos(particle_theta) * observations[j].x)-(sin(particle_theta) * observations[j].y);
      transformed.y  = particle_y + (sin(particle_theta) * observations[j].x)+(cos(particle_theta) * observations[j].y);
      transformed_observations.push_back(transformed);
    }
  
    //Data association for the predictions and transformed observations on current particle
    dataAssociation(predictions, transformed_observations);
    particles[i].weight = 1.0;

    // Multivariate_Gaussian
    double sigma_x = std_landmark[0];
    double sigma_y = std_landmark[1];
    double sigma_x_2 = pow(sigma_x,2);
    double sigma_y_2 = pow(sigma_y,2);
    double normalizer = (1.0/(2.0 * M_PI*sigma_x*sigma_y));


    for (unsigned int j = 0; j < transformed_observations.size(); j++){
      double trans_obs_x = transformed_observations[j].x;
      double trans_obs_y = transformed_observations[j].y;
      int trans_obs_id = transformed_observations[j].id;

      for (unsigned int k = 0; k < predictions.size(); k++){
        double pred_x = predictions[k].x;
        double pred_y = predictions[k].y;
        int pred_id = predictions[k].id;

        if (pred_id ==trans_obs_id){
          double observations_weight = normalizer * exp(-1.0 * ((pow((trans_obs_x - pred_x), 2)/(2.0 * sigma_x_2)) + (pow((trans_obs_y - pred_y), 2)/(2.0 * sigma_y_2))));
          particles[i].weight += observations_weight;
        }
      }
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
  // get weights
  vector<double> weights;
  for (int i =0; i < num_particles; i++){
    weights.push_back(particles[i].weight);
  }

  // get max value
  double maxWeight = *max_element(weights.begin(), weights.end());

  // initiate beta
  double beta = 0.0;

  // uniform distribution
  uniform_real_distribution<double> distDouble(0.0, maxWeight);
  uniform_int_distribution<int> distInt(0, num_particles - 1);

  // Generate index
  int index = distInt(gen);

  // Resample wheel
  vector<Particle> resampled_particles;
  for(int i=0; i < num_particles; i++){
    beta += distDouble(gen)*2.0;
    while(beta > weights[index]){
      beta -= weights[index];
      index = (index+1)%num_particles;
    }
    resampled_particles.push_back(particles[index]);
  }
  particles = resampled_particles;


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