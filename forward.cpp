#include "forward.h"
#include <math.h>
#include <iostream>
#include <vector>
#include <string>
#include "log_functions.h"

/*
 * OPERATING IN LOG SPACE
 */

typedef std::map<std::string, double> string_double_map;

std::vector<std::map<std::string, double> > forward(std::vector<std::string>& obs, 
string_double_map& start_p, 
std::map<std::string, string_double_map >& transition, 
std::map<std::string, string_double_map >& emission, 
std::map<std::string, std::vector<std::string> >& word2pos) {

	//std::vector<std::map<std::string, double> > trellis (obs.size());
	std::vector<string_double_map> trellis (obs.size());
	
	//the trellis is a vector of dictionaries trellis[time_step][tag] = probability
	
	//we only want to have the tag/states that appear	
	std::vector<std::string> initial_states = word2pos[obs[0]];
	
	//initializing the trellis, t=0 
	for (std::vector<std::string>::const_iterator state = initial_states.begin(); state != initial_states.end(); ++state) {
		double log_prob = start_p[*state] + emission[*state][obs[0]];
		trellis[0][*state] = log_prob;	
		//updated the probability for the trellis
	} 

	for (unsigned int t = 1; t < obs.size(); ++t) {
	// time step
	
		//again only want to keep the states that actually appear
		std::vector<std::string> current_states = word2pos[obs[t]];
		
		//for each state, were going to calculate the forward probability
		for (std::vector<std::string>::const_iterator y = current_states.begin(); y != current_states.end(); ++y) {
			
			/*************************************************************************
			 * calculating forward probability
			 *************************************************************************/
			//calculate the forward probability 
			int time_step = t;
			std::string current_tag = *y;
			std::string current_word = obs[time_step];
			std::string prev_word= obs[time_step-1];
			//the current tag

			//log_prob is initiaised to the log(P(w_i, t))
			double log_prob = emission[current_tag][current_word];
			
			string_double_map prev_trellis_col = trellis[time_step-1];
			//previous time step's trellis values
			
			std::vector<std::string> prev_states = word2pos[prev_word];
			//prev_states has the tags of time_step-1
			
			std::vector<double> inside_summation;
			for (std::vector<std::string>::const_iterator state = prev_states.begin(); state != prev_states.end(); ++state) {
				double log_p_transition = transition[*state][current_tag]; 
				//double log_forward = trellis[time_step-1][*state];
				double log_forward = prev_trellis_col[*state];
				double log_mult = log_p_transition + log_forward;
				inside_summation.push_back(log_mult);
			}
			//now we have the logs of the inside of the summation in inside_summation
			log_prob += logSigma(inside_summation);
			//summed them and added it the the log prob	
			/*************************************************************************
			 * calculating forward probability
			 *************************************************************************/

			trellis[t][*y] = log_prob;
		}
	}

	return trellis;
}
