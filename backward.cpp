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

std::vector<string_double_map > backward(std::vector<std::string>& obs, 
std::map<std::string, string_double_map >& transition, 
std::map<std::string, string_double_map >& emission, 
std::map<std::string, std::vector<std::string> >& word2pos) {

	//std::vector<std::map<std::string, double> > trellis (obs.size());
	std::vector<string_double_map> trellis (obs.size());
	//the trellis is a vector of dictionaries trellis[time_step][tag] = probability
	
	//we only want to have the tag/states that appear	
	unsigned int lastPos = obs.size()-1;
	std::vector<std::string> initial_states = word2pos[obs[lastPos]];
	
	//std::cout<<obs.size()<<std::endl;

	//initializing the trellis, t=end
	for (std::vector<std::string>::const_iterator state = initial_states.begin(); state != initial_states.end(); ++state) {
		trellis[lastPos][*state] = 0;
		//log(1) is 0
		//updated the probability for the trellis
	} 

	//running the viterbi algo on the trellis now
	for (int t = lastPos-1; t >= 0; --t) {
	// time step
		//again only want to keep the states that actually appear
		std::vector<std::string> current_states = word2pos[obs[t]];
		
		//for each state, were going to calculate the backward probability
		for (std::vector<std::string>::const_iterator y = current_states.begin(); y != current_states.end(); ++y) {
			
			/*************************************************************************
			 * calculating backward probability
			 *************************************************************************/
			int time_step = t;
			std::string current_tag = *y;
			std::string current_word = obs[time_step];
			std::string next_word= obs[time_step+1];
			//the current tag

			//log_prob is initiaised to the log(P(w_i, t))
			double log_prob = emission[current_tag][current_word];
			
			string_double_map prev_trellis_col = trellis[time_step+1];
			//previous time step's trellis values
			
			std::vector<std::string> prev_states = word2pos[next_word];
			//prev_states has the tags of time_step-1
			
			std::vector<double> inside_summation;
			for (std::vector<std::string>::const_iterator bw_prev_state = prev_states.begin(); bw_prev_state != prev_states.end(); ++bw_prev_state) {
				double log_p_transition = transition[*bw_prev_state][current_tag]; 
				double log_back = trellis[time_step+1][*bw_prev_state];
				double log_p_emission = emission[*bw_prev_state][next_word];
				double log_mult = log_p_transition + log_back + log_p_emission;
				inside_summation.push_back(log_mult);
			}
			//std::cout<<"pos="<<time_step<<" tag="<<*y<<"log prob="<<log_prob<<std::endl;
			//now we have the logs of the inside of the summation in inside_summation
			log_prob += logSigma(inside_summation);
			//summed them and added it the the log prob	
			/*************************************************************************
			 * calculating backward probability
			 *************************************************************************/

			trellis[t][*y] = log_prob;
		}
	}

	return trellis;
}
