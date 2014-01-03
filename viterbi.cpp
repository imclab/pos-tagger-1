#include "viterbi.h"
#include <math.h>
#include <iostream>
#include <vector>
#include <string>
#include "log_functions.h"

/*
 * OPERATING IN LOG SPACE
 */

typedef std::map<std::string, double> string_double_map;

std::vector<std::string > viterbi(std::vector<std::string>& obs, 
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
	
	std::map<std::string, std::vector<std::string> > paths;

	for (unsigned int t = 1; t < obs.size(); ++t) {
	// time step
	
		//again only want to keep the states that actually appear
		std::vector<std::string> current_states = word2pos[obs[t]];
		if (1 == t) {
			//we need to intialize the paths
			for (std::vector<std::string>::const_iterator y = current_states.begin(); y != current_states.end(); ++y) {
				std::vector<std::string> temppath;
				paths[*y] = temppath;
				paths[*y].push_back(*y);
			}
		}
		
		std::map<std::string, std::vector<std::string> > newpaths;

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
			/*
			double max_prob = log(0); //everything should be bigger than -inf
			std::string max_tag = "?"; //just some meaningless value
			*/
			double max_prob = log(0);
			std::string max_tag = prev_states[0];
			
			for (std::vector<std::string>::const_iterator state = prev_states.begin(); state != prev_states.end(); ++state) {
				double log_p_transition = transition[*state][current_tag]; 
				double log_forward = prev_trellis_col[*state];
				double log_mult = log_p_transition + log_forward;
				//if (log_mult > max_prob) {
				//std::cout<<"mult = " << log_mult << " max prob " << max_prob<<" comp "<< logGT(log_mult, max_prob)<<std::endl;
				if (logGT(log_mult, max_prob)) {
					//this is bigger, so update
					max_prob = log_mult;
					max_tag = *state;
				}
			}
			//now we have the logs of the inside of the summation in inside_summation
			log_prob += max_prob;
			//summed them and added it the the log prob	
			/*************************************************************************
			 * calculating forward probability
			 *************************************************************************/

			trellis[t][*y] = log_prob;
			newpaths[*y] = paths[max_tag];
			//std::cout<<"pushing back " << *y<<std::endl;
			newpaths[*y].push_back(*y);
		}
		paths = newpaths;
	}
	//now we have the paths and the trellis
	//we need to find the max probability in the last col and then get its sequence
	std::vector<std::string> last_states = word2pos[obs[obs.size()-1]];
	//std::string max_state = "?"; //some meanignless tag
	std::string max_state = last_states[0];
	double max_prob = log(0); //-inf
	for (std::vector<std::string>::iterator tag = last_states.begin(); tag != last_states.end(); ++tag) {
		double current_prob = trellis[obs.size()-1][*tag];	
		//if (current_prob > max_prob) {
		if (logGT(current_prob, max_prob)) {
			max_prob = current_prob;
			max_state = *tag;
		}
	}
	return paths[max_state];
	//to make it compile, if it hits the other return statement, then the code is broken
	std::cout<<"something went wrong in this code"<<std::endl;
	std::vector<std::string> temp;
	return temp;
}
