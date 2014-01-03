/*
 * pos tagger
 */
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <set>
#include <stdlib.h>

#include <cmath>

//my functions
#include "log_functions.h"
#include "viterbi.h"
#include "forward.h"
#include "backward.h"

//define negative infinity so we don't need to compute it later
double neg_inf = -1.0 * std::numeric_limits<double>::infinity(); 

using namespace std;
typedef std::vector<std::string> string_vec;
typedef std::map<string, string_vec> string_to_sentence;
typedef std::vector<std::map<std::string, double> > trellis;
typedef std::map<std::string, std::map<std::string, double> > movement_matrix;
typedef std::map<std::string, std::map<std::string, double> > count_matrix;
//movement matrix can be the emission or transition matrix

string_vec all_tags;
//all the possible tags
string_vec all_words;
//all the words that we can have
string_to_sentence word2pos;
//words corresponding to their pos tags
vector<string_vec> sentences;
//all the sentences
movement_matrix transition;
movement_matrix emission;
//the transition and emission matrices
map<std::string, double> initial_p;
//initial p is the initial probabilities

bool map_contains(map<string, double>& target_map, string key) {
	return target_map.find(key) != target_map.end();
}

void import_lexicon(string& filename) {
	std::ifstream in(filename.c_str());
	string line;
	set<string> tags;

	while(getline(in, line)) {
		stringstream iss(line);

		string word;
		bool read_word = false;
		string_vec lex;

		string temp;
		while(iss >> temp) {
			if (false == read_word) {
				word = temp;
				all_words.push_back(word);
				read_word = true;
			}
			else {
				lex.push_back(temp);
				tags.insert(temp);
			}
		}
		word2pos[word] = lex;
	}
	//move everything from sets into all_tags
	for (set<string>::iterator itr = tags.begin(); itr != tags.end(); ++itr) {
		all_tags.push_back(*itr);
	}
	return;
}

void import_sentences(string& filename) {
	ifstream in(filename.c_str());
	//open the file
	string line;
	//where line is stored

	while(getline(in, line)) {
		//read every line

		stringstream iss(line);

		string_vec sentence;
		//sentence is the representation of it as a vector
		string word;
		//the word were seeing
		while(iss >> word) {
			sentence.push_back(word);	
			//push it into the sentence
		}
		sentences.push_back(sentence);
		//store it in the big sentences vector
	}	
}

void initialize_random_transition() {
	//generate random numbers, then normalize
	for(string_vec::iterator tag_i = all_tags.begin(); tag_i != all_tags.end(); ++tag_i) {
		double sum = 0.0;
		for(string_vec::iterator tag_j = all_tags.begin(); tag_j != all_tags.end(); ++tag_j) {
			transition[*tag_i][*tag_j] = rand();
			sum += transition[*tag_i][*tag_j];
		}	
		//normalize the row
		for(string_vec::iterator tag_j = all_tags.begin(); tag_j != all_tags.end(); ++tag_j) {
			transition[*tag_i][*tag_j] /= sum;
			//convert it into log
			transition[*tag_i][*tag_j] = log(transition[*tag_i][*tag_j]);
		}
	}	
}

count_matrix initialize_transition_count() {
	//generate random numbers, then normalize
	count_matrix transition_count;
	for(string_vec::iterator tag_i = all_tags.begin(); tag_i != all_tags.end(); ++tag_i) {
		for(string_vec::iterator tag_j = all_tags.begin(); tag_j != all_tags.end(); ++tag_j) {
			transition_count[*tag_i][*tag_j] = neg_inf; 
			//neg inf because log(0)
		}	
	}	
	return transition_count;
}

void initialize_random_emission() {
	//generate random numbers then normalize
	for(string_vec::iterator tag_i = all_tags.begin(); tag_i != all_tags.end(); ++tag_i) {
		double sum = 0;
		for(string_vec::iterator word = all_words.begin(); word != all_words.end(); ++word) {
			emission[*tag_i][*word] = rand();
			sum += emission[*tag_i][*word];
		}	
		//normalize the row
		for(string_vec::iterator word = all_words.begin(); word != all_words.end(); ++word) {
			emission[*tag_i][*word] /= sum;
			emission[*tag_i][*word] = log(emission[*tag_i][*word]);
		}
	}	
	//normalize somehow
}

count_matrix initialize_emission_count() {
	count_matrix emission_count;
	//generate random numbers then normalize
	for(string_vec::iterator tag_i = all_tags.begin(); tag_i != all_tags.end(); ++tag_i) {
		for(string_vec::iterator word = all_words.begin(); word != all_words.end(); ++word) {
			emission_count[*tag_i][*word] = neg_inf; 
			//neg inf because log (0)
		}	
	}	
	//normalize somehow
	return emission_count;
}

void initialize_random_initial() {
	//generate random numbers and normalize it in this
	double sum = 0.0;
	for(string_vec::iterator tag = all_tags.begin(); tag != all_tags.end(); ++tag) {
		initial_p[*tag] = rand();	
		sum += initial_p[*tag];
	}
	//normalizing
	for(string_vec::iterator tag = all_tags.begin(); tag != all_tags.end(); ++tag) {
		initial_p[*tag] /= sum; 
		initial_p[*tag] = log(initial_p[*tag]);
	}
}

map<std::string, double> initialize_initial_count() {
	map<std::string, double> initial_count;
	for(string_vec::iterator tag = all_tags.begin(); tag != all_tags.end(); ++tag) {
		initial_count[*tag] = neg_inf; 
		//log (0)
	}
	return initial_count;
}

string_vec sentence_tags(string_vec sentence) {
	set<string> local_tags;
	for (string_vec::iterator it = sentence.begin(); it != sentence.end(); ++it) {
		string_vec word_tags = word2pos[*it];

		for (string_vec::iterator tag = word_tags.begin(); tag != word_tags.end(); ++tag) {
			local_tags.insert(*tag);
			//put the tags for each word into a set
		}
	}
	string_vec sentence_tags; 

	for (set<string>::iterator tag = local_tags.begin(); tag != local_tags.end(); ++tag) {
		//move then out of the set into the string vector
		sentence_tags.push_back(*tag);
	}
	return sentence_tags;
}

bool string_vec_contains(string_vec& vec, string search){
	return find(vec.begin(), vec.end(), search) != vec.end();
}

//P(sentence | model)
double log_p_of_sentence(trellis& log_fw_trellis, string_vec sentence) {
    //sigma t forward[t][n]
    vector<double> sum_elems;
    string_vec possible_tags = word2pos[sentence[sentence.size()-1]];
    for (string_vec::iterator tag = possible_tags.begin(); tag != possible_tags.end(); ++tag) {
        sum_elems.push_back(log_fw_trellis[sentence.size()-1][*tag]);
    }   
    return logSigma(sum_elems);
}

double log_e_count_of_tag_and_word_in_sentence(trellis& fw_trellis, trellis& bw_trellis, int time, string tag, string_vec sentence) {
	string_vec possible_tags = word2pos[sentence[time]];
        //p(t_i = t, w)
	double log_numerator = fw_trellis[time][tag] + bw_trellis[time][tag];
        //p(w)
	double prob_sentence = log_p_of_sentence(fw_trellis, sentence);
	if (isinf(log_numerator) && isinf(prob_sentence)) {
		return log_numerator;
		//-inf - -inf will cause a -NaN
	}
        //count is p(t_i = t, w) / p(w)
	return log_numerator - prob_sentence;
}

double log_e_count_of_tag_and_other_tag(trellis& fw, trellis& bw, movement_matrix& t, movement_matrix& e, int time, string tag, string other_tag, string_vec sentence) {
	double log_forward = fw[time][tag];
	double log_p_transition = t[tag][other_tag];
	double log_p_emission = e[other_tag][sentence[time+1]];
	double log_backward = bw[time+1][other_tag];
        //p(t_i = t, w_1 .. w_i) * p(t | t') * p( w_i+1 | t') * p(w_i+2 .. w_n | t_i+1 = t')
	//cout<<"fw = "<<log_forward<<" p="<<log_p_transition<<" emission="<<log_p_emission<<" bw="<<log_backward<<endl;

	double log_numerator = log_forward + log_p_transition + log_p_emission + log_backward;
	double log_denominator = log_p_of_sentence(fw, sentence); 
	if (isinf(log_numerator) && isinf(log_denominator)) {
		return log_numerator;
		//-inf - -inf will cause a -NaN
	}
	return log_numerator - log_denominator;
}

double log_e_count_of_tag_as_first(trellis& fw, trellis& bw, movement_matrix& emit, string tag, string_vec sentence) {
    //p(t_0 = t | w) = p(t_0 = t, w) / p(w)
    //p(t_0 =t, w) = pi(t) * p(w_0, t) * backward[tag][0]
    double log_pi = initial_p[tag];  
    double log_emit = emit[tag][sentence[0]]; 
    double log_back = bw[0][tag];
    
    double log_numerator = log_pi + log_emit + log_back;
    //cout<<"pi = "<<log_pi<<" log emit = "<<log_emit<<" log back = "<<log_back<<endl;
    double log_denominator = log_p_of_sentence(fw, sentence);
	if (isinf(log_numerator) && isinf(log_denominator)) {
		return log_numerator;
		//-inf - -inf will cause a -NaN
	}
    return log_numerator - log_denominator;
     
}

vector<vector<trellis> > training_iteration() {
	//train the hmm
	vector<vector<trellis> > all_trellises;
	for	(vector<string_vec>::iterator sentence = sentences.begin(); sentence != sentences.end(); ++sentence) {
		string_vec s_tags = sentence_tags(*sentence);
		//tags that appear in the sentence

		trellis fw_trellis = forward(*sentence, initial_p, transition, emission, word2pos);			
		//cout<<"calling backwards"<<endl;
		trellis bw_trellis = backward(*sentence, transition, emission, word2pos);			
		vector<trellis> current_trellis;
		//put fw and bw in then put it to the back
		current_trellis.push_back(fw_trellis);
		current_trellis.push_back(bw_trellis);
		//push it back now
		all_trellises.push_back(current_trellis);
		//at this point, we have calculated the values of the fw and bw trellises
	}

	cout<<"finished calculating trellises"<<endl;
	//setup counts now
	count_matrix log_t_count = initialize_transition_count();
	count_matrix log_e_count = initialize_emission_count();
        //init count only counts if a tag is in the first position
	map<string, double> log_init_count = initialize_initial_count();
        //log_tag_count has the expected count of each tag
	map<string, double> log_tag_count = initialize_initial_count();
	
	//count emissions	
	for (unsigned int sentence_id = 0; sentence_id < sentences.size(); ++sentence_id) {
		string_vec sentence = sentences[sentence_id]; 
		//current sentence

		//get the trellises
		vector<trellis> all_t = all_trellises[sentence_id];
		trellis fw = all_t[0];
		trellis bw = all_t[1];

		for (unsigned int time = 0; time < sentence.size(); ++time) {
			string word = sentence[time];
			string_vec tags = word2pos[word];

			for (string_vec::iterator tag_i_ptr = tags.begin(); tag_i_ptr != tags.end(); ++tag_i_ptr) {
				string tag = *tag_i_ptr;
				//e_count[tag][word] += log_e_count_of_tag_and_word_in_sentence(fw, bw, time, tag, sentence);
				double log_expected_count = log_e_count_of_tag_and_word_in_sentence(fw, bw, time, tag, sentence);
				log_e_count[tag][word] = logAdd(log_e_count[tag][word], log_expected_count);
			}
		}
	}


	//also count the transitions
	for (unsigned int sentence_id = 0; sentence_id < sentences.size(); ++sentence_id) {
		string_vec sentence = sentences[sentence_id];

		//get the trellises
		vector<trellis> all_t = all_trellises[sentence_id];
		trellis fw = all_t[0];
		trellis bw = all_t[1];
		
		for (unsigned int time = 0; time < sentence.size()-1; ++time) {
			string_vec curr_tags = word2pos[sentence[time]];
			string_vec next_tags = word2pos[sentence[time+1]];

			for (string_vec::iterator c_tag_itr = curr_tags.begin(); c_tag_itr != curr_tags.end(); ++c_tag_itr) {
				string c_tag = *c_tag_itr;
				for (string_vec::iterator n_tag_itr = next_tags.begin(); n_tag_itr != next_tags.end(); ++n_tag_itr) {
					string n_tag = *n_tag_itr;
					// += log_e_count_of_tag_and_other_tag(fw, bw, transition, emission, time, c_tag, n_tag, sentence);
					double log_expected_count = log_e_count_of_tag_and_other_tag(fw, bw, transition, emission, time, c_tag, n_tag, sentence);
					//cout<<"tag1="<<c_tag<<" tag2="<<n_tag<<endl;
					double result = logAdd(log_t_count[c_tag][n_tag], log_expected_count);
					//cout<<"t count="<<log_t_count[c_tag][n_tag]<<" ec="<<log_expected_count<<" result="<<result<<endl;
					log_t_count[c_tag][n_tag] = result;
				}
			}
		}
	}

	//count the initial values
	for (unsigned int sentence_id = 0; sentence_id < sentences.size(); ++sentence_id) {
		string_vec sentence = sentences[sentence_id];

		trellis fw = all_trellises[sentence_id][0];
                trellis bw = all_trellises[sentence_id][1];

		string_vec tags = word2pos[sentence[0]];
		for (string_vec::iterator tag_itr = tags.begin(); tag_itr != tags.end(); ++tag_itr) {
                        double val = log_e_count_of_tag_as_first(fw, bw, emission, *tag_itr, sentence);
                        //cout<<"tag "<<*tag_itr<<" val "<<val<<endl;
			//cout<<"Tag="<<*tag_itr<<endl;
			double ir = logAdd(log_init_count[*tag_itr], val);
			//cout<<"prev val="<<log_init_count[*tag_itr]<<" val = "<<val<<" ir = "<<ir<<endl;
			log_init_count[*tag_itr] = ir;
			//cout<<"new val="<<log_init_count[*tag_itr]<<endl;
		}
	}

        //count tag
        for (string_vec::iterator tag_itr = all_tags.begin(); tag_itr != all_tags.end(); ++tag_itr) {
            string tag = *tag_itr;
            for (string_vec::iterator word_itr = all_words.begin(); word_itr != all_words.end(); ++word_itr) {
                string word = *word_itr;
                //log_tag_count[tag] += e_count[tag][word]; 
		//replace with logadd or logsigma
		log_tag_count[tag] = logAdd(log_tag_count[tag], log_e_count[tag][word]);
            }
        }

	cout<<"calculated expected counts"<<endl;

        //calculate the tag transition probabilities
        for (string_vec::iterator tag_i_itr = all_tags.begin(); tag_i_itr != all_tags.end(); ++tag_i_itr) {
            string tag_i = *tag_i_itr;
            for (string_vec::iterator tag_j_itr = all_tags.begin(); tag_j_itr != all_tags.end(); ++tag_j_itr) {
                string tag_j = *tag_j_itr;
                //p(t_j | t_i) = c(t_i t_j) / c(t_i) 
		//possible point of underflow
                //transition[tag_i][tag_j] = t_count[tag_i][tag_j] / log_tag_count[tag_i];
		//cout<<"t count = "<< log_t_count[tag_i][tag_j] << " tag count = " << log_tag_count[tag_i]<<endl;
                //transition[tag_i][tag_j] = log_t_count[tag_i][tag_j] - log_tag_count[tag_i];
                transition[tag_i][tag_j] = log_t_count[tag_i][tag_j] - log_tag_count[tag_j];
            }
        }

        //calculate emission probabilities
        for (string_vec::iterator tag_itr = all_tags.begin(); tag_itr != all_tags.end(); ++tag_itr) {
            string tag = *tag_itr;
            for (string_vec::iterator word_itr = all_words.begin(); word_itr != all_words.end(); ++word_itr) {
                string word = *word_itr;
                //p(t | w) = c(t, w) / c(t) 
                //emission[tag][word] = e_count[tag][word] / log_tag_count[tag];
                emission[tag][word] = log_e_count[tag][word] - log_tag_count[tag];
		
		//log of emission can't be a positive number, that would mean the probability is > 1
            }
        }

	//divide now
        for (string_vec::iterator it = all_tags.begin(); it != all_tags.end(); ++it) {
            double numerator = log_init_count[*it];
            double denominator = log(sentences.size());
	    if (0 == log_init_count[*it]) {
		//its probability is 0
		initial_p[*it] = 0;
		}
		else {
            initial_p[*it] = numerator - denominator;
		}

        }
	
	//at this point, the new probabilities have been calculated, but they have not been normalized
	
	//normalize the initial probability
	double initial_total = log(0);
	//-inf
	for (string_vec::iterator it = all_tags.begin(); it != all_tags.end(); ++it) {
		initial_total = logAdd(initial_total, initial_p[*it]);
		//sum everything
	}
	for (string_vec::iterator it = all_tags.begin(); it != all_tags.end(); ++it) {
		initial_p[*it] -= initial_total;	
		//divide it by the total
	}

	//normalize transitions
	double transition_tag_total = log(0); //-inf 
        for (string_vec::iterator tag_i_itr = all_tags.begin(); tag_i_itr != all_tags.end(); ++tag_i_itr) {
            string tag_i = *tag_i_itr;
            for (string_vec::iterator tag_j_itr = all_tags.begin(); tag_j_itr != all_tags.end(); ++tag_j_itr) {
		string tag_j = *tag_j_itr;
		//sum everything
		transition_tag_total = logAdd(transition_tag_total, transition[tag_i][tag_j]);
	    }
            for (string_vec::iterator tag_j_itr = all_tags.begin(); tag_j_itr != all_tags.end(); ++tag_j_itr) {
		string tag_j = *tag_j_itr;
		transition[tag_i][tag_j] -= transition_tag_total;
		//divide by total
	    }
	}

	//normalize emissions
	double emission_tag_total = log(0); //-inf
        for (string_vec::iterator tag_itr = all_tags.begin(); tag_itr != all_tags.end(); ++tag_itr) {
            string tag = *tag_itr;
            for (string_vec::iterator word_itr = all_words.begin(); word_itr != all_words.end(); ++word_itr) {
                string word = *word_itr;
                emission_tag_total = logAdd(emission_tag_total, emission[tag][word]);
            }
            for (string_vec::iterator word_itr = all_words.begin(); word_itr != all_words.end(); ++word_itr) {
                	string word = *word_itr;
		    emission[tag][word] -= emission_tag_total;
            }
        }

	cout<<"calculated new probabilities"<<endl;

	return all_trellises;
}

double log_model(vector<vector<trellis> >& t) {
	vector<double> temp;
	for (unsigned int counter = 0; counter < t.size(); ++counter) {
		trellis temp_trellis = t[counter][0];
		string_vec sentence = sentences[counter];
		double t_val = log_p_of_sentence(temp_trellis, sentence);
		temp.push_back(t_val);
	}
	return logSigma(temp);
}

bool convergence_check(vector<vector<trellis> >& old, vector<vector<trellis> >& new_trellises, double epsilon) {
	//do stuff
	double old_val = log_model(old);
	double new_val = log_model(new_trellises);
	double diff = std::abs(new_val - old_val); 
	cout<<"old = "<<old_val<<endl;
	cout<<"new = "<<new_val<<endl;
	cout<<"diff = ";
	cout<<diff<<endl;
	return diff < epsilon;
}

void run_train(double epsilon) {
	cout<<"iteration: 1"<<endl;
	//trellises = training_iteration();
	vector<vector<trellis> > trellises = training_iteration();
	//run 1
	bool has_converged = false;
	unsigned int iteration_count = 2;
	while(false == has_converged) {
		cout<<"iteration: ";
		cout<<iteration_count++<<endl;
		vector<vector<trellis> > new_trellis = training_iteration();	
		has_converged = convergence_check(trellises, new_trellis, epsilon); 
		trellises = new_trellis;
		if (iteration_count > 7) {
			has_converged = true;
			//manually end it
		}
	}		
	return;
}

string_vec tag_sentences() {
	//tags sentences and also does the confusuon matrix
	string_vec tagged_sentences;

	for (vector<string_vec>::iterator sentence = sentences.begin(); sentence != sentences.end(); ++sentence) {
		//call viterbi, have the tags now print those with the sentence
		string_vec tags = viterbi( (*sentence), initial_p, transition, emission, word2pos); 		
		//make sure we aren't missing a tag
		
		stringstream ss;
		for (unsigned int pos = 0; pos < tags.size(); ++pos) {
			ss<<(*sentence)[pos];
			ss<<"_";
			ss<<tags[pos]<<" ";	
			//print out the tagged sentence
		}
		tagged_sentences.push_back(ss.str());
	}
	return tagged_sentences;
}

void print_tagged(string_vec& tagged) {
	for (string_vec::iterator sentence = tagged.begin(); sentence != tagged.end(); ++sentence) {
		cout<<(*sentence)<<endl;
	}
}

void write_tagged(string_vec&tagged, string filename) {
	freopen(filename.c_str(), "w", stdout);
	for (string_vec::iterator sentence = tagged.begin(); sentence != tagged.end(); ++sentence) {
		cout<<(*sentence)<<endl;;
	}
}

int main(int argc, char** argv) {
	(void) argc;
	(void) argv;

	string lex = "files/lexicon.txt";
	import_lexicon(lex);
	//string train_sentences = "t10.txt";
	string train_sentences = "files/train.txt";
	import_sentences(train_sentences);
	//imported sentences
	initialize_random_transition();
	initialize_random_emission();
	initialize_random_initial();
	//initialized random values
	run_train(0.01);
	cout<<"finished training, calculating paths"<<endl;
	//its converged
	string_vec tagged = tag_sentences();
	//tagged the sentences
	cout<<"writing results"<<endl;
	write_tagged(tagged, "tagged.txt");
	return 0;
}
