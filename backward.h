#ifndef _BACKWARD_H_
#define _BACKWARD_H_

#include <vector>
#include <map>
#include <string>

/*
 *	obs is a vector of strings where where each index is part of the sentece
 *	states is a vector of strings that has every possible tag
 *	start_p is a dictionary where start_p[tag] = pi of that tag
 *	transition is a dictionary of dictionaries transition[tag_i][tag_j] = probability of going from tag_i to tag_j
 *  emission is a dictionary of dictionaries where emission[tag][word] = probability of that tag emitting that word
 */

std::vector<std::map<std::string, double> > backward(std::vector<std::string>& obs,  
std::map<std::string, std::map<std::string, double> >& transition,
std::map<std::string, std::map<std::string, double> >& emission,
std::map<std::string, std::vector<std::string> >& word2pos);

#endif
