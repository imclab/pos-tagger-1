#ifndef _FORWARD_H_
#define _FORWARD_H_

#include <vector>
#include <string>
#include <map>

/*
 *	obs is a vector of strings where where each index is part of the sentece
 *	start_p is a dictionary where start_p[tag] = pi of that tag
 *	transition is a dictionary of dictionaries transition[tag_i][tag_j] = probability of going from tag_i to tag_j
 *  emission is a dictionary of dictionaries where emission[tag][word] = probability of that tag emitting that word
 */

std::vector<std::map<std::string, double> > forward(std::vector<std::string>& obs, 
std::map<std::string, double>& start_p, 
std::map<std::string, std::map<std::string, double> >& transition,
std::map<std::string, std::map<std::string, double> >& emission,
std::map<std::string, std::vector<std::string> >& word2pos);

#endif
