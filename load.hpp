#pragma once
#include <fstream>
#include <sstream>
#include <iostream>
#include "trip_def.hpp"

bool 
load_gifts(const std::string &file, gifts_t &gifts)
{
	std::ifstream fp;
	std::string line;

	gifts.clear();
	fp.open(file.c_str());
	if (!fp) {
		std::cerr << "load error:" << file << std::endl;
		return false;
	}
	getline(fp, line); // skip headeer
	while (getline(fp, line)) {
		std::istringstream is(line);
		char sep;
		int id;
		double lat, lng, weight;

		is >> id >> sep;
		is >> lat >> sep;
		is >> lng >> sep;
		is >> weight >> sep;

		gift_t gift(lat, lng, weight);
		gift.id = id;
		gifts.push_back(gift);
	}
	fp.close();
	return true;
}

bool 
load_results(const std::string &file, gifts_t &results)
{
	std::ifstream fp;
	std::string line;

	results.clear();
	fp.open(file.c_str());
	if (!fp) {
		std::cerr << "load error:" << file << std::endl;
		return false;
	}
	getline(fp, line); // skip headeer
	while (getline(fp, line)) {
		std::istringstream is(line);
		gift_t result;
		char sep;

		is >> result.id >> sep;
		is >> result.trip_id >> sep;
		results.push_back(result);
	}
	fp.close();
	return true;
}

bool 
save_gifts(const std::string &file, const gifts_t &gifts)
{
	std::ofstream fp;
	std::string line;

	fp.open(file.c_str());
	if (!fp) {
		std::cerr << "load error:" << file << std::endl;
		return false;
	}
	fp << "GiftId,TripId\n";
	for (auto gift = gifts.begin(); gift != gifts.end(); ++gift) {
		fp << gift->id << "," << gift->trip_id << std::endl;
	}
	fp.close();
	return true;
}
