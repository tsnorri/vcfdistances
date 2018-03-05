/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCFDISTANCES_CALCULATE_DISTANCE_HH
#define VCFDISTANCES_CALCULATE_DISTANCE_HH

#include <vcfdistances/types.hh>
#include <vector>


namespace vcfdistances {
	
	void calculate_distances(
		std::vector <std::string> &&inputs,
		input_format const fmt,
		char const *sample_names_dst_path,
		char const *hamming_dst_path,
		char const *jaccard_dst_path,
		char const *smd_dst_path,
		char const *intersection_sizes_dst_path,
		char const *union_sizes_dst_path
	);
}

#endif
