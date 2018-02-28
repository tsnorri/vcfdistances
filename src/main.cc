/*
 Copyright (c) 2018 Tuukka Norri
 This code is licensed under MIT license (see LICENSE for details).
 */

#include <cstdlib>
#include <dispatch/dispatch.h>
#include <iostream>
#include <libbio/assert.hh>
#include <unistd.h>
#include <string>
#include <vcfdistances/calculate_distances.hh>
#include <vcfdistances/types.hh>
#include <vector>

#ifdef __linux__
#include <pthread_workqueue.h>
#endif

#include "cmdline.h"


namespace lb	= libbio;
namespace vd	= vcfdistances;


namespace {
	
	vd::input_format input_format(enum_input_format const fmt)
	{
		switch (fmt)
		{
			case input_format_arg_uncompressed:
				return vd::input_format::UNCOMPRESSED;
				
			case input_format_arg_gzip:
				return vd::input_format::GZIP;
				
			case input_format__NULL:
			default:
				lb::fail("Unexpected value for input format.");
				return vd::input_format::UNCOMPRESSED; // Not reached.
		}
	}
}


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		exit(EXIT_FAILURE);
	
	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.
	std::cin.tie(nullptr);					// We don't require any input from the user.

#ifndef NDEBUG
	std::cerr << "Assertions have been enabled." << std::endl;
#endif
	
	// libdispatch on macOS does not need pthread_workqueue.
#ifdef __linux__
	pthread_workqueue_init_np();
#endif
	
	std::vector <std::string> inputs;
	for (std::size_t i(0); i < args_info.variants_given; ++i)
		inputs.emplace_back(args_info.variants_arg[i]);
	
	vd::calculate_distances(
		std::move(inputs),
		input_format(args_info.input_format_arg),
		args_info.output_sample_names_arg,
		args_info.output_hamming_arg,
		args_info.output_jaccard_arg,
		args_info.output_smd_arg
	);
	
	cmdline_parser_free(&args_info);
	
	dispatch_main();
	// Not reached b.c. pthread_exit() is eventually called in dispatch_main().
	return EXIT_SUCCESS;
}
