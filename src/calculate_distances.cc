/*
 Copyright (c) 2018 Tuukka Norri
 This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/function_output_iterator.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <cstdint>
//#include <execution>
#include <libbio/dispatch_fn.hh>
#include <libbio/file_handling.hh>
#include <libbio/infix_ostream_fn.hh>
#include <libbio/vcf_reader.hh>
#include <sdsl/int_vector.hpp>
#include <thread>
#include <valarray>
#include <vcfdistances/calculate_distances.hh>
#include <vector>


namespace bnu	= boost::numeric::ublas;
namespace ios	= boost::iostreams;
namespace lb	= libbio;
namespace vd	= vcfdistances;


typedef bnu::triangular_matrix <std::uint32_t, bnu::lower> triangular_uint32_matrix;
typedef bnu::triangular_matrix <float, bnu::lower> triangular_float_matrix;


namespace {
	
	class calculate_task;
	class calculate_context;
	
	// Used the idea from https://lists.boost.org/ublas/2009/05/3509.php
	// Divide a vector element-wise, store the result to lhs.
	template <typename t_lhs_el, typename t_rhs_el>
	void divide_vector(
		bnu::matrix_row <t_lhs_el> &lhs,
		bnu::matrix_row <t_rhs_el> &rhs
	)
	{
		lhs = bnu::element_div(lhs, rhs);
	}
	
	template <typename t_lhs_matrix, typename t_rhs_matrix>
	void hadamard_division(
		t_lhs_matrix &lhs,
		t_rhs_matrix &rhs
	)
	{
		// Check the dimensions.
		assert(lhs.size1() == rhs.size1()); 
		assert(lhs.size2() == rhs.size2()); 
		
		for (std::size_t i(0), count(lhs.size1()); i < count; ++i)
		{
			// Create proxies.
			auto lhs_row(bnu::matrix_row(lhs, i));
			auto rhs_row(bnu::matrix_row(rhs, i));
			divide_vector(lhs_row, rhs_row);
		}
	}
	
	
	class output_handler
	{
	protected:
		lb::file_ostream								m_stream;
		
	public:
		virtual ~output_handler() {}
		
		lb::file_ostream &stream() { return m_stream; }
		
		virtual void output(
			std::uint32_t const total_variants,
			triangular_uint32_matrix const &mutual_set_values,
			triangular_uint32_matrix const &all_set_values,
			triangular_uint32_matrix const &different_values
		) = 0;
		
	protected:
		template <typename t_matrix>
		void output_matrix(t_matrix const &matrix)
		{
			for (auto it(matrix.begin1()), end(matrix.end1()); it != end; ++it)
			{
				std::copy(
					it.begin(),
					it.end(),
					boost::make_function_output_iterator(lb::infix_ostream_fn <typename t_matrix::value_type>(m_stream, '\t'))
				);
				
				m_stream << '\n';
			}
			
			m_stream << std::flush;
		}
	};
	
	
	class hamming_output_handler final : public output_handler
	{
	public:
		virtual void output(
			std::uint32_t const total_variants,
			triangular_uint32_matrix const &mutual_set_values,
			triangular_uint32_matrix const &all_set_values,
			triangular_uint32_matrix const &different_values
		) override;
	};
	
	
	class smd_output_handler final : public output_handler
	{
	public:
		virtual void output(
			std::uint32_t const total_haplotypes,
			triangular_uint32_matrix const &mutual_set_values,
			triangular_uint32_matrix const &all_set_values,
			triangular_uint32_matrix const &different_values
		) override;
	};
	
	
	class jaccard_output_handler final : public output_handler
	{
	public:
		virtual void output(
			std::uint32_t const total_variants,
			triangular_uint32_matrix const &mutual_set_values,
			triangular_uint32_matrix const &all_set_values,
			triangular_uint32_matrix const &different_values
		) override;
	};
	
	
	class compare_two_task final
	{
	protected:
		calculate_task									*m_calculate_task{nullptr};
		std::size_t										m_i{};
		std::size_t										m_j{};
		
		std::uint32_t									m_same_value_present_in_both{};
		std::uint32_t									m_any_value_present_in_either{};
		std::uint32_t									m_different_values{};
		
	public:
		compare_two_task(
			calculate_task &calculate_task,
			std::size_t const i,
			std::size_t const j
		):
			m_calculate_task(&calculate_task),
			m_i(i),
			m_j(j)
		{
		}
		
		void execute();
		std::uint32_t same_value_present_in_both() const { return m_same_value_present_in_both; }
		std::uint32_t any_value_present_in_either() const { return m_any_value_present_in_either; }
		std::uint32_t different_values() const { return m_different_values; }
		
	protected:
		std::uint32_t count_ones(sdsl::bit_vector const &bv) const;
	};
	
	
	class calculate_task final
	{
	protected:
		calculate_context								*m_ctx{nullptr};
		lb::dispatch_ptr <dispatch_group_t>				m_group{};
		std::uint32_t									m_file_idx{0};
		std::vector <sdsl::bit_vector>					m_haplotypes;
		std::vector <std::vector <std::uint8_t>>		m_multi_alt_haplotypes;
		
	public:
		calculate_task(
			calculate_context &ctx,
			std::uint32_t const idx,
			std::vector <sdsl::bit_vector> &&haplotypes,
			std::vector <std::vector <std::uint8_t>> &&multi_alt_haplotypes
		):
			m_ctx(&ctx),
			m_file_idx(idx),
			m_haplotypes(std::move(haplotypes)),
			m_multi_alt_haplotypes(std::move(multi_alt_haplotypes))
		{
		}
		
		std::vector <sdsl::bit_vector> const &haplotypes() const { return m_haplotypes; }
		std::vector <std::vector <std::uint8_t>> const &multi_alt_haplotypes() const { return m_multi_alt_haplotypes; }
		bool operator<(calculate_task const &other) const { return m_file_idx < other.m_file_idx; }
		
		void prepare();
		void execute_comparison_tasks();
	};
	
	
	class calculate_context final
	{
	protected:
		vd::input_format								m_input_format{};
		lb::vcf_reader::sample_name_map					m_sample_names;
		std::vector <std::uint8_t>						m_ploidy;
		std::uint32_t									m_haplotype_count{0};
		std::atomic_size_t								m_variant_count{0};
		std::atomic_size_t								m_calculation_tasks{0};
		
		std::mutex										m_task_mutex;
		std::set <calculate_task>						m_tasks;
		lb::dispatch_ptr <dispatch_semaphore_t>			m_task_semaphore{};
		
		lb::dispatch_ptr <dispatch_queue_t>				m_input_queue;
		lb::dispatch_ptr <dispatch_queue_t>				m_task_queue{};
		
		std::mutex										m_distance_mutex;
		triangular_uint32_matrix						m_mutual_set_values;	// For Jaccard.
		triangular_uint32_matrix						m_all_set_values;		// For Jaccard.
		triangular_uint32_matrix						m_different_values;		// For Hamming and SMD.
		
		lb::file_ostream								m_sample_names_stream;
		std::vector <std::unique_ptr <output_handler>>	m_output_handlers;
		
	public:
		calculate_context(vd::input_format const fmt):
			m_input_format(fmt)
		{
		}
		
		calculate_context() = delete;
		calculate_context(calculate_context const &) = delete;
		calculate_context(calculate_context &&) = delete;
		
		std::uint32_t haplotype_count() const { return m_haplotype_count; }
		lb::dispatch_ptr <dispatch_queue_t> &task_queue() { return m_task_queue; }
		lb::dispatch_ptr <dispatch_semaphore_t> &task_semaphore() { return m_task_semaphore; }
		
		void cleanup();
		
		void output_sample_names(char const *dst_path) { lb::open_file_for_writing(dst_path, m_sample_names_stream, false); }
		void output_hamming_distances(char const *dst_path) { add_output_handler <hamming_output_handler>(dst_path); }
		void output_jaccard_distances(char const *dst_path) { add_output_handler <jaccard_output_handler>(dst_path); }
		void output_smd(char const *dst_path) { add_output_handler <smd_output_handler>(dst_path); }
		
		void prepare();
		void calculate_distances(std::vector <std::string> &&inputs);
		void update_matrices(
			std::size_t i,
			std::size_t j,
			std::uint32_t same_value_present_in_both,
			std::uint32_t any_value_present_in_either,
			std::uint32_t different_value
		);
			
		void report_progress(std::size_t const file_idx, std::size_t const comparisons);
		void finish_calculation_task(calculate_task &task);
		
	protected:
		void cleanup_mt() { delete this; }

		template <typename t_handler>
		void add_output_handler(char const *dst_path);
		
		void process_path(std::size_t const i, std::string const &path);
		void process_with_reader(std::size_t const i, lb::vcf_reader &reader);
		void read_headers(lb::vcf_reader &reader, bool const is_first);
		void check_ploidy(lb::vcf_reader &reader, std::vector <std::uint8_t> &ploidy);
		void read_haplotypes(
			lb::vcf_reader &reader,
			std::vector <sdsl::bit_vector> &haplotypes,
			std::vector <std::vector <std::uint8_t>> &multi_alt_haplotypes
		);
		void start_calculation_task(
			std::size_t const file_idx,
			std::vector <sdsl::bit_vector> &&haplotypes,
			std::vector <std::vector <std::uint8_t>> &&multi_alt_haplotypes
		);
		void output_sample_names_to_stream();
		void output_distances();
		void report_progress_mt(std::size_t const file_idx, std::size_t const comparisons);
	};
	
	
	void calculate_context::cleanup()
	{
		auto queue(dispatch_get_main_queue());
		lb::dispatch_async_fn(queue, [this, queue](){
			cleanup_mt();
			exit(EXIT_SUCCESS);
		});
	}
	
	
	void calculate_context::prepare()
	{
		m_input_queue.reset(dispatch_queue_create("fi.iki.tsnorri.vcfdistances.input_queue", DISPATCH_QUEUE_SERIAL));
		m_task_queue.reset(dispatch_queue_create("fi.iki.tsnorri.vcfdistances.task_queue", DISPATCH_QUEUE_CONCURRENT));
		//m_task_queue.reset(dispatch_queue_create("fi.iki.tsnorri.vcfdistances.task_queue", DISPATCH_QUEUE_SERIAL));
		m_task_semaphore.reset(dispatch_semaphore_create(2 * std::thread::hardware_concurrency())); // FIXME: some other value?
	}
	
	
	void calculate_context::calculate_distances(std::vector <std::string> &&inputs)
	{
		m_calculation_tasks = inputs.size();
		
		// Use dispatch_async with each of the input paths.
		std::size_t i(0);
		for (auto &path : inputs)
		{
			lb::dispatch_async_fn(*m_input_queue,
				[
					this,
					i,
					path{std::move(path)}
				]()
				{
					process_path(i, path);
				}
			);
			
			++i;
		}
	}
	
	
	void calculate_context::process_path(
		std::size_t const i,
		std::string const &path
	)
	{
		lb::vcf_reader reader;
		
		// Output status in this thread for now.
		lb::dispatch_async_fn(dispatch_get_main_queue(), [path](){
			std::cerr << "Processing " << path << "…" << std::endl;
		});
		
		if (vd::input_format::UNCOMPRESSED == m_input_format)
		{
			lb::vcf_stream_input <lb::file_istream> vcf_stream;
			lb::open_file_for_reading(path.c_str(), vcf_stream.input_stream());
			reader.set_input(vcf_stream);
			process_with_reader(i, reader);
		}
		else if (vd::input_format::GZIP == m_input_format)
		{
			lb::file_istream gzip_stream;
			lb::open_file_for_reading(path.c_str(), gzip_stream);
			
			lb::vcf_stream_input <ios::filtering_istream> vcf_stream;
			vcf_stream.input_stream().push(ios::gzip_decompressor());
			vcf_stream.input_stream().push(gzip_stream);
			reader.set_input(vcf_stream);
			process_with_reader(i, reader);
		}
	}
	
	
	void calculate_context::process_with_reader(std::size_t const i, lb::vcf_reader &reader)
	{
		std::vector <sdsl::bit_vector> haplotypes;
		std::vector <std::vector <std::uint8_t>> multi_alt_haplotypes;
		
		read_headers(reader, 0 == i);
		
		if (0 == i && m_sample_names_stream.is_open())
			output_sample_names_to_stream();
		
		if (m_output_handlers.size())
		{
			read_haplotypes(reader, haplotypes, multi_alt_haplotypes);
			// haplotypes now contains a bit vector for each haplotype and
			// multi_alt_haplotypes a vector for each variant that has more than one ALT.
			start_calculation_task(i, std::move(haplotypes), std::move(multi_alt_haplotypes));
		}
	}
	
	
	void calculate_context::read_headers(lb::vcf_reader &reader, bool const is_first)
	{
		reader.read_header();
		
		// Check that the sample names match.
		// Also check the ploidy.
		if (is_first)
		{
			m_sample_names = reader.sample_names();
			check_ploidy(reader, m_ploidy);
			m_haplotype_count = std::accumulate(m_ploidy.cbegin(), m_ploidy.cend(), static_cast <uint64_t>(0));
			
			// Allocate space for the distances.
			m_mutual_set_values.resize(m_haplotype_count, m_haplotype_count, false);
			m_all_set_values.resize(m_haplotype_count, m_haplotype_count, false);
			m_different_values.resize(m_haplotype_count, m_haplotype_count, false);
			
			// Zero the matrices.
			// FIXME: there is bound to be a faster way to do this but couldn't find one.
			for (std::size_t i(0); i < m_haplotype_count; ++i)
			{
				for (std::size_t j(0); j <= i; ++j)
				{
					m_mutual_set_values(i, j) = 0;
					m_all_set_values(i, j) = 0;
					m_different_values(i, j) = 0;
				}
			}
			
		}
		else
		{
			lb::always_assert(reader.sample_names() == m_sample_names);
			
			std::vector <std::uint8_t> ploidy;
			check_ploidy(reader, ploidy);
			lb::always_assert(
				std::equal(
					/* std::execution::parallel_unsequenced_policy, */ // FIXME: enable.
					std::begin(ploidy),
					std::end(ploidy),
					std::begin(m_ploidy),
					std::end(m_ploidy)
				)
			);
		}
	}
	
	
	void calculate_context::read_haplotypes(
		lb::vcf_reader &reader,
		std::vector <sdsl::bit_vector> &haplotypes,
		std::vector <std::vector <std::uint8_t>> &multi_alt_haplotypes
	)
	{
		lb::dispatch_async_fn(dispatch_get_main_queue(), [](){
			std::cerr << "Counting lines…" << std::flush;
		});
		{
			reader.reset();
			reader.set_parsed_fields(lb::vcf_field::CHROM);
			bool should_continue(false);
			do {
				reader.fill_buffer();
				should_continue = reader.parse([](lb::transient_variant const &var) -> bool { return true; });
			} while (should_continue);
		}
		
		auto const lines(reader.lineno() - reader.last_header_lineno());
		m_variant_count += lines;
		lb::dispatch_async_fn(dispatch_get_main_queue(), [lines](){
			std::cerr << " got " << lines << '.' << std::endl;
		});
		
		// Allocate bitvectors for each sample.
		// Also have a buffer to speed up reading.
		lb::dispatch_async_fn(dispatch_get_main_queue(), [](){
			std::cerr << "Allocating bitvectors…" << std::endl;
		});
		haplotypes.resize(m_haplotype_count);
		std::valarray <std::uint64_t> gt_buffer(static_cast <std::uint64_t>(0), m_haplotype_count);
		for (auto &vec : haplotypes)
		{
			sdsl::bit_vector bv(lines, 0);
			vec = std::move(bv);
		}
		
		std::uint32_t single_alt_variants(0);
		
		lb::dispatch_async_fn(dispatch_get_main_queue(), [](){
			std::cerr << "Reading variants…" << std::endl;
		});
		{
			reader.reset();
			reader.set_parsed_fields(lb::vcf_field::ALL);
			bool should_continue(false);
			std::size_t i(0); // Bit index (from 0 to 63).
			std::size_t k(0); // Word index in hapotypes.
			std::size_t lineno(0);
			do {
				reader.fill_buffer();
				should_continue = reader.parse(
					[
						this,
						&single_alt_variants,
						&haplotypes,
						&multi_alt_haplotypes,
						&gt_buffer,
						&i,
						&k,
						&lineno
					](lb::transient_variant const &var) -> bool {
					
					// If there is just one ALT, use the bit vectors. Otherwise, use the std::vector.
					if (1 < var.alts().size())
					{
						// Add a new vector.
						auto &vec(multi_alt_haplotypes.emplace_back(m_haplotype_count, 0));
						
						std::size_t j(0);
						for (auto const &sample : var.samples())
						{
							for (auto const &gt : sample.get_genotype())
								vec[j++] = gt.alt;
						}
					}
					else
					{
						// Read the value of each genotype, set the corresponding bit in buffer if needed.
						{
							std::size_t j(0);
							for (auto const &sample : var.samples())
							{
								for (auto const &gt : sample.get_genotype())
								{
									uint64_t mask(gt.alt);
									mask <<= i;
									gt_buffer[j++] |= mask;
								}
							}
						}
					
						++i;
						++single_alt_variants;
					
						// Check if we have a whole word.
						if (64 == i)
						{
							for (std::size_t j(0); j < m_haplotype_count; ++j)
							{
								auto const val(gt_buffer[j]);
								*(haplotypes[j].data() + k) = val;
							}

							i = 0;
							++k;
							gt_buffer = 0;
						}
					}
					
					++lineno;
					if (0 == lineno % 10000)
					{
						lb::dispatch_async_fn(dispatch_get_main_queue(), [lineno](){
							std::cerr << "Handled " << lineno << " lines…" << std::endl;
						});
					}
					
					return true;
				});
			} while (should_continue);
			
			// Copy the final word if needed.
			if (i)
			{
				for (std::size_t j(0); j < m_haplotype_count; ++j)
				{
					auto const val(gt_buffer[j]);
					*(haplotypes[j].data() + k) = val;
				}
			}
		}
	}
	
	
	void calculate_context::check_ploidy(lb::vcf_reader &reader, std::vector <std::uint8_t> &ploidy)
	{
		auto const sample_count(reader.sample_names().size());
		ploidy.resize(sample_count, 0);
		
		reader.reset();
		reader.set_parsed_fields(lb::vcf_field::ALL);
		
		reader.fill_buffer();
		if (!reader.parse([this, &reader, &ploidy](lb::transient_variant const &var) -> bool {
			for (auto const &kv : reader.sample_names())
			{
				auto const sample_no(kv.second);
				assert(sample_no);
				auto const &sample(var.sample(sample_no));
				auto const current_ploidy(sample.ploidy());
				ploidy[sample_no - 1] = current_ploidy;
			}
			
			return false;
		}))
		{
			lb::fail("Unable to read the first variant");
		}
	}
	
	
	void calculate_context::start_calculation_task(
		std::size_t const file_idx,
		std::vector <sdsl::bit_vector> &&haplotypes,
		std::vector <std::vector <std::uint8_t>> &&multi_alt_haplotypes
	)
	{
		calculate_task task(*this, file_idx, std::move(haplotypes), std::move(multi_alt_haplotypes));
		task.prepare();
		
		std::set <calculate_task>::iterator it;
		
		{
			std::lock_guard <std::mutex> guard(m_task_mutex);
			auto res(m_tasks.emplace(std::move(task)));
			lb::always_assert(res.second);
			it = res.first;
		}
		
		// execute_comparison_tasks is not going to affect the comparison order.
		const_cast <calculate_task &>(*it).execute_comparison_tasks();
	}
	
	
	void calculate_context::finish_calculation_task(calculate_task &task)
	{
		{
			std::lock_guard <std::mutex> guard(m_task_mutex);
			m_tasks.erase(task);
		}
		
		if (0 == --m_calculation_tasks)
		{
			output_distances();
			cleanup();
		}
	}
	
	
	void calculate_context::update_matrices(
		std::size_t i,
		std::size_t j,
		std::uint32_t same_value_present_in_both,
		std::uint32_t any_value_present_in_either,
		std::uint32_t different_value
	)
	{
		std::lock_guard <std::mutex> guard(m_distance_mutex);
		m_mutual_set_values(i, j) += same_value_present_in_both;
		m_all_set_values(i, j) += any_value_present_in_either;
		m_different_values(i, j) += different_value;
	}
	
	
	void calculate_context::output_sample_names_to_stream()
	{
		// Order guaranteed b.c. m_sample_names is an std::map.
		std::size_t sample_idx(0);
		for (auto const &kv : m_sample_names)
		{
			auto const idx(kv.second);
			assert(idx);
			auto const sample_ploidy(m_ploidy[idx - 1]);
			for (std::size_t i(0); i < sample_ploidy; ++i)
				m_sample_names_stream << kv.first << '\t' << (sample_idx++) << '\n';
		}
		m_sample_names_stream << std::flush;
	}
	
	
	void calculate_context::output_distances()
	{
		for (auto &handler : m_output_handlers)
			handler->output(m_variant_count, m_mutual_set_values, m_all_set_values, m_different_values);
	}
	
	
	template <typename t_handler>
	void calculate_context::add_output_handler(char const *dst_path)
	{
		auto &handler(m_output_handlers.emplace_back(new t_handler()));
		lb::open_file_for_writing(dst_path, handler->stream(), false);
	}
	
	
	void calculate_context::report_progress(std::size_t const file_idx, std::size_t const comparisons)
	{
		lb::dispatch_async_fn(dispatch_get_main_queue(), [this, file_idx, comparisons](){
			report_progress_mt(file_idx, comparisons);
		});
	}
	
	
	void calculate_context::report_progress_mt(std::size_t const file_idx, std::size_t const comparisons)
	{
		std::cerr << "File " << (1 + file_idx) << ": done " << comparisons << " comparisons." << std::endl;
	}
	
	
	void calculate_task::prepare()
	{
		m_group.reset(dispatch_group_create());
	}
	
	
	void calculate_task::execute_comparison_tasks()
	{
		auto queue(*(m_ctx->task_queue()));
		auto sema(*(m_ctx->task_semaphore()));
		std::size_t task_count(0);
		
		// Make sure that the group doesn't get notified immediately.
		dispatch_group_enter(*m_group);
		
		for (std::size_t i(0), haplotype_count(m_ctx->haplotype_count()); i < haplotype_count; ++i)
		{
			for (std::size_t j(0); j <= i; ++j)
			{
				++task_count;
				dispatch_semaphore_wait(sema, DISPATCH_TIME_FOREVER);
				lb::dispatch_group_async_fn(*m_group, queue,
					[this, i, j, task_count, sema](){
						// Execute a compare_two_task for (i, j).
						compare_two_task task(*this, i, j);
						task.execute();
						m_ctx->update_matrices(
							i,
							j,
							task.same_value_present_in_both(),
							task.any_value_present_in_either(),
							task.different_values()
						);
						
						if (0 == task_count % 100000)
							m_ctx->report_progress(m_file_idx, task_count);
						
						dispatch_semaphore_signal(sema);
					}
				);
			}
		}
		
		dispatch_group_leave(*m_group);
		lb::dispatch_group_notify_fn(*m_group, dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), [this](){
			m_ctx->finish_calculation_task(*this);
		});
	}
	
	
	void compare_two_task::execute()
	{
		auto const &haplotypes(m_calculate_task->haplotypes());
		auto const &multi_alt_haplotypes(m_calculate_task->multi_alt_haplotypes());
		
		// Make the comparisons.
		auto const &lhs_bv(haplotypes[m_i]);
		auto const &rhs_bv(haplotypes[m_j]);
		sdsl::bit_vector or_target_bv(lhs_bv);
		sdsl::bit_vector and_target_bv(lhs_bv);
		sdsl::bit_vector xor_target_bv(lhs_bv);
		or_target_bv |= rhs_bv;
		and_target_bv &= rhs_bv;
		xor_target_bv ^= rhs_bv;
		
		// Count the ones in the result.
		auto const or_ones(count_ones(or_target_bv));
		auto const and_ones(count_ones(and_target_bv));
		auto const xor_ones(count_ones(xor_target_bv));
		
		// Handle the multi-alt vectors.
		std::uint32_t same_nonzero_mv(0);
		std::uint32_t any_two_nonzero_mv(0);
		std::uint32_t any_different_mv(0);
		for (auto const &mv_vec : multi_alt_haplotypes)
		{
			if (mv_vec[m_i] && mv_vec[m_j])
			{
				++any_two_nonzero_mv;
				
				if (mv_vec[m_i] == mv_vec[m_j])
					++same_nonzero_mv;
			}
			
			if (mv_vec[m_i] != mv_vec[m_j])
				++any_different_mv;
		}
		
		m_same_value_present_in_both = and_ones + same_nonzero_mv;
		m_any_value_present_in_either = or_ones + any_two_nonzero_mv;
		m_different_values = xor_ones + any_different_mv;
	}
	
	
	std::uint32_t compare_two_task::count_ones(sdsl::bit_vector const &bv) const
	{
		auto const *data(bv.data());
		auto const bit_count(bv.size());
		auto const full_word_count(bit_count / 64);
		std::uint32_t retval(0);
		
		// The bit vector has 8 additional bytes if 0 == bit_count % 64.
		// In any case the padding is set to zero.
		// See https://github.com/xxsds/sdsl-lite/blob/master/include/sdsl/memory_management.hpp#L758
		// (sdsl::memory_manager::resize())
		for (std::size_t i(0); i <= full_word_count; ++i)
			retval += sdsl::bits_impl <>::cnt(data[i]);
		
		return retval;
	}
	
	
	void hamming_output_handler::output(
		std::uint32_t const total_variants,
		triangular_uint32_matrix const &mutual_set_values,
		triangular_uint32_matrix const &all_set_values,
		triangular_uint32_matrix const &different_values
	)
	{
		output_matrix(different_values);
	}
	
	
	void smd_output_handler::output(
		std::uint32_t const total_variants,
		triangular_uint32_matrix const &mutual_set_values,
		triangular_uint32_matrix const &all_set_values,
		triangular_uint32_matrix const &different_values
	)
	{
		triangular_float_matrix mat(1.0 * different_values / total_variants);
		output_matrix(mat);
	}
	
	
	void jaccard_output_handler::output(
		std::uint32_t const total_variants,
		triangular_uint32_matrix const &mutual_set_values,
		triangular_uint32_matrix const &all_set_values,
		triangular_uint32_matrix const &different_values
	)
	{
		triangular_float_matrix mat(mutual_set_values * -1.0);
		hadamard_division(mat, all_set_values);
		for (std::size_t i(0), count(mat.size1()); i < count; ++i)
		{
			for (std::size_t j(0); j <= i; ++j)
				mat(i, j) += 1.0;
		}
		
		output_matrix(mat);
	}
}


namespace vcfdistances {
	
	void calculate_distances(
		std::vector <std::string> &&inputs,
		input_format const fmt,
		char const *sample_names_dst_path,
		char const *hamming_dst_path,
		char const *jaccard_dst_path,
		char const *smd_dst_path
	)
	{
		calculate_context *ctx(new calculate_context(fmt));
		
		if (sample_names_dst_path)
			ctx->output_sample_names(sample_names_dst_path);
		
		if (hamming_dst_path)
			ctx->output_hamming_distances(hamming_dst_path);
		
		if (jaccard_dst_path)
			ctx->output_jaccard_distances(jaccard_dst_path);
		
		if (smd_dst_path)
			ctx->output_smd(smd_dst_path);
		
		ctx->prepare();
		ctx->calculate_distances(std::move(inputs));
	}
}
