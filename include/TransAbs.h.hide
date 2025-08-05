#ifndef TRANSABS_H
#define TRANSABS_H

#define _USE_MATH_DEFINES
#include <cmath>
#include "PSEvt/EventId.h"
#include "psana/Module.h"
#include <boost/shared_ptr.hpp>
#include "ImgAlgos/GlobalMethods.h"
#include "MsgLogger/MsgLogger.h"

#include "psddl_psana/camera.ddl.h"
#include "psddl_psana/evr.ddl.h"

#include <opencv2/core/core.hpp> // for types like cv::Mat
#include <fftw/fftw3.h>

#include <openmpi/mpi.h>

#include <vector>
#include <algorithm>
#include <list>
#include <fstream>	
#include <iostream>
#include <string>
#include <exception>

#include <boost/multi_array.hpp>

#include "CookieBox_pkg/dataops.h"

namespace CookieBox_pkg 
{

class TransAbs
{
	friend class CookieBox_mod;
	typedef boost::multi_array<long long,2> a2d_ll_t;
	typedef boost::multi_array_types::index_range range;

	public:

	TransAbs(void);
	~TransAbs(void);
	TransAbs & operator=( TransAbs rhs );
	TransAbs(const TransAbs & b);

		// annoyingly need to define the data getter function almost as an external function whose address can be passed to gsl
		static inline double cheb_return_data(double x,void* p) // this needs to be appropriate for the gsl-ref page 333
		{
			std::vector<double> * vecPtr = (std::vector<double> *)p;
			unsigned i = x * vecPtr->size()-1;
			return (*vecPtr)[i];
		}

	private:
		void deepcopy_data(const TransAbs & b);
		void swap(TransAbs & a, TransAbs & b);

	public:
		inline bool use(bool in){m_use = in;return m_use;}
		inline bool use(void){return m_use;}
		inline bool print(bool in){
			if (!m_use)
				return false;
			m_print = in;
			return m_print;
		}
		inline bool print(void)
		{
			if (!m_use)
				return false;
			return m_print;
		}
		
		// Dark removal operations //
		inline void dark_file(const std::string inname) {m_dark_filename = inname;}
		bool fill_dark_stage_position_list(const std::vector<std::string> in);
		bool in_dark_stage_list(void);
		inline bool remove_dark(const bool in){m_remove_dark = in;return m_remove_dark;}
		inline bool remove_dark(void){return m_remove_dark;}
		inline bool write_dark(const bool in){m_write_dark = in;return m_write_dark;}
		inline bool write_dark(void){return m_write_dark;}
		bool write_dark_file(void);
		bool fill_dark_vector(const std::string inname);
		bool fill_dark_vector(void);
		inline int npowers(void){return m_npowers;}
		inline int nderivs(void){return m_nderivs;}
		inline int npowers(const unsigned in){m_npowers = in; return m_npowers;}
		inline int nderivs(const unsigned in){
			m_nderivs = in; 
			m_records_r.resize(m_nderivs,(double*)NULL); 
			m_records_hc.resize(m_nderivs,(double*)NULL);
			return m_nderivs;}

		bool init(void);
		bool fft_init(void);
		bool fft_rows(void);


		inline unsigned refcode(const unsigned in){m_refcode = in; return m_refcode;}
		inline unsigned refcode(void){return m_refcode;}
		bool updateStage(const float in);
		bool updateStage_close(void);

		void srcStr(Source srcStr_in);
		void srcStrEvr(Source srcStr_in);

	private:bool fillavg(void);
		bool fill_fft(void);
		bool fill_accum_result(void);
		bool fft_accum_result(void);
		bool fft_result_vector_r(void);
		bool fill_result_vector_hc(void);
		bool fill_result_vector_r(void);
		bool fill_orthonormal_result_r(void);
		bool add_to_dark_vector(void);
		bool cheb_approx_result_vector_r(void);
		bool filter_result_vector_hc(void);
		bool fit_filter(std::vector<double> & s, std::vector<double> & n, std::vector<double> & f);
		bool backfft_result_vector_hc(void);
		//bool unwrap_phase(void);
		void clearFFT(void);
		template <typename T>
			bool getCentroid(boost::multi_array<T,2> & frame, T & centroid, T & width, std::vector<double> & weights);
		double process_phase_hc(double data[],const size_t sz);
			

		bool start_new_accumulation(void);
		bool close_new_accumulation(void);

		bool fill_accum(Event& evt,Env& env);
		bool isref(Event& evt, Env& env);


	public:	bool fill(Event& evt, Env& env);
		bool process_accum(void);

		bool print_out(void);
		bool print_out(const unsigned eventnum);
		bool print_out_fft(const unsigned eventnum);
		bool print_out_result(void);
		bool print_out_result_hc(void);
		bool print_out_result_r(void);

		bool print_out_accum(void);
		bool print_header(void);
		bool print_debug(void);
		bool print_header_stats(void);
		void open_file( std::string & filename );
		void close_file(void);

		Source m_srcStr;
		Source m_evr_src;

	private:
		std::ofstream * m_samplecodes;

		bool m_use,m_print,m_print_avg;
		bool m_remove_dark,m_write_dark;
		unsigned m_dark_contributions;

		unsigned m_nderivs,m_npowers;

		std::string m_filename;
		std::ofstream m_outfile;

		std::string m_stage_position;
		std::vector < std::string > m_stage_position_list;
		std::vector < std::string > m_dark_stage_position_list;
		unsigned m_refcode;

		fftw_plan m_plan_r2hc,m_plan_hc2r;
		bool m_fft_init;

		//std::vector< double > m_record_r,m_record_hc;
		double * m_record_r;
		double * m_record_hc;
		std::vector<double> m_dark_vector;
		std::string m_dark_filename;
		
		std::vector<double*> m_records_r; // for the sake of acumulating derivatives
		std::vector<double*> m_records_hc;
		std::vector< std::vector< std::vector < double > > > m_orthonorm_r; // this one taces both derivatives and powers to form Grahm Schmidt ortho basis

		std::vector< std::vector<unsigned> > m_data; // keep this for single images... just in case
		std::vector< a2d_ll_t * > m_data_accum; // accumulate and then average images here.
		std::vector< a2d_ll_t * > m_data_accum_ref; // accumulate references here.
		std::vector <long long > m_shots_accum; // count shots into accum here
		std::vector <long long > m_shots_accum_ref; // count shots into ref here
		unsigned m_accum_index;

		// change the FFT stuff to work on accumulated images //
		double * m_accum_result_r_fftwPtr;
		double * m_accum_result_hc_fftwPtr;

		boost::shared_ptr<Psana::Camera::FrameV1> m_srcPtr;
		boost::shared_ptr<Psana::EvrData::DataV3> m_evr_srcPtr;

		// 	MPI related	//
		unsigned m_root_rank,m_rank;
		unsigned m_mpi_size;
};
}


#endif
