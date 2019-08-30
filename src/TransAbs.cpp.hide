#include "CookieBox_pkg/TransAbs.h"
#include "CookieBox_pkg/CookieBox_mod.h"
#include "CookieBox_pkg/dataops.h"
#include <fstream>
#include <iostream>
#include <assert.h>
#include <iterator>
#include <algorithm>
#include <vector>
#include <complex>
#include <exception>
#include <functional> // also needed for bind1st in the std::stransform() call 
#include <boost/math/constants/constants.hpp>

#include <opencv2/core/core.hpp> // for eroding the mask image
#include <opencv2/imgproc/imgproc.hpp> // for filtering/bluring, and thresholding the mask image

#include <gsl/gsl_math.h>
#include <gsl/gsl_chebyshev.h>
#include <gsl/gsl_fit.h>



using namespace CookieBox_pkg;

namespace CookieBox_pkg{
	TransAbs::TransAbs(void)
	: m_use(false)
	, m_print(false)
	, m_plan_r2hc(NULL)
	, m_plan_hc2r(NULL)
	, m_remove_dark(false)
	, m_write_dark(false)
	, m_dark_contributions(0)
	, m_accum_result_r_fftwPtr(NULL)
	, m_accum_result_hc_fftwPtr(NULL)
	, m_nderivs(1)
	, m_npowers(1)
	{
		m_samplecodes = new std::ofstream("samplecodes.out",std::ios::out);
		(*m_samplecodes) << "# opening file initially" << std::endl;
		m_samplecodes->close();
		m_records_r.resize(m_nderivs+1,(double*)NULL);
		m_records_hc.resize(m_nderivs+1,(double*)NULL);
	}

	TransAbs::~TransAbs(void){
		if (m_outfile.is_open())
			m_outfile.close();
		if ((m_plan_r2hc != NULL) || (m_plan_hc2r != NULL) )
			clearFFT();
		if (m_data_accum.size() > 0){
			while (m_data_accum.size()>0){
				delete m_data_accum.back();
				m_data_accum.pop_back();
			}
		}
		if (m_data_accum_ref.size() > 0){
			while (m_data_accum_ref.size() > 0){
				delete m_data_accum_ref.back();
				m_data_accum_ref.pop_back();
			}
		}
		m_samplecodes->close();
	}
	void TransAbs::clearFFT(void)
	{
		if (m_plan_r2hc != NULL)
			fftw_destroy_plan(m_plan_r2hc);
		if (m_plan_hc2r != NULL)
			fftw_destroy_plan(m_plan_hc2r);
		if (m_accum_result_r_fftwPtr != NULL){
			fftw_free(m_accum_result_r_fftwPtr);
			m_accum_result_r_fftwPtr = NULL;
		}
		if (m_accum_result_hc_fftwPtr != NULL){
			fftw_free(m_accum_result_hc_fftwPtr);
			m_accum_result_hc_fftwPtr = NULL;
		}
		for (unsigned i=0;i<m_records_r.size();++i){
			if (m_records_r[i] != NULL){
				fftw_free(m_records_r[i]);
				m_records_r[i] = NULL;
			}
		}
		for (unsigned i=0;i<m_records_hc.size();++i){
			if (m_records_hc[i] != NULL){
				fftw_free(m_records_hc[i]);
				m_records_hc[i] = NULL;
			}
		}
		if ((m_plan_r2hc != NULL) || (m_plan_hc2r != NULL) )
			fftw_cleanup();
	}
	TransAbs & TransAbs::operator=( TransAbs rhs ){ // the compiler is making a pass by copy here //
		swap(*this,rhs); // using the swap makes this exception safe... 
		return *this;
	}
	TransAbs::TransAbs(const TransAbs & b)
		: m_srcPtr(b.m_srcPtr)
		  , m_use(b.m_use)
		  , m_srcStr(b.m_srcStr)
		  , m_print(b.m_print)
		  , m_filename(b.m_filename)
		  , m_root_rank(b.m_root_rank)
		  , m_rank(b.m_rank)
		  , m_mpi_size(b.m_mpi_size)
		  , m_remove_dark(b.m_remove_dark)
		  , m_write_dark(b.m_write_dark)
		  , m_dark_contributions(b.m_dark_contributions)
	{ // not so sure this is exception safe

		m_data = b.m_data;

		m_filename += ".copy";
		if (b.m_outfile.is_open())
			m_outfile.open(m_filename.c_str(),std::ios::out);
	}
	void TransAbs::deepcopy_data(const TransAbs & b)
	{
		m_data.resize(b.m_data.size());
		for ( unsigned i=0;i<b.m_data.size();++i){
			m_data[i].resize(b.m_data[i].size());
			for (unsigned j=0;j<b.m_data[i].size();++j)
				m_data[i][j] = b.m_data[i][j];
		}
	}
	void TransAbs::swap(TransAbs & a, TransAbs & b)
	{
		std::swap(a.m_use,b.m_use);
		std::swap(a.m_root_rank, b.m_root_rank);
		std::swap(a.m_rank, b.m_rank);
		std::swap(a.m_mpi_size, b.m_mpi_size);
		std::swap(a.m_srcPtr, b.m_srcPtr);
		std::swap(a.m_srcStr,b.m_srcStr);
		std::swap(a.m_remove_dark,b.m_remove_dark);
		std::swap(a.m_write_dark,b.m_write_dark);
		std::swap(a.m_dark_contributions,b.m_dark_contributions);
		
		a.m_data.swap(b.m_data);

		std::swap(a.m_print,b.m_print);
		std::swap(a.m_filename,b.m_filename);
		if (a.m_outfile.is_open() && b.m_outfile.is_open()){
			a.m_outfile.close();
			b.m_outfile.close();
			a.m_outfile.open(a.m_filename.c_str(),std::ios::app);
			b.m_outfile.open(b.m_filename.c_str(),std::ios::app);
			return;
		}
		if (a.m_outfile.is_open() && !b.m_outfile.is_open()){
			a.m_outfile.close();
			b.m_outfile.open(b.m_filename.c_str(),std::ios::app);
			return;
		}
		if (!a.m_outfile.is_open() && b.m_outfile.is_open()){
			b.m_outfile.close();
			a.m_outfile.open(a.m_filename.c_str(),std::ios::app);
			return;
		}

	}


	bool TransAbs::init(void){
		std::cerr << "entered TransAbs::init(void)" << std::endl;
		m_plan_hc2r = m_plan_r2hc = NULL;
		m_fft_init = false;
		m_rank = MPI::COMM_WORLD.Get_rank();
		return true;
	}

	bool TransAbs::fft_init(void)
	{
		//std::cerr << "entered bool TransAbs::fft_init(void)" << std::endl;
		if (!m_fft_init){
			try {
				//std::cerr << "inside try: of bool TransAbs::fft_init(void)\n\tm_data.size() = " << m_data.size();
				//std::cerr << "\n\t m_data.front().size() = " << m_data.front().size() << std::endl;
				if ( m_data.size() == 0 ){
					std::cerr << "m_data.size() == 0 in bool TransAbs::fft_init(void)" << std::endl;
					throw;
				}
				unsigned sz = m_data.front().size();
				m_records_r.clear();
				m_records_hc.clear();
				// this is for hte 0th derivative that must get computed
				m_records_r.push_back((double *) fftw_malloc(sizeof(double) * m_data.front().size()));
				m_records_hc.push_back((double *) fftw_malloc(sizeof(double) * m_data.front().size()));
				for (unsigned i=0;i<m_nderivs;++i){
					m_records_r.push_back((double *) fftw_malloc(sizeof(double) * m_data.front().size()));
					m_records_hc.push_back((double *) fftw_malloc(sizeof(double) * m_data.front().size()));
				}


				m_plan_r2hc = fftw_plan_r2r_1d(m_data.front().size(),
						m_records_r.front(),
						m_records_hc.front(), // NOTE: this cannot be changed to newdata... hc2r destroys 
						FFTW_R2HC,
						FFTW_MEASURE);
				//std::cerr << "after 1st plan def: of bool TransAbs::fft_init(void)" << std::endl;
				m_plan_hc2r = fftw_plan_r2r_1d(m_data.front().size(),
						m_records_hc.front(),
						m_records_r.front(),
						FFTW_HC2R,
						FFTW_MEASURE);
				fftw_print_plan(m_plan_r2hc);
				std::cout << "\n";
				fftw_print_plan(m_plan_hc2r);
				std::cout << "\n" << std::flush;



			}
			catch (std::exception & e) {
				std::cerr << "Failed in TransAbs::fft_init() method: " << e.what() << std::endl;
				return false;
			}
			catch (...) {
				std::cerr << "Failed in TransAbs::fft_init() method with default (m_data.size() == 0): " << std::endl;
				return false;
			}
				
			m_fft_init = true;
		}
		return m_fft_init;
	}

	void TransAbs::srcStrEvr(Source srcStr_in)
	{
		m_evr_src = srcStr_in;
	}
	void TransAbs::srcStr(Source srcStr_in){
		m_srcStr = srcStr_in;
	}

	bool TransAbs::fill(Event& evt, Env& env){
		if (write_dark() && !in_dark_stage_list()){
			return true;// returning early to avoid processing
		}

		m_srcPtr = evt.get(m_srcStr);
		if (!m_srcPtr){
			std::cerr << "couldn't get TransAbs image" << std::endl;
			return false; 
		}
		if (m_data.size() < m_srcPtr->height() || m_data.back().size() < m_srcPtr->width()){
			m_data.resize(m_srcPtr->height());
			for (unsigned i=0;i<m_data.size();++i)
				m_data.at(i).resize(m_srcPtr->width(),0);
		}
		for (unsigned i=0;i<m_srcPtr->height();i++){
			for (unsigned j=0;j<m_srcPtr->width();j++){
				m_data.at(i).at(j) =  m_srcPtr->data16()[i][j];
			}
		}

		if (!fill_accum(evt,env)){
			std::cerr << "failed to fill_accum() in TransAbs::fill() method" << std::endl;
			return false;
		}
		return true;	
	}
	bool TransAbs::fill_accum(Event& evt,Env& env)
	{
		a2d_ll_t * accumPtr = NULL;
		long long * shotsPtr = NULL;
		// std::cerr << "Made it here in TransAbs::fill_accum() method" << std::endl;
		// HERE try to do if address of m_data_accum is NULL then intialize it to the size of the image (puch back)
		if (isref(evt,env)){
			// HERE HERE is where I fail. OK here I have not yet initialized the boost::multi_array
			if (m_data_accum_ref[m_accum_index] == NULL){
				m_data_accum_ref[m_accum_index] = new a2d_ll_t(boost::extents [m_data.size()][m_data.front().size()]);
			}
			accumPtr = m_data_accum_ref[m_accum_index];
			shotsPtr = &(m_shots_accum_ref[m_accum_index]);
		} else {
			if (m_data_accum[m_accum_index] == NULL){
				m_data_accum[m_accum_index] = new a2d_ll_t(boost::extents [m_data.size()][m_data.front().size()]);
			}
			//std::cerr << "\n\n m_data_accum[m_accum_index] = " << m_data_accum[m_accum_index] << "\n\n" << std::endl;
			accumPtr = m_data_accum[m_accum_index];
			shotsPtr = &(m_shots_accum[m_accum_index]);
		}

		//std::cerr << "starting to fill the accumPtr in fill_accum() method" << std::endl;

		for (unsigned i=0;i<accumPtr->shape()[0];++i){
			for (unsigned j=0;j<accumPtr->shape()[1];++j){
				(*accumPtr)[i][j] += m_data[i][j];
			}
		}
		++(*shotsPtr);
		//std::cerr << "made it through TransAbs::fill_accum() " << std::endl;
		return true;
	}
	bool TransAbs::isref(Event& evt,Env& env)
	{
		shared_ptr<Psana::EvrData::DataV3> srcPtr = evt.get(m_evr_src);
		ndarray<const Psana::EvrData::FIFOEvent, 1> eventList = srcPtr->fifoEvents();
		try{
			//(*m_samplecodes) << "eventList[ ] = ";
			for (unsigned i=0;i<eventList.size();++i){
				//(*m_samplecodes) << eventList[i].eventCode() << "\t";
				if (eventList[i].eventCode() == m_refcode){
				//	(*m_samplecodes) << "... found, returning true" << std::endl;
					return true;
				}
			}
			//(*m_samplecodes) << std::endl;
		} catch (std::exception & e){
			std::cerr << "failed to get Evr data in TransAbs::isref() method, e = " << e.what() << std::endl;
			return false;
		}
		return false;
	}


	// File interface methods //
	void TransAbs::open_file( std::string & filename ){
		m_filename = filename;
		m_outfile.open(filename.c_str(),std::ios::out);
	}
	void TransAbs::close_file(void){
		if (m_outfile.is_open())
			m_outfile.close();
	}
	bool TransAbs::print_header(void){
		m_outfile << "#xtcav images";
		m_outfile << "\n";
		return true;
	}
	bool TransAbs::print_out(void){
		if (!m_outfile.is_open())
			return false;
		for (unsigned i=0;i<m_data.size();++i){
			for (unsigned j=0;j<m_data[i].size();++j){
				m_outfile << m_data[i][j] << "\t";
			}
			m_outfile << "\n";
		}
		m_outfile << "\n";
		

		return true;
	}

	bool TransAbs::print_out_result_hc(void)
	{
		unsigned sz = m_data_accum_ref[m_accum_index]->shape()[1];
		if (m_shots_accum[m_accum_index] == 0){
			std::cerr << "m_shots_accum[m_accum_index] == 0, failing TransAbs::print_out_result_hc() method" << std::endl;
			return false;
		}
		std::string filename(m_filename);
		filename += ".result.hc." + m_stage_position_list[m_accum_index];
		std::ofstream outfile(filename.c_str(),std::ios::out);
		if (!outfile.is_open()){
			std::cerr << "Failed to open accumout file in TransAbs::print_out_result() method" << std::endl;
			return false;
		}
		std::cerr << "Attempting to write to file " << filename << std::endl;
		outfile << "# result fft vector out, only for the 0th derivative\n";
		std::complex<double> z(m_records_hc[0][0],0.);
		outfile << std::abs(z) << "\t" << std::arg(z) << "\n";
		for (unsigned i=1;i<sz/2;++i){
			z = std::complex<double>(m_records_hc[0][i],m_records_hc[0][sz-i]);
			outfile << std::abs(z) << "\t" << std::arg(z) << "\n";
		}
		z = std::complex<double>(0.,m_records_hc[0][sz/2]);// this might be wrong, Nyquist frequency is either fully real or fully imag
		outfile << std::abs(z) << "\t" << std::arg(z) << "\n";
		outfile.close();
		return true;
	}
	bool TransAbs::print_out_result(void)
	{
		if (m_shots_accum[m_accum_index] == 0){
			std::cerr << "m_shots_accum[m_accum_index] == 0, failing TransAbs::print_out_result() method" << std::endl;
			return false;
		}
		std::string filename(m_filename);
		filename += ".result.image." + m_stage_position_list[m_accum_index];
		std::cerr << "filename(m_filename) in TransAbs::print_out_result() is " << filename << std::endl;
		std::ofstream outfile(filename.c_str(),std::ios::out);
		if (!outfile.is_open()){
			std::cerr << "Failed to open accumout file in TransAbs::print_out_result() method" << std::endl;
			return false;
		}
		outfile << "# result images out\n";
		for (unsigned i=0;i<m_data_accum[m_accum_index]->shape()[0];++i){
			for (unsigned j=0;j<m_data_accum[m_accum_index]->shape()[1];++j){
				double value = m_accum_result_r_fftwPtr[i*m_data_accum[m_accum_index]->shape()[1] + j];
				outfile << value << "\t";
			}
			outfile << "\n";
		}
		outfile.close();
		return true;
	}
	bool TransAbs::fill_accum_result(void)
	{
		try {
			boost::multi_array_ref<long long , 2> frame(*(m_data_accum[m_accum_index]));
			boost::multi_array_ref<long long , 2> reference(*(m_data_accum_ref[m_accum_index]));

			if (m_accum_result_r_fftwPtr == NULL){
				//m_accum_result_r_fftwPtr = (double*)fftw_malloc(sizeof(double) * frame.shape()[0] * frame.shape()[1]);
				m_accum_result_r_fftwPtr = fftw_alloc_real((frame.shape()[0] * frame.shape()[1]));
			}
			//boost::multi_array_ref<double,2> accum(m_accum_result_r_fftwPtr,m_data_accum[m_accum_index]->shape()); // HERE having trouble with multi_array_ref()
			for (unsigned i=0;i<frame.shape()[0];++i){
				for (unsigned j=0;j<frame.shape()[1];++j){
					double value = frame[i][j] / m_shots_accum[m_accum_index];
					double ref = reference[i][j] / m_shots_accum_ref[m_accum_index];
					if (ref != 0. && value != 0.){
						value /= ref;
						value -= 1.;
					} else {
						value = 0.;
					}
					m_accum_result_r_fftwPtr[i*frame.shape()[1] + j] = value;
				}
			}
		} catch (std::exception & e) {
			std::cerr << "failed in bool TransAbs::fill_accum_result(void) with error: " << e.what() << std::endl;
			return false;
		}
		return true;
	}

	bool TransAbs::add_to_dark_vector(void)
	{
		// m_records_r and _hc are indeed both fftw aligned 
		boost::multi_array_ref<long long , 2> frame(*(m_data_accum[m_accum_index]));
		unsigned sz = frame.shape()[1];

		try {
			if (m_records_r.front() == NULL)
				throw;
			if (m_dark_vector.size() != sz){
				m_dark_contributions = 0;
				m_dark_vector.resize(sz);
				std::fill(m_dark_vector.begin(),m_dark_vector.end(),0.);
			}
			for (unsigned i=0;i<sz;++i)
				m_dark_vector[i] += m_records_r.front()[i];
			++m_dark_contributions;
		} catch (std::exception & e) {
			std::cerr << "Failed in bool TransAbs::fill_dark_vector(void) with error: " << e.what() << std::endl;
			return false;
		}
		return true;
	}
	bool TransAbs::fill_dark_vector(const std::string inname)
	{
		m_dark_filename = inname;
		if( !fill_dark_vector()){
			std::cerr << "Failed request to fill_dark_vector() in bool TransAbs::process_accum(void)" << std::endl;
			m_remove_dark = false;
			return false;
		}
		return true;
	}
	bool TransAbs::fill_dark_vector(void)
	{
		try {
			std::cout << "in TransAbs::fill_dark_vector() method, loading dark_file " << m_dark_filename << std::endl;
			std::ifstream infile(m_dark_filename.c_str(),std::ios::in);
			if (!infile.is_open()){
				std::cerr << "\n\t\t" << m_dark_filename
					<< " didn't open in bool TransAbs::fill_dark_vector(const std::string inname)" << std::endl;
				throw;
			}
			std::vector<std::string> header;
			infile >> header;
			infile >> m_dark_vector;
			infile.close();

		} catch (std::exception & e) {
			std::cerr << "Failed in bool TransAbs::fill_dark_vector(void) with error: " << e.what() << std::endl;
			return false;
		}
		return true;
	}
	bool TransAbs::fft_result_vector_r(void)
	{
		if (m_records_r.size() == 0)
			return false;
		if (m_records_r.size() != m_records_hc.size())
			m_records_hc.resize(m_records_r.size(),(double*)NULL);
		for (unsigned i = 0;i<m_records_r.size();++i){
			if (m_records_r[i] == NULL) m_records_r[i] = (double*)fftw_malloc(sizeof(double) * m_data.front().size());
			if (m_records_hc[i] == NULL) m_records_hc[i] = (double*)fftw_malloc(sizeof(double) * m_data.front().size());
		}
		try {
			fftw_execute_r2r(m_plan_r2hc,
					m_records_r.front(), 
					m_records_hc.front()
					);
				
		} catch (std::exception &e) {
			std::cerr << "Failed in fft for bool TransAbs::fft_result_vector_r() with error: " << e.what() << std::endl;
			return false;
		}
		return true;
	}
	bool TransAbs::fft_accum_result(void)// hopefully this is legacy.. should do weighted average before FFT since not using 2DFFT
	{
		if (m_accum_result_r_fftwPtr == NULL){
			std::cerr << "Failed in bool TransAbs::fft_accum_result(void), m_accum_result_r_fftwPtr = ";
			std::cerr << m_accum_result_r_fftwPtr << std::endl;
			return false;
		}
		if (m_accum_result_hc_fftwPtr == NULL){
			m_accum_result_hc_fftwPtr = (double*)fftw_malloc(sizeof(double)
					* m_data_accum[m_accum_index]->shape()[0] 
					* m_data_accum[m_accum_index]->shape()[1]);
		}
		try {
			for (unsigned r = 0; r<m_data_accum[m_accum_index]->shape()[0] ; ++r){
				fftw_execute_r2r(m_plan_r2hc,
						m_accum_result_r_fftwPtr + r*(m_data_accum[m_accum_index]->shape()[1]) ,
						m_accum_result_hc_fftwPtr + r*(m_data_accum[m_accum_index]->shape()[1]) 
						);
			}
		} catch (std::exception & e) {
			std::cerr << "Failed TransAbs::fill_fft(void) method: " << e.what() << std::endl;
			return false;
		}
		return true;
	}
	bool TransAbs::fill_result_vector_hc(void)
	{
		if (m_accum_result_hc_fftwPtr == NULL){
			std::cerr << "Failed bool TransAbs::fill_result_vector_hc(void) bacause m_accum_result_hc_fftwPtr = ";
			std::cerr << m_accum_result_hc_fftwPtr << std::endl;
			return false;
		}
		boost::multi_array<long long,2> * imagePtr(m_data_accum_ref[m_accum_index]);

		std::vector<double> weights(imagePtr->shape()[0],(double)0); // weights for every row
		long long centroid=0,width=0;
		getCentroid(*imagePtr,centroid,width,weights);
		unsigned sz(imagePtr->shape()[1]);

		std::cout << "Fix this incorrect width calculation in bool TransAbs::fill_result_vector_hc(void)" << std::endl;
		unsigned startind = std::min(centroid-width/2,(long long)0); // HERE HERE this is not right
		unsigned stopind = std::max(centroid+width/2,(long long)weights.size());

		try {

			double sum(0.);
			// 0th element (fully real) case // 
			double val(0.);
			for (unsigned r = startind; r < stopind;++r){
				sum += weights[r];
				val += weights[r] * m_accum_result_hc_fftwPtr[r * sz + 0];
			}
			val /= sum;
			m_records_hc.front()[0] = val;
			// sz/2 element (fully imag?) case // 
			val = 0.;
			for (unsigned r = startind; r < stopind;++r){
				sum += weights[r];
				val += weights[r] * m_accum_result_hc_fftwPtr[r * sz + sz/2];
			}
			val /= sum;
			m_records_hc.front()[sz/2] = val;

			// 1st-end elements (complex) cases //
			for (unsigned c = 1;c < imagePtr->shape()[1]/2;++c){ // weighted average for every column
				sum = 0.;
				std::complex<double> z(0.,0.);
				for (unsigned r = startind; r < stopind;++r){
					sum += weights[r];
					z += std::complex<double> (
							weights[r]*m_accum_result_hc_fftwPtr[r * imagePtr->shape()[1] + c] ,
							weights[r]*m_accum_result_hc_fftwPtr[r * imagePtr->shape()[1] + imagePtr->shape()[1]-c]
							);
				}
				z /= sum;
				m_records_hc.front()[c] = z.real();
				m_records_hc.front()[imagePtr->shape()[1] - c] = z.imag();
			}
		} catch (std::exception & e) {
			std::cerr << "Failed in TransAbs::fill_fft_data(const unsigned r) method: " << e.what() << std::endl;
			return false;
		}
		return true;
	}
	bool TransAbs::fill_result_vector_r(void)
	{

		if (m_accum_result_r_fftwPtr == NULL){
			std::cerr << "Failed bool TransAbs::fill_result_vector_r(void) bacause m_accum_result_r_fftwPtr = ";
			std::cerr << m_accum_result_r_fftwPtr << std::endl;
			return false;
		}
		boost::multi_array<long long,2> * imagePtr(m_data_accum_ref[m_accum_index]);
		unsigned sz(imagePtr->shape()[1]);
		if (m_records_r.size() == 0 || m_records_r.front() == NULL){
			m_records_r.clear();
			m_records_r.push_back((double *) fftw_malloc(sizeof(double) * sz));
			for (unsigned i=0;i<m_nderivs;++i){
				m_records_r.push_back((double *) fftw_malloc(sizeof(double) * sz));
			}
		}
		std::vector<double> weights(imagePtr->shape()[0],(double)0); // weights for every row
		long long centroid=0,width=0;
		getCentroid(*imagePtr,centroid,width,weights);

		std::cout << "Fix this incorrect width calculation in bool TransAbs::fill_result_vector_hc(void)" << std::endl;
		unsigned startind = std::min(centroid-width/2,(long long)0); // HERE HERE this is not right
		unsigned stopind = std::max(centroid+width/2,(long long)weights.size());
		try {
			for (unsigned c = 0;c < imagePtr->shape()[1];++c){ // weighted average for every column
				double sum = 0.;
				double value = 0.;
				for (unsigned r = startind; r < stopind;++r){
					sum += weights[r];
					value += weights[r]*m_accum_result_r_fftwPtr[r * imagePtr->shape()[1] + c] ;
				}
				value /= sum;
				if (m_remove_dark) value -= m_dark_vector[c];
				m_records_r.front()[c] = value;
			}
		} catch (std::exception & e) {
			std::cerr << "Failed in TransAbs::fill_fft_data(const unsigned r) method: " << e.what() << std::endl;
			return false;
		}
		return true;
	}
	// Chebyshev is worse than the Weiner filter... lots of ripples that aren't real //
	bool TransAbs::cheb_approx_result_vector_r(void)
	{
		unsigned sz = m_data_accum_ref[m_accum_index]->shape()[1];
		if (sz ==0){
			std::cerr << "Failing in bool TransAbs::cheb_approx_result_vector_r(void) because m_records_r.front().size() = " << std::endl;
			return false;
		}
		try {
			if (m_records_r.size() == 0 || m_records_r.front() == NULL){
				m_records_r.clear();
				m_records_r.push_back((double *) fftw_malloc(sizeof(double) * sz));
				for (unsigned i=0;i<m_nderivs;++i){
					m_records_r.push_back((double *) fftw_malloc(sizeof(double) * sz));
				}
			}
			// Chebyshev approximation
			gsl_cheb_series * cs = gsl_cheb_alloc(50);
			gsl_function F;
			F.function = cheb_return_data;
			F.params = (void*)m_records_r.front();
			gsl_cheb_init(cs,&F,0.,1.);
			std::cerr << "\n\n gsl_cheb_eval(cs,0.5) = " << gsl_cheb_eval(cs,0.5) << std::endl;
			for (unsigned i=0;i<sz;++i){
				double x = (double) i / (double)(sz-1); 
				if (x<0.) x=0.;
				if (x>1.0) x=1.;
				m_records_r[0][i] = gsl_cheb_eval(cs,x)*sz;
			}
			gsl_cheb_free(cs);
		} catch (std::exception & e) {
			std::cerr << "Failed in bool TransAbs::cheb_approx_result_vector_r(void) wiht error: " << e.what() << std::endl;
			return false;
		}
		return true;
	}
	bool TransAbs::filter_result_vector_hc(void)
	{
		try{
			unsigned sz(m_data_accum_ref[m_accum_index]->shape()[1]);

			std::vector<double> signal(sz/2,0.);
			std::vector<double> noise(sz/2,0.);
			std::vector<double> filter(sz/2,0.);
			fit_filter(signal,noise,filter);

			for (unsigned d=0;d<m_records_hc.size();++d){
				m_records_hc[d][0] = 0.; // forces integral to be zero, removes the DC term
				m_records_hc[d][sz/2] = 0.;// imaginary component of the 0th... maybe, this could be wrong
			}


			for(unsigned i=1;i<sz/2;++i){ // later we can change this to also write the raw
				m_records_hc[0][i] *= filter[i]; 
				m_records_hc[0][sz-i] *= filter[i];
				std::complex<double> z(m_records_hc.front()[i] * filter[i], m_records_hc.front()[sz-i] * filter[i]);
				for (unsigned d = 1;d<m_records_hc.size();++d){
					// I think I don't want derivatives, but 1st only 
					// then get the square-centered for each and again and again.
					// This would make a costruction much like the cos and sin...
					// But only do this if we get better at the fititng for Wiener filtering
					// Also, think about Grahm-Schmidt orthogonalization for computing the higher modes.
					// This process will also help when applied to the simulated revival structures for projections
					// Use std::inner_product(v.begin(),v.end(),w.begin) for this
					z *= std::complex<double>(0.,(double)i);
					z /= sqrt((double)sz);
					z *= filter[i];
					m_records_hc[d][i] = z.real();
					m_records_hc[d][sz-i] = z.imag();
				}
			}
		} catch (std::exception & e){
			std::cerr << "Failed in bool TransAbs::filter_result_vector_hc(void) wiht e.what() = " << e.what() << std::endl;
			return false;
		}

		return true;
	}
	bool TransAbs::fit_filter(std::vector<double> & s, std::vector<double> & n, std::vector<double> & f)
	{
			// from hand fitting spectrum,
			// noise(x) = 10./x + .1*x**.25
			// signal(x) = (x>8?1e4*x**-3:1e4*8.**-3)
			//
			//new fitting
			//f(x)=a*0.5*( 1-erf((x-x0)/w) )+g(x)
			//g(x)=5.*(x**-.5)
		// OK, here is how we do it.  
		// First fill acomplex vector with the values
		// Take the starting value of abs() as the hard top limit
		// fit a line in log/log representation to the portion of the signal which goes from .75:.25 of that max
		// Set the signal vector to the exp(min(top limit,line(x)))
		unsigned sz(m_data_accum_ref[m_accum_index]->shape()[1]);
		s.resize(sz/2,0.);
		n.resize(sz/2,0.);
		f.resize(sz/2,0.);
		unsigned noisewin = s.size()/2;
		try{
			assert(s.size()==n.size());
			std::complex<double> z(0.,0.);
			double noisemean = 0.;
			for (unsigned i= sz/2-noisewin;i<sz/2;++i){
				z = std::complex<double>(m_records_hc.front()[i] , m_records_hc.front()[sz - i]);
				noisemean += (std::norm(z)); 
			}
			noisemean /= (double)noisewin;
			for (unsigned i = 0;i<n.size();++i)
				n[i] = (noisemean);
				//n[i] += 1./(double)(i);
			unsigned startind = 1; // we always need to start the vector from 1
			z = std::complex<double>(m_records_hc.front()[startind] , m_records_hc.front()[sz - startind]);
			double startval = std::log(std::norm(z));
			double val = startval;
			// i0 = mean of i's
			// v0 = mean of v's
			// slope = mean of vdifferences/idifferences
			// i's and v's are all in log coordinates
			// We're going to need to print the filter
			
			// HERE HERE HERE HERE //
			// OK do this differently,
			// Fit a real linear regression, two likely, one for above the 75% line and the other below it//
			// You have to use a linear regressor, these quick trics are killing you. //
			// Remember, noisemean is now in lin space, but we work with fitting in the log space.
			double noiseref = std::log(n[startind]);
			startind = 1;
			while (startind < s.size() && val>(startval*.5 + .5*noiseref)){
				++startind;
				z = std::complex<double>(m_records_hc.front()[startind] , m_records_hc.front()[sz - startind]);
				val = std::log(std::norm(z));
			}
			std::vector<double> ivals,vvals;
			ivals.clear();ivals.reserve(50);
			vvals.clear();vvals.reserve(50);
			
			// while loop for filling fitting vectors //
			unsigned stopind = startind;
			z = std::complex<double>(m_records_hc.front()[stopind] , m_records_hc.front()[sz - stopind]);
			ivals.push_back( log((double)stopind) );
			vvals.push_back( log(std::norm(z)) );
			while (stopind<s.size() && vvals.back()>noiseref){
				++stopind;
				z = std::complex<double>(m_records_hc.front()[stopind] , m_records_hc.front()[sz - stopind]);
				ivals.push_back( log((double)stopind) );
				vvals.push_back( log(std::norm(z)) );
			}
			// HERE do the gsl linear fitting... first test that there are enough points //
			//

			if (ivals.size() < 3){
				for (unsigned i=startind-1;i>std::max((int)startind-3 , 0);--i){
					z = std::complex<double>(m_records_hc.front()[i] , m_records_hc.front()[sz - i]);
					ivals.push_back( log((double)i) );
					vvals.push_back( log(std::norm(z)) );
				}
				for (unsigned i=stopind+1;i<std::min((int)stopind+3, (int)s.size());++i){
					z = std::complex<double>(m_records_hc.front()[i] , m_records_hc.front()[sz - i]);
					ivals.push_back( log((double)i) );
					vvals.push_back( log(std::norm(z)) );
				}
				std::cerr << "\n\n\t\t--- not enough points for fitting, adding more ---\n\n" << std::endl;
			}
			double c0,c1,cov00,cov01,cov11,chisq;
			gsl_fit_linear(ivals.data(),1,vvals.data(),1,ivals.size(),&c0,&c1,&cov00,&cov01,&cov11,&chisq);

			// print out the filter function when all done, at least the signal and the noise functions //
			std::string filename(m_filename);
			filename += ".wienerspecs." + m_stage_position_list[m_accum_index];
			std::ofstream outfile(filename.c_str(),std::ios::out);
			outfile << "#spectra used for the Weiner filter\n#signal\tnoise\tWiener\n";
			for (unsigned i=0;i<s.size();++i){
				z = std::complex<double>(m_records_hc.front()[startind] , m_records_hc.front()[sz - startind]);
				//double y0 = std::log(std::abs(z));
				//double i0 = std::log((double)startind);
				double y =  c1*(log((double)i)) + c0;
				s[i] = std::exp( std::min(startval,y) );
				if (s[i]>0.)
					f[i] = 1./(1.+n[i]/s[i]);
				outfile << s[i] << "\t";
				outfile << n[i] << "\t";
				outfile << f[i] << "\n";
			}
			outfile.close();

		} catch (std::exception & e) {
			std::cerr << "Failed in bool fit_signal(std::vector<double> & x, std::vector<double> & y)" << std::endl;
			return false;
		}
		return true;
	}
	bool TransAbs::backfft_result_vector_hc(void)
	{
		unsigned sz = m_data_accum_ref[m_accum_index]->shape()[1];
		try{
			for (unsigned d=0;d<m_records_hc.size();++d){ // computing m_nderivs
				fftw_execute_r2r(m_plan_hc2r,
						m_records_hc[d],
						m_records_r[d]
						);
				for (unsigned i=0;i<sz;++i) m_records_r[d][i] /= sz; // presupposes no normalization on forward transform
			}

		} catch (std::exception & e ) {
			std::cerr << "Failed in bool TransAbs::backfft_result_vector_hc(void) with error " << e.what() << std::endl;
			return false;
		}
		return true;
	}
	bool TransAbs::fill_orthonormal_result_r(void)
	{
		// HERE HERE HERE HERE //
		// Let's see if this works for grahm schmidt orthogonalization
		// Well, this isn't working as expected... need to feed the first deriv in as a 1st excited state and make sure they are all orthonormal
		unsigned sz = m_data_accum_ref[m_accum_index]->shape()[1];
		// size the 3dmatrix properly by npowers for each deriv
		// Consider using a boost::multi_array() and 3D representation [pow][dervi][vecrange]
		// For now using a 3D vec< vec< vec<> > >
		m_orthonorm_r.resize(m_npowers);
		for (unsigned p=0;p<m_npowers;++p){
			m_orthonorm_r[p].resize(m_nderivs+1);
			for (unsigned d=0;d<m_orthonorm_r[p].size();++d){
				m_orthonorm_r[p][d].resize(sz,0.);
			}
		}
		std::vector<double> u(sz,0.);
		// first treat each derivative as it's own base for GS-orthogonalization
		try{
			assert(m_records_r.size() == m_nderivs+1);
			double scale(0.);
			double proj;
			// have to assign since m_records_r[d] is double* not a std::vector
			m_orthonorm_r[0][0].assign(m_records_r[0],m_records_r[0]+sz);
			sqr_normalize(m_orthonorm_r[0][0]);
			for (unsigned d = 1;d<m_records_r.size();++d){
				u.assign(m_records_r[d],m_records_r[d]+sz);
				m_orthonorm_r[0][d] = u;
				for (unsigned pd=0;pd<d;++pd){ // pd indicates previous derivatives
					proj = projection(m_orthonorm_r[0][pd],u);
					for (unsigned c=0;c<m_orthonorm_r[0][d].size();++c){
						m_orthonorm_r[0][d][c] -= proj * m_orthonorm_r[0][pd][c];
					}
				}
				sqr_normalize(m_orthonorm_r[0][d]);
			}
			// Now move onto the powers //
			for (unsigned p=1;p<m_npowers;++p){
				// do all the previous derivatives for this power and all the previous derivatives for all previous powers
				for (unsigned d=0;d<m_records_r.size();++d){
					for (unsigned c=0;c<u.size();++c){
						u[c] = std::pow(m_orthonorm_r[0][d][c],(int)(p+1));
					}
					m_orthonorm_r[p][d] = u;
					for (unsigned pd=0;pd<d;++pd){ // this power, previous derivatives
						proj = projection(m_orthonorm_r[p][pd],u);
						for (unsigned c=0;c<m_orthonorm_r[0][d].size();++c){
							m_orthonorm_r[p][d][c] -= proj * m_orthonorm_r[p][pd][c];
						}
					}
					for (unsigned pp=0;pp<p;++pp){ // previous powers, all derivatives
						for (unsigned pd=0;pd<m_orthonorm_r[pp].size();++pd){
							proj = projection(m_orthonorm_r[pp][pd],u);
							for (unsigned c=0;c<m_orthonorm_r[p][d].size();++c){
								m_orthonorm_r[p][d][c] -= proj * m_orthonorm_r[pp][pd][c];
							}
						}
					}
					sqr_normalize(m_orthonorm_r[p][d]);
				}
			}
		} catch (std::exception & e) {
			std::cerr << "failed in filling orthonorm_r with error: " << e.what() << std::endl;
			return false;
		}


		return true;
	}
	bool TransAbs::print_out_result_r(void)
	{
		std::string filename(m_filename);
		filename += ".result.orthonorm.r." + m_stage_position_list[m_accum_index];
		std::string gnufilename(m_filename);
		gnufilename += ".result.r." + m_stage_position_list[m_accum_index] + ".gnu";
		try {
			std::ofstream gnuoutfile(gnufilename.c_str(),std::ios::out);
			std::ofstream outfile(filename.c_str(),std::ios::out);
			gnuoutfile << "# result back fft vector out\n#";
			outfile << "# orthonormal result back fft vector out (p,d) on newlines:\t";
			for (unsigned d=0;d<m_records_r.size();++d) gnuoutfile << "d" << d << "\t";
			for (unsigned d=0;d<m_records_r.size();++d) gnuoutfile << "sq_d" << d << "\t";
			for (unsigned p=0;p<m_orthonorm_r.size();++p){
				for (unsigned d=0;d<m_orthonorm_r[p].size();++d){
					outfile << "(" << p << "," << d << ") ";
				}
			}
			gnuoutfile << "\n";
			outfile << "\n";
			unsigned sz(m_data_accum_ref[m_accum_index]->shape()[1]);

			for (unsigned p=0;p<m_orthonorm_r.size();++p){
				for (unsigned d=0;d<m_orthonorm_r[p].size();++d){
					for (unsigned c=0;c<m_orthonorm_r[p][d].size();++c){
						outfile << m_orthonorm_r[p][d][c] << "\t";
					}
					outfile << "\n";
				}
			}

			for (unsigned c=0;c<sz;++c){
				for (unsigned d=0;d<m_records_r.size();++d)
					gnuoutfile << m_records_r[d][c] << "\t";
				gnuoutfile << "\n";
			}
			gnuoutfile.close();
			outfile.close();
		} catch (std::exception & e) {
			std::cerr << "Failed in bool TransAbs::print_out_result_r(void)" << std::endl;
			return false;
		}
		return true;
	}

	bool TransAbs::write_dark_file(void)
	{
		std::ofstream outfile(m_dark_filename.c_str(),std::ios::out);
		std::string gnufilename = m_dark_filename + ".gnuplot";
		std::ofstream gnuplotout(gnufilename.c_str(),std::ios::out);
		try {
			if (!outfile.is_open()){
				std::cerr << "Failed to open " << m_dark_filename << " for bool TransAbs::write_dark_file(void) call" << std::endl;
				throw;
			}
			if (!gnuplotout.is_open()){
				std::cerr << "Failed to open " << gnufilename << " for bool TransAbs::write_dark_file(void) call" << std::endl;
				throw;
			}
			if (m_dark_contributions == 0) throw;
			unsigned sz=m_dark_vector.size();
			for (unsigned i=0;i<sz;++i)
				m_dark_vector[i] /= (double)m_dark_contributions;

			// doing a one nearest neighbor smoothing of only bad data // 
			double win = 5.;
			outfile << "# dark vector file" << std::endl;
			gnuplotout <<" # dark vector file" << std::endl; 
			double value = m_dark_vector[0];
			if (!(std::isnan(m_dark_vector[1]) && std::isnan(m_dark_vector[2]))){
				if (	std::isnan(value) 
						|| (	std::abs(value-m_dark_vector[1])>win && std::abs(value-m_dark_vector[2])>win ) ){
					m_dark_vector[0] = .5*(m_dark_vector[1] + m_dark_vector[2]);
				}// else {
				//	m_dark_vector[0] = 1./3.*(value + m_dark_vector[1] + m_dark_vector[2]);
				//}
			}
			outfile << m_dark_vector[0] << "\t";
			gnuplotout << m_dark_vector[0] << "\n";
			for (unsigned i = 1 ; i<m_dark_vector.size()-1; ++i){
				value = m_dark_vector[i];
				if (!(std::isnan(m_dark_vector[i-1]) && std::isnan(m_dark_vector[i+1]))){
					if (	std::isnan(value) 
							|| (	std::abs(value-m_dark_vector[i-1])>win && std::abs(value-m_dark_vector[i+1])>win ) ){
						m_dark_vector[i] = .5*(m_dark_vector[i-1]+m_dark_vector[i+1]);
					}// else {
					//	m_dark_vector[i] = (m_dark_vector[i-1] + value + m_dark_vector[i+1]);
					//}
				}
				outfile << m_dark_vector[i] << "\t";
				gnuplotout << m_dark_vector[i] << "\n";
			}
			value = m_dark_vector.back();
			if (!(std::isnan(m_dark_vector[sz-2]) && std::isnan(m_dark_vector[sz-3]))){
				if (	std::isnan(value) 
						|| (	std::abs(value-m_dark_vector[sz-2])>win && std::abs(value-m_dark_vector[sz-3])>win ) ){
					m_dark_vector.back() = .5*(m_dark_vector[sz-2]+m_dark_vector[sz-3]);
				}// else {
				//	m_dark_vector.back() = 1./3.*(value + m_dark_vector[sz-2]+m_dark_vector[sz-3]);
				//}
			}
			outfile << m_dark_vector.back() << std::endl;
			gnuplotout << m_dark_vector.back() << std::endl;
			outfile.close();
			gnuplotout.close();
		} catch ( std::exception & e) {
			std::cerr << "Failed in bool TransAbs::write_dark(const std::string outname) with error: " << e.what() << std::endl;
			return false;
		}
		return true;
	}

	bool TransAbs::process_accum(void)
	{
		try {
			if (!fill_accum_result()){
				std::cerr << "\tFailing fill_accum() in bool TransAbs::process_accum(void)" << std::endl;
				return false;
			}
			if (!fft_init()){
				std::cerr << "\tFailing fft_init() in bool TransAbs::process_accum(void)" << std::endl;
				throw;
			}
			/*
			// Commenting out the Chebyshev stuff since it looks pretty bad compared to Weiner
			if (!cheb_approx_result_vector_r()){
				std::cerr << "\tFailed cheb_approx_result_vector_r()" ;	
				throw;
			}
			if (!print_out_result_r()){
				std::cerr << "\tFailed print_out_result_r()" ;	
				throw;
			}
			if (!print_out_result()){
				std::cerr << "Failing print_out_result() in bool TransAbs::process_accum(void)" << std::endl;
				throw;
			}
			*/
			if (m_write_dark && in_dark_stage_list()){
				if (!fill_result_vector_r()){
					std::cerr << "Failed to fill_result_vector() in bool TransAbs::process_accum(void)" << std::endl;
					m_write_dark = false;
					throw;
				}
				if (!add_to_dark_vector()){
					std::cerr << "Stage in stagelist, but couldn't fill_dark_vector()" << std::endl;
					throw;
				}
			}
			/*
			if (! (fft_accum_result() && fill_result_vector_hc()) ){
				std::cerr << "Failing fft_accum_result() && fill_result_vector_hc() in bool TransAbs::process_accum(void)" << std::endl;
				throw;
			}
			*/
			if (!m_write_dark){
				if (!( fill_result_vector_r() ))
				{
					std::cerr << "\tFailed fill_result_vector_r()" ;
					throw;
				} 
				if (! fft_result_vector_r()){
					std::cerr << "Failed to fft_result_vector_r() in bool TransAbs::process_accum(void)" << std::endl;
					throw;
				}
				if (! print_out_result_hc()){
					std::cerr << "Failing print_out_result_hc() in bool TransAbs::process_accum(void)" << std::endl;
					throw;
				}
				if (!filter_result_vector_hc()){
					std::cerr << "Failed to filter before backfft in bool TransAbs::process_accum(void)" << std::endl;
					throw;
				}
				size_t sz = m_data_accum_ref[m_accum_index]->shape()[1];
				double slope = process_phase_hc(m_records_hc.front(),sz);
				if (! backfft_result_vector_hc()){
					std::cerr << "Failing (backfft_result_vector_hc() in bool TransAbs::process_accum(void)";
					throw;
				}
				if (! fill_orthonormal_result_r()){
					std::cerr << "Failed to fill_orthonormal_result_r() in bool TransAbs::process_accum(void)";
					throw;
				}
				if ( ! print_out_result_r()){
					std::cerr << "Failing print_out_result_r() in bool TransAbs::process_accum(void)" << std::endl;
					throw;
				}
			}
		} catch (std::exception & e) {
			std::cerr << "\n\nExeption: " << e.what() << std::endl;
			return false;
		}
		return true;
	}
	bool TransAbs::updateStage(const float in)
	{
		m_stage_position.resize(16, '\0');
		size_t written = std::snprintf(&m_stage_position[0], m_stage_position.size(), "%.3f", in);
		m_stage_position.resize(written);
		std::vector<std::string>::iterator it = find(m_stage_position_list.begin(),m_stage_position_list.end(),m_stage_position);
		if (it == m_stage_position_list.end()){
			m_accum_index = m_stage_position_list.size();
			m_stage_position_list.push_back(m_stage_position);
			std::cout << "m_stage_position_list.push_back( " << m_stage_position << "); in mm" << std::endl;
			return start_new_accumulation();
		}
		m_accum_index = std::distance(m_stage_position_list.begin(),it);
		return true;
	}
	bool TransAbs::updateStage_close(void)
	{
		std::cerr << "Here in bool TransAbs::updateStage_close(void)" << std::endl;
		return close_new_accumulation();
	}
	bool TransAbs::start_new_accumulation(void)
	{
		m_samplecodes->open("samplecodes.out",std::ios::app);
		m_data_accum.push_back(NULL); 
		m_data_accum_ref.push_back(NULL); 
		m_shots_accum.push_back(0);
		m_shots_accum_ref.push_back(0);
		return true;
	}
	bool TransAbs::close_new_accumulation(void)
	{
		std::cerr << "Here in bool TransAbs::close_new_accumulation(void)" << std::endl;
		process_accum();
		m_samplecodes->close();
		return true;
	}
	bool TransAbs::fill_dark_stage_position_list(const std::vector< std::string > in)
	{
		m_dark_stage_position_list.assign(in.begin(),in.end());
		return true;
	}
	bool TransAbs::in_dark_stage_list(void)
	{
		std::vector<std::string>::iterator it = find(m_dark_stage_position_list.begin(),m_dark_stage_position_list.end(),m_stage_position);
		return it!=m_dark_stage_position_list.end();
	}

	template <typename T>
		bool TransAbs::getCentroid(boost::multi_array<T,2> & frame, T & centroid, T & width, std::vector<double> & weights)
		{
			try{
				if (weights.size() != frame.shape()[0]){
					weights.resize(frame.shape()[0],T(0));
				}
				std::fill(weights.begin(),weights.end(),(double)0);
				centroid = T(0);
				width = T(0);
				T baseline(0),sum(0);
				unsigned basewin = frame.shape()[1]/10;
				for (unsigned c=0;c<frame.shape()[1];++c){
					for (unsigned i=0;i<basewin;++i){
						baseline += frame[i][c];
					}
					baseline /= basewin;
					for (unsigned i=0;i<frame.shape()[0];++i){
						T value = frame[i][c] - baseline;
						weights[i] += (double)value;
						sum += value;
						centroid += value * T(i);
					}
				}
				centroid /= sum;
				for (unsigned c=0;c<weights.size();++c){ weights[c] /= (double)sum; }

				sum = T(0);
				for (unsigned c=0;c<frame.shape()[1];++c){
					for (unsigned i=0;i<frame.shape()[0];++i){
						T value = frame[i][c] - baseline;
						if (value > (frame[i][c] - baseline)/T(20)){
							sum += value;
							width += value*std::pow((T(i) - centroid),int(2));
						}
					}
				}
				width = T(sqrt((double)width/sum));
			} catch(std::exception & e){
				std::cerr << "Failed in templated TransAbs::getCentroid() method, error = " << e.what() << std::endl;
				return false;
			}
			std::cout << "templated getCentroid() gives centroid = " << centroid << "\t and width = " << width << std::endl;
			return true;
		}

		double TransAbs::process_phase_hc(double data[],const size_t sz)
		{
			// For the slope calculation, I need to dial the feature back to 0 slope for fitting/unwrapping subtracting slope offset 
			// This moves the 0 slope that makes unwrapping easier into the center of the feature image.
			// then add this offset back to the slope to get back to real coordinates.
			//
			//
			//Likely the only way to handle this is to do an fftshift around the index of the peak in the signal, we can rotate back either with the fft shift index or with the slope
			// Also, we could try again doing a weighted sum of differences
			//
			// Another thing is to check that the offset, si between -pi and pi before adding it back to the data.
			//
			// OK, just adding the subtraction of 2pi to offset helped immensely... 
			// I guess this somehow keeps the slope going negative.... 
			// Whatever works

			std::vector<double> inds(sz/2,0.);
			std::vector<double> amp(sz/2,0.);
			std::vector<double> phase(sz/2,0.);

			inds[0] = 0.;
			amp[0] = std::abs(std::complex<double> (data[0],0.));
			phase[0] = std::arg(std::complex<double> ((double)data[0],0.)); 

			inds[1] = 1.;
			amp[1] = std::abs(std::complex<double> ((double)data[1],(double)data[sz-1]));
			phase[1] = std::arg(std::complex<double> ((double)data[1],(double)data[sz-1]));
			int wrappings = 0;
			double offset = phase[1] - phase[0];
			if (std::abs(offset) > M_PI){
				wrappings = (int)(offset/(2.*M_PI) + 0.5); 
			}
			offset -= wrappings*2.*M_PI;
			//if (offset > 0)
			//	offset -= 2.*M_PI; // This was a magical help. keeps slopes negative.
			phase[1] -= offset;

			inds[sz/2] = (double)sz/2;
			amp[sz/2] = std::abs(std::complex<double> (data[sz/2],0.));
			phase[sz/2] = std::arg(std::complex<double> ((double)data[sz/2],0.));;
			phase[sz/2] -= offset*inds[sz/2]; 
			double delta(0.);
			delta = phase[1] - phase[0];
			wrappings = 0;
			if (std::abs(delta) > M_PI){
				wrappings = (int)(delta/(2.*M_PI) + 0.5); 
			}
			phase[1] -= wrappings*2.*M_PI;
			unsigned stopind = 3;
			for (unsigned i = 2 ; i< sz/2; ++i){
				inds[i] = (double) i;
				amp[i] = std::abs(std::complex<double> ((double)data[i],(double)data[sz-i]));
				if (i > 3 && amp[i] > .25 * amp[1]) stopind = i;
				phase[i] = std::arg(std::complex<double> ((double)data[i],(double)data[sz-i]));
				phase[i] -= offset*inds[i]; // subtracting slope offset (to effectively center feature in image on 0 phase slope
					delta = phase[i] - phase[i-1];
					int wrappings = 0;
					if (std::abs(delta) > M_PI){
						wrappings = (int)(delta/(2.*M_PI) + 0.5); 
					}
					phase[i] -= wrappings*2.*M_PI;
			}

			double c0,c1,cov00,cov01,cov11,chisq;
			//gsl_fit_linear(inds.data(),1,phase.data(),1,size_t(stopind),&c0,&c1,&cov00,&cov01,&cov11,&chisq);
			gsl_fit_mul(inds.data(),1,phase.data(),1,size_t(stopind),&c1,&cov11,&chisq);
			chisq /= (double) stopind;

			// writing files //
			std::string phasename(m_filename);
			phasename += ".phases";
			std::string ampname(m_filename);
			ampname += ".amps";
			std::string slopename(m_filename);
			slopename += ".phaseslopes";
			std::ofstream phasefile;
			std::ofstream ampfile;
			std::ofstream slopefile;

			if (m_accum_index == 0){ // open fresh to clear previous
				phasefile.open(phasename.c_str(),std::ios::out);
				ampfile.open(ampname.c_str(),std::ios::out);
				slopefile.open(slopename.c_str(),std::ios::out);
				phasefile << "#phases\n";
				ampfile << "#amps\n";
				slopefile << "#slopes, 0 = leftmost, -2pi(degenerate with pi) = rightmost, -pi = middle #stage[mm]\tslope\ty0\tchisq\n";
			} else { // open append to continue filling
				phasefile.open(phasename.c_str(),std::ios::app);
				ampfile.open(ampname.c_str(),std::ios::app);
				slopefile.open(slopename.c_str(),std::ios::app);
			}
			for (unsigned i=0;i<sz/2;++i){
				phasefile << (phase[i]+offset*inds[i]) << "\t"; // add back the pi slope
				ampfile << amp[i] << "\t";
			}
			phasefile << "\n";
			ampfile << "\n";
			slopefile << m_stage_position_list[m_accum_index] << "\t" << (c1+offset) << "\t" << chisq << "\n"; // add back the slope offset
			phasefile.close();
			ampfile.close();
			slopefile.close();
			return c1;
		}
	/*
	   bool TransAbs::unwrap_phase(void)
	   {
	// find the centroid in the vertical dimension and move outward //
	double baseline=0.,centroid=0.,sum=0.;
	unsigned basewin = m_data_fft.size()/10;
	for (unsigned i=0;i<basewin;++i){
			baseline += std::abs(m_data_fft[i][0]);
		}
		baseline /= basewin;
		for (unsigned i=0;i<m_data_fft.size();++i){
			double value = std::abs(m_data_fft[i][0]) - baseline;
			sum += value;
			centroid += value * i;
		}
		centroid /= sum;
		m_phase.resize(m_data_fft.size());
		m_amp.resize(m_data_fft.size());
		unsigned unwrapwin = 200; // ultimately pass this via config file.
		unsigned com(centroid);
		for (unsigned i=0;i<m_data_fft.size() ;++i){
			m_phase[i].resize(m_data_fft[i].size(),0.0);
			m_amp[i].resize(m_data_fft[i].size(),0.0);

			m_amp[i][0] = std::abs(m_data_fft[i][0]);
			m_phase[i][0] = std::arg(m_data_fft[i][0]);
			double delta(0.);

			for (unsigned j = 1 ; j<m_amp[i].size();++j){
				m_amp[i][j] = std::abs(m_data_fft[i][j]);
				m_phase[i][j] = std::arg(m_data_fft[i][j]);
				if (i<com+unwrapwin+1 && i > com-unwrapwin-1){
					delta = m_phase[i][j] - m_phase[i][j-1];
					int wrappings = 0;
					if (std::abs(delta) > M_PI){
						wrappings = (int)(delta/(2.*M_PI) + 0.5); // need to check if .99 of 2pi
					}
					m_phase[i][j] -= wrappings*2.*M_PI;
				}
			}
		}
		for (unsigned i=1;i<unwrapwin;++i){
			double delta(0.);
			int wrappings(0);
			for (unsigned j = 0 ; j<m_phase[i].size();++j){
				delta = m_phase[com+i][j] - m_phase[com+i-1][j];
				wrappings = 0;
				if (std::abs(delta) > M_PI){
					wrappings = (int)(delta/(2.*M_PI) + 0.5);
				}
				m_phase[com+i][j] -= wrappings*2.*M_PI;

				delta = m_phase[com-i][j] - m_phase[com-i+1][j];
				wrappings = 0;
				if (std::abs(delta) > M_PI){
					wrappings = (int)(delta/(2.*M_PI) + 0.5);
				}
				m_phase[com-i][j] -= wrappings*2.*M_PI;
			}
		}
		return true;

	}
	*/
	bool TransAbs::print_debug(void)
	{
		std::cerr << "printing debug to stderr\n";
		std::cerr << "m_data.size() = " << m_data.size() << " records\t";
		std::cerr << "m_data.front().size() = " << m_data.front().size() << " record length\n";
		std::cerr << std::flush;
		return true;
	}
	bool TransAbs::print_out(const unsigned eventnum)
	{
		std::string eventfilename(m_filename);
		eventfilename += ".event_" + boost::lexical_cast<std::string>(eventnum);
		if (m_outfile.is_open())
			m_outfile.close();
		m_outfile.open(eventfilename.c_str(),std::ios::out);
		if (!m_outfile.is_open())
			return false;
		for (unsigned i=0;i<m_data.size();++i){ 
			for (unsigned j=0;j<m_data[i].size();++j){
				m_outfile << m_data[i][j] << "\t";
			}
			m_outfile << "\n";
		}
		m_outfile << "\n";
		m_outfile.close();
		return true;

	}
}

