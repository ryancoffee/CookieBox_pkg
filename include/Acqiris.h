#ifndef ACQIRIS_H
#define ACQIRIS_H

#include "PSEvt/EventId.h"
#include "psddl_psana/bld.ddl.h"
//#include "psddl_psana/bld.h" // this is wrong and fails...
#include "psana/Module.h"
#include <boost/shared_ptr.hpp>
#include <boost/multi_array.hpp> // for the slice version of fill() method
#include "ImgAlgos/GlobalMethods.h"
#include "MsgLogger/MsgLogger.h"

#include <vector>
#include <algorithm>
//#include <map> // map didn't help my badd pass of srcPtr
#include <list>
#include <fstream>	
#include "CookieBox_pkg/dataops.h"

namespace CookieBox_pkg {

class Acqiris 
{
	typedef ndarray<const short,2> shortwf_t;
	typedef boost::multi_array_types::index_range ind_range;
	typedef boost::multi_array<long long,4>::array_view<1>::type a4d_ll_1dview_t;
	typedef boost::multi_array<long long,4>::array_view<2>::type a4d_ll_2dview_t;
	typedef boost::multi_array<long long,4>::array_view<2>::type a5d_ll_2dview_t;

	public:
		enum LimsInd {start,stop,bins};
		Acqiris (void);
		Acqiris (std::vector<unsigned>& lims_in , std::vector<unsigned>& baselims_in);
		~Acqiris(void);
		Acqiris( const Acqiris & rhs);
		Acqiris & operator=( Acqiris rhs );
	
	private:	
		void swap(Acqiris & a, Acqiris & b);
		void deepcopy_data(const Acqiris & b);


	public:
		bool init( std::vector<unsigned>& lims_in, std::vector<unsigned>& baselims_in);

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
		inline bool invert(void){return m_invert;}
		inline bool invert(bool in){m_invert = in;return m_invert;}
		bool print_out(const unsigned eventnum);
		bool print_out(void);
		void print_header(void);
		void open_file( std::string & filename );
		void close_file( void );
		std::string & printname(void);
		void evtput(Event& evt, Env& env);
		void evtput(Event& evt, Env& env, const unsigned chan);

		bool fill(Event& evt, Env& env, a5d_ll_2dview_t & slice, a4d_ll_1dview_t & shotslice);
		bool fill(Event& evt, Env& env, a4d_ll_2dview_t & slice);
		bool fill(Event& evt, Env& env);
		bool addfeatures(std::vector<float> & out);

		void srcStr(Source srcStr_in);

		inline unsigned nchannels(void){return m_nchannels;}
		inline unsigned nsamples(void){return m_lims.at(bins);}
		std::string getConfig(Env& env);

	protected:

	private:
		bool m_use,m_print,m_invert;
		std::ofstream m_outfile;
		std::string m_filename;

		bool m_getConfig;
		boost::shared_ptr<Psana::Acqiris::ConfigV1> m_ConfigPtr;
		boost::shared_ptr<Psana::Acqiris::DataDescV1> m_srcPtr; 
		Source m_srcStr;
		Pds::Src m_src;

		std::vector< std::vector<int> > m_data;

		unsigned m_nchannels;
		unsigned m_max_samples;
		std::vector<unsigned> m_lims;
		std::vector<unsigned> m_baselims;


		// 	MPI related	//
		unsigned m_root_rank,m_rank;
		unsigned m_mpi_size;
};

} // namespace CookieBox_pkg

#endif 

