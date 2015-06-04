#ifndef TIMETOOL_H
#define TIMETOOL_H

#include "PSEvt/EventId.h"
#include "psana/Module.h"
#include <boost/shared_ptr.hpp>
#include "ImgAlgos/GlobalMethods.h"
#include "MsgLogger/MsgLogger.h"

#include "psddl_psana/camera.ddl.h"


#include <vector>
#include <algorithm>
#include <map>
#include <list>

#include <fstream>

#include "CookieBox_pkg/dataops.h"

namespace CookieBox_pkg 
{
class TimeTool 
{
	public:
		//static std::vector<std::string> m_label_strings;
		enum Var {pos, amp, width};
		enum Lim {min, max, current};

		TimeTool(std::vector<unsigned>& bins_in, std::vector<double>& mins_in, std::vector<double>& maxs_in);
		TimeTool(void);
		~TimeTool(void);
		TimeTool(const TimeTool & b);
		TimeTool & operator=( TimeTool rhs );

	private:
		void deepcopy_data(const TimeTool & b);
		void swap(TimeTool & a, TimeTool & b);

	public:
		inline bool use(bool in){m_use = in; return m_use;}
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
		void srcStr(Source srcStr_in);
		bool init(std::vector<unsigned>& bins_in
				, std::vector<double>& mins_in
				, std::vector<double>& maxs_in );

		bool fill(Event& evt,Env& env);
		bool addfeatures(std::vector<float> & out);
		std::string get_features_labels(void);

		bool print_out(const unsigned eventnum);
		void print_header(void);
		void open_file( std::string & fnamein );
		void close_file(void);

		// These should somehow be templated with interface class (as virtual functions)
		// look to the old LEAP stuff
		inline unsigned bins(Var variable){ return m_bins.at(variable); }
		inline double value(Var variable, Lim limit){ return m_data.at(limit).at(variable); }
		inline unsigned index(Var variable){
			int i = (int)( m_bins.at(variable) * (m_data.at(current).at(variable) - m_data.at(min).at(variable))
					/ (m_data.at(max).at(variable) - m_data.at(min).at(variable)) );
			if (i < 0)
				return (unsigned)0;
			if (i >= (int)m_bins.at(variable) )
				return (unsigned)(m_bins.at(variable) -1);
			return (unsigned) i;
		}

		inline double pos_dbl(void){ return m_data.at(current).at(pos); }
		inline double width_dbl(void){ return m_data.at(current).at(width); }
		inline double amp_dbl(void){ return m_data.at(current).at(amp); }

		inline unsigned width_ind(void){
			unsigned i = (short)( (m_data.at(current).at(width) - m_data.at(min).at(width))
					/ (m_data.at(max).at(width) - m_data.at(max).at(width)) );
			return (i<m_bins.at(width) ? i : 0);
		}
		inline unsigned amp_ind(void){
			unsigned i = (unsigned)( (m_data.at(current).at(amp) - m_data.at(min).at(amp)) 
					/ (m_data.at(max).at(amp) - m_data.at(min).at(amp)) );
			return (i<m_bins.at(amp) ? i : 0);
		}
		inline unsigned pos_ind(void){
			unsigned i = (unsigned)( (m_data.at(current).at(pos) - m_data.at(min).at(pos))
					/ (m_data.at(max).at(pos) - m_data.at(min).at(pos)) );
			return (i<m_bins.at(pos) ? i : 0);
		}

	private:
		bool m_use,m_print;
		std::ofstream m_outfile;
		std::string m_filename;

		std::vector< std::vector<double> > m_data; 
		std::vector< std::vector<unsigned> > m_cam_data; 
		std::vector<unsigned> m_bins; // the 0th element is pos, amp, width in that order

		Source m_srcStr;
		boost::shared_ptr<Psana::Camera::FrameV1> m_srcPtr;

		// 	MPI related	//
		unsigned m_root_rank,m_rank;
		unsigned m_mpi_size;
};
}


#endif
