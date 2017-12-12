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
#include <exception>

#include <fstream>

#include "CookieBox_pkg/dataops.h"

namespace CookieBox_pkg 
{
class TimeTool 
{
	public:
		//static std::vector<std::string> m_label_strings;
		enum Var {pos, amp, width, delay};
		enum Lim {min, max, current};

		TimeTool(std::vector<unsigned>& bins_in, std::vector<double>& mins_in, std::vector<double>& maxs_in);
		TimeTool(void);
		~TimeTool(void);
		TimeTool(const TimeTool & b);
		TimeTool & operator=( TimeTool rhs );

	private:
		void deepcopy_data(const TimeTool & b);
		void swap(TimeTool & a, TimeTool & b);
		double pos2delay(void);
		bool testvalid(void);
		bool testvalid_surf(void);
		bool inslicewin(void);
		bool inslicewin(const double gdin);
		bool isref(Event& evt,Env& env);

	public:
		inline bool use(bool in){m_use = in; return m_use;}
		inline bool use(void){return m_use;}
		inline bool filled(void){return m_filled;}				
		inline bool use_filter(bool in){m_use_filter = in; return m_use_filter;}
		inline bool use_filter(void){return m_use_filter;}
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

		bool fill(Event& evt,Env& env,const double gdin);
		bool fill(Event& evt,Env& env);
		bool fill_image(Event& evt,Env& env);
		bool addfeatures(std::vector<float> & out);
		std::string get_features_labels(void);

		bool setTimeCalib(std::vector<double> & calib_in);
		void setTimeLims(void);
		std::string getTimeLims(void);
		unsigned getTimeBins(std::string & sout);
		double getTspan(std::string & sout);
		double getTspan(void);

		bool print_out(const unsigned eventnum);
		void print_header(void);
		void open_file( std::string & fnamein );
		void close_file(void);

		void setslicewin(std::vector<double> posin,
				std::vector<double> widthin,
				std::vector<double> amplin,
				std::vector<double> widthPamplin
				);

		// These should somehow be templated with interface class (as virtual functions)
		// look to the old LEAP stuff
		inline unsigned bins(Var variable){ return m_bins.at(variable); }
		inline double value(Var variable, Lim limit){ return m_data.at(limit).at(variable); }
		inline unsigned index(Var variable){
			int i;
			i = (int)( m_bins.at(variable) * (m_data.at(current).at(variable) - m_data.at(min).at(variable))
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
		inline double delay_dbl(void){ return m_data.at(current).at(delay); }

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
		inline unsigned delay_ind(void){
			unsigned i = (unsigned)( (m_data.at(current).at(delay) - m_data.at(max).at(delay))
					/ (m_data.at(min).at(delay) - m_data.at(max).at(delay)) );
			return (i<m_bins.at(delay) ? i : 0);
		}

	private:
		bool m_use,m_print,m_filled;
		bool m_use_filter;
		bool m_slicewin;
		bool m_TimeLimsSet;
		std::vector<double> m_amplwin;// = 0.1 0.8
		std::vector<double> m_fwhmwin; // = 30 90
		std::vector<double> m_poswin; // = 100 1024
		std::vector<double> m_fwhmPamplwin; // = 50 500


		std::ofstream m_outfile;
		std::string m_filename;
//		std::ofstream * m_samplecodes;
//		std::string m_samplecodesname;

		std::vector<double> m_calib;
		std::vector< std::vector<double> > m_data; 


		std::vector<unsigned> m_bins; // the 0th element is pos, amp, width in that order

		std::vector< std::vector<unsigned> > m_cam_data; 
		std::vector<double> m_signal;
		std::vector<double> m_ref;

		Source m_srcStr;
		Source m_evr_src;
		unsigned m_refcode;

		boost::shared_ptr<Psana::Camera::FrameV1> m_srcPtr;

		// 	MPI related	//
		unsigned m_root_rank,m_rank;
		unsigned m_mpi_size;
};
}


#endif
