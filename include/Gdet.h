#ifndef GDET_H
#define GDET_H

#include "PSEvt/EventId.h"
#include "psddl_psana/bld.ddl.h"
//#include "psddl_psana/bld.h" // this is wrong and fails...
#include "psana/Module.h"
#include <boost/shared_ptr.hpp>
#include "ImgAlgos/GlobalMethods.h"
#include "MsgLogger/MsgLogger.h"

#include <vector>
#include <algorithm>
//#include <map> // map didn't help my badd pass of srcPtr
#include <list>
#include <fstream>	
#include "CookieBox_pkg/dataops.h"


namespace CookieBox_pkg 
{

class Gdet 
{
	public:
		enum Var {gdet_11,gdet_12,gdet_21,gdet_22,average};
		enum Lim {max,min,current,mean};

		Gdet(std::vector<unsigned>& bins_in, std::vector<double>& mins_in, std::vector<double>& maxs_in);
		Gdet(void);
		~Gdet(void);
		Gdet(const Gdet & b);
		Gdet & operator=( Gdet rhs );

	private:
		void deepcopy_data(const Gdet & b);
		void deepcopy_accumulation(const Gdet & b);
		void swap(Gdet & a, Gdet & b);

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
		inline bool print_parallel(void){return m_print_parallel;}
		inline bool print_parallel(bool in){m_print_parallel = in;return m_print_parallel;}
		bool init( std::vector<unsigned>& bins_in
				, std::vector<double>& mins_in
				, std::vector<double>& maxs_in );

		void srcStr(Source srcStr_in);

		std::string get_features_labels(void);
		bool addfeatures(std::vector<float> & out);
		bool fill(Event& evt, Env& env);
		void stats(void);

		bool print_out(const unsigned eventnum);
		void print_header(void);
		void open_file( std::string & fnamein );
		void close_file(void);

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


		/*
		   somehow add to stats to be used in parallel coordinate representation maybe...
		   for (unsigned i=0;i<vals.size();++i){
		   m_gdets_stats[m_meanInd][i] += vals[i];
		   if (m_gdets_stats[m_minInd][i] > m_shot_gdets[i]) { m_gdets_stats[m_minInd][i] = m_shot_gdets[i];}
		   if (m_gdets_stats[m_maxInd][i] < m_shot_gdets[i]) { m_gdets_stats[m_maxInd][i] = m_shot_gdets[i];}
		   }
		   */



	private:
		bool m_use,m_print,m_print_parallel;
		std::ofstream m_outfile;
		std::string m_filename;

		std::vector< std::vector<double> > m_data;
		std::vector<unsigned> m_bins;

		std::vector< std::vector <double> > m_accumulation;

		double globalmin(void);
		double globalmax(void);

		boost::shared_ptr<Psana::Bld::BldDataFEEGasDetEnergyV1> m_srcPtr;
		Source m_srcStr;

		inline double shotAvg(void){
			return ( (m_data.at(current).at(gdet_11) + m_data.at(current).at(gdet_12)) / 2.) ;
		}

		// 	MPI related	//
		unsigned m_root_rank,m_rank;
		unsigned m_mpi_size;
};
}


#endif
