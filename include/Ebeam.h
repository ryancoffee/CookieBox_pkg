#ifndef EBEAM_H
#define EBEAM_H

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

class Ebeam
{
	public:
		//static std::vector<std::string> m_label_strings;
		enum Var {energy,charge,ltux,ltuy,ltuDx,ltuDy,pkcurr_bc2,energy_bc2,pkcurr_bc1,energy_bc1};
		enum Lim {max,min,current};

		Ebeam( std::vector<unsigned>& bins_in, std::vector<double>& mins_in, std::vector<double>& maxs_in);
		Ebeam(void);
		~Ebeam(void);
		Ebeam(const Ebeam & b);
		Ebeam & operator=(Ebeam rhs);

	private:
		void deepcopy_data(const Ebeam & b);
		void swap(Ebeam & a, Ebeam & b);

	public:
		inline bool use(void){return m_use;}
		inline bool use(bool in){m_use = in; return m_use;}
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
		bool init( std::vector<unsigned>& bins_in
				, std::vector<double>& mins_in
				, std::vector<double>& maxs_in );

		void srcStr(Source srcStr_in);

		bool fill(Event& evt, Env& env);
		bool addfeatures(std::vector<float> & out);
		std::string get_features_labels(void);

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

		inline double charge_dbl(void){ return m_data.at(current).at(charge); }
		inline double energy_dbl(void){ return m_data.at(current).at(energy); }
		inline double ltux_dbl(void){ return m_data.at(current).at(ltux); }
		inline double ltuy_dbl(void){ return m_data.at(current).at(ltuy); }
		inline double ltuDx_dbl(void){ return m_data.at(current).at(ltuDx); }
		inline double ltuDy_dbl(void){ return m_data.at(current).at(ltuDy); }
		inline double pkcurr_bc2_dbl(void){ return m_data.at(current).at(pkcurr_bc2); }
		inline double energy_bc2_dbl(void){ return m_data.at(current).at(energy_bc2); }
		inline double pkcurr_bc1_dbl(void){ return m_data.at(current).at(pkcurr_bc1); }
		inline double energy_bc1_dbl(void){ return m_data.at(current).at(energy_bc1); }


		inline unsigned charge_ind(void){
			unsigned i = (unsigned)( (m_data.at(current).at(charge) - m_data.at(min).at(charge))
					/ (m_data.at(max).at(charge) - m_data.at(min).at(charge)) );
			return (i<m_bins.at(charge) ? i : 0);
		}
		inline unsigned energy_ind(void){
			unsigned i = (unsigned)( (m_data.at(current).at(energy) - m_data.at(min).at(energy))
					/ (m_data.at(max).at(energy) - m_data.at(min).at(energy)) );
			return (i<m_bins.at(energy) ? i : 0);
		}
		inline unsigned ltux_ind(void){
			unsigned i = (unsigned)( (m_data.at(current).at(ltux) - m_data.at(min).at(ltux))
					/ (m_data.at(max).at(ltux) - m_data.at(min).at(ltux)) );
			return (i<m_bins.at(ltux) ? i : 0);
		}
		inline unsigned ltuy_ind(void){
			unsigned i = (unsigned)( (m_data.at(current).at(ltuy) - m_data.at(min).at(ltuy))
					/ (m_data.at(max).at(ltuy) - m_data.at(min).at(ltuy)) );
			return (i<m_bins.at(ltuy) ? i : 0);
		}
		inline unsigned ltuDx_ind(void){
			unsigned i = (unsigned)( (m_data.at(current).at(ltuDx) - m_data.at(min).at(ltuDx))
					/ (m_data.at(max).at(ltuDx) - m_data.at(min).at(ltuDx)) );
			return (i<m_bins.at(ltuDx) ? i : 0);
		}
		inline unsigned ltuDy_ind(void){
			unsigned i = (unsigned)( (m_data.at(current).at(ltuDy) - m_data.at(min).at(ltuDy))
					/ (m_data.at(max).at(ltuDy) - m_data.at(min).at(ltuDy)) );
			return (i<m_bins.at(ltuDy) ? i : 0);
		}
		inline unsigned pkcurr_bc2_ind(void){
			unsigned i = (unsigned)( (m_data.at(current).at(pkcurr_bc2) - m_data.at(min).at(pkcurr_bc2))
					/ (m_data.at(max).at(pkcurr_bc2) - m_data.at(min).at(pkcurr_bc2)) );
			return (i<m_bins.at(pkcurr_bc2) ? i : 0);
		}
		inline unsigned energy_bc2_ind(void){
			unsigned i = (unsigned)( (m_data.at(current).at(energy_bc2) - m_data.at(min).at(energy_bc2))
					/ (m_data.at(max).at(energy_bc2) - m_data.at(min).at(energy_bc2)) );
			return (i<m_bins.at(energy_bc2) ? i : 0);
		}
		inline unsigned pkcurr_bc1_ind(void){
			unsigned i = (unsigned)( (m_data.at(current).at(pkcurr_bc1) - m_data.at(min).at(pkcurr_bc1))
					/ (m_data.at(max).at(pkcurr_bc1) - m_data.at(min).at(pkcurr_bc1)) );
			return (i<m_bins.at(pkcurr_bc1) ? i : 0);
		}
		inline unsigned energy_bc1_ind(void){
			unsigned i = (unsigned)( (m_data.at(current).at(energy_bc1) - m_data.at(min).at(energy_bc1))
					/ (m_data.at(max).at(energy_bc1) - m_data.at(min).at(energy_bc1)) );
			return (i<m_bins.at(energy_bc1) ? i : 0);
		}
	private:
		bool m_use,m_print;
		std::ofstream m_outfile;
		std::string m_filename;

		bool testVersion(Event& evt, Env& env);

		std::vector< std::vector<double> > m_data;
		std::vector<unsigned> m_bins;

		boost::shared_ptr<Psana::Bld::BldDataEBeamV6> m_srcPtrV6;
		boost::shared_ptr<Psana::Bld::BldDataEBeamV7> m_srcPtrV7;
		Source m_srcStr;
		unsigned m_version;
		
		// 	MPI related	//
		unsigned m_root_rank,m_rank;
		unsigned m_mpi_size;
};
}


#endif
