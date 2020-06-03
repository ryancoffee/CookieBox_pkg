#ifndef COOKIEBOX_PKG_COOKIEBOX_MOD_H
#define COOKIEBOX_PKG_COOKIEBOX_MOD_H

//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id$
//
// Description:
//	Class CookieBox_mod.
//
//------------------------------------------------------------------------

//-----------------
// C/C++ Headers --
//-----------------
#include <algorithm>
#include <memory>
#include <map>
#include <fftw/fftw3.h>
#include <time.h>
#include <iterator>
#include <map>
#include <utility>
//#include <boost/cstdint.hpp>
#include <stdint.h>

//----------------------
// Base Class Headers --
//----------------------
#include "psana/Module.h"
#include "ndarray/ndarray.h"
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/multi_array.hpp>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
// A crap ton of hidden methods are buried in ImgAlgos/GlobalMethods.h
#include "ImgAlgos/GlobalMethods.h"
#include "MsgLogger/MsgLogger.h"
// to work with detector data include corresponding 
// header from psddl_psana package
#include "PSEvt/EventId.h"
#include "psddl_psana/bld.ddl.h"
#include "psddl_psana/acqiris.ddl.h"
#include "psddl_psana/camera.ddl.h"

#include <openmpi/mpi.h>

#include "CookieBox_pkg/Acqiris.h"
#include "CookieBox_pkg/Xtcav.h"
#include "CookieBox_pkg/Gdet.h"
#include "CookieBox_pkg/Ebeam.h"
#include "CookieBox_pkg/TimeTool.h"
//#include "CookieBox_pkg/TransAbs.h"
#include "CookieBox_pkg/Learning.h"
#include "CookieBox_pkg/dataops.h"

#include <fftw/fftw3.h>


//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

#include <fstream>
#include <exception>

//		---------------------
// 		-- Class Interface --
//		---------------------

// Epics PVs of interest
// AMO:LAS:DLS:05:MTR.RBV

namespace CookieBox_pkg {

enum useLogicInd {use_aq, use_eb, use_gd, use_xt, use_tt, use_learn, use_ta};
enum printLogicInd {print_aq, print_eb, print_gd, print_xt, print_tt, print_learn, print_ta};
// index order is [0=TimeTool][1=ebeam][2=gasdet][3=chan][4=integwin]
enum DimDefs {ttind,ebind,gdind,chanind,tofind};

	/// @addtogroup CookieBox_pkg

	/**
	 *  @ingroup CookieBox_pkg
	 *
	 *  @brief Example module class for psana
	 *
	 *  @note This software was developed for the LCLS project.  If you use all or 
	 *  part of it, please give an appropriate acknowledgment.
	 *
	 *  @version \$Id$
	 *
	 *  @author Ryan N. Coffee
	 */

/*
class myexeption: public exception
{
	virtual const char* what() const throw()
	{
		return "My exception happened";
	}
};

int main () {
	try
	{
		throw myex;
	}
	catch (exception& e)
	{
		cout << e.what() << '\n';
	}
	return 0;
}
*/

class CookieBox_mod 
	: public Module 
{
	typedef ndarray<const short,2> shortwf_t;

	typedef boost::multi_array<unsigned,2> a2d_u_t;
	typedef boost::multi_array<long,2> a2d_ll_t;
	typedef boost::multi_array<long double,2> a2d_ld_t;
	typedef boost::multi_array<double,2> a2d_d_t;
	typedef boost::multi_array<long,3> a3d_ll_t;
	typedef boost::multi_array<long double,3> a3d_ld_t;
	typedef boost::multi_array<double,3> a3d_d_t;
	typedef boost::multi_array<long,4> a4d_ll_t;
	typedef boost::multi_array<long double,4> a4d_ld_t;
	typedef boost::multi_array<double,4> a4d_d_t;
	typedef boost::multi_array<long,5> a5d_ll_t;
	typedef boost::multi_array<long double,5> a5d_ld_t;
	typedef boost::multi_array<double,5> a5d_d_t;

	typedef boost::multi_array<long double,5>::array_view<1>::type a5d_ld_1dview_t;
	typedef boost::multi_array<double,5>::array_view<1>::type a5d_d_1dview_t;

	// This is for slicing multi_array into views //
	typedef boost::multi_array_types::index_range range;

	typedef double (CookieBox_mod::*E2T_MemFuncPtr)(const double x, const unsigned c);
	typedef double (CookieBox_mod::*CF_MemFuncPtr)(const unsigned k, const unsigned c);
	typedef bool (CookieBox_mod::*READ_CorrFactorsMemFuncPtr)(std::string & filename);
	typedef bool (CookieBox_mod::*READ_E2T_ConversionMemFuncPtr)(std::string & filename);

	enum SpecVar {idx_kinenergy, idx_legendre};
	public:
	// Default constructor
	CookieBox_mod (const std::string& name) ;

	// Destructor
	virtual ~CookieBox_mod () ;

			/// Method which is called once at the beginning of the job
			virtual void beginJob(Event& evt, Env& env);

			/// Method which is called at the beginning of the run
			virtual void beginRun(Event& evt, Env& env);

			/// Method which is called at the beginning of the calibration cycle
			virtual void beginCalibCycle(Event& evt, Env& env);

			/// Method which is called with event data, this is the only required 
			/// method, all other methods are optional
			virtual void event(Event& evt, Env& env);

			/// Method which is called at the end of the calibration cycle
			virtual void endCalibCycle(Event& evt, Env& env);

			/// Method which is called at the end of the run
			virtual void endRun(Event& evt, Env& env);

			/// Method which is called once at the end of the job
			virtual void endJob(Event& evt, Env& env);

		protected:

		private:
			// Pointers to member function //
			E2T_MemFuncPtr e2t;
			CF_MemFuncPtr corr_factor;
			READ_CorrFactorsMemFuncPtr read_corr_factors;
			READ_E2T_ConversionMemFuncPtr read_e2t_conversion;

			fftw_plan plan_r2hc;
			fftw_plan plan_hc2r;

			// member objects //
			Ebeam m_eb;
			TimeTool m_tt;
			Gdet m_gd;
			Xtcav m_xt;
			//TransAbs m_ta;
			std::vector<Acqiris> m_aq;
			Learning m_learn;

			std::vector< std::vector< float > > m_testpca_features;

			enum writeStyle {inteb_gd,inttt_gd,intgd_eb,intgd_tt};

			// private methods //
			bool isApprovedByCounters();
			std::vector<unsigned> m_skipIDs;
		 	bool processEvent(Event& evt, Env& env); // if processing fails, return from event
			bool sampleevery(unsigned in);
			void setWriteStyle(writeStyle style);
			bool printSpectra(void);
			bool printSpectraLegendre(void);
			bool printShots(void);
			bool printIntegs(void);

			bool read_e2t_interp_conversion(std::string & filename);
			bool read_e2t_conversion_old(std::string & filename);
			bool read_corr_factors_old(std::string & filename);
			bool read_corr_factors_interp(std::string & filename);
			double corr_factor_old(const unsigned k, const unsigned c);
			double corr_factor_neon(const unsigned k, const unsigned c);
			double corr_factor_hand(const unsigned k, const unsigned c);
			double corr_factor_hand(const unsigned c);
			void e2t_setlims(const unsigned k, const unsigned c,std::vector<unsigned> & lims, std::vector<double> & weights);
			void e2t_fill_tofs(const unsigned k, std::vector<double> & tofs);
			double e2t_neon(const double x, const unsigned c);
			double e2t_hand(const double x, const unsigned c);
			double e2t_old(const double x, const unsigned c);
			bool computeLegendreCoeffs(const unsigned t, const unsigned e, const unsigned g, const unsigned k);
			bool fillLegendreVecs(void);
			//
			// Here is where we will project the time-sorted specra onto the othonormal basis from a particular stage reading of the TransAbs output.
			bool projectDiffSpectra(const std::string orthofile);
			void weinerfilter(std::vector<double *> & in,const unsigned sz); // used to smooth the projections
			bool fftDiffSpectra(void);
			bool m_makeRotorProjections;
			bool m_printRotorResiduals;


			// for exporting to Python plotter module //
			bool evtput(Event& evt, Env& env, const unsigned chan);


			// private members //
			time_t m_beginruntime,m_endruntime;

			bool m_gnuplotting;

			std::vector<unsigned> m_count_event;
			std::vector<unsigned> m_failed_event;
			unsigned m_skip_events;
			int m_last_event;
			unsigned m_print_every;

			std::string m_str_runnum,m_str_experiment;
			std::string m_datadir;

			// acq is a 10 bit data (stored in a short) 
			a4d_ll_t m_data_4d; // from typedefs, array4d boost::multi_array<long long int ,4>
			a4d_ll_t m_aqSum_4d; // from typedefs, array4d boost::multi_array<long long int,4>
			a5d_ll_t m_aqSum_5d; // from typedefs, array4d boost::multi_array<long long int ,5>
			
			// OK, now going to do a tt,gd,eb,chan,tof 5d histogram
			// Going to update the number of shots in a tt,gd,eb,chan 4d hist
			std::vector < std::vector <double> > m_legendre_vectors;
			a5d_d_t m_legendres_5d; // from typedefs, array5d boost::multi_array<double ,5>
			a4d_d_t m_avg_legendres_4d;
			a4d_d_t m_totsignal_4d;
			a2d_u_t m_avg_legendres_nsums;
			a5d_ll_t m_data_5d; // from typedefs, array5d boost::multi_array<long long int,5>
			a4d_ll_t m_shots_4d; // from typedefs, array4d boost::multi_array<long long int ,4>
			a4d_ll_t m_avgSpectra;
			a3d_ll_t m_avgSpectra_shots;

			// utility limits for defining kinetic energy spectral ranges (maybe also legendre lengths)
			std::vector<bool> m_sp_projectionMask;
			std::vector<unsigned> m_sp_binsvec;
			std::vector<double> m_sp_startsvec;
			std::vector<double> m_sp_stepsvec;
			// OK, this is for adjusting the transmission functions
			// organization is t0,c-1,c0,c1,c2,c3,... for taylor expansion including c-1 for 1/x term.
			// New organization. using taylors, so e0,c0,c1,c2,c3 f(e) = c0+c1*(e-e0)^1+c2*(e-e0)^2 ...
			// std::vector< std::vector < float > > m_transcoeffs;
			std::vector< std::vector< double > > m_corr_factors_data;
			std::vector< std::vector< double > > m_e2t_data;
			bool m_use_e2t;
			std::vector< std::map< double, double > > m_corr_factors_maps;
			//bool m_adjustTransmission;
			bool m_legendreVecs_filled;
			bool m_printLegendreFFTs;
			double m_gaussroll_center,m_gaussroll_width;
			bool m_gaussroll;
			bool m_rollvibration;
			size_t m_nrolls;

			bool m_printMarkusLegendresCompare;
			std::ofstream m_markusfileChans;
			std::ofstream m_markusfileLegs;

			bool print_patch_results(void);
			inline bool inpatch(const unsigned tin,const unsigned ein, const unsigned gin, const unsigned kin){
				return ( inwin(m_twin,tin) && inwin(m_ewin,ein) && inwin(m_gwin,gin) && inwin(m_kwin,kin) );
			}
			bool m_patch_compute;
			std::vector<double> m_patch_ys;
			std::vector<double> m_patch_legcoeffs;
			unsigned m_patch_bins;
			std::vector<unsigned> m_twin;
			std::vector<unsigned> m_ewin;
			std::vector<unsigned> m_gwin;
			std::vector<unsigned> m_kwin;

			// I think these are legacy //
			a2d_ld_t m_gdetSum_2d; // from typedefs, array4d boost::multi_array<long double,2>
			a2d_ll_t m_gdetShots_2d; // from typedefs, array4d boost::multi_array<long long int,2>
			a3d_ll_t m_gdetShots_3d; // from typedefs, array4d boost::multi_array<long long int ,4>

			// 	MPI related	//
			unsigned m_root_rank,m_rank;
			unsigned m_mpi_size;
	};

} // namespace CookieBox_pkg

#endif // COOKIEBOX_PKG_COOKIEBOX_MOD_H
