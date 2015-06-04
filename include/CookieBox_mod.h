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
#include <map>
#include <fftw/fftw3.h>
#include <time.h>

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
#include "CookieBox_pkg/Learning.h"
#include "CookieBox_pkg/dataops.h"


//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

#include <fstream>

//		---------------------
// 		-- Class Interface --
//		---------------------

// Epics PVs of interest
// AMO:LAS:DLS:05:MTR.RBV

namespace CookieBox_pkg {

enum useLogicInd {use_aq, use_eb, use_gd, use_xt, use_tt, use_learn};
enum printLogicInd {print_aq, print_eb, print_gd, print_xt, print_tt, print_learn};

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


class CookieBox_mod 
	: public Module 
{
	typedef ndarray<const short,2> shortwf_t;

	typedef boost::multi_array<long long,2> a2d_ll_t;
	typedef boost::multi_array<long double,2> a2d_ld_t;
	typedef boost::multi_array<double,2> a2d_d_t;
	typedef boost::multi_array<long long,3> a3d_ll_t;
	typedef boost::multi_array<long double,3> a3d_ld_t;
	typedef boost::multi_array<long long,4> a4d_ll_t;
	typedef boost::multi_array<long double,4> a4d_ld_t;
	typedef boost::multi_array<long long,5> a5d_ll_t;
	typedef boost::multi_array<long double,5> a5d_ld_t;

	// THis is for slicing multi_array into views //
	typedef boost::multi_array_types::index_range range;

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
			Ebeam m_eb;
			TimeTool m_tt;
			Gdet m_gd;
			Xtcav m_xt;
			std::vector<Acqiris> m_aq;
			Learning m_learn;

			data_f m_testpca_features;

			enum writeStyle {inteb_gd,inttt_gd,intgd_eb,intgd_tt};

			// private methods //
			bool isApprovedByCounters();
		 	bool processEvent(Event& evt, Env& env); // if processing fails, return from event
			bool sampleevery(unsigned in);
			void setWriteStyle(writeStyle style);
			// for exporting to Python plotter module //
			void evtput(Event& evt, Env& env, const unsigned chan);


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
			a4d_ll_t m_data_4d; // from typedefs, array4d boost::multi_array<long double,4>
			a4d_ll_t m_aqSum_4d; // from typedefs, array4d boost::multi_array<long double,4>
			a5d_ll_t m_aqSum_5d; // from typedefs, array4d boost::multi_array<long double,4>
			a2d_ld_t m_gdetSum_2d; // from typedefs, array4d boost::multi_array<long double,4>
			a2d_ll_t m_gdetShots_2d; // from typedefs, array4d boost::multi_array<long double,4>
			a3d_ll_t m_gdetShots_3d; // from typedefs, array4d boost::multi_array<long double,4>

			// 	MPI related	//
			unsigned m_root_rank,m_rank;
			unsigned m_mpi_size;
	};

} // namespace CookieBox_pkg

#endif // COOKIEBOX_PKG_COOKIEBOX_MOD_H
