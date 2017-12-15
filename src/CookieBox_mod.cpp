//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id$
//
// Description:
//	Class CookieBox_mod...
//
// Author List:
//      Ryan N. Coffee
//
//------------------------------------------------------------------------

//-----------------------
// This Class's Header --
//-----------------------
#include "CookieBox_pkg/CookieBox_mod.h"

//-----------------
// C/C++ Headers --
//-----------------
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <exception>
#include <complex>
#include <fftw/fftw3.h>
#include <cmath>
#include <functional>
#include <utility>
#include <map>

#include <boost/math/special_functions/legendre.hpp>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

// This declares this class as psana module
using namespace CookieBox_pkg;
PSANA_MODULE_FACTORY(CookieBox_mod);

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

namespace CookieBox_pkg {

	//----------------
	// Constructors --
	//----------------
	CookieBox_mod::CookieBox_mod(const std::string& name)
		: Module(name)
		  , m_skip_events()
		  , m_print_every(1)
		  , m_last_event()
		  , m_root_rank(0)
		  , m_rank(0)
		  , m_beginruntime(0)
		  , m_endruntime(0)
		  , m_makeRotorProjections(false)
	{
		// get the values from configuration or use defaults
		m_printMarkusLegendresCompare = config("printMarkusLegendresCompare",false);
		m_skipIDs = configList("skipIDs");
		m_print_every = config("print_every",10);
		if (m_print_every < 1) {m_print_every = 1;}
		std::vector<bool> use_logic = configList("use_logic");
		std::vector<bool> print_logic = configList("print_logic");
		std::vector<bool>::iterator logic_itr = use_logic.begin();
		std::vector<bool>::iterator print_itr = print_logic.begin();


		std::cout << "setting e2t pointer to e2t_neon" << std::endl;
		e2t = &CookieBox_mod::e2t_neon;
		corr_factor = &CookieBox_mod::corr_factor_neon;
		read_corr_factors = &CookieBox_mod::read_corr_factors_interp;
		read_e2t_conversion = &CookieBox_mod::read_e2t_conversion_old;


		m_datadir = configStr("datadir","tempdata/");

		logic_itr = use_logic.begin() + use_xt;
		print_itr = print_logic.begin() + print_xt;
		if (m_xt.use(*logic_itr)){
			m_xt.init( (configList("xt_downsample_steps")) , config("xt_threshold",0));
			m_xt.srcStr(configSrc("xt_source","DetInfo(:Opal1k.0)"));
			m_xt.valhist_init( (configList("xt_valhist_list") ) );
			m_xt.set_crop( (configList("xt_crop")) );
			m_xt.print(*print_itr);
			m_xt.addimages(config("add_images_as_features", false));
			m_xt.print_avg(config("avg_xt_images", false));
			if ( m_xt.remove_dark( config("xt_remove_dark",false)) )
			{
				m_xt.dark_file(configStr("xt_dark_file")) ;
				std::cerr << "Failed to load the dark image file in CookieBox_mod::CookieBox_mod() constructor" << std::endl;
			}
		}

		// HERE this is starting the TrasAbs object //
		logic_itr = use_logic.begin() + use_ta;
		print_itr = print_logic.begin() + print_ta;
		if (m_ta.use(*logic_itr)){
			m_ta.init();
			m_ta.srcStr(configSrc("ta_source","AmoEndstation-1|Opal1000-0"));
			m_ta.refcode(config("ta_refcode",46));
			m_ta.srcStrEvr(configSrc("ta_evr_src","NoDetector.0:Evr.0"));
			m_ta.print(*print_itr);
			if ( m_ta.remove_dark( config("ta_remove_dark",false)) || m_ta.write_dark(config("ta_write_dark",true)) ){
				m_ta.dark_file(configStr("ta_dark_file","amoi0314_data/transabs/darkfile.r0214")) ;
				if ( m_ta.write_dark() ){
					std::vector<std::string> darkstagelist = configList("ta_dark_stage_position_list");
					unsigned ndarks = m_ta.fill_dark_stage_position_list(darkstagelist) ;
					if (ndarks<1){
						std::cerr << "failed to pass config params to m_ta.fill_dark_stage_position_list(darkstagelist)" << std::endl;
					}
				}
				if (m_ta.remove_dark()){
					m_ta.npowers(config("ta_npowers",4));
					m_ta.nderivs(config("ta_nderivs",1));
					bool check = m_ta.fill_dark_vector();
					std::cout << "\t\tm_ta.fill_dark_vector() = " << check << std::endl;
				}
			}

		}


		logic_itr = use_logic.begin() + use_gd;
		print_itr = print_logic.begin() + print_gd;
		if (m_gd.use(*logic_itr)){
			std::vector<unsigned> gd_bins_in = configList("gd_bins_list");
			std::vector<double> gd_mins_in = configList("gd_mins_list");
			std::vector<double> gd_maxs_in = configList("gd_maxs_list");
			m_gd.init(gd_bins_in, gd_mins_in, gd_maxs_in);
			m_gd.srcStr(configSrc("gd_source","BldInfo(FEEGasDetEnergy)"));
			m_gd.print(*print_itr);
			m_gwin.resize(2);
			m_gwin = configList("patch_gd_win");

		}

		logic_itr = use_logic.begin() + use_eb;
		print_itr = print_logic.begin() + print_eb;
		if (m_eb.use(*logic_itr)){
			std::vector<unsigned> eb_bins_in = configList("eb_bins_list");
			std::vector<double> eb_mins_in = configList("eb_mins_list");
			std::vector<double> eb_maxs_in = configList("eb_maxs_list");
			m_eb.init(eb_bins_in, eb_mins_in, eb_maxs_in);
			m_eb.srcStr(configSrc("eb_source", "BldInfo(EBeam)"));
			m_eb.print(*print_itr);
			m_ewin.resize(2);
			m_ewin = configList("patch_eb_win");
		}

		logic_itr = use_logic.begin() + use_tt;
		print_itr = print_logic.begin() + print_tt;
		if (m_tt.use(*logic_itr)){
			std::vector<unsigned> tt_bins_in = (configList("tt_bins_list"));
			std::vector<double> tt_mins_in = (configList("tt_mins_list"));
			std::vector<double> tt_maxs_in = (configList("tt_maxs_list"));
			m_tt.init(tt_bins_in, tt_mins_in, tt_maxs_in);
			m_tt.use_filter(config("tt_use_filter",false));
			std::vector<double> tt_calib_in = (configList("tt_calib_list"));
			m_tt.setTimeCalib(tt_calib_in);
			m_tt.srcStr(configSrc("tt_source", "DetInfo(:Opal1k)"));
			m_tt.print(*print_itr);
			m_twin.resize(2);
			m_twin = configList("patch_tt_win");
			if (config("tt_slicewin")){
				m_tt.setslicewin(configList("tt_poswin"),
						configList("tt_fwhmwin"),
						configList("tt_amplwin"),
						configList("tt_fwhmPamplwin"));
				}
		}

		logic_itr = use_logic.begin() + use_aq;
		print_itr = print_logic.begin() + print_aq;
		if (*logic_itr){
			std::cerr << "Entered the locic_itr for initializing vector<m_aq>" << std::endl;
			unsigned totalchannels, samples;
			totalchannels = 0;
			samples = 0;
			std::vector<std::string> aq_source_list = (configList("aq_source_list"));
			std::vector<unsigned> aq_baseline_lims = (configList("aq_baseline_lims"));
			std::vector<unsigned> aq_lims = (configList("aq_lims"));
			m_aq.clear();
			m_aq.resize(aq_source_list.size());
			for (std::vector<std::string>::iterator srcItr = aq_source_list.begin(); 
				srcItr != aq_source_list.end(); 
				++srcItr){
				unsigned i = std::distance(aq_source_list.begin(), srcItr);
				m_aq.at(i).init(aq_lims,aq_baseline_lims);
				m_aq.at(i).use(*logic_itr);
				m_aq.at(i).srcStr((Source)*srcItr);
				m_aq.at(i).print(*print_itr);
				m_aq.at(i).invert(config("aq_invert",false));
				totalchannels += m_aq.at(i).nchannels();
				unsigned s = m_aq.at(i).nsamples();
				if ( samples < s)
					samples = s;
			}

			std::cerr << "Exiting the locic_itr for initializing vector<m_aq>" << std::endl;
			std::cerr << " totalchannels = " << totalchannels << std::endl;
			m_patch_compute = config("patch_compute",false);
		}
		

		logic_itr = use_logic.begin() + use_learn;
		print_itr = print_logic.begin() + print_learn;
		if (m_learn.use(*logic_itr)){ // eventually, make these switchable between library building and using
			if ((bool)config("learn_read_pca",false)){
				std::cout << "m_learn.init(" << configStr("learn_pca_read_file" ) << ");" << std::endl;
				m_learn.init(configStr("learn_pca_read_file"));
			} else {
				m_learn.init((unsigned)config("learn_pca_samples_limit",25)
						, (unsigned)config("learn_pca_features",10)
						, (float)config("learn_pca_variance", 0.9)
					    );
			}
			m_learn.print(*print_itr);
		}

		//m_gdetSum_2d.resize(boost::extents[m_eb.bins(Ebeam::energy)][m_gd.bins(Gdet::average)]);
		//m_gdetShots_2d.resize(boost::extents[m_eb.bins(Ebeam::energy)][m_gd.bins(Gdet::average)]);


		m_skip_events = config("skip_events",0);
		m_last_event = config("last_event",100);
		m_gnuplotting = config("gnuplotting",false);
		m_rank = MPI::COMM_WORLD.Get_rank();
		m_mpi_size = MPI::COMM_WORLD.Get_size();
		m_count_event.resize(m_mpi_size,0); // Watch this MPI related one when undoing MPI
		m_failed_event.resize(m_mpi_size,0);
	}

	//--------------
	// Destructor --
	//--------------
	CookieBox_mod::~CookieBox_mod ()
	{

	}

	bool CookieBox_mod::isApprovedByCounters()
	{
		if (m_count_event.at(m_rank) > m_last_event) {
			stop();
			std::cout << "Skipping shot " << m_count_event.at(m_rank) << " since is above m_last_event" << std::endl;
			return false;
		}
		std::vector<unsigned>::iterator itr;
		itr = std::find(m_skipIDs.begin(),m_skipIDs.end(),m_count_event.at(m_rank));
		if (itr != m_skipIDs.end()){ // shotID was found in the skipIDs vector //
			std::cout << "Skipping shot " << m_count_event.at(m_rank) << " since appears in skip vector" << std::endl;
			return false;
		}
		return (m_count_event.at(m_rank) > m_skip_events);
	}


	/// Method which is called once at the beginning of the job
	void CookieBox_mod::beginJob(Event& evt, Env& env)
	{
		m_rank = MPI::COMM_WORLD.Get_rank();
		std::cout << "Started My Job for rank " << m_rank << std::endl;

		if (m_rank != m_root_rank){
			std::cout << "returning from beginJob() in rank " << m_rank << std::endl;
			return;
		}
		return;
	}

	/// Method which is called at the beginning of the run
	void CookieBox_mod::beginRun(Event& evt, Env& env)
	{
		m_rank = MPI::COMM_WORLD.Get_rank();
		std::cout << "Started My Run for rank " << m_rank << std::endl;

		std::stringstream ranktail;
		ranktail << ".rank" << m_rank;

		m_str_runnum = ImgAlgos::stringRunNumber(evt);
		m_str_experiment = ImgAlgos::stringExperiment(env);

		if (m_rank == m_root_rank) std::cout << "filenames = ";
		std::string filename;

		if (m_eb.print()){
			filename = m_datadir + configStr("fname_eb_prefix",  "eb") + "-" + m_str_experiment + "-r" + m_str_runnum + ranktail.str();
			m_eb.open_file( filename );
			m_eb.print_header();
			if (m_rank == m_root_rank) std::cout << filename << "\t";
		}

		if (m_tt.print()){
			filename = m_datadir + configStr("fname_tt_prefix",  "tt") + "-" + m_str_experiment + "-r" + m_str_runnum + ranktail.str();
			m_tt.open_file( filename );
			m_tt.print_header();
			if (m_rank == m_root_rank) std::cout << filename << "\t";
		}

		if (m_gd.print()){
			filename = m_datadir + configStr("fname_gd_prefix",  "gd") + "-" + m_str_experiment + "-r" + m_str_runnum + ranktail.str();
			m_gd.open_file( filename );
			m_gd.print_header();
			if (m_rank == m_root_rank) std::cout << filename << "\t";
		}

		if (m_ta.print()){
			filename = m_datadir + configStr("fname_ta_prefix",  "ta") + "-" + m_str_experiment + "-r" + m_str_runnum + ranktail.str();
			m_ta.open_file( filename );
			m_ta.print_header();
			if (m_rank == m_root_rank) std::cout << filename << "\t";
		}

		if (m_xt.print() || m_xt.print_avg() ){
			filename = m_datadir + configStr("fname_xt_prefix",  "xt") + "-" + m_str_experiment + "-r" + m_str_runnum + ranktail.str();
			m_xt.open_file( filename );
			std::string vhfilename = filename;
			std::string statsfilename = filename;
			vhfilename += ".valhist";
			statsfilename += ".stats";
			if (m_rank == m_root_rank) std::cout << filename << "\t";
			m_xt.open_file_vh( vhfilename );
			if (m_rank == m_root_rank) std::cout << vhfilename << "\t";
			m_xt.open_file_stats( statsfilename );
			if (m_rank == m_root_rank) std::cout << statsfilename << "\t";
			m_xt.print_header();
			m_xt.print_header_vh();
			m_xt.print_header_stats();

		}

		//std::cerr << "Entering the m_aq.size() loop for initializing file" << std::endl;
		for (unsigned i=0;i<m_aq.size();++i){
			if (m_aq[i].use()){
				m_aq[i].getConfig(env);
				unsigned totalchannels = 0;
				unsigned samples = 0;
				totalchannels += m_aq.at(i).nchannels();
				if (samples < m_aq[i].nsamples())
					samples = m_aq[i].nsamples();

				if ( m_eb.use() && m_gd.use() && m_tt.use() ) {
					if (m_shots_4d.shape()[0] < m_tt.bins(TimeTool::delay)//delay) //pos)
							|| m_shots_4d.shape()[1] < m_eb.bins(Ebeam::energy)
							|| m_shots_4d.shape()[2] < m_gd.bins(Gdet::average)
							|| m_shots_4d.shape()[3] < totalchannels )
					{
						m_shots_4d.resize(boost::extents[m_tt.bins(TimeTool::delay)]//delay)]//pos)]
								[m_eb.bins(Ebeam::energy)]
								[m_gd.bins(Gdet::average)]
								[totalchannels]);
						for(long long * it = m_shots_4d.origin(); it < m_shots_4d.origin() + m_shots_4d.num_elements() ; ++it)
							*it = 0;
					}

					if (m_data_5d.shape()[0] < m_tt.bins(TimeTool::delay)//delay)//pos)
							|| m_data_5d.shape()[1] < m_eb.bins(Ebeam::energy)
							|| m_data_5d.shape()[2] < m_gd.bins(Gdet::average)
							|| m_data_5d.shape()[3] < totalchannels 
							|| m_data_5d.shape()[4] < samples )
					{
						m_data_5d.resize(boost::extents[m_tt.bins(TimeTool::delay)]//delay)]//pos)]
								[m_eb.bins(Ebeam::energy)]
								[m_gd.bins(Gdet::average)]
								[totalchannels]
								[samples]);
						for(long long * it = m_data_5d.origin(); it < m_data_5d.origin() + m_data_5d.num_elements() ; ++it)
							*it = 0;
					}
				}

				if (m_aq.at(i).print()){
					std::cerr << "setting up acqiris file" << std::endl;
					std::string acqiris = "crate" + boost::lexical_cast<std::string>(i);
					filename = m_datadir + configStr("fname_aq_prefix",  "aq") + "_" + acqiris;
					filename += "-" + m_str_experiment + "-r" + m_str_runnum + ranktail.str();
					m_aq.at(i).open_file( filename );
					if (m_rank == m_root_rank) 
						std::cout << filename << "\t";
				}
				std::cout << std::endl;
			}
		}
		if (m_aq.front().use()){
			m_makeRotorProjections = (config("sp_makeRotorProjections",false));
			m_printLegendreFFTs = (config("sp_printLegendreFFTs",false));
			m_gaussroll = (config("sp_rolloff",false));
			if (m_gaussroll){
				m_gaussroll_center = (config("sp_rolloff_center",360.));
				m_gaussroll_width = (config("sp_rolloff_width",200.));
			}
			m_use_e2t = config("aq_useEnergyToTime",false);
			if (m_use_e2t){
				m_sp_projectionMask = (configList("sp_projectionMask"));
				m_sp_binsvec = (configList("sp_bins_list"));
				m_sp_startsvec = (configList("sp_starts_list"));
				m_sp_stepsvec = (configList("sp_steps_list"));
				std::string filename = (configStr("aq_transmissions_file","amoi0314_data/corr_factors.taylor")); 
				if (! (this->*read_corr_factors)(filename)){
					std::cerr << "Failed to set the transmission correction coefficients" << std::endl;
					m_use_e2t = false;
				}
				filename = (configStr("aq_energytotime_file","amoi0314_data/energytotime.taylor")); 
				if (! (this->*read_e2t_conversion)(filename)){
					std::cerr << "Failed to set the energy to time coefficients" << std::endl;
					m_use_e2t = false;
				}
			}
		}

		if (m_learn.print()){
			filename = m_datadir + configStr("fname_learn_prefix",  "learn");
			filename += "-" + m_str_experiment + "-r" + m_str_runnum + ranktail.str();
			m_learn.open_file(filename);
			if (m_rank == m_root_rank) 
				std::cout << "set up learning file" << filename << "\n";
		}


		time(&m_beginruntime);

		return;
	}

	/// Method which is called at the beginning of the calibration cycle
	void CookieBox_mod::beginCalibCycle(Event& evt, Env& env)
	{
		m_rank = MPI::COMM_WORLD.Get_rank();
		std::cout << "Started my CalibCycle in rank " << m_rank << std::endl;
		if (m_ta.use()){
			const EpicsStore& estore = env.epicsStore();
			const std::string& value = estore.value("AMO:LAS:DLS:05:MTR.RBV", 0);
			m_ta.updateStage(boost::lexical_cast<float>(value));
			std::cout << "passing out of beginCalibCycle() in rank " << m_rank ;
			std::cout << "\t with m_ta.m_accum_index = " << m_ta.m_accum_index << std::endl;
		}
	}
	//
	/// Method which is called with event data, this is the only required 
	/// method, all other methods are optional
	void CookieBox_mod::event(Event& evt, Env& env)
	{
		m_rank = MPI::COMM_WORLD.Get_rank();
		++ m_count_event.at(m_rank);
		//std::cout << "in event number " << m_count_event.at(m_rank) << " for rank " << m_rank << std::endl;

		if (!isApprovedByCounters()) 
			return; // if out of skip and limit bounds, return from event

		if (!processEvent(evt,env)) {
			//std::cout << "processing event failed" << std::endl;
			++ m_failed_event.at(m_rank);
			return; // if processing fails, return from event
		}
		//std::cout << " finished event " << m_count_event.at(m_rank) << " for rank " << m_rank << std::endl;

		return;
	}

	

	bool CookieBox_mod::sampleevery(unsigned in = 1)
	{
		if (in < 1){in = 1;}
		return (m_count_event.at(m_rank)%in < 2);
	}

	bool CookieBox_mod::processEvent(Event& evt,Env& env)
	{
		// use evt.put(); to add variable for access later in douwnstream events //
		// std::cout << "entered CookieBox_mod::processEvent() method" << std::endl;
		/*
		if (m_count_event.at(m_rank) %m_print_every < 2 || m_count_event.at(m_rank) > 3900){
			std::cout << "beginning processing event " << m_count_event.at(m_rank) << std::endl;
		}
		*/

		
		bool gdslice = false;
		if( m_gd.use() ){
			gdslice = m_gd.fill(evt,env);
			if (!gdslice) {
				std::cerr << "Failed the Gdet for shot " << m_count_event.at(m_rank) << std::endl;
				return false;
			}
		}

		if ( m_tt.use()){
			if ( gdslice ) {
				if (!m_tt.fill(evt,env,m_gd.value()) ) { 
					//std::cerr << "Failed the TimeTool for shot " << m_count_event.at(m_rank) << std::endl;
					return false; 
				}
			} else { 
				if (!m_tt.fill(evt,env)){
					//std::cerr << "Failed the TimeTool for shot " << m_count_event.at(m_rank) << std::endl;
					return false;
				}
			}
		}

		if ( m_eb.use() && !m_eb.fill(evt,env) ) { 
			std::cerr << "Failed the Ebeam for shot " << m_count_event.at(m_rank) << std::endl;
			return false; 
		}
		
		if (m_ta.use())
		{ 
			if ((m_ta.remove_dark() || m_ta.write_dark()) && !m_ta.fill(evt,env) ){
				std::cerr << "Failed the TransAbs for shot " << m_count_event.at(m_rank) << std::endl; 
				return false;
			}
		}

		// Order matters here.  fill, then process, then runstats...
		if ( m_xt.use() && !(	m_xt.fill(evt,env) 
					&& m_xt.process_image(m_count_event.at(m_rank)) 
					&& m_xt.runstats() 
					&& m_xt.crop_filtered(m_count_event.at(m_rank))
				    )
		   ) {
			return false;
		}


		for (unsigned i=0;i< m_aq.size();++i){
			unsigned startchanind = 0;
			if (m_aq[i].use() 
				&& m_eb.use() && m_gd.use() && m_tt.use()
				&& m_eb.filled() && m_gd.filled() && m_tt.filled()){
				//a4d_ll_t::index_gen indices;
				// now switching to 5d to include the timetool data in histogram 
				a5d_ll_t::index_gen indices;
				a5d_ll_t::index_gen integindices;
				a4d_ll_t::index_gen shotindices;
				/*
				std::cout << "checking indices m_eb.index(Ebeam::energy), m_gd.index(Gdet::average) \t" 
					<< m_eb.index(Ebeam::energy) << ", "
					<< m_gd.index(Gdet::average)
					<< std::endl;
				std::cout << "m_data_4d.shape() = " 
					<< m_data_4d.shape()[0] << ", "
					<< m_data_4d.shape()[1] << ", "
					<< m_data_4d.shape()[2] << ", "
					<< m_data_4d.shape()[3]
					<< std::endl;
				*/
				a5d_ll_t::array_view<2>::type slice = 
					m_data_5d[ 
						indices	[m_tt.index(TimeTool::delay)]//delay)]//pos)]
							[m_eb.index(Ebeam::energy)]
							[m_gd.index(Gdet::average)]
							[range(startchanind , startchanind + m_aq[i].nchannels())]
							[range()]
						];
				a4d_ll_t::array_view<1>::type shotslice = 
					m_shots_4d[ 
						shotindices [m_tt.index(TimeTool::delay)]//delay)]//pos)]
							[m_eb.index(Ebeam::energy)]
							[m_gd.index(Gdet::average)]
							[range(startchanind , startchanind + m_aq[i].nchannels())]
						];

				// Ideally, make a class of accumulators (or signal accumulators) to pass all these slices and things //
				if ( m_aq[i].fill(evt,env,slice,shotslice)) { 
					startchanind += m_aq[i].nchannels();
				} else {
					std::cerr << "Failed the Acqiris, skipping the print and rest of event " 
						<< m_count_event.at(m_rank) << std::endl;
					return false;
				}
			} else {
				std::cerr << "failed to fill eb, gd or tt, so trying to call non-slice m_aq.fill()" << std::endl;
				if ( m_aq[i].use() && m_aq[i].fill(evt,env)) { 
					startchanind += m_aq[i].nchannels();
				} else {
					std::cerr << "Failed the Acqiris, skipping the print and rest of event " 
						<< m_count_event.at(m_rank) << std::endl;
					return false;
				}
			}
		}

		if ( sampleevery(m_print_every)  ){
			if (m_ta.use() && m_ta.print() ){
				if ( !( m_ta.print_out(m_count_event.at(m_rank)))){
					std::cerr << "Failed to print TransAbs for shot " << m_count_event.at(m_rank) << std::endl;	
				}
			}
			if (m_xt.use() && m_xt.print() ){
				if ( !(	m_xt.print_out( m_count_event.at(m_rank) ) 
					&& m_xt.print_filtered(m_count_event.at(m_rank)) 
					&& m_xt.print_cropped(m_count_event.at(m_rank)) 
					) )
					std::cerr << "Failed to print Xtcav for shot " << m_count_event.at(m_rank) << std::endl;
				if ( !m_xt.print_out_vh( m_count_event.at(m_rank) ))
					std::cerr << "Failed to print Xtcav valhist for shot " << m_count_event.at(m_rank) << std::endl;
				if ( !m_xt.print_out_stats( m_count_event.at(m_rank) ))
					std::cerr << "Failed to print Xtcav stats for shot " << m_count_event.at(m_rank) << std::endl;
			}
			if (m_eb.use() && m_eb.print() && !m_eb.print_out(m_count_event.at(m_rank)))
				std::cerr << "Failed to print Ebeam for shot " << m_count_event.at(m_rank) << std::endl;
			if (m_gd.use() && m_gd.print() && !m_gd.print_out(m_count_event.at(m_rank)))
				std::cerr << "Failed to print Gdet for shot " << m_count_event.at(m_rank) << std::endl;
			if (m_tt.use() && m_tt.print() && !m_tt.print_out(m_count_event.at(m_rank)))
				std::cerr << "Failed to print TimeTool for shot " << m_count_event.at(m_rank) << std::endl;
			for (unsigned i=0;i< m_aq.size();++i){
				/*
				std::cout << "in sampleevery(m_print_every) at event " << m_count_event.at(m_rank) 
					<< " and rank " << m_rank << "\n\t... trying to print m_aq.vectors" << std::endl;
				*/
				if (m_aq[i].use() && m_aq[i].print() && !m_aq[i].print_out( m_count_event.at(m_rank) ))
					std::cerr << "Failed to print Acqiris for shot " << m_count_event.at(m_rank) << std::endl;
			}

		}

		if (m_learn.use()) {
			std::vector<float> buildfeatures, projectedfeatures;
			buildfeatures.clear();
			projectedfeatures.clear();

			std::string features_labels("#");
			// This is where we need to be careful, need ebeam for nolasing, but not gas dets, ypos though might be orbit bumped //
			if (m_eb.use()){
				m_eb.addfeatures(buildfeatures);
				if (m_learn.nfeatures()==0)
					features_labels += m_eb.get_features_labels();
			}
			if (m_gd.use()){
				m_gd.addfeatures(buildfeatures);
				if (m_learn.nfeatures()==0)
					features_labels += m_gd.get_features_labels();
			}
			if (m_xt.use()){
				//m_xt.addfeatures(buildfeatures);
				m_xt.addfeatures_centroids_image(buildfeatures);
				if (m_learn.nfeatures()==0)
					features_labels += m_xt.get_features_labels();
			}
			

			if (m_learn.write_pca()){ // If writing the no_lasing PCA for use later in identifying references //
				if (m_count_event.at(m_rank) % m_print_every < 2){
					std::cout << "\t\ttWarning:\t"
						<< "we should be reading in PCA so this means it's no_lasing\n"
						<< "\t\tso no need to finish event with acquiris multi_array()\n\n"
						<< std::flush;
				}

				if (m_learn.nfeatures()==0)
					m_learn.set_features_labels(features_labels);
				if (m_learn.fill(buildfeatures)){
					if (m_count_event.at(m_rank) %m_print_every < 2){
						std::cout << "filling features at event " << m_count_event.at(m_rank) << std::endl;
					}
				} else {
					if (!m_learn.run_pca( (unsigned)config("learn_pca_features",100) ) ) {
						std::cout << "something failed in m_learn::run_pca() call" << std::endl;
					}
					if (m_learn.print()
							&& m_count_event.at(m_rank) %m_print_every < 2
							&& m_learn.project_pca(buildfeatures,projectedfeatures) )
					{
						//std::cerr << "buildfeatures.size() = " << buildfeatures.size() << std::endl;
						std::string filename(m_learn.filename());
						filename += ".projected_event." + boost::lexical_cast<std::string>(m_count_event.at(m_rank));
						std::ofstream outfile (filename.c_str(),std::ios::out);
						for (unsigned i=0;i<projectedfeatures.size();++i){
							outfile << projectedfeatures[i] << "\n";
						}
						outfile.close();

						//if (m_learn.testfeatures(buildfeatures)){
							m_testpca_features.push_back(buildfeatures);
							//std::cout << "pushing test features at event " << m_count_event.at(m_rank) << std::endl;
						//}
					}
					// HERE consider back projecting for comparisons //
				}
				return true;
			} // end write no_lasing PCA for use in identifying references

			std::cout << "m_learn.read_pca() = " << m_learn.read_pca() << std::endl;

			if (m_learn.read_pca()){
				if (m_count_event.at(m_rank) % m_print_every < 2){
					std::cout << "\tattempting to read in PCA:\t" 
						<< "this presupposes we already read in the PCA from no_lasing\n" 
						<< std::endl;
				}
				// body of with lasing actions... putting this second since it ivolves the multidim arrays eventually

				return ( m_learn.print()
					&& m_count_event.at(m_rank) %m_print_every < 2
					&& m_learn.difference_pca_images(buildfeatures
						,m_xt.eigenface_shape()
						,m_xt.nonimage_features()
						,m_count_event.at(m_rank) 
						)
					);

			}

		} // endif m_learn.use()

		/*
		std::cout << "\n\t\t --- finished event " << m_count_event.at(m_rank) 
			<< " --- \n \t\t --- at rank " << m_rank << "---" << std::endl;
		*/



		return true;
	}


	/// Method which is called at the end of the calibration cycle
	void CookieBox_mod::endCalibCycle(Event& evt, Env& env)
	{
		// Like below, collect in rank_master and print from there.
		// OK, now we need t collect with an MPI::Reduce all the ranks m_ta.m_data_accum, m_data_accum_ref, m_shots_accum, and m_shots_accum_ref
		// we will try to make m_data_accum and it's partners private, but friend them to CookieBox... if this fails, then public them.
		std::cerr << "made it into void CookieBox_mod::endCalibCycle(Event& evt, Env& env)" << std::endl;
		if (m_ta.use() && (m_ta.remove_dark() || ( m_ta.in_dark_stage_list() && m_ta.write_dark() ) ) ){
		if (m_rank == m_root_rank){
			if (m_ta.m_data_accum[m_ta.m_accum_index]->num_elements() > 0  
				&& m_ta.m_data_accum_ref[m_ta.m_accum_index]->num_elements() > 0  
				&& m_ta.m_shots_accum[m_ta.m_accum_index] > 0  
				&& m_ta.m_shots_accum_ref[m_ta.m_accum_index] > 0 
				){
				MPI::COMM_WORLD.Reduce(
					MPI::IN_PLACE,
					m_ta.m_data_accum[m_ta.m_accum_index]->origin(),// these are now vectors of pointers
					m_ta.m_data_accum[m_ta.m_accum_index]->num_elements(),
					MPI::LONG_LONG,
					MPI::SUM,
					m_root_rank); 
				MPI::COMM_WORLD.Reduce(
					MPI::IN_PLACE,
					m_ta.m_data_accum_ref[m_ta.m_accum_index]->origin(),
					m_ta.m_data_accum_ref[m_ta.m_accum_index]->num_elements(),
					MPI::LONG_LONG,
					MPI::SUM,
					m_root_rank); 
				MPI::COMM_WORLD.Reduce(
					MPI::IN_PLACE,
					&(m_ta.m_shots_accum[m_ta.m_accum_index]),
					1,
					MPI::LONG_LONG,
					MPI::SUM,
					m_root_rank); 
				MPI::COMM_WORLD.Reduce(
					MPI::IN_PLACE,
					&(m_ta.m_shots_accum_ref[m_ta.m_accum_index]),
					1,
					MPI::LONG_LONG,
					MPI::SUM,
					m_root_rank); 
			}
			std::cout << "passing further into endCalibCycle() in rank " << m_rank ;
			std::cout << "\t with m_ta.m_accum_index = " << m_ta.m_accum_index << std::endl;
		} else {
			if (m_ta.m_data_accum[m_ta.m_accum_index]->num_elements() > 0  
				&& m_ta.m_data_accum_ref[m_ta.m_accum_index]->num_elements() > 0  
				&& m_ta.m_shots_accum[m_ta.m_accum_index] > 0 
				&& m_ta.m_shots_accum_ref[m_ta.m_accum_index] > 0 
				) {
				MPI::COMM_WORLD.Reduce(
					m_ta.m_data_accum[m_ta.m_accum_index]->origin(),
					m_ta.m_data_accum[m_ta.m_accum_index]->origin(),
					m_ta.m_data_accum[m_ta.m_accum_index]->num_elements(),
					MPI::LONG_LONG,
					MPI::SUM,
					m_root_rank); 
				MPI::COMM_WORLD.Reduce(
					m_ta.m_data_accum_ref[m_ta.m_accum_index]->origin(),
					m_ta.m_data_accum_ref[m_ta.m_accum_index]->origin(),
					m_ta.m_data_accum_ref[m_ta.m_accum_index]->num_elements(),
					MPI::LONG_LONG,
					MPI::SUM,
					m_root_rank); 
				MPI::COMM_WORLD.Reduce(
					&(m_ta.m_shots_accum[m_ta.m_accum_index]),
					&(m_ta.m_shots_accum[m_ta.m_accum_index]),
					1,
					MPI::LONG_LONG,
					MPI::SUM,
					m_root_rank); 
				MPI::COMM_WORLD.Reduce(
					&(m_ta.m_shots_accum_ref[m_ta.m_accum_index]),
					&(m_ta.m_shots_accum_ref[m_ta.m_accum_index]),
					1,
					MPI::LONG_LONG,
					MPI::SUM,
					m_root_rank); 
			}
			std::cout << "returning from endCalibCycle() in rank " << m_rank;
			std::cout << "\t with m_ta.m_accum_index = " << m_ta.m_accum_index << std::endl;
			return;
		}
		std::cerr << "made it through the MPI::COMM_WORLD.Reduce() method in void CookieBox_mod::endCalibCycle(Event& evt, Env& env)" << std::endl;

		// HERE continue to print out for the root_rank
		// Ultimately, save ideally the phase unwrapped, weighted average over signal of the phase and amplitude of the FFT 
		// (and maybe back transforms of that)
		// That will then get read into the real runs for phase multiplication with the time sorted tofs.
			try{
				if ( !m_ta.updateStage_close() ){throw;}
			} catch (std::exception & e){
				std::cerr << "Something failed here in endCalibCycle() with m_ta, e.what() = " << e.what() << std::endl;
				return;
			}
		}
		return;
	}

	/// Method which is called at the end of the run
	void CookieBox_mod::endRun(Event& evt, Env& env)
	{
		time(&m_endruntime);
		double dtime = difftime(m_endruntime,m_beginruntime);
		std::cout << "processing events took " << dtime << " seconds" << std::endl;
		/*
		if ( m_ta.use() && m_ta.print_avg() ){
			std::cerr << "Failed to print TransAbs average" << std::endl;
		}
		*/
		if (m_ta.use() && m_ta.write_dark() && !m_ta.write_dark_file()){
			std::cerr << "Failing to m_ta.write_dark_file()" << std::endl;
		}
		if ( m_xt.use() && m_xt.print_avg() && ! m_xt.print_out_avg()){
			std::cerr << "Failed to print Xtcav average" << std::endl;
		}
		/*
		 *
		 * 	Here is where we would use MPI::ScatterV() or Send() and Recieve() to spawn a new parallel portion.
		 * 	We want to send a block of acq vectors to be FFTd and whatnot
		 * 	regather pattern recognize, spawn further
		 *
		 */

		m_rank = MPI::COMM_WORLD.Get_rank();

		
		std::cout << "attempting to collect MPI in rank " << m_rank << std::endl;
		
		
		if (m_rank == m_root_rank){
			if (m_data_5d.num_elements() > 0  && m_shots_4d.num_elements() > 0 ) {
				MPI::COMM_WORLD.Reduce(MPI::IN_PLACE,m_data_5d.origin(),m_data_5d.num_elements(),MPI::LONG_LONG,MPI::SUM,m_root_rank); 
				MPI::COMM_WORLD.Reduce(MPI::IN_PLACE,m_shots_4d.origin(),m_shots_4d.num_elements(),MPI::LONG_LONG,MPI::SUM,m_root_rank); 
			}
			std::cout << "passing further into endRun() in rank " << m_rank ;
			std::cout << " with failed events/total = " << m_failed_event.at(m_rank) << " / " << m_count_event.at(m_rank) << std::endl;
		} else {
			if (m_data_5d.num_elements() > 0  && m_shots_4d.num_elements() > 0 ) {
				MPI::COMM_WORLD.Reduce(m_data_5d.origin(),m_data_5d.origin(),m_data_5d.num_elements(),MPI::LONG_LONG,MPI::SUM,m_root_rank); 
				MPI::COMM_WORLD.Reduce(m_shots_4d.origin(),m_shots_4d.origin(),m_shots_4d.num_elements(),MPI::LONG_LONG,MPI::SUM,m_root_rank); 
			}
			std::cout << "returning from endRun() in rank " << m_rank ;
			std::cout << " with failed events/total = " << m_failed_event.at(m_rank) << " / " << m_count_event.at(m_rank) << std::endl;
			return;
		}
		std::cout << "finished collecting MPI" << std::endl;

		if (m_learn.use()){
			if ((bool)config("learn_write_pca",false)){ // if asked to write
				if (!m_learn.write_pca(configStr("learn_pca_write_file"))){ // if failed to write when it was asked
					std::cerr << "Seemed to have failed to m_learn.write_pca() in CookieBox_mod::endRun() method" << std::endl;
				}
			} else {
				if (!m_learn.read_pca(configStr("learn_pca_read_file"))){ // if failed to write when it was asked
					std::cerr << "Seemed to have failed to m_learn.read_pca() in CookieBox_mod::endRun() method" << std::endl;
				}
				m_learn.inspect_eigen_pca();
				m_learn.print_eigenfaces( m_xt.eigenface_shape() );
			}
			m_learn.compare_pca_images(m_testpca_features,m_xt.eigenface_shape(),m_xt.nonimage_features());

			//m_learn.compare_pca_nonimage(m_testpca_features,m_xt.nonimage_features());
			//m_learn.compare_pca(m_testpca_features);
		} // end if m_learn.use()

		// HERE 
		// do a comparison for this by constructing the pca after some number of shots, then compare for new shots as they come in.
		// When two or more features are correlated, then we should plot the data in ortho and para for those components
		// Let's also plot up the raw features as points and as parallel coords (with mean subtracted, and scaled by the variance in that dim)
		// Another thing to try is to flatten machine features that have no correlation.... 
		// check the gas dets since they are nearly perfectly correlated to make sure we're doing the right thing.

		// technically we should close our files here //
		if (m_eb.print())
			m_eb.close_file();
		if (m_tt.print())
			m_tt.close_file();
		if (m_ta.print())
			m_ta.close_file();
		if (m_xt.print())
			m_xt.close_file();
		if (m_gd.print())
			m_gd.close_file();
		// etc... except that in the object destructor I test and close the files.



		// HERE printing out the 5D histogram as ascii... this is stupid//
		std::cout << "printing out the spectra ... \t" << std::flush;
		time_t starttime,stoptime;
		float elapsed=0.;
		time(&starttime);
		if ( !(m_data_5d.num_elements() >0 && m_shots_4d.num_elements() > 0) ) {
			std::cerr << "Unfilled m_data_5d.num_elements()" << std::endl;
			return;
		}

		if (!fillLegendreVecs()){std::cerr << "Failed to fillLegendreVecs()" << std::endl;}
		// std::cerr << "Entering printSpectraLegendre() \n\n\n\n \t\t\t HERE I AM!! \n\n\n\n" << std::endl;
		if (! printSpectraLegendre()){std::cerr << "Failed to printLegendreSpectra()" << std::endl;}
		std::cerr << "Entering print_patch_results() " << std::endl;
		if (m_patch_compute && !print_patch_results()) { std::cerr << "Failed printing patch results" << std::endl;}
		//std::cout << "skipping/printing out the spectra ...\t" << std::flush;
		std::cerr << "Entering printSpectra() " << std::endl;
		//if (!printSpectra() ) {std::cerr << "Failed to printSpectra()" << std::endl;return;}
		//if (!fftDiffSpectra()){std::cerr << "Failed to fftDiffSpectra()" << std::endl;return;}			

		std::cout << "Putting back rotor projections of diff spectra for Kareem" << std::endl;
		//std::cout << "Skipping rotor projections of diff spectra until there are proper rotor bases in time domain " << std::endl;
		if (m_makeRotorProjections){
			std::string orthofilename(configStr("sim_orthofile","amoi0314_data/transabs/ta-amoi0314-r0214.rank0.result.orthonorm.r.93.055"));
			std::cerr << "Projecting spectra " << std::endl;
			if (!projectDiffSpectra(orthofilename)){std::cerr << "Failed to projectDiffSpectra(orthofilename)" << std::endl;return;}
		}

		time(&stoptime);
		elapsed = difftime(stoptime,starttime);
		std::cout << elapsed << " seconds." << std::endl;
		starttime = stoptime;
		std::cout << "skipping/printing out the integration results ...\t" << std::flush;
		/*
		if (!printIntegs()){
			std::cerr << "Failed to printIntegs()" << std::endl;
			return;
		}
		*/
		time(&stoptime);
		elapsed = difftime(stoptime,starttime);
		std::cout << elapsed << " seconds" << std::endl;

		return;
	}

	/// Method which is called once at the end of the job
	void CookieBox_mod::endJob(Event& evt, Env& env)
	{
		//m_rank = MPI::COMM_WORLD.Get_rank();
		std::cout << "endJob() called for rank " << m_rank << std::endl;
	}

	// We should do the integrations once, here, on the histogram. //
	// No need for a second 5d histogram, in fact, we should allow arbitrary hist slicing for integration on any of the indices //
	// make the arguements the index and the ranges to integrate over //
	bool CookieBox_mod::printIntegs(void)
	{

		if (m_rank != m_root_rank){
			std::cerr << "Trying to print out m_data_5d but using non-root rank, rank = " << m_rank << std::endl;
			return false;
		}
		// index order is [0=TimeTool][1=ebeam][2=gasdet][3=chan][4=integwin]
		// enum DimDefs {ttind,ebind,gdind,chanind,tofind};

		std::vector<unsigned> lowlims = (configList("aq_integlims_low"));
		std::vector<unsigned> highlims = (configList("aq_integlims_high"));
		assert(lowlims.size()==highlims.size());


		for (unsigned w=0;w<lowlims.size();++w){
			assert(lowlims[w] >= 0 && lowlims[w] < highlims[w]);
			assert(m_data_5d.shape()[tofind]>highlims[w]);
			for (unsigned e=0;e<m_data_5d.shape()[ebind];++e){
				for (unsigned g=0;g<m_data_5d.shape()[gdind];++g){
					// build filename //
					std::string filename = m_datadir + std::string("integ");
					filename+= "-" + m_str_experiment + "-r" + m_str_runnum;
					filename += "_w" + boost::lexical_cast<std::string>(w);
					filename += "_e" + boost::lexical_cast<std::string>(e);
					filename += "_g" + boost::lexical_cast<std::string>(g);
					filename += std::string(".dat");
					std::ofstream outfile(filename.c_str(),std::ios::out);
					// print header //
					outfile << "# " << lowlims[w] << "\t" << highlims[w] << "\n";
					outfile << "# ";
					for (unsigned c=0;c<m_data_5d.shape()[chanind];++c)
						outfile << "ch" << c << "\t";
					outfile << "\n";
					// foreach time compute the integral and print //
					for (unsigned t = 0; t < m_data_5d.shape()[ttind] ; ++t){
						for (unsigned c=0;c<m_data_5d.shape()[chanind];++c){
							long shots = m_shots_4d[t][e][g][c];
							double val = 0.;
							if (shots > 0){
								for (unsigned i=lowlims[w];i<highlims[w];++i){
									val += m_data_5d[t][e][g][c][i];
								}
								val /= shots;
							}
							outfile << val << "\t";
						}
						outfile << "\n";
					}
					outfile.close();
				}
			}
		}
		return true;
	}

	bool CookieBox_mod::printSpectraLegendre(void)
	{
		// index order is [0=TimeTool][1=ebeam][2=gasdet][3=chan][4=sample]
		// enum DimDefs {ttind,ebind,gdind,chanind,tofind};
		if (m_rank != m_root_rank){
			std::cerr << "Trying to print out Legendre projections of m_data_5d but using non-root rank, rank = " << m_rank << std::endl;
			return false;
		}
		// resize the avg spectra for accumulation //
		//
		/*
		   sp_bins_list = 120 10
		   sp_starts_list = 500. 0
		   sp_steps_list = 0.25 10
		   */




		m_kwin.resize(2);
		m_kwin = configList("patch_ke_win");

		m_avg_legendres_nsums.resize(boost::extents
				[m_data_5d.shape()[ebind]] // [1] ebeam bins
				[m_data_5d.shape()[gdind]] // [2] gasdet bins
		);
		m_avg_legendres_4d.resize(boost::extents
				[m_data_5d.shape()[ebind]]  // [1] ebeam bins
				[m_data_5d.shape()[gdind]]  // [2] gasdet bins
				[m_sp_binsvec[idx_kinenergy]]	// number of tof energy bins
				[m_sp_binsvec[idx_legendre]] 	// number of Legendre coefficients
				);
		m_legendres_5d.resize(boost::extents
				[m_data_5d.shape()[ttind]]  // [0] timetool bins
				[m_data_5d.shape()[ebind]]  // [1] ebeam bins
				[m_data_5d.shape()[gdind]]  // [2] gasdet bins
				[m_sp_binsvec[idx_kinenergy]]	// [0] number of tof energy bins
				[m_sp_binsvec[idx_legendre]] 	// [1] number of Legendre coefficients
				);

		m_totsignal_4d.resize(boost::extents
				[m_data_5d.shape()[ttind]]  // [0] timetool bins
				[m_data_5d.shape()[ebind]]  // [1] ebeam bins
				[m_data_5d.shape()[gdind]]  // [2] gasdet bins
				[m_sp_binsvec[idx_kinenergy]]	// [0] number of tof energy bins
				);

		for(double * it = m_legendres_5d.origin(); it < m_legendres_5d.origin() + m_legendres_5d.num_elements() ; ++it){*it = 0.;}
		for(double * it = m_avg_legendres_4d.origin(); it < m_avg_legendres_4d.origin() + m_avg_legendres_4d.num_elements() ; ++it){*it = 0.;}
		for(unsigned * it = m_avg_legendres_nsums.origin(); it < m_avg_legendres_nsums.origin() + m_avg_legendres_nsums.num_elements() ; ++it){*it = 0;}
		for(double * it = m_totsignal_4d.origin(); it < m_totsignal_4d.origin() + m_totsignal_4d.num_elements() ; ++it){*it = 0.;}

		unsigned sz = m_data_5d.shape()[0];
		unsigned ncomponents = sz/2;
		double * timeslice_r = (double *) fftw_malloc(sizeof(double) * sz);
		double * timeslice_hc = (double *) fftw_malloc(sizeof(double) * sz);

		fftw_plan plan_r2hc = fftw_plan_r2r_1d(sz,
				timeslice_r,
				timeslice_hc,
				FFTW_R2HC,
				FFTW_MEASURE
				);
		fftw_plan plan_hc2r = fftw_plan_r2r_1d(sz,
				timeslice_hc,
				timeslice_r,
				FFTW_HC2R,
				FFTW_MEASURE
				);

		if (m_patch_compute){
			m_patch_bins=0;
			m_patch_ys.clear();
			m_patch_legcoeffs.clear();
			m_patch_ys.resize(m_data_5d.shape()[chanind],0.);
			m_patch_legcoeffs.resize(m_sp_binsvec[idx_legendre],0.);
		}

		// Accumulte and print the time-averaged Legendre coefficients as well.
		// Ultimately we will project these then onto the rotor signal
		// std::vector<double> timeslice(m_legendres_5d.shape()[ttind]);

		a3d_d_t projections_3d;// needs to be in full scope
		if (m_makeRotorProjections){
			std::vector < std::vector <double> > orthobases;
			//std::string orthofilename(configStr("ta_orthofile","amoi0314_data/transabs/ta-amoi0314-r0214.rank0.result.orthonorm.r.93.055"));
			std::string orthofilename(configStr("sim_orthofile"));//,"amoi0314_data/simulation/rotorbases.dat.fouriershifingconstruct"));
			std::ifstream orthofile(orthofilename.c_str(),std::ios::in);
			try {
				std::cerr << "Pulling orthofile " << orthofilename << std::endl;
				orthofile >> orthobases;
				std::cerr << "orthobases.size() = " << orthobases.size() << std::endl;
			} catch (std::exception &e) {
				std::cerr << "Failed to read in files for projectDiffSpectra(), error: " << e.what() << std::endl;
				return false;
			}
			orthofile.close();
			if (orthobases.front().size() > sz)
				condense(orthobases,sz);
			std::cerr << "\torthobases.size() = " << orthobases.size() << "\torthobases.front().size() = " << orthobases.front().size() << std::endl;
			projections_3d.resize(boost::extents
					[orthobases.size() ] // number of orthonormal basis vectors for projections of time dependence
					[m_legendres_5d.shape()[ebind]] // number of ebeam energy bins e
					[m_legendres_5d.shape()[3]] // number of electron kinetic energy bins k
					);
			std::cerr << "\n\tprojections_3d.shape()[ "  << 0 << " ]=\t" << projections_3d.shape()[0] << "\n" << std::endl;
			std::cerr << "\n\tprojections_3d.shape()[ "  << 1 << " ]=\t" << projections_3d.shape()[1] << "\n" << std::endl;
			std::cerr << "\n\tprojections_3d.shape()[ "  << 2 << " ]=\t" << projections_3d.shape()[2] << "\n" << std::endl;
		}
		std::cerr << "\n\tm_legendres_5d.shape()[ "  << ttind << " ]=\t" << m_legendres_5d.shape()[ttind] << "\n" << std::endl;
		std::cerr << "\n\tm_legendres_5d.shape()[ "  << ebind << " ]=\t" << m_legendres_5d.shape()[ebind] << "\n" << std::endl;
		std::cerr << "\n\tm_legendres_5d.shape()[ "  << gdind << " ]=\t" << m_legendres_5d.shape()[gdind] << "\n" << std::endl;
		std::cerr << "\n\tm_legendres_5d.shape()[ "  << 3 << " ]=\t" << m_legendres_5d.shape()[3] << "\n" << std::endl;
		std::cerr << "\n\tm_legendres_5d.shape()[ "  << 4 << " ]=\t" << m_legendres_5d.shape()[4] << "\n" << std::endl;

		for (unsigned e = 0; e<m_legendres_5d.shape()[ebind] ; ++e)
		{
			//unsigned e=1;
			for ( unsigned g = 0; g < m_legendres_5d.shape()[gdind] ; ++g)
			{
				//unsigned g = 1;

				std::string totfilename = m_datadir + std::string("tot_signal");
				std::string details = "-" + m_str_experiment + "-r" + m_str_runnum;
				details += "_e" + boost::lexical_cast<std::string>(e);
				details += "_g" + boost::lexical_cast<std::string>(g);
				totfilename += details + std::string(".dat");
				std::ofstream outfile(totfilename.c_str(),std::ios::out);
				outfile << "#\tThis is the total after scaling and masking angles, time-vectors versus energy" << std::endl;

				for (unsigned t=0;t<m_legendres_5d.shape()[ttind];++t)
				{
					if (	m_printMarkusLegendresCompare 
							&& bool(g==unsigned(m_legendres_5d.shape()[gdind]/2)) 
							&& bool(e==unsigned(m_legendres_5d.shape()[ebind]/2)) 
						//	&& bool(t%4==0)
						//	&& bool(t==unsigned(m_legendres_5d.shape()[ttind]/2))
					   ) {
						details = "-" + m_str_experiment + "-r" + m_str_runnum;
						details += "_e" + boost::lexical_cast<std::string>(e);
						details += "_g" + boost::lexical_cast<std::string>(g);
						details += "_t" + boost::lexical_cast<std::string>(t) + std::string(".dat");
						std::string fn = m_datadir + std::string("MarkusCompareChannels") + details;
						m_markusfileChans.open(fn.c_str(),std::ios::out); 
						fn = m_datadir + std::string("MarkusCompareLegendres") + details;
						m_markusfileLegs.open(fn.c_str(),std::ios::out);
					}

					for ( unsigned k = 0; k< m_legendres_5d.shape()[3];++k){
						if (!computeLegendreCoeffs(t,e,g,k)) {
							std::cerr << "Failed to compute Legendre coefficients and print markusfile" << std::endl;
							return false;
						}
						outfile << m_totsignal_4d[t][e][g][k] << "\t";
					}
					if (m_markusfileLegs.is_open()) { 
						m_markusfileLegs.close(); 
					}
					if (m_markusfileChans.is_open()) { 
						m_markusfileChans.close(); 
					}
					outfile << "\n";
				}
				outfile.close();
			}
		}
		
		// need to declare these at high scope //
		std::ofstream outabsfile;
		std::ofstream outargfile;
		std::ofstream outrealfile;

		std::string tlims =  m_tt.getTimeLims();
		std::string tspanStr;
		double tspanDbl = m_tt.getTspan(tspanStr);


		for ( unsigned l=0;l<m_legendres_5d.shape()[4];++l)
		{
			//std::cerr << "\tl\t" << l << std::flush;
			for ( unsigned g = 0; g < m_legendres_5d.shape()[gdind] ; ++g)
			{
				//std::cerr << "\tg\t" << g << std::flush;
				std::string avgfilename = m_datadir + std::string("avg_legendre");
				std::string details = "-" + m_str_experiment + "-r" + m_str_runnum;
				details += "_l" + boost::lexical_cast<std::string>(l);
				details += "_g" + boost::lexical_cast<std::string>(g);
				details += std::string(".dat");
				avgfilename += details;
				std::ofstream avgoutfile(avgfilename.c_str(),std::ios::out);
				for (unsigned e = 0; e<m_legendres_5d.shape()[ebind] ; ++e)
				{
					//std::cerr << "\te\t" << e << std::flush;
					std::string filename = m_datadir + std::string("diff_legendre");
					details = "-" + m_str_experiment + "-r" + m_str_runnum;
					details += "_e" + boost::lexical_cast<std::string>(e);
					details += "_g" + boost::lexical_cast<std::string>(g);
					details += "_l" + boost::lexical_cast<std::string>(l);
					details += "_tspan" + tspanStr;
					if (m_gaussroll)
						details += std::string("_rolloff");
					details += std::string(".dat");
					filename += details;
					std::ofstream outfile(filename.c_str(),std::ios::out);
					outfile << m_tt.getTimeLims();
					if (m_printLegendreFFTs//){ 
						&& e==5
						&& g ==1
						&& l%2==0
						){
						std::string absfilename = m_datadir + std::string("fftabs_legendres");
						std::string argfilename = m_datadir + std::string("fftarg_legendres");
						std::string realfilename = m_datadir + std::string("fftreal_legendres");
						std::string filetail = "-" + m_str_experiment + "-r" + m_str_runnum;
						filetail += "_e" + boost::lexical_cast<std::string>(e);
						filetail += "_g" + boost::lexical_cast<std::string>(g);
						filetail += "_l" + boost::lexical_cast<std::string>(l);
						filetail += "_tspan" + tspanStr + "_tbins" + boost::lexical_cast<std::string>(m_legendres_5d.shape()[ttind]);
						if (m_gaussroll){
							filetail += "_rolloff_c" + boost::lexical_cast<std::string>(m_gaussroll_center);
							filetail += "_w" + boost::lexical_cast<std::string>(m_gaussroll_width);
						}
						filetail += std::string(".dat");
						absfilename += filetail;
						argfilename += filetail;
						realfilename += filetail;
						//std::cerr << "writing to: " << absfilename << "\t" << argfilename << std::endl;
						outabsfile.open(absfilename.c_str(),std::ios::out);
						outargfile.open(argfilename.c_str(),std::ios::out);
						outabsfile << "#fft abs\n";
						outargfile << "#fft arg\n";
						outrealfile.open(realfilename.c_str(),std::ios::out);
						outrealfile << "#back fft real\n";
					} // close if (m_printLegendreFFTs){
					for ( unsigned k = 0; k< m_legendres_5d.shape()[3];++k){
						//std::cerr << "\tk\t" << k << std::flush;
						double avg = m_avg_legendres_4d[e][g][k][l]/(double)m_avg_legendres_nsums[e][g];
						avgoutfile << avg << "\t";
						for (unsigned t=0;t<m_legendres_5d.shape()[ttind];++t){
							//std::cerr << "\tt\t" << t << std::endl;
							double result = 0.;
							unsigned shots = m_shots_4d[t][e][g][0];
							if (shots > 0){
								result = m_legendres_5d[t][e][g][k][l] - avg;
							}
							outfile << result << "\t";
							timeslice_r[t] = result; 
						}
						outfile << "\n";
						if (m_printLegendreFFTs
							&& outabsfile.is_open()
							&& outargfile.is_open()
							&& outrealfile.is_open())
						{ 
							detrend(timeslice_r,sz);
							if (m_gaussroll){
							//	gaussroll(timeslice_r,sz,m_gaussroll_center,m_gaussroll_width);
								sin2roll(timeslice_r,sz,m_gaussroll_center,m_gaussroll_width);
							}
							fftw_execute_r2r(plan_r2hc,
									timeslice_r,
									timeslice_hc
									);
							std::complex<double> z(timeslice_hc[0],0.0);
							outabsfile << std::abs(z) << "\t";
							outargfile << std::arg(z) << "\t";
							for(unsigned f=1;f<sz/2;++f){
								z = std::complex<double>(timeslice_hc[f],timeslice_hc[sz-f]);
								double abs = std::abs(z);
								double arg = std::arg(z);
								outabsfile << abs << "\t";
								outargfile << arg << "\t";
							}
							z = std::complex<double>(0.0,timeslice_hc[sz/2]);
							outabsfile << std::abs(z) << "\n";
							outargfile << std::arg(z) << "\n";

							if (config("sp_fixedfilter",false))
								fixedfilter(timeslice_hc,sz);
							/*
							double noise(.1);                         
							timeslice_hc[0] *= double(1);
							timeslice_hc[sz/2] *= double(1)/(double(1) + noise/10*std::pow(double(sz/2),2));
							for( size_t i=1;i<sz/2;++i){
								timeslice_hc[i] *= double(1)/(double(1) + noise/10*std::pow(double(i),2));
								timeslice_hc[sz-i] *= double(1)/(double(1) + noise/10*std::pow(double(i),2));
							}      
							*/

							fftw_execute_r2r(plan_hc2r,
									timeslice_hc,
									timeslice_r
									);
							for(unsigned i=0;i<sz;++i){
								outrealfile << timeslice_r[i] << "\t";
							}
							outrealfile << "\n";


						}// end if (m_printLegendreFFTs)

					} // end k-loop
					outfile.close();
					if (m_printLegendreFFTs){
						if (outabsfile.is_open()) {outabsfile.close();}
						if (outargfile.is_open()) {outargfile.close();}
						if (outrealfile.is_open()) {outrealfile.close();}
					}
					avgoutfile << "\n";
				} // end e-loop
				avgoutfile.close();
				// Here we want to print the fft abs and args of the timeslices
				// now we can print our projection files as maps in e and k
				if (m_makeRotorProjections)
				{
					for (unsigned b=0;b<projections_3d.shape()[0];++b)
					{
						std::string filename = m_datadir + std::string("proj_legendre");
						details = "-" + m_str_experiment + "-r" + m_str_runnum;
						details += "_g" + boost::lexical_cast<std::string>(g);
						details += "_l" + boost::lexical_cast<std::string>(l);
						details += "_b" + boost::lexical_cast<std::string>(b);
						details += "_tspan" + tspanStr;
						details += std::string(".dat");
						filename += details;
						std::ofstream projout(filename.c_str(),std::ios::out);
						for (unsigned e = 0; e<projections_3d.shape()[1] ; ++e)
						{
							for ( unsigned k = 0; k< projections_3d.shape()[2];++k)
							{
								projout << projections_3d[b][e][k] << "\t";
							}
							projout << "\n";
						}
						projout.close();
					}
				} // end if m_makeRotorProjections

			}
		}
		std::cerr << "Leaving printSpectraLegendre()" << std::endl;
		return true;

	}
	bool CookieBox_mod::print_patch_results(void)
	{
		std::string filename = m_datadir + std::string("patch_legendre") + "-" + m_str_experiment + "-r" + m_str_runnum + std::string(".dat");
		try{
			m_patch_ys /= (double)m_patch_bins;
			m_patch_legcoeffs /= (double)m_patch_bins;
			std::ofstream outfile(filename.c_str(),std::ios::out);
			outfile << m_patch_ys << m_patch_legcoeffs << std::endl;
			cerr << m_patch_ys << m_patch_legcoeffs << std::endl;
			outfile.close();
		} catch (std::exception & e) {
			std::cerr << "Failed in CookieBox_mod::print_patch_results(std::string & filename) with error: " << e.what() << std::endl;
			return false;
		}
		return true;
	}
	bool CookieBox_mod::read_corr_factors_old(std::string & filename)
	{
		try {
			std::ifstream infile(filename.c_str(),std::ios::in);
			infile >> m_corr_factors_data;
			infile.close();
			std::cout << "Calibration corr_factors:\n" << m_corr_factors_data;
		} catch (std::exception & e) {
			std::cerr << "Failed CookieBox_mod::read_corr_factors() with error: " << e.what() << std::endl;
			return false;
		}
		return true;
	}
	bool CookieBox_mod::read_corr_factors_interp(std::string & filename)
	{
		if (!m_use_e2t) return false;
		try {
			std::ifstream infile(filename.c_str(),std::ios::in);
			infile >> m_corr_factors_data; 
			infile.close();
			std::cout << "Neon calibration correction factors:\n" << m_corr_factors_data << std::endl;

			/*
			std::cerr << "with m_corr_factors_data.size(){.size()} = "
				<< m_corr_factors_data.size() << " , " << m_corr_factors_data[0].size() << std::endl;
			for (unsigned c=0;c<m_corr_factors_data[0].size();++c){
				for (unsigned i = 0; i< m_corr_factors_data.size(); ++i){
					std::cerr << "m_corr_factors_data["<< i << "]["<< c << "] = " << m_corr_factors_data[i][c] << std::endl;
				}
			}
			*/
			m_corr_factors_maps.reserve(m_corr_factors_data[0].size()-1);
			for (unsigned c=1;c<m_corr_factors_data[0].size();++c){
				m_corr_factors_maps.push_back( std::map<double,double>() );
				for (unsigned i = 0; i< m_corr_factors_data.size(); ++i){
					(m_corr_factors_maps[c-1]).insert( std::make_pair( m_corr_factors_data[i][0],m_corr_factors_data[i][c] ) );
				}
			}
		} catch (std::exception & e) {
			std::cerr << "Failed CookieBox_mod::read_corr_factors_interp() with error: " << e.what() << std::endl;
			return false;
		}
		return true;
	}
	bool CookieBox_mod::read_e2t_conversion_old(std::string & filename)
	{
		try {
			std::ifstream infile(filename.c_str(),std::ios::in);
			infile >> m_e2t_data; 
			infile.close();
			std::cout << "Calibration e2t:\n" << m_e2t_data;
			m_use_e2t = true;
		} catch (std::exception & e) {
			std::cerr << "Failed CookieBox_mod::read_e2t_conversion() with error: " << e.what() << std::endl;
			m_use_e2t = false;
		}
		return m_use_e2t;
	}
	double CookieBox_mod::corr_factor_hand(const unsigned c)
	{
		return 1./m_corr_factors_data[c][0];
	}
	double CookieBox_mod::corr_factor_neon(const unsigned k, const unsigned c)
	{
		double x = (double) k * m_sp_stepsvec[idx_kinenergy] + (double) m_sp_startsvec[idx_kinenergy];
		return interpolate<double,double>(m_corr_factors_maps[c],x);
	}

	double CookieBox_mod::corr_factor_hand(const unsigned k, const unsigned c)
	{
		return m_corr_factors_data[c][1];
	}
	double CookieBox_mod::corr_factor_old(const unsigned k, const unsigned c)
	{
		double x = (double) k * m_sp_stepsvec[idx_kinenergy] + m_sp_startsvec[idx_kinenergy];
		x -= m_corr_factors_data[c][0];
		double y = m_corr_factors_data[c][1];
		for (unsigned i=2;i<m_corr_factors_data[c].size();++i){
			y += m_corr_factors_data[c][i] * std::pow(x,(int)i-1);
		}
		return y;
	}
	void CookieBox_mod::e2t_setlims(const unsigned k, const unsigned c,std::vector<unsigned> & lims, std::vector<double> & weights)
	{
		double x1 = (double) (k * m_sp_stepsvec.front()) + m_sp_startsvec.front();
		double x2 = x1 + m_sp_stepsvec.front();
		/* for reference
		std::vector<unsigned> m_sp_binsvec = (configList("sp_bins_list"));
		std::vector<double> m_sp_startsvec = (configList("sp_starts_list"));
		std::vector<double> m_sp_stepsvec = (configList("sp_steps_list"));
		*/
		double y1 = (this->*e2t)(x1,c);
		double y2 = (this->*e2t)(x2,c);
		//if (c == 0)
		//	std::cerr << "e2t lims: " << (x1+x0) << "\t" << y1 << "\t" << (x2+x0) << "\t" << y2 << std::endl;
		lims.front() = (unsigned) std::max((int)y2,0);
		lims.back() = std::min((unsigned) y1,(unsigned)m_data_5d.shape()[4]);
		weights.front() = (double)((unsigned)y2+1) - y2;
		weights.back() = y1 - (double)((unsigned) y1); 
		if (lims.front() > lims.back()){
			std::cerr << "!! warning !!  Weights may be wrong, energy to time gave reversed limits for binning!" << std::endl;
			std::reverse(lims.begin(),lims.end());
			weights.front() = (double)((unsigned)y1+1) - y1;
			weights.back() = y2 - (double)((unsigned) y2); 
		}
		return;
	}
	void CookieBox_mod::e2t_fill_tofs(const unsigned k, std::vector<double> & tofs)
	{
		double x = (double) (k * m_sp_stepsvec.front()) + m_sp_startsvec.front();
		for (unsigned c=0;c<tofs.size();++c){
			tofs[c] = (this->*e2t)(x,c);
		}
		return;
	}
	// make functors out of these //
	double CookieBox_mod::e2t_neon(const double x, const unsigned c)
	{
		double e = (x - m_e2t_data[c][0]);
		if (e<=0) return 0.;
		return std::pow(( m_e2t_data[c][1] / e ),(int)2);
	}
	double CookieBox_mod::e2t_hand(const double x, const unsigned c)
	{
		double x0 = m_e2t_data[c][0];
		double y = m_e2t_data[c][1];
		if (x>x0)
			y += std::sqrt(m_e2t_data[c][2]/(x-x0));
		return y; 
		//t(r,a,x0,e)=sqrt(a/(e-r))+x0
	}
	double CookieBox_mod::e2t_old(const double x, const unsigned c)
	{
		double x0 = m_e2t_data[c][0];
		double y = m_e2t_data[c][1];
		for (unsigned i=2;i<m_e2t_data[c].size();++i){
			y += m_e2t_data[c][i] * std::pow((x-x0),(int)i-1);
		}
		return y;
	}

	bool CookieBox_mod::computeLegendreCoeffs(const unsigned t, const unsigned e, const unsigned g, const unsigned k)
	{
		if (m_rank != m_root_rank){
			std::cerr << "Trying to compute Legendre projections of m_data_5d but using non-root rank, rank = " << m_rank << std::endl;
			return false;
		}

		if (!fillLegendreVecs())
		{
			std::cerr << "Failed to fillLegendreVecs() in bool CookieBox_mod::computeLegendreCoeffs()" << std::endl;
			return false;
		}
		a5d_d_t::index_gen indices;
		a4d_d_t::index_gen avg_indices;
		a5d_d_t::array_view<1>::type slice = m_legendres_5d[indices [t][e][g][k][range()] ]; // slice 1D of the Legendre output
		a4d_d_t::array_view<1>::type avg_slice = m_avg_legendres_4d[indices [e][g][k][range()] ]; // slice 1D of the time average Legendre output
		unsigned nangles = m_data_5d.shape()[chanind]; // here just the number of detectors
		std::vector<double> x(nangles);
		std::vector<double> y(nangles,0.);
		try{
			unsigned shots = m_shots_4d[t][e][g][0];
			if (shots > 0 ){
				for (unsigned c=0;c<nangles; ++c){
					//x[c] = (double)c * (2.0/(double)nangles) - 1.0; // I dont think I actually need this one.
					// for Markus, somehow write this out to file //
					double sum;
					sum = 0.;
					std::vector<unsigned> slims(2,0);
					std::vector<double> weights(2,1.);
					e2t_setlims(k,c,slims,weights);
					if (slims.front() == slims.back()){
						y[c] = (1.-weights.back() + weights.front())*(double)m_data_5d[t][e][g][c][slims.front()];
					} else {
						y[c] = weights.front()*(double)m_data_5d[t][e][g][c][slims.front()];
						for (unsigned s=slims.front()+1;s<slims.back();++s){
							y[c] += (double)m_data_5d[t][e][g][c][s];
						}
						y[c] += weights.back()*(double)m_data_5d[t][e][g][c][slims.back()];
					}
					y[c] /= (double)shots; // HERE HERE HERE HERE Careful, take this back out.
					y[c] *= (this->*corr_factor)(k,c);
					if(k==0 && c == 0)
						++m_avg_legendres_nsums[e][g];
					//if (t==m_shots_4d.shape()[ttind]/2 && e==m_shots_4d.shape()[ebind]/2 && g==m_shots_4d.shape()[gdind]/2 && c==nangles-1)
					if (t==0 && e==m_shots_4d.shape()[ebind]/2 && g==m_shots_4d.shape()[gdind]/2 && c==nangles-1){
						//std::cerr << "\tsampling_e2t_calib.dat" << std::endl;
						std::string samplename("sampling_e2t_calib.dat");
						std::ofstream samplefile(samplename.c_str(),std::ios::app);	
						samplefile << k << "\t" << y;
						std::vector<double> tofs(nangles);
						e2t_fill_tofs(k,tofs);
						samplefile << tofs << std::endl;
						samplefile.close();
					}
				}
				if (m_markusfileChans.is_open()){
					//	Print out the comparison vectors to check the Legendre Porjections... only at g,e,t at 1/2 range
					m_markusfileChans << y << "\n";
				}
				slice[0] = removemean(y,m_sp_projectionMask);
				if (m_markusfileLegs.is_open()) { m_markusfileLegs << slice[0] << "\t";}
				m_totsignal_4d[t][e][g][k] = slice[0] * sum(m_sp_projectionMask);
				avg_slice[0] += slice[0];
				for (unsigned l=1;l<slice.size();++l){
					//slice[l] = std::inner_product(y.begin(),y.end(),m_legendre_vectors[l].begin(),0.);
					slice[l] = projection(y,m_legendre_vectors[l],m_sp_projectionMask);
					avg_slice[l] += slice[l];
					if (m_markusfileLegs.is_open()) { m_markusfileLegs << slice[l] << "\t";}
				}
				if (m_markusfileLegs.is_open()) { m_markusfileLegs << "\n";}
				if (m_patch_compute && inpatch(t,e,g,k)){
					m_patch_bins++;
					m_patch_ys += y;
					//m_patch_ys += slice[0];
					for (unsigned l=0;l<slice.shape()[0];++l){
						m_patch_legcoeffs[l] += slice[l];
					}
				}
			} else {
				for (unsigned l=0;l<slice.size();++l){
					slice[l] = 0.;
				}
			}
		} catch(std::exception &e){
			std::cerr << "Failed in computeLegendreCoeffs() with error: " << e.what() << std::endl;
			return false;
		}
		return true;
	}

	bool CookieBox_mod::fillLegendreVecs(void)
	{
		if (m_legendreVecs_filled == true)
			return m_legendreVecs_filled;
		try {
			std::cout << "Filling Legendre vectors" << std::endl;
			unsigned nsamples = m_data_5d.shape()[3];
			m_legendre_vectors.resize(m_sp_binsvec[idx_legendre]);
			for (unsigned l=0;l<m_legendre_vectors.size();++l){
				m_legendre_vectors[l].resize(nsamples);
			}
			for (unsigned i=0;i<nsamples;++i){
				double x= std::cos(2*M_PI*(double)i/ (double) nsamples);
				m_legendre_vectors[0][i] = boost::math::legendre_p(0,x);
				m_legendre_vectors[1][i] = boost::math::legendre_p(1,x);
				for (unsigned l=2;l<m_legendre_vectors.size();++l){
					m_legendre_vectors[l][i] = boost::math::legendre_next(l-1, x, 
							m_legendre_vectors[l-1][i], 
							m_legendre_vectors[l-2][i]);
				}
			}
			for (unsigned l=0;l<m_legendre_vectors.size();++l){ 
				sqr_normalize(m_legendre_vectors[l],m_sp_projectionMask);
			}
			std::string filename = m_datadir + "legendres-";
			filename += m_str_experiment + "-r" + m_str_runnum;
			filename += std::string(".out");
			std::ofstream outfile(filename.c_str(),std::ios::out); 
			outfile << "#legendres of increasing l" << std::endl;
			for (unsigned i=0;i<m_legendre_vectors.front().size();++i){
				for (unsigned l=0;l<m_legendre_vectors.size();++l){
					outfile << m_legendre_vectors[l][i] << "\t";
				}
				outfile << "\n";
			}
			outfile.close();
			m_legendreVecs_filled = true;
		} catch (std::exception &e) {
			std::cerr << "Failed attempt to fill vectors in CookieBox_mod::fillLegendreVecs(void) with error: " << e.what() << std::endl;
			m_legendreVecs_filled = false;
		}
		return m_legendreVecs_filled;
	}
	bool CookieBox_mod::printSpectra(void)
	{

		// index order is [0=TimeTool][1=ebeam][2=gasdet][3=chan][4=sample]
		// enum DimDefs {ttind,ebind,gdind,chanind,tofind};
		if (m_rank != m_root_rank){
			std::cerr << "Trying to print out m_data_5d but using non-root rank, rank = " << m_rank << std::endl;
			return false;
		}
		//
		// resize the avg spectra for accumulation //
		m_avgSpectra_shots.resize(boost::extents
				[m_data_5d.shape()[1]]  // ebeam bins
				[m_data_5d.shape()[2]]  // gasdet bins
				[m_data_5d.shape()[3]]  // channels
				);
		m_avgSpectra.resize(boost::extents
				[m_data_5d.shape()[1]]  // ebeam bins
				[m_data_5d.shape()[2]]  // gasdet bins
				[m_data_5d.shape()[3]]  // channels
				[m_data_5d.shape()[4]]	// samples
				);
		for(long long * it = m_avgSpectra.origin(); it < m_avgSpectra.origin() + m_avgSpectra.num_elements() ; ++it)
		{
			*it = 0;
		}
		for(long long * it = m_avgSpectra_shots.origin(); it < m_avgSpectra_shots.origin() + m_avgSpectra_shots.num_elements() ; ++it)
		{
			*it = 0;
		}
		std::cout << "skipping // print hist spectra // only computing average" << std::endl;
		for (unsigned e = 0; e<m_data_5d.shape()[1] ; ++e)
		{
			for ( unsigned g = 0; g < m_data_5d.shape()[2] ; ++g)
			{
				for ( unsigned c = 0; c < m_data_5d.shape()[3] ; ++c){
					/*
					   std::string filename = m_datadir + std::string("hist");
					   filename+= "-" + m_str_experiment + "-r" + m_str_runnum;
					   filename += "_e" + boost::lexical_cast<std::string>(e);
					   filename += "_g" + boost::lexical_cast<std::string>(g);
					   filename += "_c" + boost::lexical_cast<std::string>(c);
					   filename += std::string(".dat");
					   */
					//std::ofstream outfile(filename.c_str(),std::ios::out);
					for (unsigned sample = 0; sample < m_data_5d.shape()[4] ; ++sample){
						for (unsigned t = 0; t < m_data_5d.shape()[0] ; ++t){
							long long shots=0;
							double val = 0.;
							shots = m_shots_4d[t][e][g][c];
							if ( shots > 0){
								val = (double) m_data_5d[t][e][g][c][sample];
								if(sample==0)
									m_avgSpectra_shots[e][g][c] += shots; // integrate on timebins
								m_avgSpectra[e][g][c][sample] += m_data_5d[t][e][g][c][sample];
								val /= shots;
							}
							//outfile << val << "\t";
						}
						//outfile << "\n";
					}
					//outfile.close();
				}
			}
		}

		std::cerr << "Made it here printing" << std::endl;
		std::cout << "skipping // print difference spectra //" << std::endl;
		/*
		   for (unsigned e = 0; e<m_data_5d.shape()[1] ; ++e)
		   {
		   for ( unsigned g = 0; g < m_data_5d.shape()[2] ; ++g)
		   {
		   for ( unsigned c = 0; c < m_data_5d.shape()[3] ; ++c){
		   std::string filename = m_datadir + std::string("diffhist");
		   filename+= "-" + m_str_experiment + "-r" + m_str_runnum;
		   filename += "_e" + boost::lexical_cast<std::string>(e);
		   filename += "_g" + boost::lexical_cast<std::string>(g);
		   filename += "_c" + boost::lexical_cast<std::string>(c);
		   filename += std::string(".dat");
		   std::ofstream outfile(filename.c_str(),std::ios::out);
		   for (unsigned sample = 0; sample < m_data_5d.shape()[4] ; ++sample){
		   for (unsigned t = 0; t < m_data_5d.shape()[0] ; ++t){
		   long long shots=0;
		   double result = 0.;
		   shots = m_shots_4d[t][e][g][c];
		   if ( shots > 0){
		   result = (double)m_data_5d[t][e][g][c][sample]/(double)shots;
		   result -= (double)m_avgSpectra[e][g][c][sample]/(double)m_avgSpectra_shots[e][g][c];
		   }
		   outfile << result << "\t";
		   }
		   outfile << "\n";
		   }
		   outfile.close();
		   }
		   }
		   }
		   */



		return true;
	}
	// Here is where we will project the time-sorted specra onto the othonormal basis from a particular stage reading of the TransAbs output.
	bool CookieBox_mod::fftDiffSpectra(void) 
	{
		unsigned sz = m_data_5d.shape()[0];
		unsigned ncomponents = sz/64;
		double * timeslice_r = (double *) fftw_malloc(sizeof(double) * sz);
		double * timeslice_hc = (double *) fftw_malloc(sizeof(double) * sz);

		fftw_plan plan_r2hc = fftw_plan_r2r_1d(sz,
				timeslice_r,
				timeslice_hc, 
				FFTW_R2HC,
				FFTW_MEASURE
				);

		//for (unsigned e = 0; e<m_data_5d.shape()[1] ; ++e)
		{
			unsigned e =  5;
			//for ( unsigned g = 0; g < m_data_5d.shape()[2] ; ++g)
			{
				unsigned g = 1;
				for ( unsigned c = 0; c < m_data_5d.shape()[3] ; ++c){
					std::string filename = m_datadir + std::string("fouriers");
					filename+= "-" + m_str_experiment + "-r" + m_str_runnum;
					filename += "_e" + boost::lexical_cast<std::string>(e);
					filename += "_g" + boost::lexical_cast<std::string>(g);
					filename += "_c" + boost::lexical_cast<std::string>(c);
					filename += std::string(".dat");
					std::ofstream outfile(filename.c_str(),std::ios::out);
					outfile << "#fft real(tab)imag ";
					outfile << "#fft meanspect\tresidualmean\treal(tab)imag ";
					for(unsigned f=0;f<ncomponents;++f){
						outfile << f << "re\t" << f << "im\t";
					}
					outfile << "\n";
					for (unsigned sample = 0; sample < m_data_5d.shape()[4] ; ++sample){
						double mean = 0.;
						unsigned meancontributions = 0;
						for (unsigned t = 0; t < sz ; ++t){
							long long shots=0;
							double result = 0.;
							shots = m_shots_4d[t][e][g][c];
							if ( shots > 0){
								result = (double)m_data_5d[t][e][g][c][sample]/(double)shots;
								mean += result;
								++meancontributions;
								result -= (double)m_avgSpectra[e][g][c][sample]/(double)m_avgSpectra_shots[e][g][c];
							}
							timeslice_r[t] = result;
						}
						mean /= meancontributions;
						outfile << mean << "\t";
						mean = removemean(timeslice_r,sz);
						outfile << mean << "\t";
						fftw_execute_r2r(plan_r2hc,
								timeslice_r,
								timeslice_hc
								);
						outfile << timeslice_hc[0];
						for(unsigned f=1;f<ncomponents;++f){
							std::complex<double> z(timeslice_hc[f],timeslice_hc[sz-f]);
							double abs = std::abs(z);
							double arg = std::arg(z);
							outfile << abs << "\t" << arg << "\t";
						}
						outfile << "\n";
					}
					outfile.close();
				}
			}
		}
		if (timeslice_r != NULL){fftw_free(timeslice_r);timeslice_r = NULL;}
		if (timeslice_hc != NULL){fftw_free(timeslice_hc);timeslice_hc = NULL;}
		if (plan_r2hc != NULL)
			fftw_destroy_plan(plan_r2hc);
		return true;
	}




	bool CookieBox_mod::projectDiffSpectra(const std::string orthofilename) 
	{
		std::vector<std::string> header;
		std::vector < std::vector <double> > orthobases;
		std::ifstream orthofile(orthofilename.c_str(),std::ios::in);
		try {
			std::cerr << "Pulling orthofile " << orthofilename << std::endl;
			orthofile >> orthobases;
		} catch (std::exception &e) {
			std::cerr << "Failed to read in files for projectDiffSpectra(), error: " << e.what() << std::endl;
			return false;
		}
		orthofile.close();
		size_t nsamples = m_legendres_5d.shape()[0];
		std::cout << "size_t nsamples = m_legendres_5d.shape()[0]; \t " << nsamples << std::endl;
		std::vector<double> timeslice(nsamples,0.);

		if (orthobases.front().size() > timeslice.size()){
			std::cout << "condensing orthobases from " << orthobases.front().size() << " to ";
			condense(orthobases,timeslice.size());
			std::cout << orthobases.front().size() << " after condense(orthobases,timeslice.size()) call" << std::endl;
		}

		std::vector<double *> projections_r(orthobases.size()+1,(double*)NULL);   // add one to accept the mean as well.
		std::vector<double *> projections_hc(orthobases.size()+1,(double*)NULL);   // add one to accept the mean as well.// later add another for back projected
		//unsigned nsamples = m_data_5d.shape()[4];
		for (unsigned i=0;i<projections_r.size();++i){
			projections_r[i] = (double*)fftw_malloc(sizeof(double) * nsamples);
			projections_hc[i] = (double*)fftw_malloc(sizeof(double) * nsamples);
		}
		fftw_plan plan_r2hc = fftw_plan_r2r_1d(nsamples,
				projections_r.front(),
				projections_hc.front(), 
				FFTW_R2HC,
				FFTW_MEASURE
				);
		fftw_plan plan_hc2r = fftw_plan_r2r_1d(nsamples,
				projections_hc.front(), 
				projections_r.front(),
				FFTW_HC2R,
				FFTW_MEASURE
				);

		/*
		 * This needs to become for Legendres rather than individual channels
		 */
		std::string tspanStr;
		double tspanDbl = m_tt.getTspan(tspanStr);
		for ( unsigned l=0;l<m_legendres_5d.shape()[4];++l)
		{
			std::cout << "Writing projections for l = " << l << "\t... " << std::flush;
			for (unsigned e = m_data_5d.shape()[ebind]/4; e<m_data_5d.shape()[ebind]*3/4+1 ; ++e)
			{
				//unsigned e = 5;
				//for ( unsigned g = 0; g < m_data_5d.shape()[gdind] ; ++g)
				unsigned g = 1;
				{
					/*
					 * HERE HERE HERE HERE
					 */
					std::string filename = m_datadir + std::string("projections");
					std::string details = "-" + m_str_experiment + "-r" + m_str_runnum;
					details += "_e" + boost::lexical_cast<std::string>(e);
					details += "_g" + boost::lexical_cast<std::string>(g);
					details += "_l" + boost::lexical_cast<std::string>(l);
					details += "_tspan" + tspanStr;
					details += std::string(".dat");
					filename += details;
					//std::cerr << "printing file " << filename << std::endl;
					std::ofstream outfile(filename.c_str(),std::ios::out);
					outfile << "#rotor projections for tspan = " << tspanDbl << "\n";
					outfile << "#mean\t";
					//std::cerr << "#mean\tprojections " << "\t header written to file" << std::endl;
					for(unsigned b=0;b<orthobases.size();++b){
						outfile << b << "\t";
					}
					outfile << "\n";
					size_t numenergies = m_legendres_5d.shape()[3];
					for (size_t k = 0; k < numenergies ; ++k){
						//std::cerr << "HERE HERE HERE HERE\tk = " << k << "\t" << std::flush;
						double mean;
						unsigned meancontributions = 0;
						for (unsigned t = 0; t < nsamples ; ++t){
							mean = 0.;
							meancontributions = 0;
							long long shots=0;
							double result = 0.;
							shots = m_shots_4d[t][e][g][0];
							if ( shots > 0){
								//result = (double)m_legendres_5d[t][e][g][k][l]/(double)shots;
								result = m_legendres_5d[t][e][g][k][l] ; // slice 1D of the Legendre output
								//std::cerr << "result = m_legendres_5d[t][e][g][k][l] = " << result << std::endl;
								mean += result;
								++meancontributions;
								//result -= (double)m_avgSpectra[e][g][c][sample]/(double)m_avgSpectra_shots[e][g][c];
							}
							timeslice[t] = result;
						}
						removemean(timeslice); // don't change this... I'm doing a conditional mean removal
						mean /= meancontributions;
						outfile << mean << "\t" ;
						projections_r[0][k] = mean;
						for(unsigned b=0;b<orthobases.size();++b){
							double proj = projection(orthobases[b],timeslice);
							outfile << proj << "\t";
							projections_r[b+1][k] = proj; // filling for FFT
						}
						outfile << "\n";
						//std::cerr << "EXIT EXIT \tk = " << k << std::endl;
					}
					outfile << "\n";

					/*


					// HERE, I want the FFT of every projection, but for now only the power spectrum.
					for (unsigned i=0;i<projections_r.size();++i){
					fftw_execute_r2r(plan_r2hc,
					projections_r[i],
					projections_hc[i]
					);
					}
					// OK now the printing is silly for sake of gnuplot //
					filename += ".fft";
					std::ofstream fftoutfile(filename.c_str(),std::ios::out);
					fftoutfile << "#projections ffts ";
					//std::cerr << "printing file " << filename << std::endl;
					unsigned sample = 0;
					for (unsigned i=0;i<projections_r.size();++i){
					fftoutfile << std::norm(std::complex<double>(projections_hc[i][sample],0.));
					fftoutfile << "\t";
					}
					fftoutfile << "\n";
					for (sample=1;sample<nsamples/2;++sample){
					for (unsigned i=0;i<projections_r.size();++i){
					fftoutfile << std::norm(
					std::complex<double>(
					projections_hc[i][sample]
					,projections_hc[i][nsamples - sample]) );
					fftoutfile << "\t";
					}
					fftoutfile << "\n";
					}
					sample = nsamples/2;
					for (unsigned i=0;i<projections_r.size();++i){
					fftoutfile << std::norm(std::complex<double>(projections_hc[i][sample],0.));
					fftoutfile << "\t";
					}
					fftoutfile.close();
					// OK, now filter and back project and write out //
					// HERE.  As per conversation with nick, it may make sense to filter in the revival basis set.
					// Compute the bases to higher order... like 50 or so...
					// Then read those in and calculate the power spectrum of each overlap integral with the different pases
					// Then apply the weiner 1/(1+n/s) in this new basis also and back transform
					weinerfilter(projections_hc,nsamples);
					for (unsigned i=0;i<projections_hc.size();++i){
					fftw_execute_r2r(plan_hc2r,
					projections_hc[i],
					projections_r[i]
					);
					double scale(1./(double)nsamples);
					std::transform(projections_r[i], projections_r[i] + nsamples, projections_r[i], 
					std::bind1st(std::multiplies<double>(),scale));
					}
					// HERE add the back-projected of the weiner filtered coefficients and see if it matches the diffhist files.
					filename += ".filtered";
					//std::cerr << "printing file " << filename << std::endl;
					std::ofstream backfftoutfile(filename.c_str(),std::ios::out);
					fftoutfile << "#wiener filtered projections back ffts ";
					for (unsigned sample=0;sample<nsamples;++sample){
					for (unsigned i=0; i<projections_r.size();++i){
					backfftoutfile << projections_r[i][sample] << "\t";
					}
					backfftoutfile << "\n";
					}
					backfftoutfile.close();
					*/
					std::cout << "... written.\n" << std::flush;
					outfile.close();
				}
			}
		}
		std::cerr << "Now killing off fftw vectors" << std::endl;
		for (unsigned i=0;i<projections_r.size();++i){
			if (projections_r[i] != NULL){fftw_free(projections_r[i]);projections_r[i] = NULL;}
			if (projections_hc[i] != NULL){fftw_free(projections_hc[i]);projections_hc[i] = NULL;}
		}
		std::cerr << "Now killing off fftw plans" << std::endl;
		if (plan_r2hc != NULL)
			fftw_destroy_plan(plan_r2hc);
		if (plan_hc2r != NULL)
			fftw_destroy_plan(plan_hc2r);
		std::cerr << "Leaving the printing " << std::endl;
		return true;
	}

	void CookieBox_mod::weinerfilter(std::vector<double *> & in,const unsigned sz)
	{
		// get the limits, mean the last few in power spectrum and...
		// For now, just give it a static filter... later we can do like in TransAbs. 
		//std::vector<double> filter(sz/2+1);
		std::vector< std::vector<double> > filter(in.size());
		std::vector< std::vector<double> > powerspectra(in.size());
		std::vector<double> noiselevel(in.size(),0.);
		std::vector<double> noisestd(in.size(),0.);
		std::vector<double> signallevel(in.size(),0.);
		for (unsigned i=0;i<in.size();++i){
			filter[i].resize(sz/2+1,0.);
			powerspectra[i].resize(sz/2+1,0.);
			powerspectra[i][0] = std::norm(std::complex<double>(in[i][0],0.));
			powerspectra[i][sz/2] = std::norm(std::complex<double>(in[i][sz/2],0.));
			for (unsigned j=1;j<sz/2;++j){
				powerspectra[i][j] = std::norm(std::complex<double>(in[i][j],in[i][sz-j]));
			}
			unsigned nelems = 64;
			unsigned selems = 4;
			std::vector<double> tempvec;
			tempvec.assign(powerspectra[i].end()-nelems,powerspectra[i].end());
			meanstdlog(tempvec, noiselevel[i], noisestd[i]);
			tempvec.clear();
			tempvec.assign(powerspectra[i].begin()+1,powerspectra[i].begin()+1+selems);
				// So let's do this to get where this should cross the x=1 axis in log scale
				signallevel[i] = 0.;
				for (unsigned j=1;j<selems+1;++j){
					signallevel[i] += std::log(powerspectra[i][j] * double(j));
				}
				signallevel[i] /= selems;
				// otherwise, fit the first 10 points with 1/f and then 
				// if you multiply by x, then you get something that we can average for a 1/x crossing value
				filter[i][0] = 1.;// /(1.+exp(noiselevel[i]-signallevel[i])); // remember noise and signal levels are in log space
				for (unsigned j=1;j<sz/2+1;++j){
					if (i!=0){
						filter[i][j] = 1./(1.+ std::exp(std::log(double(j)) - std::log(10.)));
					} else {
						double logsignal = std::log( powerspectra[i][j] );
						filter[i][j] = 1.;
						filter[i][j] /= (1.+ std::exp(std::min(
										(std::log(double(j)) - std::log(10.)) ,
										(noiselevel[i] - logsignal)
										)
									) 
								);
					}
				}

			in[i][0] *= filter[i][0];
			for (unsigned j = 1;j<sz/2;++j){
				in[i][j] *= filter[i][j];
				in[i][sz-j] *= filter[i][j];
			}
			in[i][sz/2] *= filter[i][sz/2];
		}
	}

	bool CookieBox_mod::evtput(Event& evt, Env& env, const unsigned chan)
	{
		if (m_rank != m_root_rank){
			std::cerr << "Trying to print out evtput using non-root rank, rank = " << m_rank << std::endl;
			return false;
		}

		long long int sum,shots;
		long double gsum;
		std::string outkey("aq_");
		outkey += boost::lexical_cast<std::string>(chan);

		// index order is [0=ebeam][1=gasdet][2=chan][3=sample]
		// enum DimDefs {ttind,ebind,gdind,chanind,tofind};
		unsigned ngdetbins,nchans;
		ngdetbins = m_aqSum_4d.shape()[1];
		nchans = m_aqSum_4d.shape()[2];
		assert(nchans>chan);
		unsigned shape[] = {m_aqSum_4d.shape()[0], m_aqSum_4d.shape()[3]};
		boost::shared_ptr< ndarray<double,2> > outmatPtr = boost::make_shared < ndarray<double,2> > ( shape );

		for (unsigned sample = 0; sample < m_aqSum_4d.shape()[3];++sample){ 
			for (unsigned e=0;e<m_aqSum_4d.shape()[0];++e){
				gsum = 0.;
				shots = 0;
				sum = 0;
				for (unsigned g=0;g<m_gd.bins(Gdet::average);++g){
					if (m_gdetShots_2d[e][g]>0 && m_gdetSum_2d[e][g] > 0){
						sum += m_aqSum_4d[e][g][chan][sample];
						gsum += m_gdetSum_2d[e][g];
						shots += m_gdetShots_2d[e][g]; 
					}
				}
				if (shots>0){
					(*outmatPtr)[e][sample] = (double)sum/(double)shots;
				} else {
					(*outmatPtr)[e][sample] = 0.0;
				}
			}
		}
		evt.put(outmatPtr,outkey);
		return true;
	}

} // namespace CookieBox_pkg
