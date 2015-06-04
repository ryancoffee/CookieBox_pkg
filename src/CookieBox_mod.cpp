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
	{
		// get the values from configuration or use defaults
		m_print_every = config("print_every",10);
		std::vector<bool> use_logic = configList("use_logic");
		std::vector<bool> print_logic = configList("print_logic");
		std::vector<bool>::iterator logic_itr;
		std::vector<bool>::iterator print_itr;

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
			if ( m_xt.remove_dark( config("xt_remove_dark",false)) 
				&& !m_xt.dark_file(configStr("xt_dark_file")) )
			{
				std::cerr << "Failed to load the dark image file in CookieBox_mod::CookieBox_mod() constructor" << std::endl;
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
		}

		logic_itr = use_logic.begin() + use_tt;
		print_itr = print_logic.begin() + print_tt;
		if (m_tt.use(*logic_itr)){
			std::vector<unsigned> tt_bins_in = (configList("tt_bins_list"));
			std::vector<double> tt_mins_in = (configList("tt_mins_list"));
			std::vector<double> tt_maxs_in = (configList("tt_maxs_list"));
			m_tt.init(tt_bins_in, tt_mins_in, tt_maxs_in);
			m_tt.srcStr(configSrc("tt_source", "DetInfo(:Opal1k)"));
			m_tt.print(*print_itr);
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
				totalchannels += m_aq.at(i).nchannels();
				unsigned s = m_aq.at(i).nsamples();
				if ( samples < s)
					samples = s;
			}
			std::vector < std::vector<unsigned> > integlims;
			integlims.clear();
			integlims.push_back( (configList("aq_integlims_low")) );
			integlims.push_back( (configList("aq_integlims_high")) );
			if ( integlims.front().size()>0 && integlims.back().size()>0){
				for (unsigned i=0;i<m_aq.size();++i){
					m_aq[i].integlims(integlims);
				}
			}

			std::cerr << "Exiting the locic_itr for initializing vector<m_aq>" << std::endl;
			//m_data_4d.resize(boost::extents[m_eb.bins(Ebeam::energy)][m_gd.bins(Gdet::average)][totalchannels][samples]);
		}
		

		logic_itr = use_logic.begin() + use_learn;
		print_itr = print_logic.begin() + print_learn;
		if (m_learn.use(*logic_itr)){ // eventually, make these switchable between library building and using
			m_learn.init((unsigned)config("learn_pca_samples_limit",25)
				, (unsigned)config("learn_pca_features",10)
				, (float)config("learn_pca_variance", 0.9)
				);
			m_learn.print(*print_itr);
		}

		//m_gdetSum_2d.resize(boost::extents[m_eb.bins(Ebeam::energy)][m_gd.bins(Gdet::average)]);
		//m_gdetShots_2d.resize(boost::extents[m_eb.bins(Ebeam::energy)][m_gd.bins(Gdet::average)]);

		m_skip_events = config("skip_events",0);
		m_last_event = config("last_event",100);
		m_gnuplotting = config("gnuplotting",false);
		m_rank = 0;//MPI::COMM_WORLD.Get_rank();
		m_mpi_size = 1;//MPI::COMM_WORLD.Get_size();
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
			return false;
		}
		return (m_count_event.at(m_rank) > m_skip_events);
	}


	/// Method which is called once at the beginning of the job
	void CookieBox_mod::beginJob(Event& evt, Env& env)
	{
		//m_rank = MPI::COMM_WORLD.Get_rank();
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
		//m_rank = MPI::COMM_WORLD.Get_rank();
		std::cout << "Started My Run for rank " << m_rank << std::endl;

		std::stringstream ranktail;
		ranktail << ".rank" << m_rank;

		m_str_runnum = ImgAlgos::stringRunNumber(evt);
		m_str_experiment = ImgAlgos::stringExperiment(env);

		std::cout << "filenames = ";
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
		//m_rank = MPI::COMM_WORLD.Get_rank();
		std::cout << "Started my CalibCycle in rank " << m_rank << std::endl;
	}
	//
	/// Method which is called with event data, this is the only required 
	/// method, all other methods are optional
	void CookieBox_mod::event(Event& evt, Env& env)
	{
		//m_rank = MPI::COMM_WORLD.Get_rank();
		++ m_count_event.at(m_rank);
		//std::cout << "in event number " << m_count_event.at(m_rank) << " for rank " << m_rank << std::endl;

		if (!isApprovedByCounters()) {
			std::cout << "skipping\t" << std::endl;
			return; // if out of skip and limit bounds, return from event
		}
		if (!processEvent(evt,env)) {
			//std::cout << "processing event failed" << std::endl;
			++ m_failed_event.at(m_rank);
			return; // if processing fails, return from event
		}
		//std::cout << " finished event " << m_count_event.at(m_rank) << " for rank " << m_rank << std::endl;

		return;
	}

	

	bool CookieBox_mod::sampleevery(unsigned in)
	{
		return (m_count_event.at(m_rank)%in == 0);
	}

	bool CookieBox_mod::processEvent(Event& evt,Env& env)
	{
		// use evt.put(); to add variable for access later in douwnstream events //
		// std::cout << "entered CookieBox_mod::processEvent() method" << std::endl;
		
		if( m_gd.use() && !m_gd.fill(evt,env) ) {
			std::cerr << "Failed the Gdet" << std::endl;
			return false;
		}

		if ( m_tt.use() && !m_tt.fill(evt,env) ) { 
			std::cerr << "Failed the TimeTool" << std::endl;
			return false; 
		}

		if ( m_eb.use() && !m_eb.fill(evt,env) ) { 
			std::cerr << "Failed the Ebeam" << std::endl;
			return false; 
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
			if ( m_aq.at(i).use() && !m_aq.at(i).fill(evt,env)) { 
				std::cerr << "Failed the Acqiris, skipping the print and rest of event." << std::endl;
				return false;
			}
		}

		if ( sampleevery(m_print_every)  ){
			if (m_xt.print() ){
				if ( !(	m_xt.print_out( m_count_event.at(m_rank) ) 
					&& m_xt.print_cropped(m_count_event.at(m_rank)) 
					) )
					std::cerr << "Failed to print Xtcav" << std::endl;
				if ( !m_xt.print_out_vh( m_count_event.at(m_rank) ))
					std::cerr << "Failed to print Xtcav valhist" << std::endl;
				if ( !m_xt.print_out_stats( m_count_event.at(m_rank) ))
					std::cerr << "Failed to print Xtcav stats" << std::endl;
			}
			if (m_eb.print() && !m_eb.print_out(m_count_event.at(m_rank)))
				std::cerr << "Failed to print Ebeam" << std::endl;
			if (m_gd.print() && !m_gd.print_out(m_count_event.at(m_rank)))
				std::cerr << "Failed to print Gdet" << std::endl;
			if (m_tt.print() && !m_tt.print_out(m_count_event.at(m_rank)))
				std::cerr << "Failed to print TimeTool" << std::endl;
			for (unsigned i=0;i< m_aq.size();++i){
				if (m_aq.at(i).print() && !m_aq.at(i).print_out( m_count_event.at(m_rank) ))
					std::cerr << "Failed to print Acqiris" << std::endl;
			}

		}

		if (m_learn.use()) {
			std::vector<float> buildfeatures;
			std::string features_labels("#");
			/*
			if (m_eb.use())
				m_eb.addfeatures(buildfeatures);
				if (m_learn.nfeatures()==0)
					features_labels += m_eb.get_features_labels();
			if (m_gd.use())
				m_gd.addfeatures(buildfeatures);
				if (m_learn.nfeatures()==0)
					features_labels += m_gd.get_features_labels();
			*/
			if (m_xt.use()){
				//m_xt.addfeatures(buildfeatures);
				m_xt.addfeatures_centroids_image(buildfeatures);
				if (m_learn.nfeatures()==0)
					features_labels += m_xt.get_features_labels();
			}
			//std::cout << "number of features presented to m_learn.fill(buildfeatures) is " << buildfeatures.size() << std::endl;
			if (m_learn.nfeatures()==0)
				m_learn.set_features_labels(features_labels);
			if (m_learn.fill(buildfeatures)){
				if (m_count_event.at(m_rank) %m_print_every < 2)
					std::cout << "filling features at event " << m_count_event.at(m_rank) << std::endl;
			} else {
				if (!m_learn.run_pca((unsigned)config("learn_pca_features",100)) && m_count_event.at(m_rank) %m_print_every < 2){
					std::cerr << "something failed in m_pca.run_pca() call" << std::endl;
				}
				record_f projectedfeatures;
				if (m_learn.print()
				//	&& m_learn.project_pca(buildfeatures,projectedfeatures)
					&& m_count_event.at(m_rank) %m_print_every < 2)
				{
						// HERE HERE HERE HERE // 
						// this seem sto work //
						// surprisingly //
					projectedfeatures = m_learn.project_pca(buildfeatures); 
					std::string filename(m_learn.filename());
					filename += ".projected_event." + boost::lexical_cast<std::string>(m_count_event.at(m_rank));
					std::ofstream outfile (filename.c_str(),std::ios::out);
					for (unsigned i=0;i<projectedfeatures.size();++i){
						outfile << projectedfeatures[i] << "\n";
					}
					outfile.close();
				}
				// HERE consider back projecting for comparisons //

				if (m_count_event.at(m_rank) %m_print_every < 2){
					m_testpca_features.push_back(buildfeatures);
					std::cout << "pushing test features at event " << m_count_event.at(m_rank) << std::endl;
				}
			}
		} // endif m_learn.use()

		if (m_count_event.at(m_rank) %m_print_every < 2){
			std::cout << "returning from processing event " << m_count_event.at(m_rank) << std::endl;
		}

		return true; // returning early from the code;

		// HERE //
		// trying to get the slice to work //

		// std::cout << "Made it to just before slice" << std::endl;

		/*
		unsigned startchanind;
		startchanind = 0;
		for (unsigned i=0;i< m_aq.size();++i){
			if ( m_aq.at(i).use() ) { 
				a4d_ll_t::index_gen indices;
				a4d_ll_t::array_view<2>::type slice = 
					m_data_4d[ 
						indices	[m_eb.index(Ebeam::energy)]
							[m_gd.index(Gdet::average)]
							[range(startchanind,startchanind + m_aq.at(i).nchannels())]
							[range()]
						];
				if ( !m_aq.at(i).fill( evt, env, slice ) ) { 
					std::cerr << "Failed the Acqiris" << std::endl;
					return false; 
				}
				startchanind += m_aq.at(i).nchannels();
			}
		}
		// std::cout << "Made it to just after slice at event" << std::endl;
		return true;
		*/
	}


	/// Method which is called at the end of the calibration cycle
	void CookieBox_mod::endCalibCycle(Event& evt, Env& env)
	{
	}

	/// Method which is called at the end of the run
	void CookieBox_mod::endRun(Event& evt, Env& env)
	{
		time(&m_endruntime);
		float dtime = difftime(m_endruntime,m_beginruntime);
		std::cout << "processing events took " << dtime << " seconds" << std::endl;
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

		//m_rank = MPI::COMM_WORLD.Get_rank();

		
		/*
		if (m_rank == m_root_rank){
			MPI::COMM_WORLD.Reduce(MPI::IN_PLACE,m_aqSum_4d.origin(),m_aqSum_4d.num_elements(),MPI::LONG_LONG,MPI::SUM,m_root_rank); 
			MPI::COMM_WORLD.Reduce(MPI::IN_PLACE,m_gdetSum_2d.origin(),m_gdetSum_2d.num_elements(),MPI::LONG_DOUBLE,MPI::SUM,m_root_rank); 
			MPI::COMM_WORLD.Reduce(MPI::IN_PLACE,m_gdetShots_2d.origin(),m_gdetShots_2d.num_elements(),MPI::LONG_LONG,MPI::SUM,m_root_rank); 
			MPI::COMM_WORLD.Reduce(MPI::IN_PLACE,m_aqSum_baseline.data(),m_aqSum_baseline.size(),MPI::LONG_LONG,MPI::SUM,m_root_rank); 
			MPI::COMM_WORLD.Reduce(MPI::IN_PLACE,m_aqSum_baseline_nvals.data(),m_aqSum_baseline_nvals.size(),MPI::LONG_LONG,MPI::SUM,m_root_rank); 
		 } else {
			MPI::COMM_WORLD.Reduce(m_aqSum_4d.origin(),m_aqSum_4d.origin(),m_aqSum_4d.num_elements(),MPI::LONG_LONG,MPI::SUM,m_root_rank); 
			MPI::COMM_WORLD.Reduce(m_gdetSum_2d.origin(),m_gdetSum_2d.origin(),m_gdetSum_2d.num_elements(),MPI::LONG_DOUBLE,MPI::SUM,m_root_rank); 
			MPI::COMM_WORLD.Reduce(m_gdetShots_2d.origin(),m_gdetShots_2d.origin(),m_gdetShots_2d.num_elements(),MPI::LONG_LONG,MPI::SUM,m_root_rank); 
			MPI::COMM_WORLD.Reduce(m_aqSum_baseline.data(),m_aqSum_baseline.data(),m_aqSum_baseline.size(),MPI::LONG_LONG,MPI::SUM,m_root_rank); 
			MPI::COMM_WORLD.Reduce(m_aqSum_baseline_nvals.data(),m_aqSum_baseline_nvals.data(),m_aqSum_baseline_nvals.size(),MPI::LONG_LONG,MPI::SUM,m_root_rank); 
			return;
		}
		*/


		m_learn.run_pca((unsigned)config("learn_pca_features",100));
		//m_learn.run_pca((float)config("learn_pca_variance",0.95));
		m_learn.inspect_eigen_pca();
		m_learn.print_eigenfaces( m_xt.eigenface_shape() );
		// HERE
		// m_learn.compare_pca_images(m_testpca_features,m_xt.eigenface_shape(),m_xt.nonimage_features());
		m_learn.compare_pca_nonimage(m_testpca_features,m_xt.nonimage_features());
		//m_learn.compare_pca(m_testpca_features);
		
		// HERE 
		// do a comparison for this by constructing the pca after some number of shots, then compare for new shots as they come in.
		// When two or more features are correlated, then we should plot the data in ortho and para for those components
		// Let's also plot up the raw features as points and as parallel coords (with mean subtracted, and scaled by the variance in that dim)
		// Another thing to try is to flatten machine features that have no correlation.... 
		// check the gas dets since they are nearly perfectly correlated to make sure we're doing the right thing.

		std::cout << "returning from endRun() in rank " << m_rank << std::endl;
		for (unsigned i=0;i<m_failed_event.size();++i)
			std::cout << "failed events/total = " << m_failed_event.at(i) << " / " << m_count_event.at(i) << std::endl;

		// technically we should close our files here //
		if (m_gd.print())
			m_gd.close_file();
		// etc... except that in the object destructor I test and close the files.

		return;
	}

	/// Method which is called once at the end of the job
	void CookieBox_mod::endJob(Event& evt, Env& env)
	{
		//m_rank = MPI::COMM_WORLD.Get_rank();
		std::cout << "endJob() called for rank " << m_rank << std::endl;
	}

	// HERE //
	// use this evtput as a model to design the slice and integration methods for the ndarrays //
	void CookieBox_mod::evtput(Event& evt, Env& env, const unsigned chan)
	{
		long long int sum,shots;
		long double gsum;
		std::string outkey("aq_");
		outkey += boost::lexical_cast<std::string>(chan);

		// index order is [0=ebeam][1=gasdet][2=chan][3=sample]
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
					(*outmatPtr)[e][sample] *= -1.; // reverses polarity since etofs read negative // HERE use an bool m_invertAcq to toggle this //
				} else {
					(*outmatPtr)[e][sample] = 0.0;
				}
			}
		}
		evt.put(outmatPtr,outkey);
	}

} // namespace CookieBox_pkg
