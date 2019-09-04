#include "CookieBox_pkg/Acqiris.h"
#include <sstream>
#include "CookieBox_pkg/CookieBox_mod.h"
#include <assert.h>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <fstream>

using namespace CookieBox_pkg;

namespace CookieBox_pkg 
{
	Acqiris::Acqiris(std::vector<unsigned>& lims_in , std::vector<unsigned>& baselims_in)
	: m_getConfig(true)
	, m_use(false)
	, m_print(false)
	, m_nchannels(1)
	{
		init(lims_in,baselims_in);
	}

	Acqiris::Acqiris(void)
	: m_getConfig(true)
	, m_use(false)
	, m_print(false)
	, m_nchannels(1)
	, wf_y(NULL)
	, wf_ddy(NULL)
	, wf_Y_hc(NULL)
	, wf_DDY_hc(NULL)
	, plan_r2hc_Ptr(NULL)
	, plan_hc2r_Ptr(NULL)
	{
	}

	Acqiris::~Acqiris(void)
	{
		if (m_outfile.is_open())
			m_outfile.close();
		//std::cerr << "Now killing off fftw vectors" << std::endl;
		if (wf_y != NULL){ fftw_free(wf_y); }
		if (wf_ddy != NULL) { fftw_free(wf_ddy); }
		if (wf_Y_hc != NULL) { fftw_free(wf_Y_hc); }
		if (wf_DDY_hc != NULL) { fftw_free(wf_DDY_hc); }
	}
	Acqiris & Acqiris::operator= ( const Acqiris & rhs )
	{
		m_srcPtr = rhs.m_srcPtr;
		m_getConfig=rhs.m_getConfig;
		m_use=rhs.m_use;
		m_srcStr=rhs.m_srcStr;
		m_ConfigPtr=rhs.m_ConfigPtr;
		m_print=rhs.m_print;
		m_filename=rhs.m_filename;
		m_nchannels=rhs.m_nchannels;
		m_max_samples=rhs.m_max_samples;

		m_lims.resize(rhs.m_lims.size());
		m_baselims.resize(rhs.m_baselims.size());
		std::copy(rhs.m_lims.begin(),rhs.m_lims.end(),m_lims.begin());
		std::copy(rhs.m_baselims.begin(),rhs.m_baselims.end(),m_baselims.begin());
		init();
		deepcopy_data(rhs);
		return *this;
	}

	Acqiris::Acqiris(const Acqiris & rhs)
	: m_srcPtr(rhs.m_srcPtr)
	, m_getConfig(rhs.m_getConfig)
	, m_use(rhs.m_use)
	, m_srcStr(rhs.m_srcStr)
	, m_ConfigPtr(rhs.m_ConfigPtr)
	, m_print(rhs.m_print)
	, m_filename(rhs.m_filename)
	, m_nchannels(rhs.m_nchannels)
	, m_max_samples(rhs.m_max_samples)
	{ 
		m_lims.resize(rhs.m_lims.size());
		m_baselims.resize(rhs.m_baselims.size());
		std::copy(rhs.m_lims.begin(),rhs.m_lims.end(),m_lims.begin());
		std::copy(rhs.m_baselims.begin(),rhs.m_baselims.end(),m_baselims.begin());
		init();
		deepcopy_data(rhs);

	}
	void Acqiris::deepcopy_data(const Acqiris & b)
	{
		m_data.resize(b.m_data.size());
		for ( unsigned i=0;i<b.m_data.size();++i){
			m_data[i].resize(b.m_data[i].size());
			for (unsigned j=0;j<b.m_data[i].size();++j)
				m_data[i][j] = b.m_data[i][j];
		}
	}

	bool Acqiris::init( std::vector<unsigned>& lims_in
			, std::vector<unsigned>& baselims_in)
	{
		std::cerr << "lims_in.size() = " << lims_in.size() << " baselims_in.size() " << baselims_in.size() << std::endl;
		if (! ( lims_in.size() == baselims_in.size() &&  lims_in.size() == 2 ) ){
			std::cerr << "Failed sizes in Acqiris::init() method " << std::endl;
			return false;
		}
		m_lims = lims_in;
		m_baselims = baselims_in;
		m_lims.push_back(lims_in.back() - lims_in.front());
		m_baselims.push_back(baselims_in.back() - baselims_in.front());
		std::cerr << "m_lims =\t";  
		for (unsigned i=0;i<m_lims.size();++i)
			std::cerr << m_lims.at(i) << "\t";
		std::cerr << "\nm_baselims =\t";
		for (unsigned i=0;i<m_baselims.size();++i)
			std::cerr << m_baselims.at(i) << "\t";
		std::cerr << "\n" << std::flush;
		std::cerr << "m_nchannels (in Acqiris::init() ) = " << m_nchannels << std::endl;

		return init();
	}
	bool Acqiris::init()
	{

		std::cerr << "Now setting up the fftw vectors and... \n" << std::flush;

		size_t sz = m_lims.at(bins);
		wf_y = (double *) fftw_malloc(sizeof(double) * sz);
		wf_ddy = (double *) fftw_malloc(sizeof(double) * sz);
		wf_Y_hc = (double *) fftw_malloc(sizeof(double) * sz);
		wf_DDY_hc = (double *) fftw_malloc(sizeof(double) * sz);

		std::cerr << "... plans \n" << std::flush;

		setplans(rhs);

		std::cerr << "OK, fftw for Acqiris is initialized \n" << std::flush;
		return true;
	}
	void Acqiris::setplans(const Acqiris & rhs)
	{
		plan_r2hc_Ptr = rhs.plan_r2hc_Ptr;
		plan_hc2r_Ptr = rhs.plan_hc2r_dPtr;
	}

	void Acqiris::setmasterplans(fftw_plan * const forward,fftw_plan * const backward)
	{
		assert(plan_r2hc_Ptr.use_count()==0 && plan_hc2r_Ptr.use_count()==0);
		size_t sz = m_lims.at(bins);
		*forward = fftw_plan_r2hc(sz,
				wf_y,
				wf_Y,	
				FFTW_ESTIMATE);
		*backward = fftw_plan_hc2r(sz,
				wf_Y_hc,
				wf_y,
				FFTW_ESTIMATE);
		plan_r2hc_Ptr = std::make_shared<fftw_plan> (*forward);
		plan_hc2r_Ptr = std::make_shared<fftw_plan> (*backward);

	}

	void Acqiris::srcStr(Source srcStr_in){
		m_srcStr = srcStr_in;
	}
	
	void Acqiris::evtput(Event& evt, Env& env, const unsigned chan)
	{
		// use something like a 2d ndarray to pass including the chan bins to python
		std::string outkey("aq_");
		outkey += boost::lexical_cast<std::string>(chan);
		unsigned shape[] = {m_data[chan].size()};
		boost::shared_ptr< ndarray<double,1> > outmatPtr = boost::make_shared < ndarray<double,1> > ( shape );

		for (unsigned j=0;j<m_data[chan].size();++j){
			(*outmatPtr)[j] = m_data[chan][j];
		}
		evt.put(outmatPtr,outkey);
	}
	void Acqiris::evtput(Event& evt, Env& env)
	{
		// use something like this just to pass a single channel data to python

		std::string outkey("aq_allchans");

		unsigned shape[] = {m_data.size(), m_data[0].size()};
		boost::shared_ptr< ndarray<double,2> > outmatPtr = boost::make_shared < ndarray<double,2> > ( shape );

		for (unsigned i = 0; i < m_data.size();++i){ 
			for (unsigned j=0;j<m_data[i].size();++j){
				(*outmatPtr)[i][j] = m_data[i][j];
			}
		}
		evt.put(outmatPtr,outkey);
	}




	//	Eventually... //
	//	Lets use the sparsity rather than fill the entire ND array	//
	//	std::map gives a key::value pair	//

	bool Acqiris::fancyfill(Event& evt, Env& env, a5d_ll_2dview_t & slice, a4d_ll_1dview_t & shotslice)
	{
		// use this method to fill a slice of the multi_array, but incorporating the Y conv ddY filter, where we want the histogram.
		m_srcPtr = evt.get(m_srcStr);
		if (!m_srcPtr.get()){return false;}
		//std::cerr << "m_nchannels in Acqiris::fill() = " << m_nchannels << std::endl;
		if (m_getConfig){std::cout << getConfig(env);}
		//std::cerr << "m_nchannels in Acqiris::fill() after m_getConfig test = " << m_nchannels << std::endl;

		if (m_data.size() != m_nchannels)
			m_data.resize(m_nchannels);
		if (m_data.back().size() != m_lims.at(bins)){
			for (unsigned c=0;c<m_data.size();++c)
				m_data[c].resize(m_lims.at(bins),0);
		}

		for (unsigned chan=0;chan<m_nchannels;++chan){
			wf_t wf = m_srcPtr->data(chan).waveforms(); // the 2D'ness of this is for the unused segments, not the channels.
			

			const int segment = 0; // [chris ogrady] always 0 for LCLS data taken so far (a feature of the acqiris we don't use)
			long long basesum;
			basesum = 0;
			for (unsigned s = m_baselims.at(start); s < m_baselims.at(stop);++s){
				// fill baseline //
				basesum += wf[segment][s];
			}

			//HERE HERE HERE HERE working on the fftw for sake of Y conv ddY filtering


			size_t sz = m_lims.at(bins);
			for (unsigned s = 0; s < sz; ++s) {
				long int val = (wf[segment][m_lims.at(start) + s] - basesum/m_baselims.at(bins));
				if (m_invert)
					val *= -1;
				wf_y[s] = double(val);
				m_data.at(chan).at(s) = val;
			}

			fftw_execute_r2r(*plan_r2hc_Ptr.get(), wf_y, wf_Y_hc );

			unsigned bwd = sz/3;
			wf_DDY_hc[0] = 0.;
			for( unsigned s=1;s<bwd;++s){
				double lpf = (std::cos(double(s)/(double)bwd),int(2));
				wf_Y_hc[s] *= lpf;
				wf_Y_hc[sz-s] *= lpf;
				wf_DDY_hc[s] = -std::pow(double(s),int(2)) * wf_Y_hc[s];
				wf_DDY_hc[sz-s] = -std::pow(double(s),int(2)) * wf_Y_hc[sz-s];
			}	
			for( unsigned s=bwd;s<sz/2;++s){
				wf_Y_hc[s] = wf_Y_hc[sz-s] = wf_DDY_hc[s] = wf_DDY_hc[sz-s] = 0.;
			}

			fftw_execute_r2r(*plan_hc2r_Ptr.get(), wf_Y_hc, wf_y );
			fftw_execute_r2r(*plan_hc2r_Ptr.get(), wf_DDY_hc, wf_ddy );

			double thresh = 5.e-5;
			for ( unsigned s=0; s<sz; ++s){
				double res = 0.;
				if ( (wf_y[s] > 0.) && (wf_ddy[s] < 0.)){
					res = wf_y[s] * -1. * wf_ddy[s] / thresh;	
					if (res > 1){					
						slice[chan][s] += std::log(res);
					}
				}
			}
			++shotslice[chan]; 
		}
		return true;
	}

	bool Acqiris::fill(Event& evt, Env& env, a5d_ll_2dview_t & slice, a4d_ll_1dview_t & shotslice)
	{
		// use this method to fill a slice of the multi_array where we want the histogram.
		m_srcPtr = evt.get(m_srcStr);
		if (!m_srcPtr.get()){return false;}

		//std::cerr << "m_nchannels in Acqiris::fill() = " << m_nchannels << std::endl;
		if (m_getConfig){std::cout << getConfig(env);}
		//std::cerr << "m_nchannels in Acqiris::fill() after m_getConfig test = " << m_nchannels << std::endl;

		if (m_data.size() != m_nchannels)
			m_data.resize(m_nchannels);
		if (m_data.back().size() != m_lims.at(bins)){
			for (unsigned c=0;c<m_data.size();++c)
				m_data[c].resize(m_lims.at(bins),0);
		}

		for (unsigned chan=0;chan<m_nchannels;++chan){
			wf_t wf = m_srcPtr->data(chan).waveforms(); // the 2D'ness of this is for the unused segments, not the channels.

			const int segment = 0; // [chris ogrady] always 0 for LCLS data taken so far (a feature of the acqiris we don't use)
			long long basesum;
			basesum = 0;
			for (unsigned s = m_baselims.at(start); s < m_baselims.at(stop);++s){
				// fill baseline //
				basesum += wf[segment][s];
			}
			for (unsigned s = 0; s < m_lims.at(bins); ++s) {
				long int val = (wf[segment][m_lims.at(start) + s] - basesum/m_baselims.at(bins));
				if (m_invert)
					val *= -1;
				slice[chan][s] += val;
				m_data.at(chan).at(s) = val;
			}
			++shotslice[chan]; 
		}
		return true;
	}

	bool Acqiris::fill(Event& evt, Env& env, a4d_ll_2dview_t & slice)
	{
		// use this method to fill a slice of the multi_array where we want the histogram.
		m_srcPtr = evt.get(m_srcStr);
		if (!m_srcPtr.get()){return false;}

		//std::cerr << "m_nchannels in Acqiris::fill() = " << m_nchannels << std::endl;
		if (m_getConfig){std::cout << getConfig(env);}
		//std::cerr << "m_nchannels in Acqiris::fill() after m_getConfig test = " << m_nchannels << std::endl;

		if (m_data.size() != m_nchannels)
			m_data.resize(m_nchannels);
		if (m_data.back().size() != m_lims.at(bins)){
			for (unsigned c=0;c<m_data.size();++c)
				m_data[c].resize(m_lims.at(bins),0);
		}

		for (unsigned chan=0;chan<m_nchannels;++chan){
			// wf_t wf = m_srcPtr->data(chan).waveforms(); // the 2D'ness of this is for the unused segments, not the channels.
			// Also, we define it here since as an interface, 
			// it won't get used elsewhere and therefore should live on the stack.. .not heap as with new .. and delete...
			// I now also see that since I'm burying the wf type in the function, 
			// I want something in the .h file where I can tweak the types... thus the typedef ndarray<short,2> wf_t;

			const Psana::Acqiris::DataDescV1Elem& elem = m_srcPtr->data(chan);
			wf_t wf = elem.waveforms(); // the 2D'ness of this is for the unused segments, not the channels.

			const int segment = 0; // [chris ogrady] always 0 for LCLS data taken so far (a feature of the acqiris we don't use)
			long long basesum;
			basesum = 0;
			for (unsigned s = m_baselims.at(start); s < m_baselims.at(stop);++s){
				// fill baseline //
				basesum += wf[segment][s];
			}
			for (unsigned s = 0; s < m_lims.at(bins); ++s) {
				long int val = (wf[segment][m_lims.at(start) + s] - basesum/m_baselims.at(bins));
				if (m_invert)
					val *= -1;
				slice[chan][s] += val;
				m_data.at(chan).at(s) = val;
			}
		}
		return true;
	}

	bool Acqiris::fill(Event& evt, Env& env)
	{
		// HERE HERE HERE //
		// Using this one //
		// find a way to get the time-vector start value and step.
		// Use that to compute an offset
		m_srcPtr = evt.get(m_srcStr);
		if (!m_srcPtr.get()){return false;}

		std::cerr << "trying to test to getConfig() in Acqiris::fill() method" << std::endl;
		if (m_getConfig){std::cout << getConfig(env);}
		std::cerr << "m_nchannels = " << m_nchannels << std::endl;


		if (m_data.size() != m_nchannels)
			m_data.resize(m_nchannels);
		if (m_data.back().size() != m_lims.at(bins)){
			for (unsigned c=0;c<m_data.size();++c)
				m_data[c].resize(m_lims.at(bins),0);
		}
		//std::cerr << "m_data.size() = " << m_data.size() << std::endl;
		//std::cerr << "m_data.back().size() = " << m_data.back().size() << std::endl;

		for (unsigned chan=0;chan<m_nchannels;++chan){
			const unsigned segment = 0; // [chris ogrady] always 0 for LCLS data taken so far (a feature of the acqiris we don't use)
			const Psana::Acqiris::DataDescV1Elem& elem = m_srcPtr->data(chan);
			//std::cerr << "filling channel " << chan << std::flush;
			const ndarray<const Psana::Acqiris::TimestampV1, 1>& timestamps = elem.timestamp();
			double pos = timestamps[segment].pos();
			int indexFirstPoint = elem.indexFirstPoint();
			int indoffset = int(pos*2./m_sampleInterval)/2;
			wf_t wf = elem.waveforms(); 
			// the 2D'ness of this is for the unused segments, not the channels.
			// Also, we define it here since as an interface, 
			// it won't get used elsewhere and therefore should live on the stack.. .not heap as with new .. and delete...
			// I now also see that since I'm burying the wf type in the function, 
			// I want something in the .h file where I can tweak the types... thus the typedef ndarray<short,2> shortwf_t;

			// HERE HERE HERE HERE //
			// Debugging the indices and the clock start time //
			// here the 0th index can correspond to a time form 0 back to the increment size (.5ns) here 
			// This means we could sub-divide the waveforms for computing the energy simply by adjusting by this offset.
	   //const unsigned shape[] = {nbrChannels, nbrSamples};
	   //    ndarray<wform_t, 2> wf = make_ndarray<wform_t>(nbrChannels, nbrSamples);
	   //        ndarray<wtime_t, 2> wt = make_ndarray<wtime_t>(nbrChannels, nbrSamples);
	   //


			long long basesum;
			basesum = 0;
			for (unsigned s = m_baselims.at(start); s < m_baselims.at(stop);++s){
				// fill baseline //
				basesum += wf[segment][s] - m_vert_offset[chan];
			}
			//std::cerr << "\t ... got baseline channel " << chan << std::flush;
			for (unsigned s = 0; s < m_lims.at(bins); ++s) {
				//long int val = (wf[segment][m_lims.at(start) + s] - basesum/m_baselims.at(bins));
				long int val = 0;
				if (indoffset != 0){
					if ((s + indoffset > 0 ) && (s + indoffset < m_max_samples)){
						val = wf[segment][m_lims.at(start) + s + indoffset];
						val -= m_vert_offset[chan];
						val -= basesum/m_baselims.at(bins);
						val -= s*m_vert_slope[chan];
					}
				} else {
					val = wf[segment][m_lims.at(start) + s + indoffset];
					val -= m_vert_offset[chan];
					val -= basesum/m_baselims.at(bins);
					val -= s*m_vert_slope[chan];
				}
				if (m_invert)
					val *= -1;
				m_data.at(chan).at(s) = val;
			}
			//std::cerr << "\t ... filled channel " << chan << std::endl;
		}
		return true;
	}

	std::string Acqiris::getConfig(Env& env)
	{
		if (!m_getConfig)
			return "Already go tconfig";
		std::cerr << "Entered Acqiris::getConfig() method " << std::endl;
		m_ConfigPtr = env.configStore().get(m_srcStr, &m_src);
		std::stringstream ss;
		if (m_ConfigPtr != NULL) {
			const ndarray<const Psana::Acqiris::VertV1, 1>& vert = m_ConfigPtr->vert();

			ss  << "Acqiris::ConfigV1:\n"
				<< "  nbrBanks="    << m_ConfigPtr->nbrBanks()
				<< " channelMask="  << m_ConfigPtr->channelMask()
				<< " nbrChannels="  << m_ConfigPtr->nbrChannels()
				<< " nbrConvertersPerChannel=" << m_ConfigPtr->nbrConvertersPerChannel();
			if (m_rank == m_root_rank){
				const Psana::Acqiris::HorizV1& h = m_ConfigPtr->horiz();
				ss   << "\n  horiz: sampInterval=" << h.sampInterval()
					<< " delayTime="              << h.delayTime()
					<< " nbrSegments="            << h.nbrSegments()
					<< " nbrSamples="             << h.nbrSamples();

				for (unsigned ch = 0; ch < m_ConfigPtr->nbrChannels(); ++ ch) {
					const Psana::Acqiris::VertV1& v = vert[ch];
					ss  << "\n  vert(" << ch << "):"
						<< " fullScale="  << v.fullScale()
						     << " slope="      << v.slope()
						     << " offset="     << v.offset()
						     << " coupling="   << v.coupling()
						     << " bandwidth="  << v.bandwidth();
				}
			}

			// Put this where all ranks see it.
			//
			m_nchannels = m_ConfigPtr->nbrChannels();
			std::cerr << "m_nchannels = " << m_nchannels << std::endl;

			m_vert_slope.resize(m_nchannels,0.);
			m_vert_offset.resize(m_nchannels,0.);

			for (unsigned c=0 ; c<m_nchannels ; ++c){ // used for correcting even the short int data vactors
				//const Psana::Acqiris::VertV1& v = vert[c];
				m_vert_slope[c]  = vert[c].slope();
				m_vert_offset[c] = vert[c].offset();
			}
			//const Psana::Acqiris::HorizV1& h = m_ConfigPtr->horiz();
			//double sampInterval = h.sampInterval();
			//uint32_t nbrSamples = h.nbrSamples();

			m_sampleInterval = m_ConfigPtr->horiz().sampInterval();
	   
			m_max_samples = m_ConfigPtr->horiz().nbrSamples();
			if (m_max_samples < m_lims.at(stop)){
				m_lims.at(bins) = m_max_samples - m_lims.at(start); // ensure to not over-run
			}
			m_getConfig = false;
			return ss.str();
		}
		m_getConfig = true;
		return std::string("WARNING! Acqiris::ConfigV1 is not found...\n");
	}

	// File interface methods //
	void Acqiris::open_file( std::string & filename ){
		m_filename = filename;
		m_outfile.open(m_filename.c_str(),std::ios::out);
	}
	void Acqiris::close_file(void){
		if (m_outfile.is_open())
			m_outfile.close();
	}
	bool Acqiris::print_out(const unsigned eventnum)
	{
		std::string eventfilename(m_filename);
		eventfilename += ".event_" + boost::lexical_cast<std::string>(eventnum);
		if (m_outfile.is_open())
			m_outfile.close();
		m_outfile.open(eventfilename.c_str(),std::ios::out);
		if (!m_outfile.is_open())
			return false;
		for (unsigned s=0;s<m_data.at(0).size();++s){
			for (unsigned c=0;c<m_data.size();++c){
				long int val = m_data.at(c).at(s);
				m_outfile << val << "\t";
			}
			m_outfile << "\n";
		}
		m_outfile.close();
		return true;
	}
	bool Acqiris::print_out(void){
		if (!m_outfile.is_open())
			return false;
		for (unsigned c=0;c<m_data.size();++c){
			for (unsigned s=0;s<m_data.at(c).size();++s){
				int val = m_data.at(c).at(s);
				m_outfile << val << "\t";
			}
			m_outfile << "\n";
		}
		return true;
	}

	void Acqiris::print_header(void){
		m_outfile << "#AcqirisTraces along the rows, one row per channel";
		m_outfile << "\n";
	}

} // namespace CookieBox_pkg

