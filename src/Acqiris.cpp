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
	{
	}

	Acqiris::~Acqiris(void)
	{
		if (m_outfile.is_open())
			m_outfile.close();
	}
	Acqiris & Acqiris::operator=( Acqiris rhs ){ // the compiler is making a pass by copy here //
		swap(*this,rhs); // using the swap makes this exception safe... 
		return *this;
	}

	Acqiris::Acqiris(const Acqiris & b)
	: m_srcPtr(b.m_srcPtr)
	, m_getConfig(b.m_getConfig)
	, m_use(b.m_use)
	, m_srcStr(b.m_srcStr)
	, m_ConfigPtr(b.m_ConfigPtr)
	, m_print(b.m_print)
	, m_filename(b.m_filename)
	, m_nchannels(b.m_nchannels)
	, m_max_samples(b.m_max_samples)
	{ // not so sure this is exception safe
		m_lims.resize(b.m_lims.size());
		m_baselims.resize(b.m_baselims.size());
		for (unsigned i=0;i<b.m_lims.size();++i)
			m_lims[i] = b.m_lims[i];
		for (unsigned i=0;i<b.m_baselims.size();++i)
			m_baselims[i] = b.m_baselims[i];

		this->m_data = b.m_data; 

		if (m_print){
			m_filename += ".copy";
			if (b.m_outfile.is_open())
				m_outfile.open(m_filename.c_str(),std::ios::out);
		}
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
	void Acqiris::swap(Acqiris& a, Acqiris& b){
		std::swap(a.m_getConfig,b.m_getConfig);
		std::swap(a.m_use,b.m_use);
		std::swap(a.m_srcStr,b.m_srcStr);
		std::swap(a.m_ConfigPtr,b.m_ConfigPtr);
		std::swap(a.m_nchannels , b.m_nchannels);
		std::swap(a.m_max_samples , b.m_max_samples);
		std::swap(a.m_rank,b.m_rank);
		std::swap(a.m_mpi_size,b.m_mpi_size);

		m_data.swap(b.m_data);
		m_lims.swap(b.m_lims);
		m_baselims.swap(b.m_baselims);

		std::swap(a.m_print,b.m_print);
		std::swap(a.m_filename,b.m_filename);
		if (a.m_outfile.is_open() && b.m_outfile.is_open()){
			a.m_outfile.close();
			b.m_outfile.close();
			a.m_outfile.open(a.m_filename.c_str(),std::ios::app);
			b.m_outfile.open(b.m_filename.c_str(),std::ios::app);
			return;
		}
		if (a.m_outfile.is_open() && !b.m_outfile.is_open()){
			a.m_outfile.close();
			b.m_outfile.open(b.m_filename.c_str(),std::ios::app);
			return;
		}
		if (!a.m_outfile.is_open() && b.m_outfile.is_open()){
			b.m_outfile.close();
			a.m_outfile.open(a.m_filename.c_str(),std::ios::app);
			return;
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
		//std::vector<unsigned>::iterator itr = lims_in.begin();
		//std::vector<unsigned>::iterator itr_base = baselims_in.begin();
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
		return true;
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
			shortwf_t wf = m_srcPtr->data(chan).waveforms(); // the 2D'ness of this is for the unused segments, not the channels.

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
				m_data.at(chan).at(s) += val;
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
			shortwf_t wf = m_srcPtr->data(chan).waveforms(); // the 2D'ness of this is for the unused segments, not the channels.
			// Also, we define it here since as an interface, 
			// it won't get used elsewhere and therefore should live on the stack.. .not heap as with new .. and delete...
			// I now also see that since I'm burying the wf type in the function, 
			// I want something in the .h file where I can tweak the types... thus the typedef ndarray<short,2> shortwf_t;

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
				m_data.at(chan).at(s) += val;
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
			//std::cerr << "filling channel " << chan << std::flush;
			shortwf_t wf = m_srcPtr->data(chan).waveforms(); // the 2D'ness of this is for the unused segments, not the channels.
			// Also, we define it here since as an interface, 
			// it won't get used elsewhere and therefore should live on the stack.. .not heap as with new .. and delete...
			// I now also see that since I'm burying the wf type in the function, 
			// I want something in the .h file where I can tweak the types... thus the typedef ndarray<short,2> shortwf_t;

			// HERE HERE HERE HERE //
			// Debugging the indices and the clock start time //
			// here the 0th index can correspond to a time form 0 back to the increment size (.5ns) here 
			// This means we could sub-divide the waveforms for computing the energy simply by adjusting by this offset.
	//		doublewt_t wt = m_srcPtr->data(chan).wftime();  // this will fail... 
	   //const unsigned shape[] = {nbrChannels, nbrSamples};
	   //    ndarray<wform_t, 2> wf = make_ndarray<wform_t>(nbrChannels, nbrSamples);
	   //        ndarray<wtime_t, 2> wt = make_ndarray<wtime_t>(nbrChannels, nbrSamples);
	   //


			const int segment = 0; // [chris ogrady] always 0 for LCLS data taken so far (a feature of the acqiris we don't use)
			long long basesum;
			basesum = 0;
	//		double step = std::abs(wt[segment][1] - wt[segment][0]);
	//		int offset = int(2.*wt[segment][0]/step);
			for (unsigned s = m_baselims.at(start); s < m_baselims.at(stop);++s){
				// fill baseline //
				basesum += wf[segment][s];
			}
			//std::cerr << "\t ... got baseline channel " << chan << std::flush;
			for (unsigned s = 0; s < m_lims.at(bins); ++s) {
				long int val = (wf[segment][m_lims.at(start) + s] - basesum/m_baselims.at(bins));
				//long int val = (wf[segment][m_lims.at(start) + s + offset] - basesum/m_baselims.at(bins));
				if (m_invert)
					val *= -1;
				m_data.at(chan).at(s) += val;
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

				const ndarray<const Psana::Acqiris::VertV1, 1>& vert = m_ConfigPtr->vert();
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
			m_nchannels = m_ConfigPtr->nbrChannels();
			std::cerr << "m_nchannels = " << m_nchannels << std::endl;
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
				int val = m_data.at(c).at(s);
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

