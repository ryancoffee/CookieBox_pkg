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
	{
		init(lims_in,baselims_in);
	}

	Acqiris::Acqiris(void)
	: m_getConfig(true)
	, m_use(false)
	, m_print(false)
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
		this->m_integlims = b.m_integlims;
		//deepcopy_data(b);
		//deepcopy_integlims(b);

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
	void Acqiris::deepcopy_integlims(const Acqiris & b)
	{
		m_integlims.resize(b.m_integlims.size());
		for ( unsigned i=0;i<b.m_integlims.size();++i){
			m_integlims[i].resize(b.m_integlims[i].size());
			for (unsigned j=0;j<b.m_integlims[i].size();++j)
				m_integlims[i][j] = b.m_integlims[i][j];
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
		m_integlims.swap(b.m_integlims);

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
			(*outmatPtr)[j] *= -1.; // reverses polarity since etofs read negative // HERE use an bool m_invertAcq to toggle this //
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
				(*outmatPtr)[i][j] *= -1.; // reverses polarity since etofs read negative // HERE use an bool m_invertAcq to toggle this //
			}
		}
		evt.put(outmatPtr,outkey);
	}




	//	Eventually... //
	//	Lets use the sparsity rather than fill the entire ND array	//
	//	std::map gives a key::value pair	//

	bool Acqiris::fill(Event& evt, Env& env, a4d_ll_2dview_t & slice)
	{
		// use this method to fill a slice of the multi_array where we want the histogram.
		if (m_getConfig){std::cout << getConfig(env);}

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
				slice[chan][s] += (wf[segment][m_lims.at(start) + s] - basesum/m_baselims.at(bins));
			}
		}
		return true;
	}
	bool Acqiris::fill(Event& evt, Env& env)
	{
		m_srcPtr = evt.get(m_srcStr);
		if (!m_srcPtr.get()){return false;}

		//std::cerr << "trying to test to getConfig() in Acqiris::fill() method" << std::endl;
		if (m_getConfig){std::cout << getConfig(env);}
		//std::cerr << "m_nchannels = " << m_nchannels << std::endl;


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

			const int segment = 0; // [chris ogrady] always 0 for LCLS data taken so far (a feature of the acqiris we don't use)
			long long basesum;
			basesum = 0;
			for (unsigned s = m_baselims.at(start); s < m_baselims.at(stop);++s){
				// fill baseline //
				basesum += wf[segment][s];
			}
			//std::cerr << "\t ... got baseline channel " << chan << std::flush;
			for (unsigned s = 0; s < m_lims.at(bins); ++s) {
				m_data.at(chan).at(s) += (wf[segment][m_lims.at(start) + s] - basesum/m_baselims.at(bins));
			}
			//std::cerr << "\t ... filled channel " << chan << std::endl;
		}
		return true;
	}

	void Acqiris::integlims(std::vector < std::vector <unsigned> > & in)
	{
		std::vector <unsigned>::iterator itr_lows = in.front().begin();
		std::vector <unsigned>::iterator itr_highs = in.back().begin(); 

		m_integlims.resize(in.size());
		m_integlims.front().assign(itr_lows,in.front().end());
		m_integlims.back().assign(itr_highs,in.back().end());
	}
	
	std::string Acqiris::getConfig(Env& env)
	{
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
				m_outfile << m_data.at(c).at(s) << "\t";
			}
			m_outfile << "\n";
		}
		return true;
	}
	bool Acqiris::print_out(void){
		if (!m_outfile.is_open())
			return false;
		for (unsigned c=0;c<m_data.size();++c){
			for (unsigned s=0;s<m_data.at(c).size();++s){
				m_outfile << m_data[c][s] << "\t";
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

