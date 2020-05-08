#include "CookieBox_pkg/Acqiris.h"
#include <sstream>
#include "CookieBox_pkg/CookieBox_mod.h"
#include <assert.h>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <memory>
#include <fftw/fftw3.h>
#include <cmath>
//#include "CookieBox_pkg/Constants.hpp"

using namespace CookieBox_pkg;

namespace CookieBox_pkg 
{
	Acqiris::Acqiris(std::vector<unsigned>& lims_in , std::vector<unsigned>& baselims_in)
	: m_getConfig(true)
	, m_use(false)
	, m_print(false)
	, m_nchannels(1)
	{
		// this is actually about 1/2 of the single counts, but the log2 function will send this to below 1 and get ignored
		m_thresh = std::vector<double>(16,1.e5);
		init(lims_in,baselims_in);
	}

	Acqiris::Acqiris(void)
	: m_getConfig(true)
	, m_use(false)
	, m_print(false)
	, m_nchannels(1)
	, wf_y(NULL)
	, wf_ddy(NULL)
	, wf_dy(NULL)
	, wf_Y_hc(NULL)
	, wf_DDY_hc(NULL)
	, wf_DY_hc(NULL)
	, plan_r2hc_Ptr(NULL)
	, plan_hc2r_Ptr(NULL)
	{
		// this is actually about 1/2 of the single counts, but the log2 function will send this to below 1 and get ignored
		m_thresh = std::vector<double>(16,1.e5);
	}

	Acqiris::~Acqiris(void)
	{
		if (m_outfile.is_open())
			m_outfile.close();
		if (wf_y != NULL){ fftw_free(wf_y); }
		if (wf_ddy != NULL) { fftw_free(wf_ddy); }
		if (wf_dy != NULL) { fftw_free(wf_dy); }
		if (wf_Y_hc != NULL) { fftw_free(wf_Y_hc); }
		if (wf_DDY_hc != NULL) { fftw_free(wf_DDY_hc); }
		if (wf_DY_hc != NULL) { fftw_free(wf_DY_hc); }
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
		m_thresh.resize(rhs.m_thresh.size());
		std::copy(rhs.m_thresh.begin(),rhs.m_thresh.end(),m_thresh.begin());
		pointplans(rhs);
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
		pointplans(rhs);
		m_thresh.resize(rhs.m_thresh.size());
		std::copy(rhs.m_thresh.begin(),rhs.m_thresh.end(),m_thresh.begin());
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
	void Acqiris::setthresh(std::vector<double> in)
	{
		m_thresh.resize(in.size(),1e5);
		std::copy(in.begin(),in.end(),m_thresh.begin());
	}
	void Acqiris::setbwdlim(std::vector<size_t> in)
	{
		m_bwd_lim.resize(in.size(),200);
		std::copy(in.begin(),in.end(),m_bwd_lim.begin());
	}
	void Acqiris::setbwdnr(std::vector<float> in)
	{
		m_bwd_nr.resize(in.size(),0.1);
		std::copy(in.begin(),in.end(),m_bwd_nr.begin());
	}

	bool Acqiris::init()
	{
		size_t sz = m_lims.at(bins);
		wf_y = (double *) fftw_malloc(sizeof(double) * sz);
		wf_ddy = (double *) fftw_malloc(sizeof(double) * sz);
		wf_dy = (double *) fftw_malloc(sizeof(double) * sz);
		wf_Y_hc = (double *) fftw_malloc(sizeof(double) * sz);
		wf_DDY_hc = (double *) fftw_malloc(sizeof(double) * sz);
		wf_DY_hc = (double *) fftw_malloc(sizeof(double) * sz);
		plan_r2hc_Ptr = NULL;
		plan_hc2r_Ptr = NULL;
		return true;
	}
	void Acqiris::pointplans(const Acqiris & rhs)
	{
		plan_r2hc_Ptr = rhs.plan_r2hc_Ptr;
		plan_hc2r_Ptr = rhs.plan_hc2r_Ptr;
	}
	void Acqiris::pointplans(fftw_plan * const forward,fftw_plan * const backward)
	{
		plan_r2hc_Ptr = forward;
		plan_hc2r_Ptr = backward;
	}
	void Acqiris::detachplans()
	{
		plan_r2hc_Ptr = NULL;
		plan_hc2r_Ptr = NULL;
	}

	void Acqiris::setmasterplans(fftw_plan * const forward,fftw_plan * const backward)
	{
		//assert(plan_r2hc_Ptr.use_count()==0 && plan_hc2r_Ptr.use_count()==0);
		size_t sz = m_lims.at(bins);
		*forward = fftw_plan_r2r_1d(sz,
				wf_y,
				wf_Y_hc,	
				FFTW_R2HC,
				FFTW_ESTIMATE);
		*backward = fftw_plan_r2r_1d(sz,
				wf_Y_hc,
				wf_y,
				FFTW_HC2R,
				FFTW_ESTIMATE);
		plan_r2hc_Ptr = forward;
		plan_hc2r_Ptr = backward;

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

	bool Acqiris::filtered_fill(Event& evt, Env& env, a5d_ll_2dview_t & slice, a4d_ll_1dview_t & shotslice)
	{
		// use this method to fill a slice of the multi_array, but incorporating the Y conv ddY filter, where we want the histogram.
		m_srcPtr = evt.get(m_srcStr);
		if (!m_srcPtr.get()){return false;}
		if (m_getConfig){std::cout << getConfig(env);}

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

			size_t sz = m_lims.at(bins);
			for (unsigned s = 0; s < sz; ++s) {
				long int val = (wf[segment][m_lims.at(start) + s] - basesum/m_baselims.at(bins));
				if (m_invert)
					val *= -1;
				wf_y[s] = double(val);
			}
			// fftw makes un-normalized transforms, so we need to divid by root(n) for each transform (or n for each back and forth
			fftw_execute_r2r(*plan_r2hc_Ptr, wf_y, wf_Y_hc );
			// fftw makes un-normalized transforms, so we need to divid by root(n) for each transform (or n for each back and forth
			scaleVec(sz,wf_Y_hc,1.0/std::sqrt(double(sz)));
			size_t bwd = m_bwd_lim[chan];
			float nr = m_bwd_nr[chan];
			for( size_t s=1;s<bwd;++s){
				double lpf = 0.5*(1.0 + std::cos(double(s)/(double)bwd *  M_PI));
				double f = lpf/(nr+lpf); // poor man's weiner filter
				wf_Y_hc[s] *= f;
				wf_Y_hc[sz-s] *= f;
			}	
			for( unsigned s=bwd;s<sz/2;++s){
				wf_Y_hc[s] = wf_Y_hc[sz-s] = wf_DDY_hc[s] = wf_DDY_hc[sz-s] = 0.;
			}
			wf_Y_hc[sz/2] = 0.;
			fftw_execute_r2r(*plan_hc2r_Ptr, wf_Y_hc, wf_y);
			scaleVec(sz,wf_y,1.0/std::sqrt(double(sz)));
			for (size_t s=0;s<sz;++s){
				slice[chan][s] += (int16_t)std::max(wf_y[s]/m_thresh[chan],0.);
				m_data[chan][s] = (uint16_t)std::max(wf_y[s]/m_thresh[chan],0.);
			}

			
			++shotslice[chan]; 
		}
		return true;

	}
	bool Acqiris::ydy_fill(Event& evt, Env& env, a5d_ll_2dview_t & slice, a4d_ll_1dview_t & shotslice)
	{
		// use this method to fill a slice of the multi_array, but incorporating the Y conv ddY filter, where we want the histogram.
		m_srcPtr = evt.get(m_srcStr);
		if (!m_srcPtr.get()){return false;}
		if (m_getConfig){std::cout << getConfig(env);}

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

			size_t sz = m_lims.at(bins);
			for (unsigned s = 0; s < sz; ++s) {
				long int val = (wf[segment][m_lims.at(start) + s] - basesum/m_baselims.at(bins));
				if (m_invert)
					val *= -1;
				wf_y[s] = double(val);
			}
			// fftw makes un-normalized transforms, so we need to divid by root(n) for each transform (or n for each back and forth
			fftw_execute_r2r(*plan_r2hc_Ptr, wf_y, wf_Y_hc );

			unsigned bwd = sz/3;
			wf_DY_hc[0] = 0.;
			for( unsigned s=1;s<bwd;++s){
				double lpf = (std::cos(double(s)/(double)bwd *  M_PI),int(2));
				wf_Y_hc[s] *= lpf;
				wf_Y_hc[sz-s] *= lpf;
				wf_DY_hc[s] = -std::pow(double(s),int(2))/std::pow(double(bwd),int(2)) * wf_Y_hc[s];
				wf_DY_hc[sz-s] = -std::pow(double(s),int(2))/std::pow(double(bwd),int(2)) * wf_Y_hc[sz-s];
			}	

			for( unsigned s=bwd;s<sz/2;++s){
				wf_Y_hc[s] = wf_Y_hc[sz-s] = wf_DY_hc[s] = wf_DY_hc[sz-s] = 0.;
			}

			fftw_execute_r2r(*plan_hc2r_Ptr, wf_Y_hc, wf_y);
			fftw_execute_r2r(*plan_hc2r_Ptr, wf_DY_hc, wf_dy );

			for ( unsigned s=0; s<sz; ++s){
				double res = 0.;
				m_data.at(chan).at(s) = int16_t(res);
				if ( (wf_y[s] > 0.) && (wf_dy[s] < 0.)){
					res = wf_y[s] * -1. * wf_dy[s] * std::pow(sz,int(-2)) / m_thresh[chan]; // this normalizes	
					res = log2(res); // using the to supress the pileup
					if (res >= 1.){					
						slice[chan][s] += (int16_t)std::max((int)res,0);
						m_data.at(chan).at(s) = (int16_t)std::max((int)res,0);
					}
				}
			}
			++shotslice[chan]; 
		}
		return true;
	}
	bool Acqiris::yddy_fill(Event& evt, Env& env, a5d_ll_2dview_t & slice, a4d_ll_1dview_t & shotslice)
	{
		// use this method to fill a slice of the multi_array, but incorporating the Y conv ddY filter, where we want the histogram.
		m_srcPtr = evt.get(m_srcStr);
		if (!m_srcPtr.get()){return false;}
		if (m_getConfig){std::cout << getConfig(env);}

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

			size_t sz = m_lims.at(bins);
			for (unsigned s = 0; s < sz; ++s) {
				long int val = (wf[segment][m_lims.at(start) + s] - basesum/m_baselims.at(bins));
				if (m_invert)
					val *= -1;
				wf_y[s] = double(val);
			}

			// fftw makes un-normalized transforms, so we need to divid by root(n) for each transform (or n for each back and forth
			fftw_execute_r2r(*plan_r2hc_Ptr, wf_y, wf_Y_hc );

			unsigned bwd = sz/8; // HERE HERE HERE HERE manually playing with the low pass filter bandwidth
			wf_DDY_hc[0] = 0.;
			for( unsigned s=1;s<bwd;++s){
				double lpf = (std::cos(double(s)/(double)bwd *  M_PI),int(2));
				wf_Y_hc[s] *= lpf;
				wf_Y_hc[sz-s] *= lpf;
				wf_DDY_hc[s] = -std::pow(double(s),int(2))/std::pow(double(bwd),int(2)) * wf_Y_hc[s];
				wf_DDY_hc[sz-s] = -std::pow(double(s),int(2))/std::pow(double(bwd),int(2)) * wf_Y_hc[sz-s];
			}	
			for( unsigned s=bwd;s<sz/2;++s){
				wf_Y_hc[s] = wf_Y_hc[sz-s] = wf_DDY_hc[s] = wf_DDY_hc[sz-s] = 0.;
			}

			fftw_execute_r2r(*plan_hc2r_Ptr, wf_Y_hc, wf_y);
			fftw_execute_r2r(*plan_hc2r_Ptr, wf_DDY_hc, wf_ddy );

			for ( unsigned s=0; s<sz; ++s){
				double res = 0.;
				m_data.at(chan).at(s) = int16_t(res);
				if ( (wf_y[s] > 0.) && (wf_ddy[s] < 0.)){
					res = wf_y[s] * -1. * wf_ddy[s] * std::pow(sz,int(-2)) / m_thresh[chan]; // this normalizes	
					res = log2(res); // using the to supress the pileup
					if (res >= 1.){					
						slice[chan][s] += (uint16_t)std::max((int)res,0);
						m_data.at(chan).at(s) = (uint16_t)std::max((int)res,0);
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
				int16_t val = (wf[segment][m_lims.at(start) + s] - basesum/m_baselims.at(bins));
				if (m_invert)
					val *= -1;
				val /= m_thresh[chan]; // here I'm using m_thresh[chan] to scale down so integer values are low... really this should fit in 8 bits
				slice[chan][s] += (uint16_t)std::max(val,(int16_t)0);
				m_data.at(chan).at(s) = (uint16_t)std::max(val,(int16_t)0);
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

		if (m_getConfig){std::cout << getConfig(env);}

		if (m_data.size() != m_nchannels)
			m_data.resize(m_nchannels);
		if (m_data.back().size() != m_lims.at(bins)){
			for (unsigned c=0;c<m_data.size();++c)
				m_data[c].resize(m_lims.at(bins),0);
		}

		for (unsigned chan=0;chan<m_nchannels;++chan){
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
				int16_t val = (wf[segment][m_lims.at(start) + s] - basesum/m_baselims.at(bins));
				if (m_invert)
					val *= -1;
				val /= m_thresh[chan]; // here I'm using m_thresh[chan] to scale down so integer values are low... really this should fit in 8 bits
				slice[chan][s] += (uint16_t)std::max(val,(int16_t)0);
				m_data.at(chan).at(s) = (uint16_t)std::max(val,(int16_t)0);
			}
		}
		return true;
	}

	bool Acqiris::fill(Event& evt, Env& env)
	{
		m_srcPtr = evt.get(m_srcStr);
		if (!m_srcPtr.get()){return false;}

		if (m_getConfig){std::cout << getConfig(env);}


		if (m_data.size() != m_nchannels)
			m_data.resize(m_nchannels);
		if (m_data.back().size() != m_lims.at(bins)){
			for (unsigned c=0;c<m_data.size();++c)
				m_data[c].resize(m_lims.at(bins),0);
		}

		for (unsigned chan=0;chan<m_nchannels;++chan){
			const unsigned segment = 0; // [chris ogrady] always 0 for LCLS data taken so far (a feature of the acqiris we don't use)
			const Psana::Acqiris::DataDescV1Elem& elem = m_srcPtr->data(chan);
			const ndarray<const Psana::Acqiris::TimestampV1, 1>& timestamps = elem.timestamp();
			double pos = timestamps[segment].pos();
			int indexFirstPoint = elem.indexFirstPoint();
			int indoffset = int(pos*2./m_sampleInterval)/2;
			wf_t wf = elem.waveforms(); 

			long long basesum;
			basesum = 0;
			for (unsigned s = m_baselims.at(start); s < m_baselims.at(stop);++s){
				// fill baseline //
				basesum += wf[segment][s] - m_vert_offset[chan];
			}
			for (unsigned s = 0; s < m_lims.at(bins); ++s) {
				//long int val = (wf[segment][m_lims.at(start) + s] - basesum/m_baselims.at(bins));
				int16_t val = 0;
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
				val /= m_thresh[chan]; // here I'm using m_thresh[chan] to scale down so integer values are low... really this should fit in 8 bits
				m_data.at(chan).at(s) = (uint16_t)std::max(val,(int16_t)0);
			}
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
	bool Acqiris::fillFancyFiltered(std::vector< std::vector< int16_t > > & filtered_data)
	{
		// HERE HERE HERE HERE manually playing with the low pass filter bandwidth
		if (m_data.size() != m_nchannels){
			std::cerr << "Trying to filter an empty m_data\n" << std::flush;
			return false;
		}
		size_t sz = m_data[0].size();
		filtered_data.resize(m_data.size());
		for (size_t c=0;c<m_data.size();++c){
			filtered_data[c].resize(sz);
			for (size_t s=0;s<sz;++s){
				wf_y[s] = double(m_data[c][s]);
			}
			fftw_execute_r2r(*plan_r2hc_Ptr, wf_y, wf_Y_hc );
			// fftw makes un-normalized transforms, so we need to divid by root(n) for each transform (or n for each back and forth
			scaleVec(sz,wf_Y_hc,1.0/std::sqrt(double(sz)));
			size_t bwd = m_bwd_lim[c];
			float nr = m_bwd_nr[c];
			for( size_t s=1;s<bwd;++s){
				double lpf = 0.5*(1.0 + std::cos(double(s)/(double)bwd *  M_PI));
				double f = lpf/(nr+lpf); // poor man's weiner filter
				wf_Y_hc[s] *= f;
				wf_Y_hc[sz-s] *= f;
			}	
			for( unsigned s=bwd;s<sz/2;++s){
				wf_Y_hc[s] = wf_Y_hc[sz-s] = wf_DDY_hc[s] = wf_DDY_hc[sz-s] = 0.;
			}
			wf_Y_hc[sz/2] = 0.;
			fftw_execute_r2r(*plan_hc2r_Ptr, wf_Y_hc, wf_y);
			scaleVec(sz,wf_y,1.0/std::sqrt(double(sz)));
			for (size_t s=0;s<sz;++s){
				filtered_data[c][s] = int16_t(wf_y[s]);
			}
		}

		
		return true;
	}
	bool Acqiris::print_out(const unsigned eventnum)
	{
		std::string eventfilename(m_filename);
		eventfilename += ".event_" + boost::lexical_cast<std::string>(eventnum);
		if (m_outfile.is_open())
			m_outfile.close();
		m_outfile.open(eventfilename.c_str(),std::ios::out);
		std::vector< std::vector< int16_t > > filtered_data;
		bool filteredexists = fillFancyFiltered(filtered_data);
		std::ofstream filtered_outfile;
		if (filteredexists){
			std::string filteredeventfilename(m_filename);
			filteredeventfilename += ".filtered.event_" + boost::lexical_cast<std::string>(eventnum);
			filtered_outfile.open(filteredeventfilename.c_str(),std::ios::out);
		}
		if (!m_outfile.is_open())
			return false;
		if (filteredexists && !filtered_outfile.is_open())
			return false;
		if (m_data.size() != m_nchannels) 
			return false;
		for (unsigned s=0;s<m_data.at(0).size();++s){
			for (unsigned c=0;c<m_data.size();++c){
				m_outfile << m_data.at(c).at(s) << "\t";
				if(filteredexists) filtered_outfile << filtered_data.at(c).at(s) << "\t";
			}
			m_outfile << "\n";
			if (filteredexists) filtered_outfile << "\n";
		}
		m_outfile.close();
		if (filteredexists && filtered_outfile.is_open()) filtered_outfile.close();
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

