#include "CookieBox_pkg/TimeTool.h"
#include "CookieBox_pkg/CookieBox_mod.h"
#include <fstream>
#include <iostream>
#include <assert.h>
#include <iterator>
#include <exception>
#include <functional> 

//#include <initializer_list>
// This declares this class as psana module
using namespace CookieBox_pkg;
//PSANA_MODULE_FACTORY(CookieBox_mod)

using namespace CookieBox_pkg;

namespace CookieBox_pkg{
	TimeTool::TimeTool(std::vector<unsigned>& bins_in
			, std::vector<double>& mins_in
			, std::vector<double>& maxs_in)
	: m_use(false)
	, m_print(false)
	, m_filled(false)
	, m_slicewin(false)
	, m_TimeLimsSet(false)
	{
		std::cout << "\n\n\t\t\t----PAY ATTENTION TO THIS----\n\n"
			<< "Lesson for TimeTool:: from TransAbs:: \t (here 'spectrum' means FFT of the images)\n"
			<< " - Implement the same weighting based on the reference amp and use to weight the averaging of FFT values\n" 
			<< " - Work with the difference images in the end and set the DC 0th term in FFT to zero\n"
			<< " - Ultimately we will only work with the slope of the unwrapped phase\n"
			<< " - Compute the noise spectrum from the single ref / accum ref and update running average of noise spectrum\n"
			<< " - Compute the signal spectrum from the single shot / accum ref and update running avarage of signal spectrum\n"
			<< " - Use the true single shot signal spectrum, but the running average of single shot reference spectra for Weiner\n"
			<< std::endl;
		init(bins_in,mins_in,maxs_in);
		//m_label_strings.assign({"position","amp","width"});
		/*
		std::string m_samplecodesname("samplecodes.out");
		m_samplecodesname += boost::lexical_cast<std::string>(m_rank);
		m_samplecodes = new std::ofstream(m_samplecodesname.c_str(),std::ios::out);
		(*m_samplecodes) << "# opeining file initially" << std::endl;
		*/
	}

	TimeTool::TimeTool(void)
	: m_use(false)
	, m_print(false)
	, m_filled(false)
	, m_TimeLimsSet(false)
	{
		//m_label_strings.assign({"position","amp","width"});
	}
	TimeTool::~TimeTool(void){
		if (m_outfile.is_open())
			m_outfile.close();
			/*
		if (m_samplecodes->is_open())
			m_samplecodes->close();
			*/
	}
	TimeTool & TimeTool::operator=( TimeTool rhs ){ // the compiler is making a pass by copy here //
		swap(*this,rhs); // using the swap makes this exception safe... 
		return *this;
	}
	TimeTool::TimeTool(const TimeTool & b)
	: m_use(b.m_use)
	, m_print(b.m_print)
	, m_filled(b.m_filled)
	, m_srcStr(b.m_srcStr)
	, m_srcPtr(b.m_srcPtr)
	, m_root_rank(b.m_root_rank)
	, m_rank(b.m_rank)
	, m_filename(b.m_filename)
	, m_TimeLimsSet(b.m_TimeLimsSet)
	{ // not so sure this is exception safe
		m_data = b.m_data;
		m_cam_data = b.m_cam_data;

		m_filename += ".copy";
		if (b.m_outfile.is_open())
			m_outfile.open(m_filename.c_str(),std::ios::out);
	}
	void TimeTool::deepcopy_data(const TimeTool & b)
	{
		m_data.resize(b.m_data.size());
		for ( unsigned i=0;i<b.m_data.size();++i){
			m_data[i].resize(b.m_data[i].size());
			for (unsigned j=0;j<b.m_data[i].size();++j)
				m_data[i][j] = b.m_data[i][j];
		}
	}
	void TimeTool::swap(TimeTool & a, TimeTool & b)
	{
		std::swap(a.m_use,b.m_use);
		std::swap(a.m_print,b.m_print);
		std::swap(a.m_filled,b.m_filled);
		std::swap(a.m_srcStr,b.m_srcStr);
		std::swap(a.m_srcPtr,b.m_srcPtr);
		std::swap(a.m_root_rank,b.m_root_rank);
		std::swap(a.m_rank,b.m_rank);
		std::swap(a.m_filename,b.m_filename);
		m_data.swap(b.m_data);
		m_cam_data.swap(b.m_cam_data);
	}


	// By adding an unused sourceStr into init, then we could virtual void these all
	bool TimeTool::init(std::vector<unsigned>& bins_in
			, std::vector<double>& mins_in
			, std::vector<double>& maxs_in)
	{
		// HERE //
		// Not sure what will happen when we don't give this a time-bins //
		if (! (bins_in.size() == mins_in.size() && bins_in.size() == maxs_in.size()) ){
			std::cerr << "Failed size matching in TimeTool::init() method" << std::endl;
			return false;
		}
		std::vector<unsigned>::iterator itr_u;
		itr_u = bins_in.begin();
		m_bins.assign(itr_u,bins_in.end());
		if (mins_in.size()>3)
			std::cout << "we need to make sure the calibration is read into TimeTool since using the time distribution" << std::endl;

		m_data.resize(m_bins.size());

		std::vector<double>::iterator itr_d;
		itr_d = mins_in.begin();
		m_data.at(min).assign(itr_d,mins_in.end());
		itr_d = maxs_in.begin();
		m_data.at(max).assign(itr_d,maxs_in.end());
		m_data.at(current).assign(m_data.at(min).size(),0.);
		m_filled = false;
		return true;
	}
	bool TimeTool::setTimeCalib(std::vector<double> & calib_in){
		// input vector should be (a,b,c,x0) values that turn x into time in picoseconds
		// change this, let it be (x0,a,b,c,d,...) and the length of the vector gives the number of powers for the fitting
		//
		// quart(x) = a + b*(x-x0) + c*(x-x0)**2 + d*(x-x0)**3 + e*(x-x0)**4 
		// These give phase slope in radians versus stage setting... phase =0 means left (blue) end of spectrum, phase = -2pi means right (red) side
		// In this case, stage drove the continuum stage
		// Large number is shorter whitelight path (I believe)
		//
		// f(x)=(x<-2.*pi ? (x+2.*pi) : x) * -1024./2./pi
		// plot file u (f($2)):(($1-93)/.3*2.),quart(x)
		// $1 is the stage setting (gets converted to delay)
		// $2 is the slope of the phase in radians, -2pi means right (red) side
		//
		// Final set of parameters            Asymptotic Standard Error
		// =======================            ==========================
		//
		// a               = 0.722173         +/- 0.002769     (0.3834%)
		// b               = -0.00200986      +/- 1.499e-05    (0.7457%)
		// c               = 1.5119e-06       +/- 2.942e-08    (1.946%)
		// d               = -1.11271e-09     +/- 9.956e-11    (8.948%)
		//
		// x0 = 512
		// so, delay in ps is T = a + b*(x-x0)**1 + c*(x-x0)**2 ...
		//
		m_calib.assign(calib_in.begin(),calib_in.end());
		setTimeLims();
	}
	double TimeTool::getTspan(void){
		return double(m_data.at(max).at(delay) - m_data.at(min).at(delay));
	}
	double TimeTool::getTspan(std::string & s){
		double d = getTspan();
		char buffer[10];
		sprintf(buffer,"%.4g",d);
		s = buffer;
		//s = boost::lexical_cast<std::string>(buffer);
		return d; 
	}
	unsigned TimeTool::getTimeBins(std::string & s){
		unsigned n = m_bins.at(delay);
		s = boost::lexical_cast<std::string>(m_bins.at(delay));
		return n;
	}
	std::string TimeTool::getTimeLims(void){
		std::string s("#");
		s += "\tbins = " + boost::lexical_cast<std::string>(m_bins.at(delay));
		s += "\tmin = " + boost::lexical_cast<std::string>(m_data.at(min).at(delay));
		s += "\tmax = " + boost::lexical_cast<std::string>(m_data.at(max).at(delay));
		s += "\n";
		return s;
	}
	void TimeTool::setTimeLims(void){
		std::vector<double> y(2,0.);
		y[0] = polynomial<double,double>(m_data.at(min).at(pos),m_calib);
		y[1] = polynomial<double,double>(m_data.at(max).at(pos),m_calib);
		std::sort(y.begin(),y.end());
		m_data.at(min).at(delay) = y[0];
		m_data.at(max).at(delay) = y[1];
		std::cout << "Time limits = " << y << std::endl;
		m_TimeLimsSet = true;
		return;
	}
	double TimeTool::pos2delay(void){
		m_data.at(current).at(delay) = polynomial<double,double>(m_data.at(current).at(pos),m_calib);
		return m_data.at(current).at(delay);
	}
	void TimeTool::srcStr(Source srcStr_in){
		m_srcStr = srcStr_in;
	}

	std::string TimeTool::get_features_labels(void)
	{
		return "pos\tamp\twidth\tdelay\t";
	}
	bool TimeTool::addfeatures(std::vector<float> & out)
	{
		for (unsigned i = 0; i< m_data.at(current).size(); ++i){
			out.push_back((float)m_data.at(current)[i]);
		}
		return true;
	}
	bool TimeTool::fill_image(Event& evt,Env& env)
	{
		m_srcPtr = evt.get(m_srcStr);
		if (!m_srcPtr){
			std::cerr << "couldn't get TransAbs image" << std::endl;
			return false; 
		}
		if (m_cam_data.size() != m_srcPtr->height() || m_cam_data.back().size() != m_srcPtr->width()){
			m_cam_data.resize(m_srcPtr->height());
			for (unsigned i=0;i<m_cam_data.size();++i)
				m_cam_data.at(i).resize(m_srcPtr->width(),0);
		}
		for (unsigned i=0;i<m_srcPtr->height();i++){
			for (unsigned j=0;j<m_srcPtr->width();j++){
				m_cam_data.at(i).at(j) =  m_srcPtr->data16()[i][j];
			}
		}
	}
	bool TimeTool::fill(Event& evt,Env& env){
		m_filled = false;
		try {
			//m_data.at(current).at(ps) = (double)(env.epicsStore().value("TTSPEC:FLTPOS_PS"));
			m_data.at(current).at(pos) = (double) env.epicsStore().value("TTSPEC:FLTPOS");
			m_data.at(current).at(width) = (double)(env.epicsStore().value("TTSPEC:FLTPOSFWHM"));
			m_data.at(current).at(amp) = (double)(env.epicsStore().value("TTSPEC:AMPL"));
			if (m_data.at(current).size()>3 && m_TimeLimsSet){
				pos2delay();
			}
		} catch (std::exception & e) {
			std::cerr << "Failed TimeTool::fill\t" << e.what() << std::endl;
			m_data.at(current).at(pos) = 0.0;
			m_data.at(current).at(width) = 0.0;
			m_data.at(current).at(amp) = 0.0;
			if (m_data.at(current).size()>3){
				pos2delay();
			}
			m_filled = true;
			return false;
		}
		m_filled = true;
		if (m_use_filter)
			return testvalid_surf();	
		return true;
	}
	bool TimeTool::fill(Event& evt,Env& env,const double gdin){
		fill(evt,env);
		return inslicewin(gdin);
	}

	bool TimeTool::testvalid(void){
		float a(0),b(0);
		a = 500.;
		b = -5.;
		return ((a*m_data.at(current).at(amp) + b) < m_data.at(current).at(width));
	}
	bool TimeTool::testvalid_surf(void){
		// HERE these are all based on hand fitting the surfaces for runs 170-172 TimeTool scatter plots
		//below this one
		//g(x,y)=(aa+bb*(x-x0))*y + ss*sin(ww*(1+cc*(x-x0))*(x-x0)) + qq*(x-x0)**2 + oo
		//above this one
		//h(x,y)=(aa2+bb2*(x-x0))*y + ss2*sin(ww*(1+cc*(x-x0))*(x-x0)) + qq2*(x-x0)**2 + oo2 + oo3/x
		double *x = &m_data.at(current).at(pos);
		double *y = &m_data.at(current).at(amp);
		const float aa(140.),bb(0.1),x0(540.),ss(5.),ww(0.042),cc(-6e-4),qq(1.0e-4),oo(30.);
		const float aa2(70.),bb2(0.05),ss2(7.5),qq2(3.33e-5),oo2(20.),oo3(2.e3);
		return (
				m_data.at(current).at(width) 
				< (aa+bb*(*x-x0))* (*y) + ss*sin(ww*(1+cc*(*x-x0))*(*x-x0)) + qq*pow( (*x-x0) , int(2)) + oo
				&&
				m_data.at(current).at(width) 
				> (aa2+bb2*(*x-x0))* (*y) + ss2*sin(ww*(1+cc*(*x-x0))*(*x-x0)) + qq2*pow((*x-x0),(int)2) + oo2 + oo3/(*x) 
			);

	}
	void TimeTool::setslicewin(std::vector<double> posin,
				std::vector<double> widthin,
				std::vector<double> amplin,
				std::vector<double> widthPamplin){
		m_slicewin = true;
		m_poswin = posin;
		m_fwhmwin = widthin;
		m_amplwin = amplin;
		m_fwhmPamplwin = widthPamplin;
	}
	bool TimeTool::inslicewin(const double gdin){
		bool result = (m_data.at(current).at(amp)/gdin > (1.-1./400*m_data.at(current).at(pos)));
		return (inslicewin() && result);
	}

	bool TimeTool::inslicewin(void){
		bool result = 
			(inwin(m_data.at(current).at(pos),m_poswin) &&
			inwin(m_data.at(current).at(width),m_fwhmwin) &&
			inwin(m_data.at(current).at(amp),m_amplwin) &&
			inwin(double(m_data.at(current).at(width))/m_data.at(current).at(amp),m_fwhmPamplwin)
			);
		return result;
	}

	void TimeTool::open_file(std::string & filename){
		m_filename = filename;
		m_outfile.open(m_filename.c_str(),std::ios::out);
	}
	void TimeTool::close_file(void){
		if (m_outfile.is_open())
			m_outfile.close();
	}
	bool TimeTool::print_out(const unsigned eventnum)
	{
		if (!m_outfile.is_open())
			m_outfile.open(m_filename.c_str(),std::ios::out);
		if (!m_outfile.is_open()){
			std::cerr << "TimeTool::print_out() failed because m_outfile.is_open() failed" << std::endl;
			return false;
		}
		m_outfile << eventnum << "\t";
		m_outfile << index(pos) << "\t" << index(delay) << "\t";
		for (unsigned i=0;i<m_data.at(current).size();++i){
			m_outfile << m_data.at(current).at(i) << "\t";
		}
		m_outfile << "\n";
		return true;
	}
	void TimeTool::print_header(void){
		//enum Var {pos, amp, width, delay};
		m_outfile << "#eventnum\tindex(pos)\tindex(delay)\tpos\tamp\twidth\tdelay";
		m_outfile << "\n";
	}
	bool TimeTool::isref(Event& evt,Env& env)
	{
		shared_ptr<Psana::EvrData::DataV3> srcPtr = evt.get(m_evr_src);
		ndarray<const Psana::EvrData::FIFOEvent, 1> eventList = srcPtr->fifoEvents();
		try{
			//(*m_samplecodes) << "eventList[ ] = ";
			for (unsigned i=0;i<eventList.size();++i){
			//	(*m_samplecodes) << eventList[i].eventCode() << "\t";
				if (eventList[i].eventCode() == m_refcode){
					//(*m_samplecodes) << "... found, returning true" << std::endl;
					return true;
				}
			}
			//(*m_samplecodes) << std::endl;
		} catch (std::exception & e){
			std::cerr << "failed to get Evr data in TransAbs::isref() method, e = " << e.what() << std::endl;
			return false;
		}
		return true;
	}
}

