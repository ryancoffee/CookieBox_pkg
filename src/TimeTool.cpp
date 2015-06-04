#include "CookieBox_pkg/TimeTool.h"
#include "CookieBox_pkg/CookieBox_mod.h"
#include <fstream>
#include <iostream>
#include <assert.h>
#include <iterator>

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
	{
		init(bins_in,mins_in,maxs_in);
		//m_label_strings.assign({"position","amp","width"});
	}
	TimeTool::TimeTool(void)
	: m_use(false)
	, m_print(false)
	{
		//m_label_strings.assign({"position","amp","width"});
	}
	TimeTool::~TimeTool(void){
		if (m_outfile.is_open())
			m_outfile.close();
	}
	TimeTool & TimeTool::operator=( TimeTool rhs ){ // the compiler is making a pass by copy here //
		swap(*this,rhs); // using the swap makes this exception safe... 
		return *this;
	}
	TimeTool::TimeTool(const TimeTool & b)
	: m_use(b.m_use)
	, m_print(b.m_print)
	, m_srcStr(b.m_srcStr)
	, m_srcPtr(b.m_srcPtr)
	, m_root_rank(b.m_root_rank)
	, m_rank(b.m_rank)
	, m_filename(b.m_filename)
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
		if (! (bins_in.size() == mins_in.size() && bins_in.size() == maxs_in.size()) ){
			std::cerr << "Failed size matching in TimeTool::init() method" << std::endl;
			return false;
		}
		std::vector<unsigned>::iterator itr_u;
		itr_u = bins_in.begin();
		m_bins.assign(itr_u,bins_in.end());

		m_data.resize(m_bins.size());

		std::vector<double>::iterator itr_d;
		itr_d = mins_in.begin();
		m_data.at(min).assign(itr_d,mins_in.end());
		itr_d = maxs_in.begin();
		m_data.at(max).assign(itr_d,maxs_in.end());

		m_data.at(current).assign(m_data.at(min).size(),0.);
		return true;
	}
	void TimeTool::srcStr(Source srcStr_in){
		m_srcStr = srcStr_in;
	}

	std::string TimeTool::get_features_labels(void)
	{
		return "pos\tamp\twidth\t";
	}
	bool TimeTool::addfeatures(std::vector<float> & out)
	{
		for (unsigned i = 0; i< m_data.at(current).size(); ++i){
			out.push_back((float)m_data.at(current)[i]);
		}
		return true;
	}
	bool TimeTool::fill(Event& evt,Env& env){
		//m_data.at(current).at(ps) = (double)(env.epicsStore().value("TTSPEC:FLTPOS_PS"));
		m_data.at(current).at(pos) = (double) env.epicsStore().value("TTSPEC:FLTPOS");
		m_data.at(current).at(width) = (double)(env.epicsStore().value("TTSPEC:FLTPOSFWHM"));
		m_data.at(current).at(amp) = (double)(env.epicsStore().value("TTSPEC:AMPL"));
		return true;	
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
		for (unsigned i=0;i<m_data.at(current).size();++i){
			m_outfile << m_data.at(current).at(i) << "\t";
		}
		m_outfile << "\n";
		return true;
	}
	void TimeTool::print_header(void){
		//enum Var {pos, amp, width};
		m_outfile << "#pos\tamp\twidth";
		m_outfile << "\n";
	}
}

