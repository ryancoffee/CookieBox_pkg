#include "CookieBox_pkg/Gdet.h"
#include "CookieBox_pkg/CookieBox_mod.h"
#include <fstream>
#include <iostream>
#include <assert.h>
#include <iterator>

using namespace CookieBox_pkg;

namespace CookieBox_pkg{
	Gdet::Gdet(std::vector<unsigned>& bins_in
			, std::vector<double>& mins_in
			, std::vector<double>& maxs_in)
	: m_use(false)
	, m_print(false)
	  , m_print_parallel(false)
	{
		init(bins_in, mins_in, maxs_in);
	}
	Gdet::Gdet(void)
		: m_use(false)
		, m_print(false)
	{
	}

	Gdet::~Gdet(void){
		if (m_outfile.is_open())
			m_outfile.close();
	}
	Gdet & Gdet::operator=( Gdet rhs ){ // the compiler is making a pass by copy here //
		swap(*this,rhs); // using the swap makes this exception safe... 
		return *this;
	}
	Gdet::Gdet(const Gdet & b)
	: m_use(b.m_use)
	, m_print(b.m_print)
	, m_print_parallel(b.m_print_parallel)
	, m_srcPtr(b.m_srcPtr)
	, m_srcStr(b.m_srcStr)
	, m_root_rank(b.m_root_rank)
	, m_rank(b.m_rank)
	, m_mpi_size(b.m_mpi_size)
	, m_filename(b.m_filename)
	{
		m_bins.resize(b.m_bins.size());
		for (unsigned i=0;i<m_bins.size();++i)
			m_bins[i] = b.m_bins[i];
		deepcopy_data(b);
		deepcopy_accumulation(b);
		if (m_print){
			m_filename += ".copy";
			if (b.m_outfile.is_open())
				m_outfile.open(m_filename.c_str(),std::ios::out);
		}
	}
	void Gdet::deepcopy_data(const Gdet & b)
	{
		m_data.resize(b.m_data.size());
		for ( unsigned i=0;i<b.m_data.size();++i){
			m_data[i].resize(b.m_data[i].size());
			for (unsigned j=0;j<b.m_data[i].size();++j)
				m_data[i][j] = b.m_data[i][j];
		}
	}
	void Gdet::deepcopy_accumulation(const Gdet & b)
	{
		m_accumulation.resize(b.m_accumulation.size());
		for ( unsigned i=0;i<b.m_accumulation.size();++i){
			m_accumulation[i].resize(b.m_accumulation[i].size());
			for (unsigned j=0;j<b.m_accumulation[i].size();++j)
				m_accumulation[i][j] = b.m_accumulation[i][j];
		}
	}
	void Gdet::swap(Gdet & a, Gdet & b)
	{
		std::swap(a.m_use,b.m_use);
		std::swap(a.m_print,b.m_print);
		std::swap(a.m_print_parallel,b.m_print_parallel);
		std::swap(a.m_srcPtr,b.m_srcPtr);
		std::swap(a.m_srcStr,b.m_srcStr);
		std::swap(a.m_root_rank,b.m_root_rank);
		std::swap(a.m_rank,b.m_rank);
		std::swap(a.m_mpi_size,b.m_mpi_size);

		m_bins.swap(b.m_bins);
		m_data.swap(b.m_data);
		m_accumulation.swap(b.m_accumulation);

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



	bool Gdet::init( std::vector<unsigned>& bins_in
			, std::vector<double>& mins_in
			, std::vector<double>& maxs_in)
	{
		// Add one move element to hold the shot average //
		assert( (bins_in.size() == mins_in.size()) & (bins_in.size() == maxs_in.size()) );
		m_bins.resize(bins_in.size());
		std::vector<unsigned>::iterator itr_u;
		itr_u = bins_in.begin();
		m_bins.assign(itr_u,bins_in.end());

		m_data.resize(m_bins.size()); // this is always 3 since it's for current, min and max values
		for (unsigned i=0;i<m_data.size();++i){
			m_data[i].resize(m_bins.size());
		}

		std::vector<double>::iterator itr_d;
		itr_d = mins_in.begin();
		m_data.at(min).assign(itr_d,mins_in.end());
		itr_d = maxs_in.begin();
		m_data.at(max).assign(itr_d,maxs_in.end());
		m_data.at(current).assign(m_data.at(min).size(),0.);

		m_data.at(min).push_back(globalmin());
		m_data.at(max).push_back(globalmax());
		m_data.at(current).push_back(0.);
		m_bins.push_back(m_bins.front()); // this is silly... maybe take the max number of bins.

		m_accumulation.clear();

		return true;
	}

	void Gdet::srcStr(Source srcStr_in){
		m_srcStr = srcStr_in;

	}

	std::string Gdet::get_features_labels(void)
	{
		return "gdet_11\tgdet_12\tgdet_21\tgdet_22\taverage\t";
	}
	bool Gdet::addfeatures(std::vector<float> & out)
	{
		for (unsigned i = 0; i< m_data.at(current).size(); ++i){
			out.push_back((float)m_data.at(current)[i]);
		}
		return true;
	}
	bool Gdet::fill(Event& evt, Env& env){
		m_srcPtr = evt.get(m_srcStr);
		if ( !m_srcPtr.get()) { 
			std::cerr << "Couldn't get m_srcPtr = evt.get(m_srcStr); in  Gdet::fill() method" << std::endl;
			return false;	
		}
		m_data.at(current).at(gdet_11) = m_srcPtr->f_11_ENRC();
		m_data.at(current).at(gdet_12) = m_srcPtr->f_12_ENRC();
		m_data.at(current).at(gdet_21) = m_srcPtr->f_21_ENRC();
		m_data.at(current).at(gdet_22) = m_srcPtr->f_22_ENRC();
		m_data.at(current).at(average) = shotAvg();

		/*
		std::cerr << "Gdet::m_data = ";
		for (unsigned i=0;i<m_data.at(current).size();++i){
			std::cerr << m_data.at(current)[i] << " ";
		}
		std::cerr << std::endl;
		*/

		return true;	
	}
	void Gdet::stats(void)
	{
		//compute and fill mins, maxes, means
		long unsigned shots = m_accumulation.size();
		for (unsigned i=0;i<m_data.at(current).size();++i){
			m_data[mean][i] *= ((double)(shots)/(double)(shots + 1));
			m_data[mean][i] += m_data[current][i]/(double)(shots + 1);
			if (m_data[min][i] > m_data[current][i]) { m_data[min][i] = m_data[current][i];}
			if (m_data[max][i] < m_data[current][i]) { m_data[max][i] = m_data[current][i];}
		}
		if (m_print){ // append to accumulation
			m_accumulation.push_back(m_data.at(current));
		}
	}


	double Gdet::globalmin(void){
		double val;
		val = m_data.at(min).front();
		for (unsigned i = 0;i<average;++i){
			if (val > m_data.at(min).at(i))
				val = m_data.at(min).at(i);
		}
		return val;
	}
	double Gdet::globalmax(void){
		double val;
		val = m_data.at(max).front();
		for (unsigned i = 0;i<average;++i){
			if (val < m_data.at(max).at(i))
				val = m_data.at(max).at(i);
		}
		return val;
	}
	void Gdet::open_file(std::string & filename){
		m_filename = filename;
		m_outfile.open(m_filename.c_str(),std::ios::out);
	}
	void Gdet::close_file(void){
		if (m_outfile.is_open())
			m_outfile.close();
	}
	bool Gdet::print_out(const unsigned eventnum)
	{
		if (!m_outfile.is_open())
			m_outfile.open(m_filename.c_str(),std::ios::out);
		if (!m_outfile.is_open()){
			std::cerr << "failing to open Gdet file" << std::endl;
			return false;
		}
		m_outfile << eventnum << "\t";
		for (unsigned i=0;i<m_data.at(current).size();++i){
			m_outfile << m_data.at(current).at(i) << "\t";
		}
		m_outfile << "\n";
		return true;
	}
	void Gdet::print_header(void){
		//enum Var {gdet_preatten1,gdet_preatten2,gdet_postatten1,gdet_postatten2,average};
		m_outfile << "#gdet_preatten1\tgdet_preatten2\tgdet_postatten1\tgdet_postatten2\taverage";
		m_outfile << "\n";
	}
}


