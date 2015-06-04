#include "CookieBox_pkg/Ebeam.h"
#include "CookieBox_pkg/CookieBox_mod.h"
#include <fstream>
#include <iostream>
#include <assert.h>
#include <iterator>

using namespace CookieBox_pkg;

namespace CookieBox_pkg{
	Ebeam::Ebeam( std::vector<unsigned>& bins_in
			, std::vector<double>& mins_in
			, std::vector<double>& maxs_in)
	: m_version(0)
	, m_use(false)
	, m_print(false)
	{
		init(bins_in,mins_in,maxs_in);
	}

	Ebeam::Ebeam(void)
	: m_version(0)
	, m_use(false)
	, m_print(false)
	{
//		m_label_strings = {"charge","energy","ltux","ltuy","ltuDx","ltuDy","pkcurr_bc2","energy_bc2","pkcurr_bc1","energy_bc1"};
	}
	// By adding an unused sourceStr into TimeTool::init(), then we could virtual these all
	//
	Ebeam::~Ebeam(void){
		if (m_outfile.is_open())
			m_outfile.close();
	}
	Ebeam & Ebeam::operator=( Ebeam rhs ){ // the compiler is making a pass by copy here //
		swap(*this,rhs); // using the swap makes this exception safe... 
		return *this;
	}
	Ebeam::Ebeam(const Ebeam & b)
	: m_use(b.m_use)
	, m_print(b.m_print)
	, m_srcPtrV6(b.m_srcPtrV6)
	, m_srcPtrV7(b.m_srcPtrV7)
	, m_srcStr(b.m_srcStr)
	, m_version(b.m_version)
	, m_root_rank(b.m_root_rank)
	, m_rank(b.m_rank)
	, m_mpi_size(b.m_mpi_size)
	{ // not so sure this is exception safe
		m_bins.resize(b.m_bins.size());
		for (unsigned i=0;i<b.m_bins.size();++i)
			m_bins[i] = b.m_bins[i];

		//deepcopy_data(b);
		this->m_data = b.m_data; 

		if (m_print){
			m_filename += ".copy";
			if (b.m_outfile.is_open())
				m_outfile.open(m_filename.c_str(),std::ios::out);
		}
	}
	void Ebeam::deepcopy_data(const Ebeam & b)
	{
		m_data.resize(b.m_data.size());
		for ( unsigned i=0;i<b.m_data.size();++i){
			m_data[i].resize(b.m_data[i].size());
			for (unsigned j=0;j<b.m_data[i].size();++j)
				m_data[i][j] = b.m_data[i][j];
		}
	}
	void Ebeam::swap(Ebeam & a, Ebeam & b)
	{
		std::swap(a.m_use,b.m_use);
		std::swap(a.m_srcPtrV6,b.m_srcPtrV6);
		std::swap(a.m_srcPtrV7,b.m_srcPtrV7);
		std::swap(a.m_srcStr,b.m_srcStr);
		std::swap(a.m_version,b.m_version);
		std::swap(a.m_root_rank,b.m_root_rank);
		std::swap(a.m_rank,b.m_rank);
		std::swap(a.m_mpi_size,b.m_mpi_size);

		m_data.swap(b.m_data);

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


	bool Ebeam::init( std::vector<unsigned>& bins_in
			, std::vector<double>& mins_in
			, std::vector<double>& maxs_in)
	{

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
		return true;
	}

	void Ebeam::srcStr(Source srcStr_in){
		m_srcStr = srcStr_in;

	}

	bool Ebeam::testVersion(Event& evt, Env& env){
		m_srcPtrV6 = evt.get(m_srcStr);
		if ( m_srcPtrV6.get() ){
			m_version = 6;
			return true;
		}
		
		m_srcPtrV7 = evt.get(m_srcStr);
		if ( m_srcPtrV7.get() ){
			m_version = 7;
			return true;
		}

		m_version = 0;
		return false;
	}
	
	std::string Ebeam::get_features_labels(void)
	{
		return "energy\tcharge\tltux\tltuy\tltuDx\tltuDy\tpkcurr_bc2\tenergy_bc2\tpkcurr_bc1\tenergy_bc1\t";
	}
	bool Ebeam::addfeatures(std::vector<float> & out)
	{
		for (unsigned i = 0; i< m_data.at(current).size(); ++i){
			out.push_back((float)m_data.at(current)[i]);
		}
	}
	bool Ebeam::fill(Event& evt, Env& env){
		if (m_version == 0) {
			if ( ! testVersion( evt, env) ) // this checks if untested, tests if not
				return false;
		}
		
		
		switch( m_version ){
			case 6:
				m_srcPtrV6 = evt.get(m_srcStr);
				m_data.at(current).at(energy)=m_srcPtrV6->ebeamL3Energy();
				m_data.at(current).at(charge)=m_srcPtrV6->ebeamCharge();
				m_data.at(current).at(ltux)=m_srcPtrV6->ebeamLTUPosX();
				m_data.at(current).at(ltuy)=m_srcPtrV6->ebeamLTUPosY();
				m_data.at(current).at(ltuDx)=m_srcPtrV6->ebeamLTUAngX();
				m_data.at(current).at(ltuDy)=m_srcPtrV6->ebeamLTUAngY();
				m_data.at(current).at(pkcurr_bc2)=m_srcPtrV6->ebeamPkCurrBC2();
				m_data.at(current).at(energy_bc2)=m_srcPtrV6->ebeamEnergyBC2();
				m_data.at(current).at(pkcurr_bc1)=m_srcPtrV6->ebeamPkCurrBC1();
				m_data.at(current).at(energy_bc1)=m_srcPtrV6->ebeamEnergyBC1();
				return true;	
			
			case 7:
				m_srcPtrV7 = evt.get(m_srcStr);
				m_data.at(current).at(energy)=m_srcPtrV7->ebeamL3Energy();
				m_data.at(current).at(charge)=m_srcPtrV7->ebeamCharge();
				m_data.at(current).at(ltux)=m_srcPtrV7->ebeamLTUPosX(); 
				m_data.at(current).at(ltuy)=m_srcPtrV7->ebeamLTUPosY();
				m_data.at(current).at(ltuDx)=m_srcPtrV7->ebeamLTUAngX();
				m_data.at(current).at(ltuDy)=m_srcPtrV7->ebeamLTUAngY();
				m_data.at(current).at(pkcurr_bc2)=m_srcPtrV7->ebeamPkCurrBC2();
				m_data.at(current).at(energy_bc2)=m_srcPtrV7->ebeamEnergyBC2();
				m_data.at(current).at(pkcurr_bc1)=m_srcPtrV7->ebeamPkCurrBC1();
				m_data.at(current).at(energy_bc1)=m_srcPtrV7->ebeamEnergyBC1();
				m_version = 7;
				return true;
			
			default:
				std::cerr << "Couldn't get Ebeam::fill() method" << std::endl;
				
				return testVersion(evt, env); // this checks if untested, tests if not
		}
		return false; // shouldn't actually ever fall here.

	}

	// File interface methods //
	void Ebeam::open_file( std::string & filename ){
		m_outfile.open(filename.c_str(),std::ios::out);
	}
	void Ebeam::close_file(void){
		if (m_outfile.is_open())
			m_outfile.close();
	}
	bool Ebeam::print_out(const unsigned eventnum)
	{
		if (!m_outfile.is_open())
			m_outfile.open(m_filename.c_str(),std::ios::out);
		if (!m_outfile.is_open())
			return false;
		m_outfile << eventnum << "\t";
		for (unsigned i=0;i<m_data.at(current).size();++i){
			m_outfile << m_data.at(current).at(i) << "\t";
		}
		m_outfile << "\n";
		return true;
	}
	void Ebeam::print_header(void){
		//enum Var {energy,charge,ltux,ltuy,ltuDx,ltuDy,pkcurr_bc2,energy_bc2,pkcurr_bc1,energy_bc1};
		m_outfile << "#energy\tcharge\tltux\tltuy\tltuDx\tltuDy\tpkcurr_bc2\tenergy_bc2\tpkcurr_bc1\tenergy_bc1";
		m_outfile << "\n";
	}
}

