#include "CookieBox_pkg/Xtcav.h"
#include "CookieBox_pkg/CookieBox_mod.h"
#include <fstream>
#include <iostream>
#include <assert.h>
#include <iterator>
#include <algorithm>

#include <opencv2/core/core.hpp> // for eroding the mask image
#include <opencv2/imgproc/imgproc.hpp> // for filtering/bluring, and thresholding the mask image


using namespace CookieBox_pkg;

namespace CookieBox_pkg{
	Xtcav::Xtcav(void)
	: m_use(false)
	, m_print(false)
	, m_addimagesasfeatures(false)
	, m_compute_average(false)
	, m_average_nimages(0)
	, m_remove_dark(false)
	, m_threshold(10)
	{
		m_downsample_steps.resize(2,0);
	}

	Xtcav::~Xtcav(void){
		if (m_outfile.is_open())
			m_outfile.close();
		if (m_outfile_vh.is_open())
			m_outfile_vh.close();
		if (m_outfile_stats.is_open())
			m_outfile_stats.close();
	}
	Xtcav & Xtcav::operator=( Xtcav rhs ){ // the compiler is making a pass by copy here //
		swap(*this,rhs); // using the swap makes this exception safe... 
		return *this;
	}
	Xtcav::Xtcav(const Xtcav & b)
	: m_srcPtr(b.m_srcPtr)
	, m_use(b.m_use)
	, m_srcStr(b.m_srcStr)
	, m_print(b.m_print)
	, m_addimagesasfeatures(b.m_addimagesasfeatures)
	, m_filename(b.m_filename)
	, m_filename_vh(b.m_filename_vh)
	, m_filename_stats(b.m_filename_stats)
	, m_root_rank(b.m_root_rank)
	, m_rank(b.m_rank)
	, m_mpi_size(b.m_mpi_size)
	, m_compute_average(b.m_compute_average)
	, m_average_nimages(b.m_average_nimages)
	, m_remove_dark(b.m_remove_dark)
	, m_threshold(b.m_threshold)
	{ // not so sure this is exception safe

		m_valhist = b.m_valhist;
		m_valhist_lims = b.m_valhist_lims;
		m_data = b.m_data;
		m_data_stats = b.m_data_stats;

		m_filename += ".copy";
		m_filename_vh += ".copy";
		m_filename_stats += ".copy";
		if (b.m_outfile.is_open())
			m_outfile.open(m_filename.c_str(),std::ios::out);
		if (b.m_outfile_vh.is_open())
			m_outfile_vh.open(m_filename_vh.c_str(),std::ios::out);
		m_average = b.m_average.clone();
	}
	void Xtcav::deepcopy_data(const Xtcav & b)
	{
		m_data.resize(b.m_data.size());
		for ( unsigned i=0;i<b.m_data.size();++i){
			m_data[i].resize(b.m_data[i].size());
			for (unsigned j=0;j<b.m_data[i].size();++j)
				m_data[i][j] = b.m_data[i][j];
		}
	}
	void Xtcav::swap(Xtcav & a, Xtcav & b)
	{
		std::swap(a.m_use,b.m_use);
		std::swap(a.m_root_rank, b.m_root_rank);
		std::swap(a.m_rank, b.m_rank);
		std::swap(a.m_mpi_size, b.m_mpi_size);
		std::swap(a.m_srcPtr, b.m_srcPtr);
		std::swap(a.m_srcStr,b.m_srcStr);
		std::swap(a.m_addimagesasfeatures,b.m_addimagesasfeatures);
		std::swap(a.m_compute_average,b.m_compute_average);
		std::swap(a.m_average_nimages,b.m_average_nimages);
		std::swap(a.m_remove_dark,b.m_remove_dark);
		std::swap(a.m_threshold,b.m_threshold);
		
		a.m_valhist.swap(b.m_valhist);
		a.m_valhist_lims.swap(b.m_valhist_lims);
		a.m_data.swap(b.m_data);
		a.m_data_stats.swap(b.m_data_stats);

		cv::Mat tmp = a.m_average.clone();
		a.m_average = b.m_average.clone();
		b.m_average = tmp.clone();

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


		std::swap(a.m_filename_vh,b.m_filename_vh);
		if (a.m_outfile_vh.is_open() && b.m_outfile_vh.is_open()){
			a.m_outfile_vh.close();
			b.m_outfile_vh.close();
			a.m_outfile_vh.open(a.m_filename_vh.c_str(),std::ios::app);
			b.m_outfile_vh.open(b.m_filename_vh.c_str(),std::ios::app);
			return;
		}
		if (a.m_outfile_vh.is_open() && !b.m_outfile_vh.is_open()){
			a.m_outfile_vh.close();
			b.m_outfile_vh.open(b.m_filename_vh.c_str(),std::ios::app);
			return;
		}
		if (!a.m_outfile_vh.is_open() && b.m_outfile_vh.is_open()){
			b.m_outfile_vh.close();
			a.m_outfile_vh.open(a.m_filename_vh.c_str(),std::ios::app);
			return;
		}
		
		std::swap(a.m_filename_stats,b.m_filename_stats);
		if (a.m_outfile_stats.is_open() && b.m_outfile_stats.is_open()){
			a.m_outfile_stats.close();
			b.m_outfile_stats.close();
			a.m_outfile_stats.open(a.m_filename_stats.c_str(),std::ios::app);
			b.m_outfile_stats.open(b.m_filename_stats.c_str(),std::ios::app);
			return;
		}
		if (a.m_outfile_stats.is_open() && !b.m_outfile_stats.is_open()){
			a.m_outfile_stats.close();
			b.m_outfile_stats.open(b.m_filename_stats.c_str(),std::ios::app);
			return;
		}
		if (!a.m_outfile_stats.is_open() && b.m_outfile_stats.is_open()){
			b.m_outfile_stats.close();
			a.m_outfile_stats.open(a.m_filename_stats.c_str(),std::ios::app);
			return;
		}

	}


	bool Xtcav::init(std::vector<unsigned> downsamplesteps,unsigned in = 10)
	{
		m_threshold = in;
		m_downsample_steps.resize(2);
		std::vector<unsigned>::iterator it = downsamplesteps.begin();
		m_downsample_steps.assign(it,downsamplesteps.end());
		return true;
	}

	bool Xtcav::init(unsigned in = 10)
	{
		m_threshold = in;
		m_downsample_steps.resize(2);
		m_downsample_steps.at(0) = 8;
		m_downsample_steps.at(1) = 2;
		return true;
	}

	void Xtcav::srcStr(Source srcStr_in){
		m_srcStr = srcStr_in;
	}

	std::string Xtcav::get_features_labels(void)
	{
		std::string s("centroid_x\tcentroid_y\twidth_x\twidth_y\tdark_npix\tmid_npix\thigh_npix\t");
		if (m_addimagesasfeatures)
			s += "flat_images\t";
		return s;
	}

	bool Xtcav::addfeatures_centroids_image(std::vector<float> & out)
	{
		out.clear();
		out.push_back(centroid_x);
		out.push_back(centroid_y);
		m_nonimage_features = out.size();
		if (m_addimagesasfeatures && m_cropped_downsampled.total() > 0){
			out.reserve(m_nonimage_features + m_cropped_downsampled.total());
			for (unsigned i=0;i<m_cropped_downsampled.rows;++i){
				for (unsigned j=0;j<m_cropped_downsampled.cols;++j){
					out.push_back(m_cropped_downsampled.at<float>(i,j));
				}
			}
		}
		return true;
	}
	bool Xtcav::addfeatures(std::vector<float> & out)
	{
		if (m_data_stats.size() >0 ){
			for (unsigned i = 0; i< m_data_stats.size(); ++i){
				out.push_back((float)m_data_stats[i]);
			}
		}
		m_nonimage_features = out.size();
		if (m_addimagesasfeatures && m_cropped_downsampled.total() > 0){
			out.reserve(m_nonimage_features + m_cropped_downsampled.total());
			for (unsigned i=0;i<m_cropped_downsampled.rows;++i){
				for (unsigned j=0;j<m_cropped_downsampled.cols;++j){
					out.push_back(m_cropped_downsampled.at<float>(i,j));
				}
			}
		}
		return true;
	}

	bool Xtcav::fill(Event& evt, Env& env){
		m_srcPtr = evt.get(m_srcStr);
		if (!m_srcPtr){
			//std::cerr << "couldn't get xtcav image" << std::endl;
			return false; 
		}
		if (m_data.size() < m_srcPtr->height() || m_data.back().size() < m_srcPtr->width()){
			m_data.resize(m_srcPtr->height());
			for (unsigned i=0;i<m_data.size();++i)
				m_data.at(i).resize(m_srcPtr->width(),0);
		}
		//std::cerr << "m_dark_image.size() = " << m_dark_image.size() << " and m_data.size() = " << m_data.size() << std::endl;
		if (m_remove_dark && m_dark_image.size()==m_data.size() && m_dark_image.front().size() == m_data.front().size()){
			//std::cerr << "attempting to remove dark image in Xtcav::fill() method " << std::endl;
			for (unsigned i=0;i<m_srcPtr->height();i++){
				for (unsigned j=0;j<m_srcPtr->width();j++){
					int val = (int)m_srcPtr->data16()[i][j] - (int)m_dark_image.at(i).at(j);
					m_data.at(i).at(j) = ( val>0 ? val : 0);
				}
			}
		} else {
			for (unsigned i=0;i<m_srcPtr->height();i++){
				for (unsigned j=0;j<m_srcPtr->width();j++){
					m_data.at(i).at(j) =  m_srcPtr->data16()[i][j];
				}
			}
		}
		if (m_compute_average && !fillavg()){
			std::cerr << "failed to fill average in Xtcav::fill() method" << std::endl;
			return false;
		}
		return true;	
	}
	bool Xtcav::fillavg(void)
	{
		++m_average_nimages;
		if (m_average.rows != m_data.size() || m_average.cols != m_data[0].size()){
			m_average = cv::Mat::zeros(m_data.size(),m_data[0].size(),CV_32F);
		}
		for (unsigned i=0;i<m_average.rows;++i){
			for (unsigned j=0;j<m_average.cols;++j){
				m_average.at<float>(i,j) += (float)m_data[i][j];
			}
		}
		return true;
	}

	int Xtcav::identify(Event& evt, Env& env){
		/*
		   if (!m_library.exists())
		   return -1;
		   */
		return 0;
	}
	bool Xtcav::loadlibrary(std::ifstream& file)
	{
		return false;
	}


	void Xtcav::valhist_init(std::vector<unsigned> in ) //(configList"xt_valhist_list"))
	{
		std::vector<unsigned>::iterator itr = in.begin();
		m_valhist_lims.assign(itr,in.end());
	}
	bool Xtcav::runstats(void)
	{
		//enum StatsParam {centroid_x, centroid_y,width_x,width_y,dark_npix,mid_npix,high_npix};
		//
		// Be sure to do erosion and thresholding and background before calculating the centroid 
		
		if (m_data_stats.size() < 7)
			m_data_stats.assign(7,0);
		if (m_valhist.size() != m_valhist_lims[vh_bins])
			m_valhist.assign(m_valhist_lims[vh_bins],0);
		std::fill(m_data_stats.begin(),m_data_stats.end(),0);
		std::fill(m_valhist.begin(),m_valhist.end(),0);
		if (m_filtered_image.rows == m_data.size() && m_filtered_image.cols == m_data.front().size()){
			long double sum;
			sum = 0.;
			for (unsigned i=0;i<m_filtered_image.rows;++i){
				for (unsigned j=0;j<m_filtered_image.cols;++j){
					float val =  m_filtered_image.at<float>(i,j);
					if (val > 0.){
						sum += val;
						m_data_stats.at(centroid_x) += (long unsigned)((float)j * val);
						m_data_stats.at(centroid_y) += (long unsigned)((float)i * val);
						++m_data_stats.at(brightness(m_data[i][j]));
						vh_incr(m_data[i][j]);
					}
				}
			}
			if (sum == 0.){
				std::cerr << "in Xtcav::runstats() method, (long double) sum == 0.0" << std::endl;
				return false;
			}
			m_data_stats.at(centroid_x) /= (long unsigned)sum;
			m_data_stats.at(centroid_y) /= (long unsigned)sum;
			return true;
		}
		long long sum;
		sum = 0;
		for (unsigned i=0;i<m_data.size();++i){
			for (unsigned j=0;j<m_data[i].size();++j){
				unsigned val =  m_data[i][j];
				sum += val;
				m_data_stats.at(centroid_x) += j * val;
				m_data_stats.at(centroid_y) += i * val;
				++m_data_stats.at(brightness(val));
				vh_incr(val);
			}
		}
		if (sum == 0){
			std::cerr << "in Xtcav::runstats() method, (long long int)sum == 0" << std::endl;
			return false;
		}
		m_data_stats.at(centroid_x) /= sum;
		m_data_stats.at(centroid_y) /= sum;
		return true;
	}
	bool Xtcav::set_crop(std::vector<int> in)
	{
		// HERE crop is useful for dealing with the particular shape as an identifier in an eigen (or fisher) faces style algorithm
		m_crop_win = in;
	}
	bool Xtcav::crop_filtered(const unsigned eventnum)
	{
		// HERE crop is useful for dealing with the particular shape as an identifier in an eigen (or fisher) faces style algorithm
		if (m_crop_win.at(0) <=0 || m_crop_win.at(0) > m_filtered_image.cols)
			m_crop_win.at(0) = m_filtered_image.cols;
		if (m_crop_win.at(1) <=0 || m_crop_win.at(1) > m_filtered_image.rows)
			m_crop_win.at(1) = m_filtered_image.rows;
		int left = 0;
		int right = m_filtered_image.cols;
		int top = 0;
		int bottom = m_filtered_image.rows;
		left = m_data_stats.at(centroid_x) - m_crop_win.at(0)/2;
		right = m_data_stats.at(centroid_x) + m_crop_win.at(0)/2;
		top = m_data_stats.at(centroid_y) - m_crop_win.at(1)/2;
		bottom = m_data_stats.at(centroid_y) + m_crop_win.at(1)/2;
		try{
			m_cropped_downsampled = cv::Mat::zeros(
					m_crop_win.at(1)/m_downsample_steps.at(0)+1
					, m_crop_win.at(0)/m_downsample_steps.at(1)+1
					, CV_32F);
			for (int r = std::max(0,top); 
					r < std::min(bottom,(int)(m_cropped_downsampled.rows - m_downsample_steps.at(0))); 
					r += m_downsample_steps.at(0)){
				for (int c = std::max(0,left); 
						c < std::min(right,(int)(m_filtered_image.cols - m_downsample_steps.at(1))); 
						c += m_downsample_steps.at(1)){
					unsigned nelems;
					nelems = 0;
					for (unsigned i=r;i<std::min((int)(r+m_downsample_steps.at(0)),m_filtered_image.rows);++i){
						for (unsigned j=c;j<std::min((int)(c+m_downsample_steps.at(1)),m_filtered_image.cols);++j){
							m_cropped_downsampled.at<float>(r/m_downsample_steps.at(0),c/m_downsample_steps.at(1)) += (m_filtered_image.at<float>(i,j));
							++nelems;
						}
					}
					if (nelems > 1)
						m_cropped_downsampled.at<float>(r/m_downsample_steps.at(0),c/m_downsample_steps.at(1)) /= nelems;
				}
			}

		} catch (cv::Exception & e) {
			std::cerr << "Failed to crop image in Xtcav::crop_filtered() method for event " << eventnum ;
			std::cerr << " m_crop_win = " << m_crop_win.at(0) << ", " << m_crop_win.at(1) << std::endl;
			std::cerr << " with cv::Exception: " << e.msg << std::endl;
			return false;
		}

		return true;
	}

	bool Xtcav::print_cropped(const unsigned eventnum)
	{
		std::string filename(m_filename);
		filename += ".event_" + boost::lexical_cast<std::string>(eventnum);
		filename += ".filtered.cropped.downsampled";
		std::ofstream outfile(filename.c_str(),std::ios::out);
		if (!outfile.is_open()){
			std::cerr << "Failed to open file in Xtcav::print_cropped() method" << std::endl;
			return false;
		}
		outfile << "# cropped filtered downsampled image\n";
		for (unsigned i=0;i<m_cropped_downsampled.rows;++i){
			for (unsigned j=0; j< m_cropped_downsampled.cols; ++j){
				outfile << m_cropped_downsampled.at<float>(i,j) << "\t";
			}
			outfile << "\n";
		}
		outfile.close();
		return true;
	}
	void Xtcav::valhist(std::vector<unsigned> & out)
	{
		std::vector<unsigned>::iterator itr = m_valhist.begin();
		out.assign(itr,m_valhist.end());
	}

	// File interface methods //
	void Xtcav::open_file( std::string & filename ){
		m_filename = filename;
		m_outfile.open(filename.c_str(),std::ios::out);
	}
	void Xtcav::close_file(void){
		if (m_outfile.is_open())
			m_outfile.close();
	}
	void Xtcav::open_file_vh( std::string & filename ){
		m_filename_vh = filename;
		m_outfile_vh.open(m_filename_vh.c_str(),std::ios::out);
	}
	void Xtcav::open_file_stats( std::string & filename ){
		m_filename_stats = filename;
		m_outfile_stats.open(m_filename_stats.c_str(),std::ios::out);
	}
	void Xtcav::close_file_vh(void){
		if (m_outfile_vh.is_open())
			m_outfile_vh.close();
	}
	void Xtcav::close_file_stats(void){
		if (m_outfile_stats.is_open())
			m_outfile_stats.close();
	}
	bool Xtcav::print_out_vh(const unsigned eventnum)
	{
		std::string eventfilename(m_filename_vh);
		eventfilename += ".event_" + boost::lexical_cast<std::string>(eventnum);
		if (m_outfile_vh.is_open())
			m_outfile_vh.close();
		m_outfile_vh.open(eventfilename.c_str(),std::ios::out);
		if (!m_outfile_vh.is_open())
			return false;
		double x;
		for (unsigned i=0;i<m_valhist.size();++i){
			x = (double)i/m_valhist_lims[vh_bins] * (m_valhist_lims[vh_high] - m_valhist_lims[vh_low]);
			x += m_valhist_lims[vh_low];
			m_outfile_vh << x << "\t" << m_valhist[i] << "\n";
		}
		m_outfile_vh.close();
		return true;
	}
	bool Xtcav::print_out_vh(void)
	{
		if (!m_outfile_vh.is_open())
			m_outfile_vh.open(m_filename_vh.c_str(),std::ios::app);
		double x;
		if (!m_outfile_vh.is_open())
			return false;
		for (unsigned i=0;i<m_valhist.size();++i){
			x = m_valhist_lims[vh_low] + (m_valhist_lims[vh_high] - m_valhist_lims[vh_low]) / m_valhist_lims[vh_bins];
			m_outfile_vh << x << "\t" << m_valhist[i] << "\n";
		}
		return true;
	}
	bool Xtcav::print_header_stats(void){
		//enum StatsParam {centroid_x, centroid_y,width_x,width_y,dark_npix,mid_npix,high_npix};
		if (!m_outfile_stats.is_open())
			m_outfile_stats.open(m_filename_stats.c_str(),std::ios::out);
		if (!m_outfile_stats.is_open())
			return false;
	//	m_outfile_stats << "#xtcav stats\n";
		m_outfile_stats << "#{centroid_x,centroid_y,width_x,width_y,dark_npix,mid_npix,high_npix}\n";
		return true;
	}
	bool Xtcav::print_header_vh(void){
		m_outfile_vh << "xtcav valhist\n";
		m_outfile_vh << "x\tsum\n";
	}
	bool Xtcav::print_header(void){
		m_outfile << "#xtcav images";
		m_outfile << "\n";
	}
	bool Xtcav::print_out(void){
		if (!m_outfile.is_open())
			return false;
		for (unsigned i=0;i<m_data.size();++i){
			for (unsigned j=0;j<m_data[i].size();++j){
				m_outfile << m_data[i][j] << "\t";
			}
			m_outfile << "\n";
		}
		m_outfile << "\n";
		

		return true;
	}
	// HERE
	// need a method to read in the average of dark images.
	// Also, we should think about also writing and reading binary

	// HERE use this for loading the face recognizer once we have it.
	// void FaceRecognizer::load(const string& filename) (C++ function, in FaceRecognizer)
	bool Xtcav::dark_file(const std::string inname)
	{
		std::cerr << "in Xtcav::dark_file() method, loading dark_file " << inname << std::endl;
		std::ifstream infile(inname.c_str(),std::ios::in);
		if (!infile.is_open()){
			std::cerr << "load dark image failed in v::dark_file() method" << std::endl;
			return false;
		}
		std::vector<std::string> header;
		infile >> header;
		infile >> m_dark_image;
		infile.close();
		return true;

		/* // HERE this is a test to read directly from some image format... but this doesn't read as floats. 
		try {
			m_dark_image = imread(inname, 0); // 0== grayscale cv::CV_LOAD_IMAGE_GRAYSCALE );
		} catch (cv::Exception& err) {
			std::cerr << "load dark image failed in v::dark_file() method, err = " << err.msg << std::endl;
			return false;
		}
		return true;
		*/
	}
	bool Xtcav::print_out_avg(void)
	{
		if (m_average_nimages == 0){
			std::cerr << "m_average_nimages == 0, failing Xtcav::print_out_avg() method" << std::endl;
			return false;
		}
		std::string filename(m_filename);
		filename += ".average";
		std::cerr << "filename(m_filename) in Xtcav::print_out_avg() is " << filename << std::endl;
		std::ofstream averageout(filename.c_str(),std::ios::out);
		if (!averageout.is_open()){
			std::cerr << "Failed to open averageout file in Xtcav::print_avg() method" << std::endl;
			return false;
		}
		averageout << "# average image out\n";
		for (unsigned i=0;i<m_average.rows;++i){
			for (unsigned j=0;j<m_average.cols;++j){
				m_average.at<float>(i,j) /= (float) m_average_nimages;
				averageout << m_average.at<float>(i,j) << "\t";
			}
			averageout << "\n";
		}
		averageout.close();
		return true;
	}

	bool Xtcav::process_image(const unsigned eventnum)
	{
		// hard coding a reduction in pixels
		if (m_data.size() == 0){
			std::cerr << "Data did not get loaded before call to Xtcav::process_image() method" << std::endl;
			return false;
		}
		m_filtered_image = cv::Mat::zeros(m_data.size(),m_data[0].size(),CV_32F);
		for (unsigned i=0;i<m_filtered_image.rows;++i){
			for (unsigned j=0;j<m_filtered_image.cols;++j){
				m_filtered_image.at<float>(i,j) += m_data[i][j]; // explicit copy for sake of erosion operations
			}
		}
		cv::Mat kernel = cv::Mat::ones(3,3,CV_32F);
		cv::Point anchor(-1,-1);
		unsigned niterations = 1;
		double sigx = 3.;
		double sigy = 5.;
		cv::Size ksize(3,5);

		try {
		cv::erode(m_filtered_image,m_filtered_image, kernel, anchor, niterations) ;
		cv::threshold(m_filtered_image,m_filtered_image,(double)m_threshold,1000.,3); //3 = threshold to zero, not binary
		cv::dilate(m_filtered_image,m_filtered_image, kernel, anchor, niterations) ;
		// HERE the winning combination it seems so far is erode, threshold, dilate to get a sparse-ish representation of the original
		// without all the hot pixels
		//
		//cv::GaussianBlur(m_filtered_image,m_filtered_image,ksize,sigx,sigy,cv::BORDER_DEFAULT);
		// try thresholding lower since this is after the gaussian filter...
		// The lasing signal is not much above the dark_image subtracted background.
		// HERE is where we do the assignments and then dialate until we have a boundary, 
		// Still HERE...
		// then we use that boundary to make the current projections
		// remember, first we are after indentifying the correct reference
		// cv::dilate(m_filtered_image, m_filtered_image, kernel, anchor, niterations) ;
		} catch (cv::Exception & e) {
			std::cerr << "Failed to filter image in Xtcav::process_image() method, error: " << e.msg << std::endl;
			return false;
		}
		/*
		if (m_print && !print_filtered(eventnum)){
			return false;
		}
		*/
		return true;
	}
	bool Xtcav::print_debug(void)
	{
		std::cerr << "printing debug to stderr\n";
		std::cerr << "m_data.size() = " << m_data.size() << " records\t";
		std::cerr << "m_data.front().size() = " << m_data.front().size() << " record length\n";
		std::cerr << std::flush;
		return true;
	}
	bool Xtcav::print_filtered(const unsigned eventnum)
	{
		std::string filename(m_filename);
		filename += ".event_" + boost::lexical_cast<std::string>(eventnum) + ".filtered";
		std::ofstream filteredfile(filename.c_str(),std::ios::out);
		if (!filteredfile.is_open()){
			std::cerr << "failed to open file for filtered image in Xtcav::process_image() method" << std::endl;
			return false;
		}
		for (unsigned i=0;i<m_filtered_image.rows;++i){
			for (unsigned j=0;j<m_filtered_image.cols;++j){
				filteredfile << m_filtered_image.at<float>(i,j) << "\t";
			}
			filteredfile << "\n";
		}
		filteredfile << "\n";
		filteredfile.close();
		return true;
	}
	bool Xtcav::print_out(const unsigned eventnum)
	{
		std::string eventfilename(m_filename);
		eventfilename += ".event_" + boost::lexical_cast<std::string>(eventnum);
		if (m_outfile.is_open())
			m_outfile.close();
		m_outfile.open(eventfilename.c_str(),std::ios::out);
		if (!m_outfile.is_open())
			return false;
		for (unsigned i=0;i<m_data.size();++i){ 
			for (unsigned j=0;j<m_data[i].size();++j){
				m_outfile << m_data[i][j] << "\t";
			}
			m_outfile << "\n";
		}
		m_outfile << "\n";
		m_outfile.close();
		return true;

	}
	bool Xtcav::print_out_stats(const unsigned eventnum)
	{
		if (!m_outfile_stats.is_open())
			m_outfile_stats.open(m_filename_stats.c_str(),std::ios::app);
		if (!m_outfile_stats.is_open())
			return false;

		m_outfile_stats << eventnum << "\t";
		for (unsigned i=0;i<m_data_stats.size();++i)
			m_outfile_stats << m_data_stats[i] << "\t";
		m_outfile_stats << "\n";
		return true;
	}
}

