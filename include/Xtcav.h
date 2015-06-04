#ifndef XTCAV_H
#define XTCAV_H

#include "PSEvt/EventId.h"
#include "psana/Module.h"
#include <boost/shared_ptr.hpp>
#include "ImgAlgos/GlobalMethods.h"
#include "MsgLogger/MsgLogger.h"

#include "psddl_psana/camera.ddl.h"

#include <opencv2/core/core.hpp> // for types like cv::Mat

#include <vector>
#include <algorithm>
#include <list>
#include <fstream>	
#include <iostream>
#include <string>

#include "CookieBox_pkg/dataops.h"

namespace CookieBox_pkg 
{

class Xtcav
{

	public:
	enum ValHistLims {vh_low,vh_high,vh_bins};
	enum StatsParam {centroid_x, centroid_y,width_x,width_y,dark_npix,mid_npix,high_npix};

	Xtcav(void);
	~Xtcav(void);
	Xtcav & operator=( Xtcav rhs );
	Xtcav(const Xtcav & b);

	private:
		void deepcopy_data(const Xtcav & b);
		void swap(Xtcav & a, Xtcav & b);

	public:
		inline bool use(bool in){m_use = in;return m_use;}
		inline bool use(void){return m_use;}
		inline bool compute_avg(void){ return m_compute_average;}
		inline bool compute_avg(bool in){ 
			if (!m_use){
				m_print_avg = false;
				m_compute_average = false;
				return false;
			}
			m_compute_average = in;
			return m_compute_average;
		}
		inline bool print_avg(void){ return m_print_avg; }
		inline bool print_avg(bool in){
			if (!m_use){
				m_print_avg = false;
				m_compute_average = false;
				return false;
			}
			m_print_avg = in;
			m_compute_average = in;
			return in;
		}
		inline bool print(bool in){
			if (!m_use)
				return false;
			m_print = in;
			return m_print;
		}
		inline bool print(void)
		{
			if (!m_use)
				return false;
			return m_print;
		}

		bool init(std::vector<unsigned> downsamplesteps,unsigned in);
		bool init(unsigned in);
		bool set_crop(std::vector<int> in);

		void srcStr(Source srcStr_in);

	private:bool fillavg(void);
	public:	bool fill(Event& evt, Env& env);
		bool addfeatures(std::vector<float> & out);
		bool addfeatures_centroids_image(std::vector<float> & out);
		std::string get_features_labels(void);
		
		// anticipating learning //
		int identify(Event& evt, Env& env);
		bool loadlibrary(std::ifstream& file);
		bool crop_filtered(const unsigned eventnum);
		std::vector<int> m_crop_win;


		bool runstats(void);
		void valhist(std::vector<unsigned> & out);
		void valhist_init(std::vector<unsigned> in ); //(configList"xt_valhist_list"))
		inline void addimages(bool in = false){m_addimagesasfeatures = in;}

		inline unsigned threshold(void){return m_threshold;}
		inline unsigned threshold(unsigned in){m_threshold = in;return m_threshold;}
		inline bool remove_dark(void){return m_remove_dark;}
		inline bool remove_dark(bool in){m_remove_dark = in;return m_remove_dark;}
		bool dark_file(const std::string inname);
		data_f m_dark_image;
		inline cv::Size eigenface_shape(void){return m_cropped_downsampled.size();}
		inline unsigned nonimage_features(void){return m_nonimage_features;}

		bool process_image(const unsigned eventnum);
		bool print_out(void);
		bool print_out(const unsigned eventnum);
		bool print_filtered(const unsigned eventnum);
		bool print_cropped(const unsigned eventnum);
		bool print_out_avg(void);
		bool print_header(void);
		bool print_out_vh(void);
		bool print_out_vh(const unsigned eventnum);
		bool print_out_stats(const unsigned eventnum);
		bool print_header_vh(void);
		bool print_header_stats(void);
		void open_file( std::string & filename );
		void close_file(void);
		void open_file_vh( std::string & filename );
		void open_file_stats( std::string & filename );
		void close_file_stats(void);
		void close_file_vh(void);

		Source m_srcStr;

	private:
		bool m_use,m_print,m_print_avg;
		std::string m_filename,m_filename_vh,m_filename_stats;
		std::ofstream m_outfile,m_outfile_vh,m_outfile_stats;

		std::vector<unsigned> m_valhist;
		std::vector<unsigned> m_valhist_lims; // from configList() this gets low high nbins

		bool m_remove_dark;
		unsigned m_threshold;

		std::vector< std::vector<unsigned> > m_data;
		std::vector< long long > m_data_stats;
		cv::Mat m_average;
		cv::Mat m_filtered_image;
		cv::Mat m_cropped_image;
		cv::Mat m_cropped_downsampled;
		std::vector<unsigned> m_downsample_steps;
		unsigned m_nonimage_features;

		unsigned m_average_nimages;
		bool m_compute_average;

		bool m_addimagesasfeatures;

		boost::shared_ptr<Psana::Camera::FrameV1> m_srcPtr;

		// 	MPI related	//
		unsigned m_root_rank,m_rank;
		unsigned m_mpi_size;

		inline StatsParam brightness(long long in)
		{
			if (in<300)
				return dark_npix;
			if (in>1500)
				return high_npix;
			return mid_npix;
		}

		inline unsigned vh_incr(unsigned in) // inlined for speed since occurs for every pixel
		{
			int ind;
			ind = (int)( (double) m_valhist_lims[vh_bins] * (in - m_valhist_lims[vh_low])
					/ (double)(m_valhist_lims[vh_high] - m_valhist_lims[vh_low]) );
			if (ind < 0){
				++m_valhist.front();
				return 0;
			}
			if (ind > (int)(m_valhist_lims[vh_bins]-1) ){
				++m_valhist.back();
				return (m_valhist_lims[vh_bins]-1);
			}
			++m_valhist[ind];
			return (unsigned)ind;
		}
};
}


#endif
