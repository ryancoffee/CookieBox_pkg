#ifndef LEARNING_H
#define LEARNING_H

#include "PSEvt/EventId.h"
#include "psddl_psana/bld.ddl.h"
//#include "psddl_psana/bld.h" // this is wrong and fails...
#include "psana/Module.h"
#include <boost/shared_ptr.hpp>
#include <boost/multi_array.hpp> // for the slice version of fill() method
#include "ImgAlgos/GlobalMethods.h"
#include "MsgLogger/MsgLogger.h"

#include <opencv2/opencv.hpp>
#include <opencv2/opencv_modules.hpp> // these two somewhere define CV_PCA_DATA_AS_ROW
#include <opencv2/core/core.hpp> // this include the PCA algorithm


#include <vector>
#include <algorithm>
//#include <map> // map didn't help my badd pass of srcPtr
#include <list>
#include <fstream>	
#include "CookieBox_pkg/dataops.h"

namespace CookieBox_pkg {

class Learning 
{
	public:
		Learning(void);
		~Learning(void);
		Learning( const Learning & rhs);
		Learning & operator=( Learning rhs );
	
	private:	
		void swap(Learning & a, Learning & b);
		void deepcopy_data(const Learning & b);

	public:
		bool init(const unsigned nsamplesin, const unsigned nfeaturesin, const float variancein);
		bool init( const unsigned nfeaturesin, const float variancein);

		inline bool use(bool in){m_use = in;return m_use;}
		inline bool use(void){return m_use;}
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

		bool print_eigenfaces(cv::Size shapein);
		bool print_out(void);
		void print_header(void);
		void open_file( std::string & filename );
		inline std::string filename(void) {return m_filename;}
		void close_file( void );
		void evtput(Event& evt, Env& env);

		bool fill(std::vector <float> & in);
		bool fill_image(std::vector< std::vector < float > > & in);
		void set_features_labels(std::string & in_labels);

		inline unsigned nfeatures(void){return m_features.cols;}
		inline unsigned nsamples(void){return m_classifiers.rows;}

		bool run_pca(const float retainvariance ); // for setting that the variance is not significantly reduced
		bool run_pca(const unsigned ncomponents ); // for setting the number of principle components
		bool project_pca(record_f & in,record_f & out);
		record_f & project_pca(record_f & in);
		bool compare_pca(data_f & in);
		bool compare_pca_nonimage(data_f in,const unsigned nonimagefeatures);
		bool compare_pca_images(data_f & in,cv::Size shape,const unsigned nonimagefeatures);
		bool inspect_eigen_pca(void);

	protected:


	private:

		bool m_use,m_print;
		std::ofstream m_outfile;
		std::string m_filename;

		cv::Mat m_classifiers; 
		//std::vector< std::vector< float > > m_features;
		cv::Mat m_features, m_features_means;
		std::string m_features_labels;
		std::vector< cv::Mat > m_images; // later, consider flattening these images when reading in.
		unsigned m_pca_samples_limit;
		cv::Mat m_pca_features;
		cv::Mat m_eigenfeatures; // this will collect the results of the PCA eigen basis projection
		// eventually, we will use m_eigenfeatures to predict via KNN 
		// also, we could make a switch to flatten the images and add them to the featurs.

		cv::PCA m_pca;
		//cv::PCA* m_pcaPtr; // this didn't help... somehow there is an error when using CV_PCA_DATA_AS_ROW
		unsigned m_pca_nfeatures;
		float m_pca_retainvariance;

		// 	MPI related	//
		unsigned m_root_rank,m_rank;
		unsigned m_mpi_size;
};

} // namespace CookieBox_pkg

#endif 

