#include "CookieBox_pkg/Learning.h"
#include <sstream>
#include "CookieBox_pkg/CookieBox_mod.h"
#include <assert.h>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <math.h>
#include <exception>

using namespace CookieBox_pkg;

namespace CookieBox_pkg 
{

	Learning::Learning(void)
	: m_use(false)
	, m_print(false)
	, m_read_pca(false)
	, m_write_pca(false)
	{
	}

	Learning::~Learning(void)
	{
		if (m_outfile.is_open())
			m_outfile.close();
	}
	Learning & Learning::operator=( Learning rhs ){ // the compiler is making a pass by copy here //
		swap(*this,rhs); // using the swap makes this exception safe... 
		return *this;
	}

	Learning::Learning(const Learning & b)
	: m_use(b.m_use)
	, m_print(b.m_print)
	, m_filename(b.m_filename)
	, m_read_pca(b.m_read_pca)
	, m_write_pca(b.m_write_pca)
	{ // not so sure this is exception safe
		m_features = b.m_features.clone(); 
		m_classifiers = b.m_classifiers.clone(); 

		m_images = b.m_images;
		m_pca_features = b.m_pca_features;
		m_eigenfeatures = b.m_eigenfeatures; 
		m_pca = b.m_pca; // I'm sure this will fail.
		m_pca= b.m_pca; // I'm sure this will fail.

		if (m_print){
			m_filename += ".copy";
			if (b.m_outfile.is_open())
				m_outfile.open(m_filename.c_str(),std::ios::out);
		}
	}
	void Learning::deepcopy_data(const Learning & b)
	{
		m_features = b.m_features.clone();
		/*
		m_features.resize(b.m_features.size());
		for ( unsigned i=0;i<b.m_features.size();++i){
			m_features[i].resize(b.m_features[i].size());
			for (unsigned j=0;j<b.m_features[i].size();++j)
				m_features[i][j] = b.m_features[i][j];
		}
		*/
	}
	void Learning::swap(Learning & a, Learning & b){
		std::swap(a.m_use,b.m_use);
		std::swap(a.m_rank,b.m_rank);
		std::swap(a.m_mpi_size,b.m_mpi_size);



		a.m_images.swap(b.m_images);

		cv::Mat tmp;
		tmp = a.m_features.clone();
		b.m_features.copyTo(a.m_features);
		tmp.copyTo(b.m_features);
		tmp.release();
		tmp = a.m_pca_features.clone();
		b.m_pca_features.copyTo(a.m_pca_features);
		tmp.copyTo(b.m_pca_features);
		tmp.release();
		tmp = a.m_eigenfeatures.clone();
		b.m_eigenfeatures.copyTo(a.m_eigenfeatures);
		tmp.copyTo(b.m_eigenfeatures); 
		tmp.release();
		tmp = a.m_classifiers.clone();
		b.m_classifiers.copyTo(a.m_classifiers);
		tmp.copyTo(b.m_classifiers); 
		tmp.release();


		cv::PCA tmppca(a.m_pca);
		a.m_pca = b.m_pca; // I'm sure this will fail.
		b.m_pca = tmppca;

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


	bool Learning::init(const std::string & filename)
	{
		m_read_pca = true;
		m_write_pca = false;
		if (read_pca( filename )){
			m_pca_nfeatures = m_pca.eigenvectors.cols; 
		} else {
			std::cerr << "Failed to read file in Learning::init(const std::string & filename) method" << std::endl;
			m_read_pca = false;
			return false;
		}
		return true;
	}
	bool Learning::init(const unsigned nsamplesin, const unsigned nfeaturesin = 25, const float variancein = 0.9)
	{
		m_read_pca = false;
		m_write_pca = true;
		m_pca_samples_limit = nsamplesin;
		m_pca_nfeatures = nfeaturesin;
		m_pca_retainvariance = variancein;
		return true;
	}
	bool Learning::init( const unsigned nfeaturesin = 25,  const float variancein = 0.9)
	{
		m_pca_nfeatures = nfeaturesin;
		m_pca_retainvariance = variancein;
		return true;
	}

	void Learning::evtput(Event& evt, Env& env)
	{
		// use something like this just to pass data to python

		std::string outkey("features");

		unsigned shape[] = {m_features.rows, m_features.cols};
		boost::shared_ptr< ndarray<double,2> > outmatPtr = boost::make_shared < ndarray<double,2> > ( shape );

		for (unsigned i = 0; i < m_features.rows;++i){ 
			for (unsigned j=0;j<m_features.cols;++j){
				(*outmatPtr)[i][j] = m_features.at<float>(i,j);
			}
		}
		shape[1] = 1;
		boost::shared_ptr< ndarray<double,1> > outvecPtr = boost::make_shared < ndarray<double,1> > ( shape );
		evt.put(outmatPtr,outkey);
		std::string outkey_classifiers("classifiers");
		for (unsigned i=0;i<m_classifiers.rows;++i){
			(*outvecPtr)[i] = m_classifiers.at<float>(i);
		}
	}

	bool Learning::write_pca(const std::string & filename )
	{       
		if (!m_write_pca)
			return false;
		try{
			cv::FileStorage fs(filename,cv::FileStorage::WRITE);
			fs << "mean" << m_pca.mean;
			fs << "eigenvectors" << m_pca.eigenvectors;
			fs << "eigenvalues" << m_pca.eigenvalues;
			fs.release();
		} catch (cv::Exception &e){
			std::cerr << "Failed in PCA::write() method, error: " << e.msg << std::endl;
			return false;
		}
		return true;
	}
	bool Learning::read_pca(const std::string & filename )
	{       
		if (!m_read_pca)
			return false;
		try{
			cv::FileStorage fs(filename,cv::FileStorage::READ);
			fs["mean"] >> m_pca.mean ;
			fs["eigenvectors"] >> m_pca.eigenvectors ;
			fs["eigenvalues"] >> m_pca.eigenvalues ;
			fs.release();
		} catch (cv::Exception & e){
			std::cerr << "Failed in PCA::write() method, error: " << e.msg << std::endl;
			return false;
		}
		return true;

	}
	
	// need a rerun_pca() method since now we check to see if eigenvecs are filled to test if we need to compute //
	bool Learning::run_pca(const unsigned ncomponents )
	{
		if (m_pca.eigenvectors.rows != 0){
			//std::cout << "skipping Learning::run_pca() since m_pca.eigenvectors.rows != 0 indicates already computed" << std::endl;
			return true;
		}
		std::cout << "entered Learning::run_pca() method, first time (I think)..." << std::endl;
		time_t timer,start;
		time(&start);

		std::cerr << "m_features.rows = " << m_features.rows << "\n";
		std::cerr << "m_features.cols = " << m_features.cols << std::endl;
		std::cerr << "m_features_means.rows = " << m_features_means.rows << std::endl;
		std::cerr << "m_features_means.cols = " << m_features_means.cols << std::endl;
		try {
			//std::cerr << "CV_PCA_DATA_AS_ROW = " << CV_PCA_DATA_AS_ROW << "\nCV_PCA_DATA_AS_COL = " << CV_PCA_DATA_AS_COL << std::endl;
			//m_pca(m_pca_features, cv::Mat(), CV_PCA_DATA_AS_ROW, ncomponents); // this is weird  why no cv::PCA::DATA_AS_ROW
			//m_pcaPtr = new cv::PCA(m_pca_features, cv::Mat(), CV_PCA_DATA_AS_ROW, ncomponents); // this is weird  why no cv::PCA::DATA_AS_ROW
			m_pca(m_features, m_features_means,CV_PCA_DATA_AS_ROW ,ncomponents);
			std::cout << "Just after m_pca() ctor call in Learning::run_pca() method"
				<< "\n\t m_features_means.rows = " << m_features_means.rows 
				<< "\n\t m_features_means.cols = " << m_features_means.cols
				<< "\n\t m_pca.eigenvectors.rows = " << m_pca.eigenvectors.rows 
				<< "\n\t m_pca.eigenvectors.cols = " << m_pca.eigenvectors.cols
				<< "\n\t m_pca.mean.rows = " << m_pca.mean.rows 
				<< "\n\t m_pca.mean.cols = " << m_pca.mean.cols
				<< std::endl;
			std::cout << "m_features_means.rows = " << m_features_means.rows << std::endl;
			std::cout << "m_features_means.cols = " << m_features_means.cols << std::endl;
		} catch (cv::Exception& err) {
			std::cerr << "build m_pca failed in Learning::run_pca() method, err = " << err.msg << std::endl;
			return false;
		}
		std::cout << "leaving Learning::run_pca() method" << std::endl;
		time(&timer);
		float pca_time = difftime(timer,start);
		std::cout << "leaving Learning::run_pca() method after pca_time = " << pca_time << "seconds" << std::endl;
		return true;
	}
	bool Learning::run_pca(float retainvariance = 0.95 )
	{
		time_t timer,start;
		time(&start);
		std::cout << "\nentered Learning::run_pca() method" << std::endl;
		if (retainvariance < 0.){
			std::cerr << "bad retain variance in Learning::run_pca() method" << std::endl;
			return false;
		}
		if (retainvariance > 1.)
			retainvariance = 1/retainvariance;
		if (retainvariance < 0.5)
			retainvariance = 1. - retainvariance;

		/*
		if( !setfeatures_pca()){
			std::cerr << "failed ot Learning::setfeatures_pca() in Learning::run_pca() method" << std::endl;
			return false;
		}
		*/

		try {
			//m_pca(m_pca_features, cv::Mat(), CV_PCA_DATA_AS_ROW, retainvariance); // this is weird  why no cv::PCA::DATA_AS_ROW
			//m_pcaPtr = new cv::PCA(m_pca_features, cv::Mat(), CV_PCA_DATA_AS_ROW, retainvariance); // this is weird  why no cv::PCA::DATA_AS_ROW
			m_pca(m_features, m_features_means,CV_PCA_DATA_AS_ROW, retainvariance); 
			std::cerr << "Just after m_pca() ctor call in Learning::run_pca() method"
				<< "\n\t m_features_means.rows = " << m_features_means.rows 
				<< "\n\t m_features_means.cols = " << m_features_means.cols
				<< "\n\t m_pca.eigenvectors.rows = " << m_pca.eigenvectors.rows 
				<< "\n\t m_pca.eigenvectors.cols = " << m_pca.eigenvectors.cols
				<< "\n\t m_pca.mean.rows = " << m_pca.mean.rows 
				<< "\n\t m_pca.mean.cols = " << m_pca.mean.cols
				<< std::endl;
		} catch (cv::Exception& err) {
			std::cerr << "build m_pca failed in Learning::run_pca() method, err = " << err.msg << std::endl;
			return false;
		}
		std::cout << "ran PCA with Learning::run_pca(float retainvariance = 0.95 ), m_pcaPtr->eigenvalues.rows = " << m_pca.eigenvalues.rows << std::endl;

		time(&timer);
		float pca_time = difftime(timer,start);
		std::cout << "leaving Learning::run_pca() method after pca_time = " << pca_time << "seconds" << std::endl;
		return true;
	}

	std::vector<float> & Learning::project_pca(std::vector<float> & in)
	{
		std::vector<float> out;
		cv::Mat point;  
		try {           
			point = m_pca.project(cv::Mat(in).t());
			//point = m_pca.mean;
			point.copyTo(out);
		} catch (cv::Exception& err){
			std::cerr << "projection failed in Learning::project_pca() method, err: " << err.msg 
				<< "\n\t m_features_means.rows = " << m_features_means.rows 
				<< "\n\t m_features_means.cols = " << m_features_means.cols
				<< "\n\t m_pca.eigenvectors.rows = " << m_pca.eigenvectors.rows 
				<< "\n\t m_pca.eigenvectors.cols = " << m_pca.eigenvectors.cols
				<<  std::endl;
			out.clear();
			return out;
		}       

		return out;
	}
	bool Learning::project_pca(std::vector<float> & in,std::vector<float> & out)
	{
		cv::Mat point;
		try {
			// HERE is the error and it never works nor writes the pca.mean...
			point = m_pca.project(cv::Mat(in).t());
			//point = m_pca.mean;
			point.copyTo(out);
		} catch (cv::Exception& err){
			std::cerr << "projection failed in Learning::project_pca() method, err: " << err.msg 
				<< "\n\t m_features_means.rows = " << m_features_means.rows 
				<< "\n\t m_features_means.cols = " << m_features_means.cols
				<< "\n\t m_pca.eigenvectors.rows = " << m_pca.eigenvectors.rows 
				<< "\n\t m_pca.eigenvectors.cols = " << m_pca.eigenvectors.cols
				<<  std::endl;
			return false;
		}
		return true;
	}

	bool Learning::compare_pca_nonimage(std::vector< std::vector< float > > in,const unsigned nonimagefeatures)
	{
		cv::Mat resultmat;
		cv::Mat inproj;
			std::string infilename(m_filename);
			std::string outfilename(m_filename);
			infilename += ".testpcanonimage";
			outfilename += ".testpcanonimage.result";
			try {
				std::ofstream innonimage(infilename.c_str(),std::ios::out);
				std::ofstream outnonimage(outfilename.c_str(),std::ios::out);
				innonimage << "#pca input nonimage, PCA on the image features too, cols correspond to imagename numbering\n";
				outnonimage << "#pca result nonimage, PCA on the image features too, cols correspond to imagename numbering\n";
				inproj = m_pca.project(cv::Mat(in));
				resultmat = m_pca.backProject(inproj);
				for (unsigned j=0;j<inproj.cols;++j){
					for (unsigned i=0;i<nonimagefeatures;++i){
						innonimage << inproj.at<float>(i,j) << "\t";
						outnonimage << resultmat.at<float>(i,j) << "\t";
					}
					innonimage << "\n";
					outnonimage << "\n";
				}
				innonimage.close();
				outnonimage.close();
			} catch (cv::Exception& err){
				std::cerr << "projection failed in Learning::compare_pca_nonimage() method, err: " << err.msg <<  std::endl;
				return false;
			} catch (std::exception & e) {
				std::cerr << "projection failed in Learning::compare_pca_nonimage() method, err: " << e.what() << std::endl;
				return false;
			}
		return true;
	}
	// HERE HERE HERE HERE //
	// OK make a new method since we are doing this anyway row by row, might as well do it for each event //
	// and name the images according to MPI rank //
	bool Learning::difference_pca_images(std::vector<float> & in,cv::Size shape,const unsigned nonimagefeatures,const unsigned event){
		cv::Mat result;
		cv::Mat point;
		std::string outfilename(m_filename);
		outfilename += ".diffimage.result." + boost::lexical_cast<std::string>(event);
		try {
			std::ofstream outimage(outfilename.c_str(),std::ios::out);
			outimage << "#pca difference result image, note: still doing PCA on the non-image features too\n";
			point = m_pca.project(cv::Mat(in).t());
			result = m_pca.backProject(point);
			for (unsigned i=0;i<shape.height;++i){
				for (unsigned j=0;j<shape.width;++j){
					float value = in.at(nonimagefeatures + i*shape.width + j);
					if (value > 0.0)
						value -= result.at<float>(nonimagefeatures + i*shape.width + j);
					outimage << value << "\t";
				}
				outimage << "\n";
			}
			outimage.close();
		} catch (cv::Exception& err){
			std::cerr << "projection failed in Learning::difference_pca_images() method, err: " << err.msg <<  std::endl;
			return false;
		} catch (std::exception & e) {
			std::cerr << "file handeling failed in Learning::difference_pca_images() method, err: " << e.what() << std::endl;
			return false;
		}
		return true;
	}
	bool Learning::compare_pca_images(std::vector<float> & in,cv::Size shape,const unsigned nonimagefeatures,const unsigned event){
		cv::Mat result;
		cv::Mat point;
		std::string infilename(m_filename);
		std::string outfilename(m_filename);
		infilename += ".testpcaimage." + boost::lexical_cast<std::string>(event);
		outfilename += ".testpcaimage.result." + boost::lexical_cast<std::string>(event);
		try {
			std::ofstream inimage(infilename.c_str(),std::ios::out);
			std::ofstream outimage(outfilename.c_str(),std::ios::out);
			inimage << "#pca input image, note: still doing PCA on the non-image features too\n";
			outimage << "#pca result image, note: still doing PCA on the non-image features too\n";
			point = m_pca.project(cv::Mat(in).t());
			result = m_pca.backProject(point);
			for (unsigned i=0;i<shape.height;++i){
				for (unsigned j=0;j<shape.width;++j){
					inimage << in.at(nonimagefeatures + i*shape.width + j) << "\t";
					outimage << result.at<float>(nonimagefeatures + i*shape.width + j) << "\t";
				}
				inimage << "\n";
				outimage << "\n";
			}
			inimage.close();
			outimage.close();
		} catch (cv::Exception& err){
			std::cerr << "projection failed in Learning::compare_pca_images() method, err: " << err.msg <<  std::endl;
			return false;
		} catch (std::exception & e) {
			std::cerr << "file handeling failed in Learning::compare_pca_images() method, err: " << e.what() << std::endl;
			return false;
		}
		return true;
	}
	bool Learning::compare_pca_images(std::vector< std::vector< float > > & in,cv::Size shape,const unsigned nonimagefeatures)
	{
		cv::Mat result;
		cv::Mat point;
		for (unsigned r=0;r<in.size();++r){
			std::string infilename(m_filename);
			std::string outfilename(m_filename);
			infilename += ".testpcaimages." + boost::lexical_cast<std::string>(r);
			outfilename += ".testpcaimages.result." + boost::lexical_cast<std::string>(r);
			// HERE //
			// we could figure a way to project and back project the entire matrix...
			// But then we have to seperately unflatten the rows for the images...
			try {
				std::ofstream inimage(infilename.c_str(),std::ios::out);
				std::ofstream outimage(outfilename.c_str(),std::ios::out);
				inimage << "#pca input image, note: still doing PCA on the non-image features too\n";
				outimage << "#pca result image, note: still doing PCA on the non-image features too\n";
				point = m_pca.project(cv::Mat(in.at(r)).t());
				result = m_pca.backProject(point);
				for (unsigned i=0;i<shape.height;++i){
					for (unsigned j=0;j<shape.width;++j){
						inimage << in[r].at(nonimagefeatures + i*shape.width + j) << "\t";
						outimage << result.at<float>(nonimagefeatures + i*shape.width + j) << "\t";
					}
					inimage << "\n";
					outimage << "\n";
				}
				inimage.close();
				outimage.close();
			} catch (cv::Exception& err){
				std::cerr << "projection failed in Learning::compare_pca_images() method, err: " << err.msg <<  std::endl;
				return false;
			} catch (std::exception & e) {
				std::cerr << "filehandeling failed in Learning::compare_pca_images() method, err: " << e.what() << std::endl;
				return false;
			}
		}
		return true;
	}
	bool Learning::compare_pca(std::vector< std::vector< float > > & in)
	{
		cv::Mat result;
		cv::Mat point;
		for (unsigned i=0;i<in.size();++i){
			std::string filename(m_filename);
			filename += ".testpca." + boost::lexical_cast<std::string>(i);
			// HERE //
			// we could figure a way to project and back project the entire matrix...
			// But then we have to seperately unflatten the rows for the images...
			try {
				std::ofstream outfile(filename.c_str(),std::ios::out);
				outfile << "#in\tpca_result\n";
				point = m_pca.project(cv::Mat(in.at(i)).t());
				result = m_pca.backProject(point);
				for (unsigned j=0;j<in[i].size();++j){
					outfile << in[i][j] << "\t" << result.at<float>(j) << "\n";
				}
				outfile.close();
			} catch (cv::Exception& err){
				std::cerr << "projection failed in Learning::compare_pca() method, err: " << err.msg <<  std::endl;
				return false;
			}
		}
		return true;
	}

	bool Learning::inspect_eigen_pca(void)
	{
		std::string eigenfilename(m_filename);
		eigenfilename += ".eigen";
		std::ofstream eigenout(eigenfilename.c_str(),std::ios::out);
		std::cerr << "attempting to write to " << eigenfilename << std::endl;
		if (!eigenout.is_open()){
			std::cerr << "failing to open file for Learning::inspect_eigen_pca() method" << std::endl;
			return false;
		}
		
		eigenout << "#should work out to be nfeatures wide and pca_features long matrix\n";
		eigenout << "#" << m_features_labels << "\n";
		for (unsigned i=0;i<m_pca.eigenvectors.rows;++i){
			for (unsigned j=0;j<m_pca.eigenvectors.cols;++j){
				// This is outputting into column vectors for ultimately a parallel representation;
				eigenout << std::pow(m_pca.eigenvectors.at<float>(i,j),2) << "\t";
			}
			eigenout << "\n";
		}
		eigenout << "\n";
		eigenout.close();
		return true;
	}

	bool Learning::print_eigenfaces(cv::Size shapein)
	{
		std::ofstream outfile;
		unsigned nonimage_features = m_pca.eigenvectors.cols - (shapein.area());
		for (unsigned f=0;f<m_pca.eigenvectors.rows;++f){
			std::string filename(m_filename);
			filename += ".eigenface.";
			filename += boost::lexical_cast<std::string>(f);
			outfile.open(filename.c_str(),std::ios::out);
			cv::Mat facevectors = m_pca.eigenvectors(cv::Range::all(),cv::Range(nonimage_features , m_pca.eigenvectors.cols));
			cv::Mat faceout = facevectors.row(f).reshape(1,shapein.height); // reshape(1 channel, n rows)
			for (unsigned i=0;i<faceout.rows;++i){
				for (unsigned j=0;j<faceout.cols;++j){
					outfile << faceout.at<float>(i,j) << "\t";
				}
				outfile << "\n";
			}
			outfile.close();
		}
		
	}


			// HERE // 
			// example of filestorage reading //
			// primarily for after the learning has happened and been written to xml format //
			/*
			 *
			// READ THE STORED LEARN DATA
			cv::FileStorage fs("test.xml", cv::FileStorage::READ);
			cv::Mat stored_mean, stored_eigenvalues, stored_vectors, stored_projections, stored_labels;
			int stored_num_componants = fs["num_components"];
			fs["mean"] >> stored_mean;
			fs["eigenvalues"] >> stored_eigenvalues;
			fs["eigenvectors"] >> stored_vectors;
			fs["projections"] >> stored_projections;
			fs["labels"] >> stored_labels;
			 */
	void Learning::set_features_labels(std::string & in_labels)
	{
		m_features_labels = in_labels;
	}
	// HERE
	// PCA 50 images for the eigen images
	// repeate for n groups of 50
	// 
	// OK, the pipeline should be something like ... remove the average dark image, then threshold, then calculate centroids(keep for features, 
	// then crop image, then feed training set (no lasing) to eigenfaces to reduce the dimensions (PCA) and k-means classify those,
	// Now use the classified images and feed images and classifications into Fisherfaces training
	// Now use the fisherfaces to give coefficients for the rest of the (no lasing) set
	// Now cluster the fisher coefficients and calssify them
	//
	// Now run with lasing and use the fisher model to project image and get classifier from the coefficients
	// subtract the no lasing reference that is associated with that fisher class
	// identify the islands that remain for centroids, widths, and separations
	// Classify the with lasing Neon spectra and use these clustered to define classifiers for the fisher output 
	// Add these as features to a kmeans and cluster to find classifiers for the x-ray Neon spectra
	// Train the regressor model on the spectra-as-clasifiers and the fisher output to make a spectrum prediction based on fisher output.
	bool Learning::fill(std::vector<float> & in)
	{
		if (m_read_pca && m_pca.eigenvalues.rows != 0)
			return false;
		if (m_features.rows >= m_pca_samples_limit)
			return false;
		cv::Mat inMat(cv::Mat(in,true).t());// try making this false later
		if (m_features.rows == 0){
			m_features = inMat;
			m_features.reserve(1000); // this reserves for at least 1000 future m_features.push_back(inMat) operations
			return true;
		}
		if (m_features.cols != in.size()){
			std::cerr << "unequal features length in Learning::fill() method" << std::endl;
			return false;
		}
		m_features.push_back(inMat); // HERE, try making this cpydata to false and see if it fails
		return true;
	}

	// File interface methods //
	void Learning::open_file( std::string & filename ){
		m_filename = filename;
		m_outfile.open(m_filename.c_str(),std::ios::out);
	}
	void Learning::close_file(void){
		if (m_outfile.is_open())
			m_outfile.close();
	}
	// Eventually write out a parallel representation of the feature list and color based on classifier //
	bool Learning::print_out(void)
	{
		if (!m_outfile.is_open())
			m_outfile.open(m_filename.c_str(),std::ios::app);
		if (!m_outfile.is_open()){
			std::cerr << "Somehow failed to open the features output file in Learning::print_out() meathod" << std::endl;
			return false;
		}
		if (m_classifiers.rows != m_features.rows){
			std::cerr << "mis-matched number of samples between classifiers and features" << std::endl;
			return false;
		}
		// feels like we need something for event ID matching here //

		for (unsigned i=0;i<m_classifiers.rows;++i){
			m_outfile << m_classifiers.at<int>(i) << "\t";
			for (unsigned f=0;f<m_features.cols;++f){
				m_outfile << m_features.at<float>(i,f) << "\t";
			}
			m_outfile << "\n";
		}
		return true;
	}

	void Learning::print_header(void){
		if (!m_outfile.is_open())
			m_outfile.open(m_filename.c_str(),std::ios::out);
		m_outfile << "# classifier ID\tfeature list\n";
	}

} // namespace CookieBox_pkg

