#ifndef DATAOPS_H
#define DATAOPS_H

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

namespace CookieBox_pkg {

	typedef std::vector<float> record_f;
	typedef  std::vector< record_f > data_f;
	typedef std::vector<std::string> record_s;
	typedef  std::vector< record_s > data_s;

	std::istream& operator >> ( std::istream & ins, record_f & record );
	std::istream& operator >> ( std::istream& ins, data_f & data );
	std::istream& operator >> ( std::istream & ins, record_s & record );
	std::istream& operator >> ( std::istream& ins, data_s & data );

}

#endif
