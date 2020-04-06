#include "lungModelling.h"

#ifndef _HELPER_
#define _HELPER_

class helper{
public:
	double uniformRandom(void);
	double normalRandom(void);
	std::vector<std::vector<float>> generate3D(Vector3f start, Vector3f end, float radiusMultiplier, float radius, int numOfPOints, double multiplier, int mode);
	float round2d(float input);
	double round2d(double input);
	float roundatdecimal(float input, unsigned int decimal);
	double roundatdecimal(double input, unsigned int decimal);
	void reset(std::vector<int>& target);
	float getSpecificMean(std::vector<float> input);
	void absValueVector(std::vector<float> &input);
	void printVectorToFile(std::vector<float> in, std::string filename);
	void printVectorToFile(std::vector<int> in, std::string filename);
	void printVectorToFile(std::vector<std::vector<int>> in, std::string filename);
	void printVectorToFile(std::vector<std::vector<float>> in, std::string filename);
	Vector3f projPointOnVector(Vector3f point, Vector3f Pv1, Vector3f Pv2);
	Vector3f projPointOnPlane(Vector3f randomPoint, Vector3f pointUponPlane, Vector3f normalOfPlane);
	bool valueExistsInVector(std::vector<int> & in, int control);
	bool valueExistsInVector(std::vector<float> & in, float control);
	bool valueExistsInVector(std::vector<double> & in, double control);
	template<typename T> std::vector<T> multelementwize(std::vector<T> inputVector1, std::vector<T> inputVector2);
	Vector3f getNormalFromTriangle(Vector3f a, Vector3f b, Vector3f c);
	float cosineAngleBetweenVectors(Vector3f a, Vector3f b);
	Vector3f stdtoeigenvec(std::vector<float> in);
	std::vector<float> eigentostdvec(Vector3f in);
	std::vector<std::vector<float>> eigentostdvec(std::vector <Vector3f> in, std::vector<std::vector<float>> & out);


	Vector3d vectorRotation(Vector3d r, double theta, Vector3d v);
	void gen_random(char *s, const int len);

	bool write_file_binary(std::string const & filename,
		char const * data, size_t const bytes)
	{
		std::ofstream b_stream(filename.c_str(),
			std::fstream::out | std::fstream::binary);
		if (b_stream)
		{
			b_stream.write(data, bytes);
			return (b_stream.good());
		}
		return false;
	}

};
#endif _HELPER_