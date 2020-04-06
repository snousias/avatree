#include "lungModelling.h"

std::vector<std::vector<float>> helper::eigentostdvec(std::vector <Vector3f> in, std::vector<std::vector<float>> & out) {
	
	helper functions;
	out.clear();
	for (int j = 0; j <in.size(); j++) {
		out.push_back(functions.eigentostdvec(in[j]));
	}
	return out;
}


Vector3f helper::stdtoeigenvec(std::vector<float> in) {
	Vector3f p;

	p[0] = in[0];
	p[1] = in[1];
	p[2] = in[2];

	return p;
}

std::vector<float> helper::eigentostdvec(Vector3f in) {
	std::vector<float> p;

	p.resize(3);

	p[0] = in[0];
	p[1] = in[1];
	p[2] = in[2];

	return p;
}

Vector3d helper::vectorRotation(Vector3d r, double theta, Vector3d v){
	v.normalize();
	r.normalize();

	double c = cos(2 * Pi*theta / 360);
	double s = sin(2 * Pi*theta / 360);
	double t = 1 - c;
	double ux = r[0];
	double uy = r[1];
	double uz = r[2];

	MatrixXd T(3, 3);

	T << (t*ux*ux) + c, (t*ux*uy) - (s*uz), (t*ux*uz) + (s*uy),
		(t*uy*ux) + (s*uz), (t*uy*uy) + c, (t*uy*uz) - (s*ux),
		(t*uz*ux) - (s*uy), (t*uz*uy) + (s*ux), (t*uz*uz) + c;

	return T*v;
}

void helper::printVectorToFile(std::vector<int> in, std::string filename){
	std::ofstream file;
	file.open(filename);
	for (int i = 0; i < in.size(); i++){
		file << in.at(i) << std::endl;
	}
	file.close();
}

void helper::printVectorToFile(std::vector<std::vector<int>> in, std::string filename){
	std::ofstream file;
	file.open(filename);
	for (int i = 0; i < in.size(); i++){
		for (int j = 0; j < in[i].size(); j++){
			file << in.at(i).at(j)<<"-";
			
		}
		file << std::endl;
	}
	file.close();
}

void helper::printVectorToFile(std::vector<std::vector<float>> in, std::string filename){
	std::ofstream file;
	file.open(filename);
	for (int i = 0; i < in.size(); i++){
		for (int j = 0; j < in[i].size(); j++){
			file << in.at(i).at(j) << "-";

		}
		file << std::endl;
	}
	file.close();
}



void helper::absValueVector(std::vector<float> &input){
	for (int i = 0; i < input.size(); i++){ if (input.at(i) < 0){ input.at(i) = -input.at(i); } }
}








float helper::getSpecificMean(std::vector<float> inputVector){
	float mean = 1.0;
	int numOfSamples = 0;
	//sort(inputVector.begin(), inputVector.end());
	//inputVector.erase(unique(inputVector.begin(), inputVector.end()), inputVector.end());
	//int trim = (int)(inputVector.size()*0.2);
	//inputVector.erase(inputVector.begin(), inputVector.begin()+trim);
	//inputVector.erase(inputVector.end()-trim, inputVector.end());

	if (inputVector.size() > 0){
		for (int i = 0; i < inputVector.size(); i++){
			if (inputVector.at(i) > 0){
				mean += inputVector.at(i);
				numOfSamples++;
			}
		}
		//mean = mean / inputVector.size();
		mean = mean / numOfSamples;
	}
	//For Median
	//mean = inputVector.at((int)((float)inputVector.size() / 2.0));
	return mean;
}




void helper::printVectorToFile(std::vector<float> in, std::string filename){
	std::ofstream file;
	file.open(filename);
	for (int i = 0; i < in.size(); i++){
		file << in.at(i) << std::endl;
	}
	file.close();
}


void helper::reset(std::vector<int>& target) {
	while (target.size() != 0)  target.pop_back();
	//std::cout << "Reset complete.New size is "<< target.size() << std::endl;
}


float helper::round2d(float input){
	//Rounding function
	return roundf(input * 100) / 100;
}

double helper::round2d(double input){
	//Rounding function
	return (double)(roundf((float)input * 100) / 100);
}

float helper::roundatdecimal(float input, unsigned int decimal){
	//Rounding function
	return roundf(input * (10 * decimal)) / (10 * decimal);
}

double helper::roundatdecimal(double input, unsigned int decimal){
	//Rounding function
	return (double)roundf(input * (10 * decimal)) / (10 * decimal);
}

double helper::uniformRandom(void)
{
	return ((double)(rand()) + 1.) / ((double)(RAND_MAX)+1.);
}
// return a normally distributed random number
double helper::normalRandom(void)
{
	double u1 = uniformRandom();
	double u2 = uniformRandom();
	return cos(8.*atan(1.)*u2)*sqrt(-2.*log(u1));
}


std::vector<std::vector<float>> helper::generate3D(Vector3f start, Vector3f end, float radiusMultiplier,float radius, int numOfPOints, double multiplier, int mode){

	Vector3f vecOfRef;
	Vector3f  vec = end - start;
	Vector3f point;

	Vector3f zpositive = { 0.0, 0.0, 1.0 };

	Vector3f ypositive = { 0.0, 1.0, 1.0 };

	Vector3f xpositive = { 1.0, 0.0, 0.0 };

	if (abs(vec.dot(zpositive)) != vec.norm()*zpositive.norm()){
		vecOfRef = vec.cross(zpositive);
	}
	else if (abs(vec.dot(ypositive)) != vec.norm()*ypositive.norm()){
		vecOfRef = vec.cross(ypositive);
	}
	else if (abs(vec.dot(xpositive)) != vec.norm()*xpositive.norm()){
		vecOfRef = vec.cross(xpositive);
	}

	Vector3f axis;
	std::vector<std::vector<float>> geom;
	axis = vec;
	axis.normalize();
	vecOfRef.normalize();



	//Predefinened radius for each generation
	if (mode == 1){
		for (int i = 0; i < numOfPOints; i++){

			//**Uniform sampling in angle---------------------------------------------------------
			//float angle = (Pi/3.0)*normalRandom()+Pi;
			float angle = (float)(2 * Pi)*uniformRandom();  //Angle in rad
			//**Rotation---------------------------------------------------------
			Quaternion<float> q;
			q = AngleAxis<float>(angle, axis);
			point = q*vecOfRef;
			point.normalize();
			//**Select sistance distribution-----------------------------------------
			//float dist = 0.3* normalRandom() + 0.5;	
			//float dist = 0.08 + 0.84*uniformRandom();
			float dist = (float)1.5*uniformRandom() - 0.25;
			//float dist = uniformRandom();
			
			
			//**Place point------------------------------------------------------------
			Vector3f result = start + dist*vec + (radiusMultiplier*radius)*point;     // + (0.0833* normalRandom()*point); //Randomize mesh ,simulate irregularities ??
			//Vector3f result = start + dist*vec + (0.1666* normalRandom() * (10.0 / multiplier)*point);     
			//Vector3f result = start + dist*vec + (0.1666* normalRandom() * (10.0)*point);     
			//**Push back to geometry
			geom.push_back({ result[0], result[1], result[2] });


		}
	}

	//Different density for each generation
	if (mode == 2){
		for (int i = 0; i < (int)((double)numOfPOints / (14.72/radius)); i++){   //Density mode //powf(2.0,(float)multiplier)
			//**Uniform sampling in angle---------------------------------------------------------
			//float angle = (Pi/3.0)*normalRandom()+Pi;
			float angle = (float)(2 * Pi)*uniformRandom();
			//**Rotation---------------------------------------------------------
			Quaternion<float> q;
			q = AngleAxis<float>(angle, axis);
			point = q*vecOfRef;
			point.normalize();
			//**Select sistance distribution-----------------------------------------
			//float dist = 0.1666* normalRandom() + 0.5;	
			//float dist = 0.05 + 0.9*uniformRandom();
			float dist = uniformRandom();
			//Place point------------------------------------------------------------
			
			Vector3f result = start + dist*vec + (2 * 0.166* normalRandom() *(radiusMultiplier*14.72)*point);
			
			geom.push_back({ result[0], result[1], result[2] });
		}
	}

	//Diffenent sigma for each generation
	if (mode == 3){
		for (int i = 0; i < numOfPOints; i++){   //Density mode
			//**Uniform sampling in angle---------------------------------------------------------
			//float angle = (Pi/3.0)*normalRandom()+Pi;
			float angle = (float)(2 * Pi)*uniformRandom();
			//**Rotation---------------------------------------------------------
			Quaternion<float> q;
			q = AngleAxis<float>(angle, axis);
			point = q*vecOfRef;
			point.normalize();
			//**Select sistance distribution-----------------------------------------
			//float dist = 0.1666* normalRandom() + 0.5;	
			//float dist = 0.05 + 0.9*uniformRandom();
			float dist = uniformRandom();
			//Place point------------------------------------------------------------
			//Vector3f result = start + dist*vec + (radius *point);     // + (0.0833* normalRandom()*point); //Randomize mesh ,simulate irregularities ??
			Vector3f result = start + dist*vec + (0.1666* normalRandom() * (radius *radiusMultiplier)*point);


			geom.push_back({ result[0], result[1], result[2] });


		}
	}


	return geom;
}


Vector3f helper::projPointOnVector(Vector3f point, Vector3f P0, Vector3f P1){


	Vector3f an = point - P0;
	Vector3f bn = P1 - P0;
	float projab = an.dot(bn) / (bn.norm());
	Vector3f projection = P0 + projab*bn / (bn.norm());
	return projection;
}

Vector3f helper::projPointOnPlane(Vector3f randomPoint, Vector3f pointUponPlane, Vector3f normalOfPlane){


	Vector3f projectionOnNormal;
	normalOfPlane.normalize();
	projectionOnNormal = projPointOnVector(randomPoint, pointUponPlane, pointUponPlane + normalOfPlane);



	Vector3f projection = randomPoint - projectionOnNormal;


	return projection;
}


bool helper::valueExistsInVector(std::vector<int> & in, int control){
	bool isThere = false;
	for (int j = 0; j < in.size(); j++){
		if (in[j] == control){ isThere = true; }
	}
	return isThere;
}

bool helper::valueExistsInVector(std::vector<float> & in, float control){
	bool isThere = false;
	for (int j = 0; j < in.size(); j++){
		if (in[j] == control){ isThere = true; }
	}
	return isThere;
}

bool helper::valueExistsInVector(std::vector<double> & in, double control){
	bool isThere = false;
	for (int j = 0; j < in.size(); j++){
		if (in[j] == control){ isThere = true; }
	}
	return isThere;
}

void helper::gen_random(char *s, const int len) {
	static const char alphanum[] =
		"0123456789"
		"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		"abcdefghijklmnopqrstuvwxyz";

	for (int i = 0; i < len; ++i) {
		s[i] = alphanum[rand() % (sizeof(alphanum) - 1)];
	}

	s[len] = 0;
}


//Element wize vector multiplication
template<typename T> std::vector<T> helper::multelementwize(std::vector<T> inputVector1, std::vector<T> inputVector2){

	std::vector<T> result;

	if (inputVector1.size() == inputVector2.size()){

		for (int i = 0; i < inputVector1.size(); i++)
		{
			result.at(i) = inputVector1.at(i)*inputVector2.at(2);
		}

	}

	return result;
}


Vector3f helper::getNormalFromTriangle(Vector3f a, Vector3f b, Vector3f c){

	Vector3f e[2], n;

	e[0] = a - b;
	e[1] = a - c;
	n = e[0].cross(e[1]);
	n.normalize();
	n = -n;
	if (isnan(n[0])){ n[0] = 0.0; }
	if (isnan(n[1])){ n[1] = 0.0; }
	if (isnan(n[2])){ n[2] = 0.0; }

	return n;
}

float helper::cosineAngleBetweenVectors(Vector3f a, Vector3f b){
	float metric;

	metric = (a.dot(b))/(a.norm()*b.norm());


	return metric;
}