#include "lungModelling.h"

unsigned int split(const std::string &txt, std::vector<std::string> &strs, char ch)
{
	unsigned int pos = txt.find(ch);
	unsigned int initialPos = 0;
	strs.clear();

	// Decompose statement
	while (pos != std::string::npos) {
		strs.push_back(txt.substr(initialPos, pos - initialPos + 1));
		initialPos = pos + 1;

		pos = txt.find(ch, initialPos);
	}

	// Add the last one
	strs.push_back(txt.substr(initialPos, (std::min)((int)pos, (int)txt.size()) - initialPos + 1));

	return strs.size();
}

std::string extract_ints(std::ctype_base::mask category, std::string str, std::ctype<char> const& facet)
{
	using std::strlen;

	char const *begin = &str.front(),
		*end = &str.back();

	auto res = facet.scan_is(category, begin, end);

	begin = &res[0];
	end = &res[strlen(res)];

	return std::string(begin, end);
}

std::string extract_ints(std::string str)
{
	return extract_ints(std::ctype_base::digit, str,
		std::use_facet<std::ctype<char>>(std::locale("")));
}

void dotObj::exportSelectedIndicesToFile(std::string name) {
	std::ofstream outfile;
	std::string output = name;
	outfile.open(output);

	for (unsigned int i = 0; i < this->selectedVertices.size(); i++)
	{
		outfile << this->selectedVertices.at(i) << std::endl;
	}
	outfile.close();

	return;
}

std::stringstream dotObj::toOFF(void) {
	std::stringstream outfile;

	if (this->edges.size() == 0) {
		this->getEdges();
	}

	outfile << "OFF" << std::endl;
	outfile << this->vertices.size() << " " << this->faces.size() << " " << "0" << std::endl;
	for (unsigned int i = 0; i < this->vertices.size(); i++)
	{
		outfile << std::setprecision(PRECISION) << this->vertices.at(i).at(0) << " " << std::setprecision(PRECISION) << this->vertices.at(i).at(1) << " " << std::setprecision(PRECISION) << this->vertices.at(i).at(2) << std::endl;
	}

	for (unsigned int i = 0; i < this->faces.size(); i++)
	{
		outfile << "3 " << this->faces.at(i).at(0) - 1 << " " << this->faces.at(i).at(3) - 1 << " " << this->faces.at(i).at(6) - 1 << std::endl;
	}

	return outfile;
}

void dotObj::toOFF(std::string name) {
	if (this->edges.size() == 0) { this->getEdges(); }
	std::ofstream outfile;
	std::string output = name;
	outfile.open(output);
	outfile << "OFF" << std::endl;
	outfile << this->vertices.size() << " " << this->faces.size() << " " << "0" << std::endl;
	for (unsigned int i = 0; i < this->vertices.size(); i++)
	{
		outfile << std::setprecision(PRECISION) << this->vertices.at(i).at(0) << std::setprecision(PRECISION) << " " << this->vertices.at(i).at(1) << " " << std::setprecision(PRECISION) << this->vertices.at(i).at(2) << std::endl;
	}

	for (unsigned int i = 0; i < this->faces.size(); i++)
	{
		outfile << "3 " << this->faces.at(i).at(0) - 1 << " " << this->faces.at(i).at(3) - 1 << " " << this->faces.at(i).at(6) - 1 << std::endl;
	}

	outfile.close();

	return;
}

std::ostringstream dotObj::initializeFromFileE(std::string thefile, bool verbose, bool readvertices, bool readnormals, bool readfaces, bool readedges) {
	std::ifstream myfile(thefile);
	std::ostringstream exportedStream;
	float output[9];
	std::string line;
	std::string lineVals;
	std::string lineVals2;
	std::string val;
	std::string val0, val1, val2, val3, val4, val5, val6, val7, val8;
	int totlength = 0;
	for (int k = 0; k < 9; k++) {
		output[k] = 0.0f;
	}

	int group = 0;

	this->segment_property_map.clear();
	this->vertices.clear();
	this->normals.clear();
	this->faces.clear();

	if (myfile.is_open())
	{
		if (verbose) { std::cout << "Building obj" << std::endl; }

		while (!myfile.eof())
		{
			getline(myfile, line);
			if (line[0] == 'v' && line[1] == ' ') {
				exportedStream << line << "\r\n";

				lineVals = line.substr(2);
				val0 = lineVals.substr(0, lineVals.find(' '));
				output[0] = (float)atof(val0.c_str());

				lineVals = lineVals.substr(val0.length() + 1);
				val1 = lineVals.substr(0, lineVals.find(' '));
				output[1] = (float)atof(val1.c_str());

				lineVals = lineVals.substr(val1.length() + 1);
				val2 = lineVals;
				output[2] = (float)atof(val2.c_str());

				this->vertices.push_back({ output[0], output[1], output[2], 0.0f });

				//std::cout  << " " << "v" << " " << output[0] << " " << output[1] << " " << output[2] << std::endl;
			}

			if (line[0] == 'v' && line[1] == 'n') {
				exportedStream << line << "\r\n";
				lineVals = line.substr(3);
				val0 = lineVals.substr(0, lineVals.find(' '));
				output[0] = (float)atof(val0.c_str());

				val1 = lineVals.substr(val0.length() + 1, lineVals.find(' '));
				output[1] = (float)atof(val1.c_str());

				val2 = lineVals.substr(lineVals.find_last_of(' ') + 1);
				output[2] = (float)atof(val2.c_str());

				this->normals.push_back({ output[0], output[1], output[2] });

				//std::cout  << " " << "vn" << " " << output[0] << " " << output[1] << " " << output[2] << std::endl;
			}

			if (line[0] == 'g' && line[1] == ' ') {
				exportedStream << line << "\r\n";
				group = group + 1;
				lineVals = line.substr(2);
				segment_property_name.push_back(lineVals);
			}

			if (line[0] == 'f' && line[1] == ' ') {
				exportedStream << line << "\r\n";
				for (int i = 0; i < 9; i++) { output[i] = 0; }

				//First Set
				lineVals = line.substr(2);
				//Term #0
				val0 = lineVals.substr(0, lineVals.find('/', 0));
				output[0] = (float)atof(val0.c_str());
				totlength = val0.length();
				lineVals = lineVals.substr(totlength + 1);
				//Term #1
				val1 = lineVals.substr(0, lineVals.find('/', 0));
				output[1] = (float)atof(val1.c_str());
				totlength = val1.length();
				lineVals = lineVals.substr(totlength + 1);
				//Term #2
				val2 = lineVals.substr(0, lineVals.find(' '));
				output[2] = (float)atof(val2.c_str());
				totlength = val2.length();
				lineVals = lineVals.substr(totlength + 1);

				//Term #3
				val3 = lineVals.substr(0, lineVals.find('/', 0));
				output[3] = (float)atof(val3.c_str());
				totlength = val3.length();
				lineVals = lineVals.substr(totlength + 1);
				//Term #4
				val4 = lineVals.substr(0, lineVals.find('/', 0));
				output[4] = (float)atof(val4.c_str());
				totlength = val4.length();
				lineVals = lineVals.substr(totlength + 1);
				//Term #5
				val5 = lineVals.substr(0, lineVals.find(' ', 0));
				output[5] = (float)atof(val5.c_str());
				totlength = val5.length();
				lineVals = lineVals.substr(totlength + 1);
				//Term #6
				val6 = lineVals.substr(0, lineVals.find('/', 0));
				output[6] = (float)atof(val6.c_str());
				totlength = val6.length();
				lineVals = lineVals.substr(totlength + 1);
				//Term #7
				val7 = lineVals.substr(0, lineVals.find('/', 0));
				output[7] = (float)atof(val7.c_str());
				totlength = val7.length();
				lineVals = lineVals.substr(totlength + 1);
				//Term #8
				val8 = lineVals.substr(0, lineVals.find(' ', 0));
				output[8] = (float)atof(val8.c_str());

				this->faces.push_back({ (int)output[0], (int)output[1], (int)output[2], (int)output[3], (int)output[4], (int)output[5], (int)output[6], (int)output[7], (int)output[8] });

				this->segment_property_map.push_back(group);

				//if (group >= 0){ this->segment_property_map.push_back(group); }

				//std::cout <<  " " << "f" << " " << output[0] << "/" << output[1] << "/" << output[2] << " " << output[3] << "/" << output[4] << "/" << output[5] << " " << output[6] << "/" << output[7] << "/" << output[8] << std::endl;
			}

			if (line[0] == 'l' && line[1] == ' ') {
				exportedStream << line << "\r\n";
				lineVals = line.substr(2);
				val0 = lineVals.substr(0, lineVals.find(" "));
				output[0] = (float)atof(val0.c_str());
				val1 = lineVals.substr(val0.length() + 1);
				//val1 = lineVals.substr(val0.length() + 1, lineVals.find(' '));
				output[1] = (float)atof(val1.c_str());

				this->lines.push_back({ (int)output[0], (int)output[1] });
			}
		}
	}

	if (this->segment_property_map.size() > 0) {
		//WritePerVertex
		for (int k = 0; k < this->vertices.size(); k++) {
			segment_property_map_per_vertex.push_back(0);
		}

		for (int k = 0; k < this->faces.size(); k++) {
			int v0, v1, v2, group;
			group = this->segment_property_map.at(k);
			v0 = this->faces.at(k).at(0) - 1;
			v1 = this->faces.at(k).at(3) - 1;
			v2 = this->faces.at(k).at(6) - 1;
			segment_property_map_per_vertex.at(v0) = group;
			segment_property_map_per_vertex.at(v1) = group;
			segment_property_map_per_vertex.at(v2) = group;
		}
	}

	if (readedges) { if (this->edges.size() == 0) { this->getEdges(); } }

	if (verbose) {
		std::cout << "vertices loaded : " << this->vertices.size() << std::endl;
		std::cout << "normals loaded : " << this->normals.size() << std::endl;
		std::cout << "faces loaded : " << this->faces.size() << std::endl;
		std::cout << "Load OBJ Complete" << std::endl;
	}

	return exportedStream;
}

void dotObj::initializeFromFile(std::stringstream & filestream, std::string fileType, bool verbose, bool readvertices, bool readnormals, bool readfaces, bool readedges) {
	std::string line;
	std::vector<int> positions;
	std::vector<std::string> nums;
	std::string s;
	bool v = false;
	bool f = false;

	this->segment_property_map.clear();
	this->vertices.clear();
	this->normals.clear();
	this->faces.clear();

	int numberOfLines = 0;
	while (!filestream.eof()) {
		getline(filestream, line);
		numberOfLines++;
	}
	filestream.clear(); // clear bad state after eof
	filestream.seekg(0);

	if (fileType == "OFF") {
		if (verbose) { std::cout << "Reading off file" << std::endl; }
		getline(filestream, line);
		if (line != "OFF") {
			std::cout << "Reading OFF file error" << std::endl;
			return;
		}

		getline(filestream, line);
		positions.clear();
		for (int i = 0; i < line.length(); ++i)
		{
			if (line[i] == ' ') {
				positions.push_back(i);
			}
		}

		while (!filestream.eof()) {
			getline(filestream, line);
			positions.clear();
			nums.clear();
			f = false;
			v = false;
			positions.push_back(0);
			for (int i = 0; i < line.length(); ++i)
			{
				if (line[i] == ' ') {
					positions.push_back(i);
				}
			}
			positions.push_back(line.length());

			for (int i = 0; i < positions.size() - 1; i++) {
				int p1 = positions[i];
				int p2 = positions[i + 1];
				s = line.substr(p1, p2 - p1);
				nums.push_back(s);
			}

			for (int i = 0; i < nums.size(); i++) {
				if (nums[i].find('.') != nums[i].npos) {
					v = true;
					f = false;
				}
			}

			if (!v && line[0] == '3' && line[1] == ' ') {
				f = true;
				v = false;
			}

			if (v) {
				this->vertices.push_back({ (float)atof(nums[0].c_str()), (float)atof(nums[1].c_str()), (float)atof(nums[2].c_str()) });
			}

			if (f) {
				this->faces.push_back({ 1 + (int)atof(nums[1].c_str()), 0, 0, 1 + (int)atof(nums[2].c_str()), 0, 0, 1 + (int)atof(nums[3].c_str()), 0, 0 });
			}
		}

		std::cout << "Reading OFF file complete" << std::endl;
	}

	if (fileType == "OBJ") {
	}
	if (fileType == "PLY") {
	}
	return;
}

void dotObj::initializeFromVTK(std::string thefile, bool verbose, bool readvertices, bool readnormals, bool readfaces, bool readedges) {
	std::ifstream myfile(thefile);
	std::ostringstream exportedStream;
	float output[9];
	std::string line;
	std::string lineVals;
	std::string lineVals2;
	std::string val;
	std::string val0, val1, val2, val3, val4, val5, val6, val7, val8;
	int totlength = 0;
	for (int k = 0; k < 9; k++) {
		output[k] = 0.0f;
	}

	int group = 0;

	this->segment_property_map.clear();
	this->vertices.clear();
	this->normals.clear();
	this->faces.clear();

	if (myfile.is_open())
	{
		if (verbose) { std::cout << "Building vtk model" << std::endl; }

		while (!myfile.eof())
		{
			getline(myfile, line);
			//std::cout << line << std::endl;

			if (isdigit(line[0]) && line[1] != ' ') {
				lineVals = line;
				val0 = lineVals.substr(0, lineVals.find(' '));
				output[0] = (float)atof(val0.c_str());
				lineVals = lineVals.substr(val0.length() + 1);
				val1 = lineVals.substr(0, lineVals.find(' '));
				output[1] = (float)atof(val1.c_str());
				lineVals = lineVals.substr(val1.length() + 1);
				val2 = lineVals;
				output[2] = (float)atof(val2.c_str());
				this->vertices.push_back({ output[0], output[1], output[2] });
			}

			if (isdigit(line[0]) && line[1] == ' ') {
				lineVals = line.substr(2);
				val0 = lineVals.substr(0, lineVals.find(" "));
				output[0] = (float)atof(val0.c_str());
				val1 = lineVals.substr(val0.length() + 1);
				output[1] = (float)atof(val1.c_str());
				this->lines.push_back({ (int)output[0] + 1, (int)output[1] + 1 });
			}
		}
	}

	if (readedges) { if (this->edges.size() == 0) { this->getEdges(); } }

	if (verbose) {
		std::cout << "Vertices loaded : " << this->vertices.size() << std::endl;
		std::cout << "Lines loaded : " << this->lines.size() << std::endl;
		std::cout << "Load OBJ Complete" << std::endl;
	}

	std::vector<std::vector<float>> Vs;
	for (int i = 0; i < this->vertices.size(); i++) {
		bool exists = false;
		for (int j = 0; j < Vs.size(); j++) {
			if (
				(this->functions.roundatdecimal(this->vertices[i][0], 3) ==
					this->functions.roundatdecimal(Vs[j][0], 3)) &&
					(this->functions.roundatdecimal(this->vertices[i][1], 3) ==
						this->functions.roundatdecimal(Vs[j][1], 3)) &&
						(this->functions.roundatdecimal(this->vertices[i][2], 3) ==
							this->functions.roundatdecimal(Vs[j][2], 3))
				) {
				exists = true;
			}
		}
		if (!exists) {
			Vs.push_back(this->vertices[i]);
		}
	}

	std::vector<std::vector<int>> Ls;
	for (int i = 0; i < this->lines.size(); i++) {
		int L0 = this->lines[i][0] - 1;
		int L1 = this->lines[i][1] - 1;
		int L0new = -1;
		int L1new = -1;
		std::vector<float> V0 = this->vertices[L0];
		std::vector<float> V1 = this->vertices[L1];
		for (int k = 0; k < Vs.size(); k++) {
			if (
				(this->functions.roundatdecimal(V0[0], 3) ==
					this->functions.roundatdecimal(Vs[k][0], 3)) &&
					(this->functions.roundatdecimal(V0[1], 3) ==
						this->functions.roundatdecimal(Vs[k][1], 3)) &&
						(this->functions.roundatdecimal(V0[2], 3) ==
							this->functions.roundatdecimal(Vs[k][2], 3))
				) {
				L0new = k;
			}
			if (
				(this->functions.roundatdecimal(V1[0], 3) ==
					this->functions.roundatdecimal(Vs[k][0], 3)) &&
					(this->functions.roundatdecimal(V1[1], 3) ==
						this->functions.roundatdecimal(Vs[k][1], 3)) &&
						(this->functions.roundatdecimal(V1[2], 3) ==
							this->functions.roundatdecimal(Vs[k][2], 3))
				) {
				L1new = k;
			}
		}

		if ((L1new == -1) || (L0new == -1)) {
			std::cout << "error" << std::endl;
		}

		Ls.push_back({ L0new + 1, L1new + 1 });
	}
	std::cout << "Read VTK Complete" << std::endl;
	this->vertices = Vs;
	this->lines = Ls;
	return;
}

void dotObj::initializeFromFile(std::string thefile, std::string fileType, bool verbose, bool readvertices, bool readnormals, bool readfaces, bool readedges) {
	std::string line;
	std::vector<int> positions;
	std::vector<std::string> nums;
	std::string s;
	bool v = false;
	bool f = false;

	std::ifstream fileInput(thefile);

	this->segment_property_map.clear();
	this->vertices.clear();
	this->normals.clear();
	this->faces.clear();

	int numberOfLines = 0;
	while (!fileInput.eof()) {
		getline(fileInput, line);
		numberOfLines++;
	}
	fileInput.clear(); // clear bad state after eof
	fileInput.seekg(0);

	if (fileType == "OFF") {
		if (verbose) { std::cout << "Reading off file" << std::endl; }
		getline(fileInput, line);
		if (line != "OFF") {
			std::cout << "Reading OFF file error" << std::endl;
			return;
		}

		getline(fileInput, line);
		positions.clear();
		for (int i = 0; i < line.length(); ++i)
		{
			if (line[i] == ' ') {
				positions.push_back(i);
			}
		}

		while (!fileInput.eof()) {
			getline(fileInput, line);
			positions.clear();
			nums.clear();
			f = false;
			v = false;
			positions.push_back(0);
			for (int i = 0; i < line.length(); ++i)
			{
				if (line[i] == ' ') {
					positions.push_back(i);
				}
			}
			positions.push_back(line.length());

			for (int i = 0; i < positions.size() - 1; i++) {
				int p1 = positions[i];
				int p2 = positions[i + 1];
				s = line.substr(p1, p2 - p1);
				nums.push_back(s);
			}

			for (int i = 0; i < nums.size(); i++) {
				if (nums[i].find('.') != nums[i].npos) {
					v = true;
					f = false;
				}
			}

			if (!v && line[0] == '3' && line[1] == ' ') {
				f = true;
				v = false;
			}

			if (v) {
				this->vertices.push_back({ (float)atof(nums[0].c_str()), (float)atof(nums[1].c_str()), (float)atof(nums[2].c_str()) });
			}

			if (f) {
				this->faces.push_back({ 1 + (int)atof(nums[1].c_str()), 0, 0, 1 + (int)atof(nums[2].c_str()), 0, 0, 1 + (int)atof(nums[3].c_str()), 0, 0 });
			}
		}

		std::cout << "Reading OFF file complete" << std::endl;
	}

	if (fileType == "OBJ") {
	}
	if (fileType == "PLY") {
	}
	return;
}

void dotObj::initializeFromFile(std::string thefile, bool verbose, bool readvertices, bool readnormals, bool readfaces, bool readedges) {
	std::ifstream myfile(thefile);
	float output[9];
	std::ostringstream exportedStream;
	std::string line;
	std::string lineVals;
	std::string lineVals2;
	std::string val;
	std::string val0, val1, val2, val3, val4, val5, val6, val7, val8;
	int totlength = 0;
	for (int k = 0; k < 9; k++) {
		output[k] = 0.0f;
	}

	int group = -1;
	this->segment_property_map.clear();
	this->vertices.clear();
	this->normals.clear();
	this->faces.clear();

	if (myfile.is_open())
	{
		if (verbose) { std::cout << "Building obj" << std::endl; }

		while (!myfile.eof())
		{
			getline(myfile, line);
			if (line[0] == 'v' && line[1] == ' ') {
				exportedStream << line;

				lineVals = line.substr(2);
				val0 = lineVals.substr(0, lineVals.find(' '));
				output[0] = (float)atof(val0.c_str());

				lineVals = lineVals.substr(val0.length() + 1);
				val1 = lineVals.substr(0, lineVals.find(' '));
				output[1] = (float)atof(val1.c_str());

				lineVals = lineVals.substr(val1.length() + 1);
				val2 = lineVals;
				output[2] = (float)atof(val2.c_str());

				this->vertices.push_back({ output[0], output[1], output[2], 0.0f });

				//std::cout  << " " << "v" << " " << output[0] << " " << output[1] << " " << output[2] << std::endl;
			}

			if (line[0] == 'v' && line[1] == 'n') {
				exportedStream << line;
				lineVals = line.substr(3);
				val0 = lineVals.substr(0, lineVals.find(' '));
				output[0] = (float)atof(val0.c_str());

				val1 = lineVals.substr(val0.length() + 1, lineVals.find(' '));
				output[1] = (float)atof(val1.c_str());

				val2 = lineVals.substr(lineVals.find_last_of(' ') + 1);
				output[2] = (float)atof(val2.c_str());

				this->normals.push_back({ output[0], output[1], output[2] });

				//std::cout  << " " << "vn" << " " << output[0] << " " << output[1] << " " << output[2] << std::endl;
			}

			if (line[0] == 'g' && line[1] == ' ') {
				exportedStream << line;
				group = group + 1;
				lineVals = line.substr(2);
				segment_property_name.push_back(lineVals);
			}

			if (line[0] == 'f' && line[1] == ' ') {
				exportedStream << line;
				for (int i = 0; i < 9; i++) { output[i] = 0; }

				//First Set
				lineVals = line.substr(2);
				//Term #0
				val0 = lineVals.substr(0, lineVals.find('/', 0));
				output[0] = (float)atof(val0.c_str());
				totlength = val0.length();
				lineVals = lineVals.substr(totlength + 1);
				//Term #1
				val1 = lineVals.substr(0, lineVals.find('/', 0));
				output[1] = (float)atof(val1.c_str());
				totlength = val1.length();
				lineVals = lineVals.substr(totlength + 1);
				//Term #2
				val2 = lineVals.substr(0, lineVals.find(' '));
				output[2] = (float)atof(val2.c_str());
				totlength = val2.length();
				lineVals = lineVals.substr(totlength + 1);

				//Term #3
				val3 = lineVals.substr(0, lineVals.find('/', 0));
				output[3] = (float)atof(val3.c_str());
				totlength = val3.length();
				lineVals = lineVals.substr(totlength + 1);
				//Term #4
				val4 = lineVals.substr(0, lineVals.find('/', 0));
				output[4] = (float)atof(val4.c_str());
				totlength = val4.length();
				lineVals = lineVals.substr(totlength + 1);
				//Term #5
				val5 = lineVals.substr(0, lineVals.find(' ', 0));
				output[5] = (float)atof(val5.c_str());
				totlength = val5.length();
				lineVals = lineVals.substr(totlength + 1);
				//Term #6
				val6 = lineVals.substr(0, lineVals.find('/', 0));
				output[6] = (float)atof(val6.c_str());
				totlength = val6.length();
				lineVals = lineVals.substr(totlength + 1);
				//Term #7
				val7 = lineVals.substr(0, lineVals.find('/', 0));
				output[7] = (float)atof(val7.c_str());
				totlength = val7.length();
				lineVals = lineVals.substr(totlength + 1);
				//Term #8
				val8 = lineVals.substr(0, lineVals.find(' ', 0));
				output[8] = (float)atof(val8.c_str());

				this->faces.push_back({ (int)output[0], (int)output[1], (int)output[2], (int)output[3], (int)output[4], (int)output[5], (int)output[6], (int)output[7], (int)output[8] });

				//this->segment_property_map.push_back(group);

				if (group >= 0) { this->segment_property_map.push_back(group); }

				//std::cout <<  " " << "f" << " " << output[0] << "/" << output[1] << "/" << output[2] << " " << output[3] << "/" << output[4] << "/" << output[5] << " " << output[6] << "/" << output[7] << "/" << output[8] << std::endl;
			}

			if (line[0] == 'l' && line[1] == ' ') {
				exportedStream << line;
				lineVals = line.substr(2);
				val0 = lineVals.substr(0, lineVals.find(" "));
				output[0] = (float)atof(val0.c_str());
				val1 = lineVals.substr(val0.length() + 1);
				//val1 = lineVals.substr(val0.length() + 1, lineVals.find(' '));
				output[1] = (float)atof(val1.c_str());

				this->lines.push_back({ (int)output[0], (int)output[1] });
			}
		}
	}

	if (this->segment_property_map.size() > 0) {
		//WritePerVertex
		this->segment_property_map_per_vertex.clear();
		for (int k = 0; k < this->vertices.size(); k++) {
			this->segment_property_map_per_vertex.push_back(0);
		}
		for (int k = 0; k < this->faces.size(); k++) {
			int v0, v1, v2, group;
			group = this->segment_property_map.at(k);
			v0 = this->faces.at(k).at(0) - 1;
			v1 = this->faces.at(k).at(3) - 1;
			v2 = this->faces.at(k).at(6) - 1;
			this->segment_property_map_per_vertex.at(v0) = group;
			this->segment_property_map_per_vertex.at(v1) = group;
			this->segment_property_map_per_vertex.at(v2) = group;
		}
	}

	if (readedges) { if (this->edges.size() == 0) { this->getEdges(); } }

	if (verbose) {
		std::cout << "vertices loaded : " << this->vertices.size() << std::endl;
		std::cout << "normals loaded : " << this->normals.size() << std::endl;
		std::cout << "faces loaded : " << this->faces.size() << std::endl;
		std::cout << "Load OBJ Complete" << std::endl;
	}
	return;
}

void dotObj::initializeFromPLY(std::string fnamein, std::string fnameout, bool verbose, bool readvertices, bool readnormals, bool readfaces, bool readedges) {
	std::ifstream myfile(fnamein);
	float output[9];
	std::string line;
	std::string lineVals;
	std::string val;
	std::string val0, val1, val2, val3, val4, val5, val6, val7, val8;
	int nr = 0;
	int fr = 0;
	int totlength = 0;
	for (int k = 0; k < 9; k++) {
		output[k] = 0.0f;
	}
	bool isvert = false;

	std::ofstream exported;
	if (fnameout != "")
	{
		exported.open(fnameout);
	}

	std::cout << "Open file" << std::endl;
	if (myfile.is_open())
	{
		if (verbose) { std::cout << "Building obj" << std::endl; }

		while (!myfile.eof())
		{
			getline(myfile, line);

			if (((line[0] >= '0' && line[0] <= '9') || line[0] == '-')) {
				isvert = false;
				for (int i = 0; i < strlen(line.c_str()); i++) {
					if (line[i] == '.') { isvert = true; }
				}
				if (isvert) {
					lineVals = line;// .substr(0);
					//std::cout << "line is"<<lineVals << std::endl;
					val0 = lineVals.substr(0, lineVals.find(' '));
					output[0] = (float)atof(val0.c_str());

					lineVals = lineVals.substr(val0.length() + 1);
					val1 = lineVals.substr(0, lineVals.find(' '));
					output[1] = (float)atof(val1.c_str());

					lineVals = lineVals.substr(val1.length() + 1);
					val2 = lineVals;
					output[2] = (float)atof(val2.c_str());

					this->vertices.push_back({ output[0], output[1], output[2], 0.0f });
					nr++;
					if (fnameout != "")
					{
						exported << "v" << " " << output[0] << " " << output[1] << " " << output[2] << "\n";
					}
				}

				//Faces
				else if (!isvert) {
					for (int i = 0; i < 9; i++) { output[i] = 0; }

					//First Set
					lineVals = line.substr(2);
					//std::cout << "line is"<<lineVals << std::endl;
					val0 = lineVals.substr(0, lineVals.find(' '));
					output[0] = (float)atof(val0.c_str()) + 1;
					//std::cout << val0 << std::endl;
					lineVals = lineVals.substr(val0.length() + 1);
					val1 = lineVals.substr(0, lineVals.find(' '));
					output[3] = (float)atof(val1.c_str()) + 1;
					//std::cout << val1 << std::endl;
					lineVals = lineVals.substr(val1.length() + 1);
					val2 = lineVals;
					output[6] = (float)atof(val2.c_str()) + 1;
					//std::cout << val2 << std::endl;

					fr++;
					if (fnameout != "")
					{
						exported << "f" << " " << (int)output[0] << " " << (int)output[3] << " " << (int)output[6] << "\n";
					}
					this->faces.push_back({ (int)output[0], (int)output[1], (int)output[2], (int)output[3], (int)output[4], (int)output[5], (int)output[6], (int)output[7], (int)output[8] });
				}
			}
		}
	}
	if (fnameout != "")
	{
		exported.close();
	}
	if (verbose) {
		std::cout << "vertices loaded : " << nr << std::endl;
		std::cout << "normals loaded : " << this->normals.size() << std::endl;
		std::cout << "faces loaded : " << fr << std::endl;
		std::cout << "Load OBJ Complete" << std::endl;
	}
	return;
}

std::string  dotObj::exportToFile(std::string name, bool overwrite)
{
	std::ofstream outfile;
	//Current time
	//time_t t = time(0);
	//struct tm*  now = localtime(&t);
	//std::string output = "model" + std::to_string(now->tm_year) + "-" + std::to_string(now->tm_mon) + "-" + std::to_string(now->tm_mday) + "-" + std::to_string(now->tm_hour) + "-" + std::to_string(now->tm_min) + "-" + std::to_string(now->tm_sec) + ".obj";
	//Produce new OBJ
	std::string output = name + ".obj";
	if (overwrite) { outfile.open(output); }
	else { outfile.open(output, std::ofstream::app); }
	for (int i = 0; i < this->vertices.size(); i++)
	{
		outfile << "v" << " " << std::setprecision(PRECISION) << this->vertices.at(i).at(0) << " " << std::setprecision(PRECISION) << this->vertices.at(i).at(1) << " " << std::setprecision(PRECISION) << this->vertices.at(i).at(2) << "\n";
	}
	for (int i = 0; i < this->normals.size(); i++)
	{
		outfile << "vn" << " " << std::setprecision(PRECISION) << this->normals.at(i).at(0) << " " << std::setprecision(PRECISION) << this->normals.at(i).at(1) << " " << std::setprecision(PRECISION) << this->normals.at(i).at(2) << "\n";
	}
	for (int i = 0; i < this->faces.size(); i++)
	{
		//outfile << "f" << " " << faces.at(i).at(0) << "/" << faces.at(i).at(1) << "/" << faces.at(i).at(2) << " " << faces.at(i).at(3) << "/" << faces.at(i).at(4) << "/" << faces.at(i).at(5) << " " << faces.at(i).at(6) << "/" << faces.at(i).at(7) << "/" << faces.at(i).at(8) << std::endl;
		if (this->normals.size() > 0) {
			outfile << "f" << " " << this->faces.at(i).at(0) << "//" << this->faces.at(i).at(2) << " " << this->faces.at(i).at(3) << "//" << this->faces.at(i).at(5) << " " << this->faces.at(i).at(6) << "//" << this->faces.at(i).at(8) << "\n";
		}
		else {
			outfile << "f" << " " << this->faces.at(i).at(0) << " " << this->faces.at(i).at(3) << " " << this->faces.at(i).at(6) << "\n";
		}
	}
	for (int i = 0; i < this->lines.size(); i++)
	{
		outfile << "l" << " " << this->lines.at(i).at(0) << " " << this->lines.at(i).at(1) << " " << "\n";
	}

	outfile << "\n";
	outfile.close();

	std::cout << "OBJ Export Complete" << std::endl;

	return output;
};

std::string  dotObj::exportToXYZ(std::string name, bool overwrite)
{
	std::ofstream outfile;
	std::string output = name + ".xyz";

	if (overwrite) { outfile.open(output); }
	else { outfile.open(output, std::ofstream::app); }
	for (int i = 0; i < this->vertices.size(); i++)
	{
		outfile << std::setprecision(PRECISION) << this->vertices.at(i).at(0) << " " << std::setprecision(PRECISION) << this->vertices.at(i).at(1) << " " << std::setprecision(PRECISION) << this->vertices.at(i).at(2) << " " << "\n";
	}

	outfile.close();

	std::cout << "XYZ Export Complete" << std::endl;

	return output;
};

std::string  dotObj::exportToFile(std::string name, std::string type, bool overwrite)
{
	std::ofstream outfile;
	//Current time
	//time_t t = time(0);
	//struct tm*  now = localtime(&t);
	//std::string output = "model" + std::to_string(now->tm_year) + "-" + std::to_string(now->tm_mon) + "-" + std::to_string(now->tm_mday) + "-" + std::to_string(now->tm_hour) + "-" + std::to_string(now->tm_min) + "-" + std::to_string(now->tm_sec) + ".obj";
	//Produce new OBJ
	std::string output = name + "." + type;

	if (overwrite) { outfile.open(output); }
	else { outfile.open(output, std::ofstream::app); }
	for (unsigned int i = 0; i < this->vertices.size(); i++)
	{
		outfile << "v" << " " << std::setprecision(PRECISION) << this->vertices.at(i).at(0) << " " << std::setprecision(PRECISION) << this->vertices.at(i).at(1) << " " << std::setprecision(PRECISION) << this->vertices.at(i).at(2) << "\n";
	}
	for (unsigned int i = 0; i < this->normals.size(); i++)
	{
		outfile << "vn" << " " << std::setprecision(PRECISION) << this->normals.at(i).at(0) << " " << std::setprecision(PRECISION) << this->normals.at(i).at(1) << " " << std::setprecision(PRECISION) << this->normals.at(i).at(2) << "\n";
	}
	for (unsigned int i = 0; i < this->faces.size(); i++)
	{
		//outfile << "f" << " " << faces.at(i).at(0) << "/" << faces.at(i).at(1) << "/" << faces.at(i).at(2) << " " << faces.at(i).at(3) << "/" << faces.at(i).at(4) << "/" << faces.at(i).at(5) << " " << faces.at(i).at(6) << "/" << faces.at(i).at(7) << "/" << faces.at(i).at(8) << std::endl;

		if (this->normals.size() > 0) {
			outfile << "f" << " " << this->faces.at(i).at(0) << "//" << this->faces.at(i).at(2) << " " << this->faces.at(i).at(3) << "//" << this->faces.at(i).at(5) << " " << this->faces.at(i).at(6) << "//" << this->faces.at(i).at(8) << "\n";
		}
		else {
			outfile << "f" << " " << this->faces.at(i).at(0) << " " << this->faces.at(i).at(3) << " " << this->faces.at(i).at(6) << "\n";
		}
	}
	for (unsigned int i = 0; i < this->lines.size(); i++)
	{
		outfile << "l" << " " << this->lines.at(i).at(0) << " " << this->lines.at(i).at(1) << " " << "\n";
	}
	outfile.close();

	std::cout << "OBJ Export Complete" << std::endl;

	return output;
};

void dotObj::exportToFileSegmentedPerFace(std::string output) {
	std::vector<int> seg, v, n;
	std::vector<int>::iterator it;

	seg = this->segment_property_map;
	sort(seg.begin(), seg.end());
	seg.erase(unique(seg.begin(), seg.end()), seg.end());

	std::ofstream outfile;

	outfile.open(output);

	for (int i = 0; i < this->vertices.size(); i++) {
		outfile << "v" << " " << std::setprecision(PRECISION) << this->vertices.at(i).at(0) << " " << std::setprecision(PRECISION) << this->vertices.at(i).at(1) << " " << std::setprecision(PRECISION) << this->vertices.at(i).at(2) << std::endl;
	}

	for (int i = 0; i < this->normals.size(); i++) {
		outfile << "vn" << " " << std::setprecision(PRECISION) << this->normals.at(i).at(0) << " " << std::setprecision(PRECISION) << this->normals.at(i).at(1) << " " << std::setprecision(PRECISION) << this->normals.at(i).at(2) << std::endl;
	}

	std::string name;

	if (seg.size() > 0) {
		for (int j = 0; j < seg.size(); j++) {
			for (int i = 0; i < this->segment_property_map.size(); i++) {
				if (this->segment_property_map.at(i) == seg.at(j)) {
					this->selectedFaces.push_back(i);
					if (segment_property_name.size() > 0) {
						name = segment_property_name[i];
					}
				}
			}

			if (segment_property_name.size() > 0) {
				//outfile << "g " << segment_property_name[seg.at(j)] << std::endl;
				outfile << "g " << name << std::endl;
			}
			else {
				outfile << "g " << "obj_" << std::to_string(seg.at(j)) << std::endl;
			}

			for (int q = 0; q < this->selectedFaces.size(); q++) {
				int k = this->selectedFaces.at(q);

				int v1 = this->faces.at(k).at(0);
				int v2 = this->faces.at(k).at(3);
				int v3 = this->faces.at(k).at(6);

				int n1 = this->faces.at(k).at(2);
				int n2 = this->faces.at(k).at(5);
				int n3 = this->faces.at(k).at(8);

				outfile << "f" << " " << v1 << "//" << n1 << " " << v2 << "//" << n2 << " " << v3 << "//" << n3 << std::endl;
			}

			this->selectedFaces.clear();
		}
	}

	outfile.close();
	std::cout << "Segmented Model Exported" << std::endl;

	return;
}

void dotObj::exportToFileSegmented(std::string output) {
	std::vector<int> seg, v, n;
	std::vector<int>::iterator it;

	seg = this->segment_property_map;
	sort(seg.begin(), seg.end());
	seg.erase(unique(seg.begin(), seg.end()), seg.end());

	std::ofstream outfile;

	outfile.open(output);

	for (int i = 0; i < this->vertices.size(); i++) {
		outfile << "v" << " " << std::setprecision(PRECISION) << this->vertices.at(i).at(0) << " " << std::setprecision(PRECISION) << this->vertices.at(i).at(1) << " " << std::setprecision(PRECISION) << this->vertices.at(i).at(2) << std::endl;
	}

	for (int i = 0; i < this->normals.size(); i++) {
		outfile << "vn" << " " << std::setprecision(PRECISION) << this->normals.at(i).at(0) << " " << std::setprecision(PRECISION) << this->normals.at(i).at(1) << " " << std::setprecision(PRECISION) << this->normals.at(i).at(2) << std::endl;
	}

	std::string name;

	if (seg.size() > 0) {
		for (int j = 0; j < seg.size(); j++) {
			for (int i = 0; i < this->segment_property_map.size(); i++) {
				if (this->segment_property_map.at(i) == seg.at(j)) {
					this->selectedFaces.push_back(i);
					if (segment_property_name.size() > 0) {
						name = segment_property_name[j];
					}
				}
			}

			if (segment_property_name.size() > 0) {
				//outfile << "g " << segment_property_name[seg.at(j)] << std::endl;
				outfile << "g " << name << std::endl;
			}
			else {
				outfile << "g " << "obj_" << std::to_string(seg.at(j)) << std::endl;
			}

			for (int q = 0; q < this->selectedFaces.size(); q++) {
				int k = this->selectedFaces.at(q);

				int v1 = this->faces.at(k).at(0);
				int v2 = this->faces.at(k).at(3);
				int v3 = this->faces.at(k).at(6);

				int n1 = this->faces.at(k).at(2);
				int n2 = this->faces.at(k).at(5);
				int n3 = this->faces.at(k).at(8);

				outfile << "f" << " " << v1 << "//" << n1 << " " << v2 << "//" << n2 << " " << v3 << "//" << n3 << std::endl;
			}

			this->selectedFaces.clear();
		}
	}

	outfile.close();
	std::cout << "Segmented Model Exported" << std::endl;

	return;
}

void dotObj::exportToFileSegmented(std::string output, std::string keyword) {
	std::vector<int> seg, v, n;
	std::vector<int>::iterator it;

	seg = this->segment_property_map;
	sort(seg.begin(), seg.end());
	seg.erase(unique(seg.begin(), seg.end()), seg.end());

	std::ofstream outfile;
	int outletEnum = 0;

	outfile.open(output + ".obj");

	for (int i = 0; i < this->vertices.size(); i++) {
		outfile << "v" << " " << std::setprecision(PRECISION) << this->vertices.at(i).at(0) << " " << std::setprecision(PRECISION) << this->vertices.at(i).at(1) << " " << std::setprecision(PRECISION) << this->vertices.at(i).at(2) << std::endl;
	}

	for (int i = 0; i < this->normals.size(); i++) {
		outfile << "vn" << " " << std::setprecision(PRECISION) << this->normals.at(i).at(0) << " " << std::setprecision(PRECISION) << this->normals.at(i).at(1) << " " << std::setprecision(PRECISION) << this->normals.at(i).at(2) << std::endl;
	}

	if (seg.size() > 0) {
		for (int j = 0; j < seg.size(); j++) {
			for (int i = 0; i < this->segment_property_map.size(); i++) {
				if (this->segment_property_map.at(i) == seg.at(j)) {
					this->selectedFaces.push_back(i);
				}
			}

			std::string leftOrRight;
			int generationId = (int)((float)seg.at(j) / EXPMAXNUMOFBRANCH);
			int sidenumber = (int)((float)(seg.at(j) - generationId * EXPMAXNUMOFBRANCH) / (EXPMAXNUMOFBRANCH / 10));
			if (sidenumber == 1) { leftOrRight = "L"; }
			if (sidenumber == 2) { leftOrRight = "R"; }
			int branchId = seg.at(j) - generationId * EXPMAXNUMOFBRANCH;

			//branchId = branchId;

			std::string identStr;
			char enumc[4];
			if (generationId == 99) {
				identStr = "outlet";
				//outfile << "g " << identStr << "_" << branchId << "_" << keyword << std::endl;
				sprintf(enumc, "%04d", outletEnum);
				outfile << "g " << identStr << "_" << enumc << std::endl;
				outletEnum++;
			}
			else {
				identStr = "Generation";
				outfile << "g " << identStr << "_" << std::setfill('0') << std::setw(1) << generationId + 1 - 1 << "_" << std::setfill('0') << std::setw(4) << branchId << "_" << leftOrRight << std::endl;
			}

			//outfile << "g " << identStr << "_" << std::to_string(seg.at(j)) << "_" << keyword << std::endl;

			for (int q = 0; q < this->selectedFaces.size(); q++) {
				int k = this->selectedFaces.at(q);

				int v1 = this->faces.at(k).at(0);
				int v2 = this->faces.at(k).at(3);
				int v3 = this->faces.at(k).at(6);

				int n1 = this->faces.at(k).at(2);
				int n2 = this->faces.at(k).at(5);
				int n3 = this->faces.at(k).at(8);

				outfile << "f" << " " << v1 << "//" << n1 << " " << v2 << "//" << n2 << " " << v3 << "//" << n3 << std::endl;
			}

			this->selectedFaces.clear();
		}
	}

	outfile.close();
	std::cout << "Segmented Model Export Complete" << std::endl;

	return;
}

void dotObj::getFaceIndicesFromVertexList(std::vector<int> &v, std::vector<int> &f) {
	f.clear();
	for (int i = 0; i < this->faces.size(); i++) {
		for (int j = 0; j < v.size(); j++) {
			if ((this->faces[i][0] == v[j] + 1) ||
				(this->faces[i][3] == v[j] + 1) ||
				(this->faces[i][6] == v[j] + 1))
			{
				f.push_back(i);
			}
		}
	}
	f.erase(unique(f.begin(), f.end()), f.end());

	return;

}

void dotObj::exportToFileSegmentedPerLine(std::string fileType, std::string fileName) {
	if (fileType == "obj") {
		std::ofstream outfile;
		outfile.open(fileName+"."+fileType);
		for (int i = 0; i < this->vertices.size(); i++) {
			outfile << "v" << " " << std::setprecision(PRECISION) << this->vertices.at(i).at(0) << " " << std::setprecision(PRECISION) << this->vertices.at(i).at(1) << " " << std::setprecision(PRECISION) << this->vertices.at(i).at(2) << std::endl;
		}
		for (int i = 0; i < this->normals.size(); i++) {
			outfile << "vn" << " " << std::setprecision(PRECISION) << this->normals.at(i).at(0) << " " << std::setprecision(PRECISION) << this->normals.at(i).at(1) << " " << std::setprecision(PRECISION) << this->normals.at(i).at(2) << std::endl;
		}

		std::vector<int> seg = this->segment_property_map_per_line;
		sort(seg.begin(), seg.end());
		seg.erase(unique(seg.begin(), seg.end()), seg.end());

		this->selectedSubset.clear();

		if (seg.size() > 0) {
			for (int j = 0; j < seg.size(); j++) {
				std::cout<<std::to_string(seg.at(j)) << std::endl;

				for (int i = 0; i < this->segment_property_map_per_line.size(); i++) {
					if (this->segment_property_map_per_line.at(i) == seg.at(j)) {
						this->selectedSubset.push_back(i);
					}
				}

				outfile << "g " << "obj_" << std::to_string(seg.at(j)) << std::endl;
				for (int q = 0; q < this->selectedSubset.size(); q++) {
					int k = this->selectedSubset.at(q);

					outfile << "l" << " " << this->lines.at(k).at(0) << " " << this->lines.at(k).at(1) << " " << "\n";
				}

				this->selectedSubset.clear();
			}
		}

		outfile.close();
		std::cout << ".obj model exported successfully" << std::endl;
	}
	return;
}


void dotObj::initializeDotObj(std::string filename, bool verbose) {
	

	std::string basename, extname, tmpname;
	const std::string suffix("tmp");
	std::string::size_type idx = filename.find('.');



	if (idx == std::string::npos) {
		// file name does not contain any period
		tmpname = filename + '.' + suffix;
	}
	else {
		/* split file name into base name and extension
		 * - base name contains all characters before the period
		 * - extension contains all characters after the period
		 */
		basename = filename.substr(0, idx);
		extname = filename.substr(idx + 1);
		if (extname.empty()) {
			// contains period but no extension: append tmp
			tmpname = filename;
			tmpname += suffix;
		}
		else if (extname == suffix) {
			// replace extension tmp with xxx
			tmpname = filename;
			tmpname.replace(idx + 1, extname.size(), "xxx");
		}
		else {
			// replace any extension with tmp
			tmpname = filename;
			tmpname.replace(idx + 1, std::string::npos, suffix);
		}
	}

	// print file name and temporary name
	std::cout << basename << " => " << extname << std::endl;




	

	if (extname == "off") {
		std::string line;
		std::vector<int> positions;
		std::vector<std::string> nums;
		std::string s;
		bool v = false;
		bool f = false;

		std::ifstream fileInput(filename);

		this->segment_property_map.clear();
		this->vertices.clear();
		this->normals.clear();
		this->faces.clear();

		int numberOfLines = 0;
		while (!fileInput.eof()) {
			getline(fileInput, line);
			numberOfLines++;
		}
		fileInput.clear(); // clear bad state after eof
		fileInput.seekg(0);
		if (verbose) { std::cout << "Reading off file" << std::endl; }
		getline(fileInput, line);
		if (line != "OFF") {
			std::cout << "Reading OFF file error" << std::endl;
			return;
		}

		getline(fileInput, line);
		positions.clear();
		for (int i = 0; i < line.length(); ++i)
		{
			if (line[i] == ' ') {
				positions.push_back(i);
			}
		}

		while (!fileInput.eof()) {
			getline(fileInput, line);
			positions.clear();
			nums.clear();
			f = false;
			v = false;
			positions.push_back(0);
			for (int i = 0; i < line.length(); ++i)
			{
				if (line[i] == ' ') {
					positions.push_back(i);
				}
			}
			positions.push_back(line.length());

			for (int i = 0; i < positions.size() - 1; i++) {
				int p1 = positions[i];
				int p2 = positions[i + 1];
				s = line.substr(p1, p2 - p1);
				nums.push_back(s);
			}

			for (int i = 0; i < nums.size(); i++) {
				if (nums[i].find('.') != nums[i].npos) {
					v = true;
					f = false;
				}
			}

			if (!v && line[0] == '3' && line[1] == ' ') {
				f = true;
				v = false;
			}

			if (v) {
				this->vertices.push_back({ (float)atof(nums[0].c_str()), (float)atof(nums[1].c_str()), (float)atof(nums[2].c_str()) });
			}

			if (f) {
				this->faces.push_back({ 1 + (int)atof(nums[1].c_str()), 0, 0, 1 + (int)atof(nums[2].c_str()), 0, 0, 1 + (int)atof(nums[3].c_str()), 0, 0 });
			}
		}

		std::cout << "Reading OFF file complete" << std::endl;
		return;
	}


	if (extname == "obj") {
		std::cout << "Reading .obj file" << std::endl;
		std::ifstream myfile(filename);
		float output[9];
		std::ostringstream exportedStream;
		std::string line;
		std::string lineVals;
		std::string lineVals2;
		std::string val;
		std::string val0, val1, val2, val3, val4, val5, val6, val7, val8;
		int totlength = 0;
		for (int k = 0; k < 9; k++) {
			output[k] = 0.0f;
		}

		int group = -1;
		this->segment_property_map.clear();
		this->vertices.clear();
		this->normals.clear();
		this->faces.clear();

		if (myfile.is_open())
		{
			if (verbose) { std::cout << "Building obj" << std::endl; }

			while (!myfile.eof())
			{
				getline(myfile, line);
				if (line[0] == 'v' && line[1] == ' ') {
					exportedStream << line;

					lineVals = line.substr(2);
					val0 = lineVals.substr(0, lineVals.find(' '));
					output[0] = (float)atof(val0.c_str());

					lineVals = lineVals.substr(val0.length() + 1);
					val1 = lineVals.substr(0, lineVals.find(' '));
					output[1] = (float)atof(val1.c_str());

					lineVals = lineVals.substr(val1.length() + 1);
					val2 = lineVals;
					output[2] = (float)atof(val2.c_str());

					this->vertices.push_back({ output[0], output[1], output[2], 0.0f });

					mVertex vert;
					vert.position = Vector3f(output[0], output[1], output[2]);
					this->mVertices.push_back(vert);

					//std::cout  << " " << "v" << " " << output[0] << " " << output[1] << " " << output[2] << std::endl;
				}

				/*if (line[0] == 'v' && line[1] == 'n') {
					exportedStream << line;
					lineVals = line.substr(3);
					val0 = lineVals.substr(0, lineVals.find(' '));
					output[0] = (float)atof(val0.c_str());
					val1 = lineVals.substr(val0.length() + 1, lineVals.find(' '));
					output[1] = (float)atof(val1.c_str());
					val2 = lineVals.substr(lineVals.find_last_of(' ') + 1);
					output[2] = (float)atof(val2.c_str());
					this->normals.push_back({ output[0], output[1], output[2] });
					//std::cout  << " " << "vn" << " " << output[0] << " " << output[1] << " " << output[2] << std::endl;
				}*/
				/*if (line[0] == 'g' && line[1] == ' ') {
					exportedStream << line;
					group = group + 1;
					lineVals = line.substr(2);
					segment_property_name.push_back(lineVals);
				}*/



				if (line[0] == 'f' && line[1] == ' ') {
					exportedStream << line;
					for (int i = 0; i < 9; i++) { output[i] = 0; }

					//First Set
					lineVals = line.substr(2);
					//Term #0
					val0 = lineVals.substr(0, lineVals.find('/', 0));
					output[0] = (float)atof(val0.c_str());
					totlength = val0.length();
					lineVals = lineVals.substr(totlength + 1);
					//Term #1
					val1 = lineVals.substr(0, lineVals.find('/', 0));
					output[1] = (float)atof(val1.c_str());
					totlength = val1.length();
					lineVals = lineVals.substr(totlength + 1);
					//Term #2
					val2 = lineVals.substr(0, lineVals.find(' '));
					output[2] = (float)atof(val2.c_str());
					totlength = val2.length();
					lineVals = lineVals.substr(totlength + 1);

					//Term #3
					val3 = lineVals.substr(0, lineVals.find('/', 0));
					output[3] = (float)atof(val3.c_str());
					totlength = val3.length();
					lineVals = lineVals.substr(totlength + 1);
					//Term #4
					val4 = lineVals.substr(0, lineVals.find('/', 0));
					output[4] = (float)atof(val4.c_str());
					totlength = val4.length();
					lineVals = lineVals.substr(totlength + 1);
					//Term #5
					val5 = lineVals.substr(0, lineVals.find(' ', 0));
					output[5] = (float)atof(val5.c_str());
					totlength = val5.length();
					lineVals = lineVals.substr(totlength + 1);
					//Term #6
					val6 = lineVals.substr(0, lineVals.find('/', 0));
					output[6] = (float)atof(val6.c_str());
					totlength = val6.length();
					lineVals = lineVals.substr(totlength + 1);
					//Term #7
					val7 = lineVals.substr(0, lineVals.find('/', 0));
					output[7] = (float)atof(val7.c_str());
					totlength = val7.length();
					lineVals = lineVals.substr(totlength + 1);
					//Term #8
					val8 = lineVals.substr(0, lineVals.find(' ', 0));
					output[8] = (float)atof(val8.c_str());

					this->faces.push_back({ (int)output[0], (int)output[1], (int)output[2], (int)output[3], (int)output[4], (int)output[5], (int)output[6], (int)output[7], (int)output[8] });

					//this->segment_property_map.push_back(group);
					//if (group >= 0) { this->segment_property_map.push_back(group); }
					//std::cout <<  " " << "f" << " " << output[0] << "/" << output[1] << "/" << output[2] << " " << output[3] << "/" << output[4] << "/" << output[5] << " " << output[6] << "/" << output[7] << "/" << output[8] << std::endl;
				}

				/*if (line[0] == 'l' && line[1] == ' ') {
					exportedStream << line;
					lineVals = line.substr(2);
					val0 = lineVals.substr(0, lineVals.find(" "));
					output[0] = (float)atof(val0.c_str());
					val1 = lineVals.substr(val0.length() + 1);
					//val1 = lineVals.substr(val0.length() + 1, lineVals.find(' '));
					output[1] = (float)atof(val1.c_str());
					this->lines.push_back({ (int)output[0], (int)output[1] });
				}*/
			}
		}

		/*if (this->segment_property_map.size() > 0) {
			//WritePerVertex
			this->segment_property_map_per_vertex.clear();
			for (int k = 0; k < this->vertices.size(); k++) {
				this->segment_property_map_per_vertex.push_back(0);
			}
			for (int k = 0; k < this->faces.size(); k++) {
				int v0, v1, v2, group;
				group = this->segment_property_map.at(k);
				v0 = this->faces.at(k).at(0) - 1;
				v1 = this->faces.at(k).at(3) - 1;
				v2 = this->faces.at(k).at(6) - 1;
				this->segment_property_map_per_vertex.at(v0) = group;
				this->segment_property_map_per_vertex.at(v1) = group;
				this->segment_property_map_per_vertex.at(v2) = group;
			}
		}*/

		

		if (verbose) {
			std::cout << "vertices loaded : " << this->vertices.size() << std::endl;
			std::cout << "normals loaded : " << this->normals.size() << std::endl;
			std::cout << "faces loaded : " << this->faces.size() << std::endl;
			std::cout << "Load OBJ Complete" << std::endl;
		}
		return;
	}
	if (extname == "ply") {

		std::ifstream myfile(filename);
		float output[9];
		std::string line;
		std::string lineVals;
		std::string val;
		std::string val0, val1, val2, val3, val4, val5, val6, val7, val8;
		int nr = 0;
		int fr = 0;
		int totlength = 0;
		for (int k = 0; k < 9; k++) {
			output[k] = 0.0f;
		}
		bool isvert = false;

		std::ofstream exported;

		std::cout << "Open file" << std::endl;
		if (myfile.is_open())
		{
			if (verbose) { std::cout << "Building obj" << std::endl; }

			while (!myfile.eof())
			{
				getline(myfile, line);

				if (((line[0] >= '0' && line[0] <= '9') || line[0] == '-')) {
					isvert = false;
					for (int i = 0; i < strlen(line.c_str()); i++) {
						if (line[i] == '.') { isvert = true; }
					}
					if (isvert) {
						lineVals = line;// .substr(0);
						//std::cout << "line is"<<lineVals << std::endl;
						val0 = lineVals.substr(0, lineVals.find(' '));
						output[0] = (float)atof(val0.c_str());

						lineVals = lineVals.substr(val0.length() + 1);
						val1 = lineVals.substr(0, lineVals.find(' '));
						output[1] = (float)atof(val1.c_str());

						lineVals = lineVals.substr(val1.length() + 1);
						val2 = lineVals;
						output[2] = (float)atof(val2.c_str());

						this->vertices.push_back({ output[0], output[1], output[2], 0.0f });
						nr++;
						
					}

					//Faces
					else if (!isvert) {
						for (int i = 0; i < 9; i++) { output[i] = 0; }

						//First Set
						lineVals = line.substr(2);
						//std::cout << "line is"<<lineVals << std::endl;
						val0 = lineVals.substr(0, lineVals.find(' '));
						output[0] = (float)atof(val0.c_str()) + 1;
						//std::cout << val0 << std::endl;
						lineVals = lineVals.substr(val0.length() + 1);
						val1 = lineVals.substr(0, lineVals.find(' '));
						output[3] = (float)atof(val1.c_str()) + 1;
						//std::cout << val1 << std::endl;
						lineVals = lineVals.substr(val1.length() + 1);
						val2 = lineVals;
						output[6] = (float)atof(val2.c_str()) + 1;
						//std::cout << val2 << std::endl;

						fr++;
					
						this->faces.push_back({ (int)output[0], (int)output[1], (int)output[2], (int)output[3], (int)output[4], (int)output[5], (int)output[6], (int)output[7], (int)output[8] });
					}
				}
			}
		}
		
		if (verbose) {
			std::cout << "vertices loaded : " << nr << std::endl;
			std::cout << "normals loaded : " << this->normals.size() << std::endl;
			std::cout << "faces loaded : " << fr << std::endl;
			std::cout << "Load OBJ Complete" << std::endl;
		}
		return;
	}
	







	return;
}
