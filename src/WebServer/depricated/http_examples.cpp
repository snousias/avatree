#include "modellerToViewer.h"

using namespace std;
//Added for the json-example:
using namespace boost::property_tree;

//Added for the default_resource example

int main_server() {
	//HTTP-server at port 8080 using 1 thread
	//Unless you do more heavy non-threaded processing in the resources,
	//1 thread is usually faster than several threads

	HttpServer server;
	std::string str;
	dotObj buffer;
	currentProjectUniversal = new project();
	currentProjectUniversal->lungmodel.push_back(buffer);
	server.config.port = DEFPORT;
	currentProjectUniversal->tempModel = new dotObj();

	//GET-example for the path /info
	//Responds with request-information
	server.resource["^/info$"]["GET"] = [](shared_ptr<HttpServer::Response> response, shared_ptr<HttpServer::Request> request) {
		stringstream content_stream;
		content_stream << "<h1>Request from " << request->remote_endpoint_address << " (" << request->remote_endpoint_port << ")</h1>";
		content_stream << request->method << " " << request->path << " HTTP/" << request->http_version << "<br>";
		for (auto& header : request->header) {
			content_stream << header.first << ": " << header.second << "<br>";
		}

		//find length of content_stream (length received using content_stream.tellp())
		content_stream.seekp(0, ios::end);

		*response << "HTTP/1.1 200 OK\r\nContent-Length: " << content_stream.tellp() << "\r\n\r\n" << content_stream.rdbuf();
	};

	/*server.resource["^/uid/([0-9]+)/pid/([0-9]+)$"]["GET"] = [&server](shared_ptr<HttpServer::Response> response, shared_ptr<HttpServer::Request> request){
		try {
			string uidstr = request->path_match[1];
			std::cout << "User ID is :" << uid << std::endl;
			uid = stoi(uidstr);

			string pidstr = request->path_match[2];
			std::cout << "Patient ID is :" << pid << std::endl;
			pid = stoi(pidstr);
			auto ifs = make_shared<ifstream>();
			if (*ifs) {
				auto length = ifs->tellg();
				ifs->seekg(0, ios::beg);
				*response << "HTTP/1.1 200 OK\r\n" << "Server Online!" << "Content-Length: " << 14 << "\r\n\r\n";
				default_resource_send(server, response, ifs);
			}
			else
				throw invalid_argument("could not read file");
		}
		catch (const exception &e) {
			string content = "Could not open path " + request->path + ": " + e.what();
			*response << "HTTP/1.1 400 Bad Request\r\nContent-Length: " << content.length() << "\r\n\r\n" << content;
		}
	};*/

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	server.resource["^/([a-zA-Z0-9?=&.]+)/modelUpload"]["POST"] = [](shared_ptr<HttpServer::Response> response, shared_ptr<HttpServer::Request> request) {
		thread work_thread([response, request] {
			this_thread::sleep_for(chrono::seconds(1));
			string sessionID = request->path_match[0];
			sessionID = sessionID.substr(sessionID.find("/", 0) + 1, sessionID.find("/", 1) - 1);
			std::cout << "Session ID is " << sessionID << std::endl; //SHOW SESSION ID
			project * currentProject = new project();
			std::string filepath = "uploaded.obj";
			std::ofstream outFile;
			outFile.open(filepath);
			stringstream ss;
			outFile << request->content.rdbuf();
			outFile.close();
			dotObj buffer;
			buffer.initializeFromFile(filepath);

			currentProject->lungmodel.push_back(buffer);
			//lungmodel.at(lungmodel.size() - 1)= buffer;
			currentProject->tempModel = new dotObj();
			*currentProject->tempModel = buffer;
			string content = "Model Uploaded";
			currentProject->token = sessionID;

			//projectsList.push_back(*currentProject);

			if (projectsList.size() == 0){
				projectsList.push_back(*currentProject);
			}
			else if (projectsList.size() < 10){
				projectsList.push_back(*currentProject);
			}
			else{
				projectsList.erase(projectsList.begin());
				projectsList.push_back(*currentProject);
			}

			std::cout << "Current List Size : " << projectsList.size() << std::endl;

			*response << "HTTP/1.1 200 OK\r\nContent-Length: " << content.length() << "\r\n\r\n" << content;
		});
		work_thread.detach();
	};

	server.resource["^/([a-zA-Z0-9?=&.]+)/extendSelection$"]["GET"] = [&server](shared_ptr<HttpServer::Response> response, shared_ptr<HttpServer::Request> request) {
		thread work_thread([response, request] {
			this_thread::sleep_for(chrono::seconds(1));

			string sessionID = request->path_match[0];
			sessionID = sessionID.substr(sessionID.find("/", 0) + 1, sessionID.find("/", 1) - 1);
			std::cout << "Session ID is " << sessionID << std::endl; //SHOW SESSION ID
			project * currentProject = new project();
			bool initProcess = false;
			/*for (int i = 0; i < projects.size(); i++){
				if (projects.at(i).token == sessionID){
				currentProject = &projects.at(i);
				std::cout << "Accessing list item # " << i << std::endl;
				}
				}*/

			for (std::list<project>::iterator it = projectsList.begin(); it != projectsList.end(); ++it){
				if (it->token == sessionID){
					currentProject = &*it;
					initProcess = true;
					//std::cout << "Accessing list item # " << i << std::endl;
				}
			}

			if ((currentProject->lungmodel.size() > 0) && (initProcess)){
				if (currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices.size() > 0){
					currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).extendSelection();
				}

				std::vector<int> selectedFaces;
				for (int i = 0; i < currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).faces.size(); i++){
					for (int j = 0; j < currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices.size(); j++){
						if ((currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).faces.at(i).at(0) == currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices.at(j) + 1) ||
							(currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).faces.at(i).at(3) == currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices.at(j) + 1) ||
							(currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).faces.at(i).at(6) == currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices.at(j) + 1))
						{
							selectedFaces.push_back(i);
						}
					}
				}
				selectedFaces.erase(unique(selectedFaces.begin(), selectedFaces.end()), selectedFaces.end());
				std::ostringstream out;
				for (int i = 0; i < selectedFaces.size(); i++){
					out << selectedFaces[i] << "\r\n";
				}
				std::string str = out.str();
				*response << "HTTP/1.1 200 OK\r\nContent-Length: " << str.length() << "\r\n\r\n" << str;
			}
			else{
				*response << "HTTP/1.1 200 OK\r\nContent-Length: " << 4 << "\r\n\r\n" << "Done";
			}
		});
		work_thread.detach();
	};

	//GET-example for the path /match/[number], responds with the matched string in path (number)
	//For instance a request GET /match/123 will receive: 123
	server.resource["^/([a-zA-Z0-9?=&.]+)/fetchModel$"]["GET"] = [&server](shared_ptr<HttpServer::Response> response, shared_ptr<HttpServer::Request> request) {
		thread work_thread([response, request] {
			this_thread::sleep_for(chrono::seconds(1));

			string sessionID = request->path_match[0];
			sessionID = sessionID.substr(sessionID.find("/", 0) + 1, sessionID.find("/", 1) - 1);
			std::cout << "Session ID is " << sessionID << std::endl; //SHOW SESSION ID
			project * currentProject = new project();
			bool initProcess = false;
			/*for (int i = 0; i < projects.size(); i++){
			if (projects.at(i).token == sessionID){
			currentProject = &projects.at(i);
			std::cout << "Accessing list item # " << i << std::endl;
			}
			}*/
			for (std::list<project>::iterator it = projectsList.begin(); it != projectsList.end(); ++it){
				if (it->token == sessionID){
					currentProject = &*it;
					initProcess = true;
					//std::cout << "Accessing list item # " << i << std::endl;
				}
			}

			if ((currentProject->lungmodel.size() > 0) && (initProcess)){
				std::ostringstream output;
				//dotObj m = currentProject->lungmodel.at(currentProject->lungmodel.size() - 1);
				dotObj m = *currentProject->tempModel;
				for (int i = 0; i < m.vertices.size(); i++){
					output << fixed << "v " << std::setprecision(9) << m.vertices.at(i).at(0) << " " << std::setprecision(9) << m.vertices.at(i).at(1) << " " << std::setprecision(9) << m.vertices.at(i).at(2) << "\r\n";
				}

				for (int i = 0; i < m.normals.size(); i++){
					output << fixed << "vn " << std::setprecision(9) << m.normals.at(i).at(0) << " " << std::setprecision(9) << m.normals.at(i).at(1) << " " << std::setprecision(9) << m.normals.at(i).at(2) << "\r\n";
				}

				for (int i = 0; i < m.faces.size(); i++){
					output << "f " << m.faces.at(i).at(0) << "//" << m.faces.at(i).at(2) << " " << m.faces.at(i).at(3) << "//" << m.faces.at(i).at(5) << " " << m.faces.at(i).at(6) << "//" << m.faces.at(i).at(8) << "\r\n";
				}
				std::string str = output.str();
				*response << "HTTP/1.1 200 OK\r\nContent-Length: " << str.length() << "\r\n\r\n" << str;
			}
			else{
				std::string str = "No model";
				*response << "HTTP/1.1 200 OK\r\nContent-Length: " << str.length() << "\r\n\r\n" << str;
			}
		});
		work_thread.detach();
	};

	//GET-example for the path /match/[number], responds with the matched string in path (number)
	//For instance a request GET /match/123 will receive: 123
	server.resource["^/([a-zA-Z0-9?=&.]+)/fetchModelSegmented$"]["GET"] = [&server](shared_ptr<HttpServer::Response> response, shared_ptr<HttpServer::Request> request) {
		thread work_thread([response, request] {
			this_thread::sleep_for(chrono::seconds(1));

			string sessionID = request->path_match[0];
			sessionID = sessionID.substr(sessionID.find("/", 0) + 1, sessionID.find("/", 1) - 1);
			std::cout << "Session ID is " << sessionID << std::endl; //SHOW SESSION ID
			project * currentProject = new project();
			bool initProcess = false;
			for (std::list<project>::iterator it = projectsList.begin(); it != projectsList.end(); ++it){
				if (it->token == sessionID){
					currentProject = &*it;
					initProcess = true;
					//std::cout << "Accessing list item # " << i << std::endl;
				}
			}

			if ((currentProject->lungmodel.size() > 0) && (initProcess)){
				std::ostringstream outfile;
				//dotObj m = currentProject->lungmodel.at(currentProject->lungmodel.size() - 1);
				dotObj m = *currentProject->tempModel;

				std::vector<int> seg, v, n;
				std::vector<int>::iterator it;

				seg = m.segment_property_map;
				sort(seg.begin(), seg.end());
				seg.erase(unique(seg.begin(), seg.end()), seg.end());

				for (int i = 0; i < m.vertices.size(); i++){
					outfile << fixed << "v" << " " << std::setprecision(9) << m.vertices.at(i).at(0) << " " << std::setprecision(9) << m.vertices.at(i).at(1) << " " << std::setprecision(9) << m.vertices.at(i).at(2) << "\r\n";
				}

				for (int i = 0; i < m.normals.size(); i++){
					outfile << fixed << "vn" << " " << std::setprecision(9) << m.normals.at(i).at(0) << " " << std::setprecision(9) << m.normals.at(i).at(1) << " " << std::setprecision(9) << m.normals.at(i).at(2) << "\r\n";
				}

				if (seg.size() > 0){
					m.selectedFaces.clear();
					for (int j = 0; j < seg.size(); j++){
						for (int i = 0; i < m.segment_property_map.size(); i++){
							if (m.segment_property_map.at(i) == seg.at(j)){
								m.selectedFaces.push_back(i);
							}
						}

						if ((m.segment_property_name.size() > 0) && (m.segment_property_name.size() == seg.size())){
							outfile << "g " << m.segment_property_name[seg.at(j)] << "\r\n";
						}
						else{
							outfile << "g " << "obj_" << std::to_string(seg.at(j)) << "\r\n";
						}


						for (int q = 0; q < m.selectedFaces.size(); q++){
							int k = m.selectedFaces.at(q);

							int v1 = m.faces.at(k).at(0);
							int v2 = m.faces.at(k).at(3);
							int v3 = m.faces.at(k).at(6);

							int n1 = m.faces.at(k).at(2);
							int n2 = m.faces.at(k).at(5);
							int n3 = m.faces.at(k).at(8);

							outfile << "f" << " " << v1 << "//" << n1 << " " << v2 << "//" << n2 << " " << v3 << "//" << n3 << "\r\n";
						}

						m.selectedFaces.clear();
					}
				}

				std::string str = outfile.str();
				*response << "HTTP/1.1 200 OK\r\nContent-Length: " << str.length() << "\r\n\r\n" << str;
			}
			else{
				std::string str = "No model";
				*response << "HTTP/1.1 200 OK\r\nContent-Length: " << str.length() << "\r\n\r\n" << str;
			}
		});
		work_thread.detach();
	};

	server.resource["^/([a-zA-Z0-9?=&.]+)/reFetchModel$"]["GET"] = [&server](shared_ptr<HttpServer::Response> response, shared_ptr<HttpServer::Request> request) {
		thread work_thread([response, request] {
			this_thread::sleep_for(chrono::seconds(2));

			string sessionID = request->path_match[0];
			sessionID = sessionID.substr(sessionID.find("/", 0) + 1, sessionID.find("/", 1) - 1);
			std::cout << "Session ID is " << sessionID << std::endl; //SHOW SESSION ID
			project * currentProject = new project();
			bool initProcess = false;
			/*for (int i = 0; i < projects.size(); i++){
			if (projects.at(i).token == sessionID){
			currentProject = &projects.at(i);
			std::cout << "Accessing list item # " << i << std::endl;
			}
			}*/
			for (std::list<project>::iterator it = projectsList.begin(); it != projectsList.end(); ++it){
				if (it->token == sessionID){
					currentProject = &*it;
					initProcess = true;
					//std::cout << "Accessing list item # " << i << std::endl;
				}
			}

			if (initProcess){
				std::ostringstream output;
				dotObj m = currentProject->lungmodel.at(currentProject->lungmodel.size() - 1);
				*currentProject->tempModel = m;
				//dotObj m = *tempModel;
				for (int i = 0; i < m.vertices.size(); i++){
					output << fixed << "v " << std::setprecision(9) << m.vertices.at(i).at(0) << " " << std::setprecision(9) << m.vertices.at(i).at(1) << " " << std::setprecision(9) << m.vertices.at(i).at(2) << "\r\n";
				}

				for (int i = 0; i < m.normals.size(); i++){
					output << fixed << "vn " << std::setprecision(9) << m.normals.at(i).at(0) << " " << std::setprecision(9) << m.normals.at(i).at(1) << " " << std::setprecision(9) << m.normals.at(i).at(2) << "\r\n";
				}

				for (int i = 0; i < m.faces.size(); i++){
					output << "f " << m.faces.at(i).at(0) << "//" << m.faces.at(i).at(2) << " " << m.faces.at(i).at(3) << "//" << m.faces.at(i).at(5) << " " << m.faces.at(i).at(6) << "//" << m.faces.at(i).at(8) << "\r\n";
				}
				std::string str = output.str();
				*response << "HTTP/1.1 200 OK\r\nContent-Length: " << str.length() << "\r\n\r\n" << str;
			}
			else{
				std::string str = "No model loaded";
				*response << "HTTP/1.1 200 OK\r\nContent-Length: " << str.length() << "\r\n\r\n" << str;
			}
		});
		work_thread.detach();
	};

	server.resource["^/([a-zA-Z0-9?=&.]+)/saveChanges$"]["GET"] = [&server](shared_ptr<HttpServer::Response> response, shared_ptr<HttpServer::Request> request) {
		thread work_thread([response, request] {
			this_thread::sleep_for(chrono::seconds(1));

			string sessionID = request->path_match[0];
			sessionID = sessionID.substr(sessionID.find("/", 0) + 1, sessionID.find("/", 1) - 1);
			std::cout << "Session ID is " << sessionID << std::endl; //SHOW SESSION ID
			project * currentProject = new project();
			bool initProcess = false;

			for (std::list<project>::iterator it = projectsList.begin(); it != projectsList.end(); ++it){
				if (it->token == sessionID){
					currentProject = &*it;
					initProcess = true;
					//std::cout << "Accessing list item # " << i << std::endl;
				}
			}

			if (initProcess){
				currentProject->lungmodel.push_back(*currentProject->tempModel);
				std::string str = "Saved";
				*response << "HTTP/1.1 200 OK\r\nContent-Length: " << str.length() << "\r\n\r\n" << str;
			}
			else{
				std::string str = "No model loaded";
				*response << "HTTP/1.1 200 OK\r\nContent-Length: " << str.length() << "\r\n\r\n" << str;
			}
		});
		work_thread.detach();
	};

	server.resource["^/([a-zA-Z0-9?=&.]+)/select([a-zA-Z0-9?=&.]+)$"]["GET"] = [&server](shared_ptr<HttpServer::Response> response, shared_ptr<HttpServer::Request> request) {
		thread work_thread([response, request] {
			this_thread::sleep_for(chrono::seconds(1));

			string sessionID = request->path_match[0];
			sessionID = sessionID.substr(sessionID.find("/", 0) + 1, sessionID.find("/", 1) - 1);
			std::cout << "Session ID is " << sessionID << std::endl; //SHOW SESSION ID
			project * currentProject = new project();
			bool initProcess = false;
			for (std::list<project>::iterator it = projectsList.begin(); it != projectsList.end(); ++it){
				if (it->token == sessionID){
					currentProject = &*it;
					initProcess = true;
				}
			}

			if (initProcess){
				string result = request->path_match[2];
				std::cout << "Query url is " << result << std::endl;
				result = result.substr(result.find("?") + 1);
				std::stringstream test(result);
				std::string segment;
				std::vector<std::string> seglist;
				while (std::getline(test, segment, '&'))
				{
					seglist.push_back(segment);
					std::cout << segment << "\n";
				}
				for (int i = 0; i < seglist.size(); i++){
					string varname = seglist[i].substr(0, seglist[i].find("="));
					string varvalue = seglist[i].substr(seglist[i].find("=") + 1);
					std::cout << varname << "|" << varvalue << "\n";

					if (varname == "index"){
						currentProject->seedPoint = atoi(varvalue.c_str());
						std::cout << "Face index:" << varvalue << "\n";
					}

					if (varname == "size")
					{
						currentProject->_brushSizeBox = atoi(varvalue.c_str());
						std::cout << "Selection size:" << varvalue << "\r\n";
					}

					if (varname == "brush")
					{
						if (varvalue == "true"){
							currentProject->brushSelectionFunctionality = true;
							currentProject->partSelectionFunctionality = false;
						}
						else{
							currentProject->brushSelectionFunctionality = false;
							currentProject->partSelectionFunctionality = true;
						}
						std::cout << "Use brush" << varvalue << "\n";
					}
				}

				if (currentProject->brushSelectionFunctionality){
					if (currentProject->lungmodel.size() > 0){
						if (currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).vertices.size() > 0){
							if (currentProject->_brushSizeBox > 0){
								currentProject->seedPoint = currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).faces[currentProject->seedPoint][0] - 1;

								currentProject->currentSelection = currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices;
								currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices.clear();
								currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selector(currentProject->seedPoint, currentProject->_brushSizeBox);
								currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices.insert(currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices.end(), currentProject->currentSelection.begin(), currentProject->currentSelection.end());
								sort(currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices.begin(), currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices.end());
								currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices.erase(unique(currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices.begin(), currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices.end()), currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices.end());

								std::vector<int> selectedFaces;
								for (int i = 0; i < currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).faces.size(); i++){
									for (int j = 0; j < currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices.size(); j++){
										if ((currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).faces.at(i).at(0) == currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices.at(j) + 1) ||
											(currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).faces.at(i).at(3) == currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices.at(j) + 1) ||
											(currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).faces.at(i).at(6) == currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices.at(j) + 1))
										{
											selectedFaces.push_back(i);
										}
									}
								}
								selectedFaces.erase(unique(selectedFaces.begin(), selectedFaces.end()), selectedFaces.end());

								std::ostringstream out;
								for (int i = 0; i < selectedFaces.size(); i++){
									out << selectedFaces[i] << "\r\n";
								}
								std::string str = out.str();
								*response << "HTTP/1.1 200 OK\r\nContent-Length: " << str.length() << "\r\n\r\n" << str;
							}
							else
							{
								currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices.clear();
								std::string str = "Clear Selection";
								*response << "HTTP/1.1 200 OK\r\nContent-Length: " << str.length() << "\r\n\r\n" << str;
							}
						}
						else{
							std::string str = "No model loaded";
							*response << "HTTP/1.1 200 OK\r\nContent-Length: " << str.length() << "\r\n\r\n" << str;
						}
					}
				}

				if (currentProject->partSelectionFunctionality){
					if (currentProject->lungmodel.size() > 0){
						if (currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).vertices.size() > 0){
							std::vector<int> v = currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).segment_property_map_per_vertex;
							if (v.size() > 0){
								int group = v.at(currentProject->seedPoint);
								currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices.clear();
								for (int i = 0; i < v.size(); i++){
									if (v.at(i) == group){
										currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices.push_back(i);
									}
								}
							}
							std::vector<int> selectedFaces;
							for (int i = 0; i < currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).faces.size(); i++){
								for (int j = 0; j < currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices.size(); j++){
									if ((currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).faces.at(i).at(0) == currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices.at(j) + 1) ||
										(currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).faces.at(i).at(3) == currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices.at(j) + 1) ||
										(currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).faces.at(i).at(6) == currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices.at(j) + 1))
									{
										selectedFaces.push_back(i);
									}
								}
							}
							selectedFaces.erase(unique(selectedFaces.begin(), selectedFaces.end()), selectedFaces.end());
							std::ostringstream out;
							for (int i = 0; i < selectedFaces.size(); i++){
								out << selectedFaces[i] << "\r\n";
							}
							std::string str = out.str();
							*response << "HTTP/1.1 200 OK\r\nContent-Length: " << str.length() << "\r\n\r\n" << str;
						}
						else{
							std::string str = "No model loaded";
							*response << "HTTP/1.1 200 OK\r\nContent-Length: " << str.length() << "\r\n\r\n" << str;
						}
					}
				}
			}
			else{
				std::string str = "No model loaded";
				*response << "HTTP/1.1 200 OK\r\nContent-Length: " << str.length() << "\r\n\r\n" << str;
			}
		});
		work_thread.detach();
	};
	//
	//
	//
	server.resource["^/([a-zA-Z0-9?=&.]+)/narrow/stats$"]["GET"] = [&server](shared_ptr<HttpServer::Response> response, shared_ptr<HttpServer::Request> request) {
		thread work_thread([response, request] {
			this_thread::sleep_for(chrono::seconds(1));

			string sessionID = request->path_match[0];
			sessionID = sessionID.substr(sessionID.find("/", 0) + 1, sessionID.find("/", 1) - 1);
			std::cout << "Session ID is " << sessionID << std::endl; //SHOW SESSION ID
			project * currentProject = new project();
			bool initProcess = false;
			/*for (int i = 0; i < projects.size(); i++){
			if (projects.at(i).token == sessionID){
			currentProject = &projects.at(i);
			std::cout << "Accessing list item # " << i << std::endl;
			}
			}*/
			for (std::list<project>::iterator it = projectsList.begin(); it != projectsList.end(); ++it){
				if (it->token == sessionID){
					currentProject = &*it;
					initProcess = true;
					//std::cout << "Accessing list item # " << i << std::endl;
				}
			}

			if (initProcess){
				try {
					std::string result = "";
					result = "{\"sessionID\":\"\",\"stats\":{\"localDIameterBeforeNarrowing\":\"" + std::to_string(currentProject->localDIameterBeforeNarrowing) + "\"" + "," +
						"\"localDIameterAfterNarrowing\":\"" + std::to_string(currentProject->localDIameterAfterNarrowing) + "\"" + "," +
						"\"narrowingRatio\":\"" + std::to_string(currentProject->narrowingRatio) + "\"" + "}}";

					*response << "HTTP/1.1 200 OK\r\nContent-Length: " << result.length() << "\r\n\r\n" << result;
				}
				catch (const exception &e) {
					string content = "Could not open path " + request->path + ": " + e.what();
					*response << "HTTP/1.1 400 Bad Request\r\nContent-Length: " << content.length() << "\r\n\r\n" << content;
				}
			}
			else{
				std::string str = "No model loaded";
				*response << "HTTP/1.1 200 OK\r\nContent-Length: " << str.length() << "\r\n\r\n" << str;
			}
		});
		work_thread.detach();
	};

	//GET-example for the path /match/[number], responds with the matched string in path (number)
	//For instance a request GET /match/123 will receive: 123
	server.resource["^/([a-zA-Z0-9?=&.]+)/narrow([a-zA-Z0-9?=&.]+)$"]["GET"] = [&server](shared_ptr<HttpServer::Response> response, shared_ptr<HttpServer::Request> request) {
		thread work_thread([response, request] {
			this_thread::sleep_for(chrono::seconds(1));

			string sessionID = request->path_match[0];
			sessionID = sessionID.substr(sessionID.find("/", 0) + 1, sessionID.find("/", 1) - 1);
			std::cout << "Session ID is " << sessionID << std::endl; //SHOW SESSION ID
			project * currentProject = new project();
			bool initProcess = false;

			for (std::list<project>::iterator it = projectsList.begin(); it != projectsList.end(); ++it){
				if (it->token == sessionID){
					currentProject = &*it;
					initProcess = true;
					//std::cout << "Accessing list item # " << i << std::endl;
				}
			}
			if (initProcess){
				string result = request->path_match[2];
				std::cout << "Query url is " << result << std::endl;
				try {
					double contractionPercentageBox = 0.8;
					int contractionStrengthBox = 12;
					bool useCustomFunction = false;
					double frequencyBox = 0.2;
					bool useInterpolation = false;
					double interpolationPercentage = 0.5;
					bool useExtendedSmoothing = true;
					int  numberOfRaysSDF = 25;
					int numberOfClustersSDF = 4;
					double lamdaSDF = 0.0001;
					result = result.substr(result.find("?") + 1);
					std::stringstream test(result);
					std::string segment;
					std::vector<std::string> seglist;
					while (std::getline(test, segment, '&'))
					{
						seglist.push_back(segment);
					}
					for (int i = 0; i < seglist.size(); i++){
						string varname = seglist[i].substr(0, result.find("="));
						string varvalue = seglist[i].substr(result.find("=") + 1);

						if (varname == "iterations"){
							contractionStrengthBox = stoi(varvalue);
							std::cout << "iterations" << varvalue << "\n";
						}
						if (varname == "percentage")
						{
							contractionPercentageBox = stof(varvalue);
							std::cout << "percentage" << varvalue << "\n";
						}
						if (varname == "smoothing")
						{
							if (varvalue == "true"){
								useExtendedSmoothing = true;
							}
							else{
								useExtendedSmoothing = false;
							}
							std::cout << "smoothing" << varvalue << "\n";
						}
					}

					std::cout << "Narrowing commencing" << std::endl;
					if (currentProject->lungmodel.size() > 0){
						if (currentProject->lungmodel.at(currentProject->lungmodel.size() - 1).selectedVertices.size() > 0){
							simulation *sim;
							sim = new simulation();
							sim->lungmodel = currentProject->lungmodel;
							sim->narrow(contractionPercentageBox, contractionStrengthBox,
								useCustomFunction, frequencyBox, useInterpolation,
								interpolationPercentage, useExtendedSmoothing,
								numberOfRaysSDF, numberOfClustersSDF, lamdaSDF);

							currentProject->localDIameterBeforeNarrowing = sim->localDIameterBeforeNarrowing;
							currentProject->localDIameterAfterNarrowing = sim->localDIameterAfterNarrowing;
							currentProject->narrowingRatio = sim->narrowingRatio;

							currentProject->tempModel = sim->tempModel;
							//lungmodel.push_back(*tempModel);//Caution
							//Reload new model
						}
						else{
							std::cout << "\n";
							std::cout << "Error no area selected" << "\n";
						}
					}

					*response << "HTTP/1.1 200 OK\r\nContent-Length: " << result.length() << "\r\n\r\n" << result;
				}
				catch (const exception &e) {
					string content = "Could not open path " + request->path + ": " + e.what();
					*response << "HTTP/1.1 400 Bad Request\r\nContent-Length: " << content.length() << "\r\n\r\n" << content;
				}
			}
			else{
				std::string str = "No model loaded";
				*response << "HTTP/1.1 200 OK\r\nContent-Length: " << str.length() << "\r\n\r\n" << str;
			}
		});
		work_thread.detach();
	};

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Add resources using path-regex and method-string, and an anonymous function
	//POST-example for the path /string, responds the posted string
	server.resource["^/string$"]["POST"] = [](shared_ptr<HttpServer::Response> response, shared_ptr<HttpServer::Request> request) {
		//Retrieve string:
		auto content = request->content.string();
		//request->content.string() is a convenience function for:
		//stringstream ss;
		//ss << request->content.rdbuf();
		//string content=ss.str();
		std::cout << content << std::endl;
		*response << "HTTP/1.1 200 OK\r\nContent-Length: " << content.length() << "\r\n\r\n" << content;
	};

	//Default GET-example. If no other matches, this anonymous function will be called.
	//Will respond with content in the web/-directory, and its subdirectories.
	//Default file: index.html
	//Can for instance be used to retrieve an HTML 5 client that uses REST-resources on this server
	server.default_resource["GET"] = [&server](shared_ptr<HttpServer::Response> response, shared_ptr<HttpServer::Request> request) {
		try {
			auto ifs = make_shared<ifstream>();
			if (*ifs) {
				auto length = ifs->tellg();
				ifs->seekg(0, ios::beg);
				*response << "HTTP/1.1 200 OK\r\n" << "Server Online!" << "Content-Length: " << 14 << "\r\n\r\n";
				//*response << "HTTP/1.1 200 OK\r\n" << cache_control << etag << "Content-Length: " << length << "\r\n\r\n";
				default_resource_send(server, response, ifs);
			}
			else
				throw invalid_argument("could not read file");
		}
		catch (const exception &e) {
			string content = "Could not open path " + request->path + ": " + e.what();
			*response << "HTTP/1.1 400 Bad Request\r\nContent-Length: " << content.length() << "\r\n\r\n" << content;
		}
	};

	thread server_thread([&server](){
		server.start();
	});

	this_thread::sleep_for(chrono::seconds(1));
	int thePort = DEFPORT;
	HttpClient client("localhost:" + std::to_string(thePort));
	auto r2 = client.request("POST", "/string", "Server online");
	server_thread.join();

	return 0;
}

void default_resource_send(const HttpServer &server, const shared_ptr<HttpServer::Response> &response,
	const shared_ptr<ifstream> &ifs) {
	//read and send 128 KB at a time
	static vector<char> buffer(131072); // Safe when server is running on one thread
	streamsize read_length;
	if ((read_length = ifs->read(&buffer[0], buffer.size()).gcount()) > 0) {
		response->write(&buffer[0], read_length);
		if (read_length == static_cast<streamsize>(buffer.size())) {
			server.send(response, [&server, response, ifs](const boost::system::error_code &ec) {
				if (!ec)
					default_resource_send(server, response, ifs);
				else
					cerr << "Connection interrupted" << endl;
			});
		}
	}
}

int main(void){
	main_server();
	//main_server_start();

	return 0;
}