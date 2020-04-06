#ifndef _MODELLERTOVIEWER_H_
#define _MODELLERTOVIEWER_H_

#include "server_http.hpp"
#include "client_http.hpp"
#define BOOST_SPIRIT_THREADSAFE
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
//Added for the default_resource example
#include <fstream>
#include <boost/filesystem.hpp>
#include <vector>
#include <lungModelling.h>
#include <algorithm>
#ifdef HAVE_OPENSSL
#include "crypto.hpp"
#endif


typedef SimpleWeb::Server<SimpleWeb::HTTP> HttpServer;
typedef SimpleWeb::Client<SimpleWeb::HTTP> HttpClient;



#include "lungModelling.h"

//extern int seedPoint;
//extern int _brushSizeBox;
//extern bool	brushSelectionFunctionality;
//extern bool	partSelectionFunctionality;
//extern dotObj *tempModel;
//extern std::vector<dotObj> lungmodel;
//extern std::vector<int> currentSelection;

extern int uid, pid;
extern std::vector<project> projects;
extern std::list<project> projectsList;
extern project *currentProjectUniversal;


using namespace std;
//Added for the json-example:
using namespace boost::property_tree;
void default_resource_send(const HttpServer &server, const shared_ptr<HttpServer::Response> &response, const shared_ptr<ifstream> &ifs);
int modelling_rest_api();

#endif