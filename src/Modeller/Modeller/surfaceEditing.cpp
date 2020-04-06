#include "lungModelling.h"




void dotObj::simplificationEdgeCollapse(int numberOfEdges){

	PolyhedronSimpleCartesian surface_mesh;

	std::stringstream is=this->toOFF();
	is >> surface_mesh;

	if (!CGAL::is_triangle_mesh(surface_mesh)){
		std::cerr << "Input geometry is not triangulated." << std::endl;
		return;
	}
	// This is a stop predicate (defines when the algorithm terminates).
	// In this example, the simplification stops when the number of undirected edges
	// left in the surface mesh drops below the specified number (1000)
	CGAL::Surface_mesh_simplification::Count_stop_predicate<PolyhedronSimpleCartesian> stop(numberOfEdges);

	// This the actual call to the simplification algorithm.
	// The surface mesh and stop conditions are mandatory arguments.
	// The index maps are needed because the vertices and edges
	// of this surface mesh lack an "id()" field.
	int r = CGAL::Surface_mesh_simplification::edge_collapse
		(surface_mesh
		, stop
		, CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index, surface_mesh))
		.halfedge_index_map(get(CGAL::halfedge_external_index, surface_mesh))
		.get_cost(CGAL::Surface_mesh_simplification::Edge_length_cost <PolyhedronSimpleCartesian>())
		.get_placement(CGAL::Surface_mesh_simplification::Midpoint_placement<PolyhedronSimpleCartesian>())
		);

	std::cout << "\nFinished...\n" << r << " edges removed.\n"
		<< (surface_mesh.size_of_halfedges() / 2) << " final edges.\n";

	std::ofstream os("out.off"); 
	os << surface_mesh;

	dotObj subsampled;
	for (PolyhedronSimpleCartesian::Vertex_const_iterator vi = surface_mesh.vertices_begin(); vi != surface_mesh.vertices_end(); ++vi) {
		subsampled.vertices.push_back({ (float)vi->point().x(), (float)vi->point().y(), (float)vi->point().z() });
        /*writer.write_vertex( ::CGAL::to_double( vi->point().x()),
                             ::CGAL::to_double( vi->point().y()),
                             ::CGAL::to_double( vi->point().z()));*/
    }
	typedef CGAL::Inverse_index< PolyhedronSimpleCartesian::Vertex_const_iterator> Index;
	Index index(surface_mesh.vertices_begin(), surface_mesh.vertices_end());
	for (PolyhedronSimpleCartesian::Facet_const_iterator fi = surface_mesh.facets_begin(); fi != surface_mesh.facets_end(); ++fi) {
		PolyhedronSimpleCartesian::Halfedge_around_facet_const_circulator hc = fi->facet_begin();
		PolyhedronSimpleCartesian::Halfedge_around_facet_const_circulator hc_end = hc;
		std::size_t n = circulator_size(hc);
		CGAL_assertion(n >= 3);
		//writer.write_facet_begin(n);
		subsampled.faces.push_back({ 0, 0, 0, 0, 0, 0, 0, 0, 0 });
		int pp = 0;
		do {
			subsampled.faces.back().at(pp) = 1+index[PolyhedronSimpleCartesian::Vertex_const_iterator(hc->vertex())];
			pp = pp + 3;

			//writer.write_facet_vertex_index(index[PolyhedronSimpleCartesian::Vertex_const_iterator(hc->vertex())]);

			++hc;
		} while (hc != hc_end);
		//writer.write_facet_end();
	}
	(*this) = subsampled;




	/*BOOST_FOREACH(boost::graph_traits<PolyhedronSimpleCartesian>::vertex_descriptor v, CGAL::vertices(surface_mesh))
	{
		//std::cout << v->point().x() << v->point().y() << v->point().z() << std::endl;
		subsampled.vertices.push_back({ (float)v->point().x(), (float)v->point().y(), (float)v->point().z() });
	}*/
	



	/*for (PolyhedronSimpleCartesian::Facet_iterator i = surface_mesh.facets_begin(); i != surface_mesh.facets_end(); ++i) {
		PolyhedronSimpleCartesian::Halfedge_around_facet_circulator j = i->facet_begin();
		// Facets in polyhedral surfaces are at least triangles.
		CGAL_assertion(CGAL::circulator_size(j) >= 3);

		//std::cout << CGAL::circulator_size(j) << ' ';

		subsampled.faces.push_back({0,0,0,0,0,0,0,0,0});
		int pp = 0;
		do {
			subsampled.faces.back().at(pp) = std::distance(surface_mesh.vertices_begin(), j->vertex())+1;
			pp = pp + 3;
			//std::cout << ' ' << std::distance(surface_mesh.vertices_begin(), j->vertex());

		} while (++j != i->facet_begin());

		//std::cout << std::endl;

	}*/

	//this->vertices = subsampled.vertices;
	//this->faces = subsampled.faces;
	//this->recalculateNormals();

	//subsampled.exportToFile("simplifiedModel");


	


	return ;
}

void extract_k_ring(Vertex_handle v,
	int k,
	std::vector<Vertex_handle>& qv)
{
	std::map<Vertex_handle, int>  D;
	qv.push_back(v);
	D[v] = 0;
	std::size_t current_index = 0;
	int dist_v;
	while (current_index < qv.size() && (dist_v = D[qv[current_index]]) < k)
	{
		v = qv[current_index++];
		Polyhedron::Halfedge_around_vertex_circulator e(v->vertex_begin()), e_end(e);
		do {
			Vertex_handle new_v = e->opposite()->vertex();
			if (D.insert(std::make_pair(new_v, dist_v + 1)).second) {
				qv.push_back(new_v);
			}
		} while (++e != e_end);
	}
}




void dotObj::refine_fair(void){
	std::stringstream input = this->toOFF();
	Polyhedron poly;
	if (!input || !(input >> poly) || poly.empty()
		|| !CGAL::is_triangle_mesh(poly)) {
		std::cerr << "Not a valid input file." << std::endl;
		return;
	}
	std::vector<Polyhedron::Facet_handle>  new_facets;
	std::vector<Vertex_handle> new_vertices;
	CGAL::Polygon_mesh_processing::refine(poly,
		CGAL::faces(poly),
		std::back_inserter(new_facets),
		std::back_inserter(new_vertices),
		CGAL::Polygon_mesh_processing::parameters::density_control_factor(2.));
	std::ofstream refined_off("refined.off");
	refined_off << poly;
	refined_off.close();
	std::cout << "Refinement added " << new_vertices.size() << " vertices." << std::endl;
	Polyhedron::Vertex_iterator v = poly.vertices_begin();
	std::advance(v, 82/*e.g.*/);
	std::vector<Vertex_handle> region;
	extract_k_ring(v, 12/*e.g.*/, region);
	bool success = CGAL::Polygon_mesh_processing::fair(poly, region);
	std::cout << "Fairing : " << (success ? "succeeded" : "failed") << std::endl;
	std::ofstream faired_off("faired.off");
	faired_off << poly;
	faired_off.close();



	dotObj refined;
	for (Polyhedron::Vertex_const_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi) {
		refined.vertices.push_back({ (float)vi->point().x(), (float)vi->point().y(), (float)vi->point().z() });
		/*writer.write_vertex( ::CGAL::to_double( vi->point().x()),
		::CGAL::to_double( vi->point().y()),
		::CGAL::to_double( vi->point().z()));*/
	}
	typedef CGAL::Inverse_index< Polyhedron::Vertex_const_iterator> Index;
	Index index(poly.vertices_begin(), poly.vertices_end());
	for (Polyhedron::Facet_const_iterator fi = poly.facets_begin(); fi != poly.facets_end(); ++fi) {
		Polyhedron::Halfedge_around_facet_const_circulator hc = fi->facet_begin();
		Polyhedron::Halfedge_around_facet_const_circulator hc_end = hc;
		std::size_t n = circulator_size(hc);
		CGAL_assertion(n >= 3);
		//writer.write_facet_begin(n);
		refined.faces.push_back({ 0, 0, 0, 0, 0, 0, 0, 0, 0 });
		int pp = 0;
		do {
			refined.faces.back().at(pp) = 1 + index[Polyhedron::Vertex_const_iterator(hc->vertex())];
			pp = pp + 3;

			//writer.write_facet_vertex_index(index[PolyhedronSimpleCartesian::Vertex_const_iterator(hc->vertex())]);

			++hc;
		} while (hc != hc_end);
		//writer.write_facet_end();
	}
	(*this) = refined;


	return;
}



typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor     edge_descriptor;
struct halfedge2edge
{
	halfedge2edge(const Mesh& m, std::vector<edge_descriptor>& edges)
		: m_mesh(m), m_edges(edges)
	{}
	void operator()(const halfedge_descriptor& h) const
	{
		m_edges.push_back(edge(h, m_mesh));
	}
	const Mesh& m_mesh;
	std::vector<edge_descriptor>& m_edges;
};

namespace PMP = CGAL::Polygon_mesh_processing;

void dotObj::isotropicRemeshing(void){
	std::stringstream input = this->toOFF();
	Mesh mesh;
	if (!input || !(input >> mesh) || !CGAL::is_triangle_mesh(mesh)) {
		std::cerr << "Not a valid input file." << std::endl;
		return ;
	}
	double target_edge_length = 0.04;
	unsigned int nb_iter = 3;
	std::cout << "Split border...";
	std::vector<edge_descriptor> border;
	PMP::border_halfedges(CGAL::faces(mesh),
		mesh,
		boost::make_function_output_iterator(halfedge2edge(mesh, border)));
	PMP::split_long_edges(border, target_edge_length, mesh);
	std::cout << "done." << std::endl;
	std::cout << "Start remeshing of model" 
		<< " (" << num_faces(mesh) << " faces)..." << std::endl;
	PMP::isotropic_remeshing(
		CGAL::faces(mesh),
		target_edge_length,
		mesh,
		PMP::parameters::number_of_iterations(nb_iter)
		.protect_constraints(true));
	std::ofstream refined_off("isomeshed.off");
	refined_off << mesh;
	refined_off.close();


	return;
}