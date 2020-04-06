#ifndef _CGALDEPENDENCIES_
#define _CGALDEPENDENCIES_


#include <CGAL/bilateral_smooth_point_set.h>

#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/properties_Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

#include <CGAL/compute_average_spacing.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <boost/function_output_iterator.hpp>

#include <CGAL/grid_simplify_point_set.h>

#include <CGAL/hierarchy_simplify_point_set.h>

#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
#include <CGAL/IO/read_xyz_points.h>

#include <CGAL/IO/write_xyz_points.h>

#include <CGAL/mesh_segmentation.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Mean_curvature_flow_skeletonization.h>

#include <CGAL/pca_estimate_normals.h>
#include <CGAL/property_map.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/property_map.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>


#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/tags.h>
#include <CGAL/wlop_simplify_and_regularize_point_set.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
#include <CGAL/Polygon_mesh_processing/refine.h>
#include <CGAL/Polygon_mesh_processing/fair.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/boost/graph/selection.h>



#include <utility> // defines std::pair
#include <list>
#include <fstream>

//#define CGAL_LINKED_WITH_TBB

// Concurrency
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif


//NORMAL ORIENTATION
typedef CGAL::Exact_predicates_inexact_constructions_kernel KernelEx;
//typedef CGAL::Simple_cartesian<double> KernelEx;
typedef KernelEx::Point_3 PointEx;
typedef KernelEx::Vector_3 VectorEx;
typedef std::pair<PointEx, VectorEx> PointVectorPair;
//POISSON
typedef KernelEx::FT FT;
typedef CGAL::Point_with_normal_3<KernelEx> Point_with_normal;
typedef KernelEx::Sphere_3 Sphere;
typedef std::vector<Point_with_normal> PointList;
typedef CGAL::Polyhedron_3<KernelEx> PolyhedronEx;
typedef CGAL::Poisson_reconstruction_function<KernelEx> Poisson_reconstruction_function;
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
typedef CGAL::Implicit_surface_3<KernelEx, Poisson_reconstruction_function> Surface_3;
typedef CGAL::Polyhedron_3<KernelEx, CGAL::Polyhedron_items_with_id_3>  Polyhedron;
typedef CGAL::Surface_mesh<KernelEx::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::edge_descriptor edge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef Polyhedron::Vertex_handle Vertex_handle;



////SIMPLE_CARTESIAN
//SKELETONIZATION
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Point_3 Vector;
typedef CGAL::Surface_mesh<Point> Triangle_mesh;
typedef boost::graph_traits<Triangle_mesh>::vertex_descriptor vertex_descriptor;
typedef CGAL::Mean_curvature_flow_skeletonization<Triangle_mesh> Skeletonization;
typedef Skeletonization::Skeleton Skeleton;
typedef Skeleton::vertex_descriptor Skeleton_vertex;
typedef Skeleton::edge_descriptor Skeleton_edge;
typedef std::pair<Point, Vector> PointVectorPairSimpleCartecian;
//SKELETONIZATION BASED SEGMENTATION
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> PolyhedronSimpleCartesianItemsWithID3;
typedef CGAL::Polyhedron_3<Kernel> PolyhedronSimpleCartesian;
typedef boost::graph_traits<PolyhedronSimpleCartesianItemsWithID3>::vertex_descriptor vertex_descriptor_PS;
typedef boost::graph_traits<PolyhedronSimpleCartesianItemsWithID3>::halfedge_descriptor halfedge_descriptor_PS;
typedef boost::graph_traits<PolyhedronSimpleCartesianItemsWithID3>::face_descriptor face_descriptor_PS;
typedef CGAL::Mean_curvature_flow_skeletonization<PolyhedronSimpleCartesianItemsWithID3> SkeletonizationPolyhedronSimpleCartesian;
typedef SkeletonizationPolyhedronSimpleCartesian::Skeleton Skeleton_PS;
typedef Skeleton_PS::vertex_descriptor Skeleton_vertex_PS;



template<class ValueType> struct Facet_with_id_pmap
	: public boost::put_get_helper<ValueType&,
	Facet_with_id_pmap<ValueType> >
{
	typedef Polyhedron::Facet_const_handle key_type;
	typedef ValueType value_type;
	typedef value_type& reference;
	typedef boost::lvalue_property_map_tag category;

	Facet_with_id_pmap(
		std::vector<ValueType>& internal_vector
		) : internal_vector(internal_vector) { }

	reference operator[](key_type key) const
	{
		return internal_vector[key->id()];
	}
private:
	std::vector<ValueType>& internal_vector;
};


// Property map associating a facet with an integer as id to an
// element in a vector stored internally
template<class ValueType> struct Facet_with_id_pmap_s : public boost::put_get_helper<ValueType&, Facet_with_id_pmap_s<ValueType> >
{
	typedef face_descriptor_PS key_type;
	typedef ValueType value_type;
	typedef value_type& reference;
	typedef boost::lvalue_property_map_tag category;
	Facet_with_id_pmap_s(std::vector<ValueType>& internal_vector) : internal_vector(internal_vector) { }

	reference operator[](key_type key) const
	{
		return internal_vector[key->id()];
	}
private:
	std::vector<ValueType>& internal_vector;
};






struct Vector_pmap_wrapper{
	std::vector<bool>& vect; Vector_pmap_wrapper(std::vector<bool>& v) : vect(v) {}
	friend bool get(const Vector_pmap_wrapper& m, face_descriptor f)
	{
		return m.vect[f];
	}
	friend void put(const Vector_pmap_wrapper& m, face_descriptor f, bool b)
	{
		m.vect[f] = b;
	}
};






#endif