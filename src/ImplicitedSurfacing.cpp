#include "ImplicitedSurfacing.h"
#include "mymesh/readers.h"
#include "mymesh/utility.h"
#include "mymesh/my_mesh.h"
#include "mymesh/YaPLY.hpp"

#include "libplyxx/libplyxx/libplyxx.h"
#include <ctime>
#include <chrono>
#include <iomanip>
#include<unordered_map>
#include <list>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/spectral_surface_reconstruction.h>
#include <CGAL/poisson_surface_reconstruction.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/property_map.h>
#include <CGAL/Implicit_reconstruction_function.h>
//#include <CGAL/Vector_3.h>

typedef std::chrono::high_resolution_clock Clock;


static Surfacer *p_ImplicitSurfacer;

static int TriProc(int in_i1, int in_i2, int in_i3, VERTICES vs) {
    const R3Pt pt = vs.ptr[in_i1].position;

    //    bool bOutside = false;
    //    for ( int j = 0; j < 3; j++ ) {
    //        const double dWidth = s_ptMax[j] - s_ptMin[j];
    //        if ( pt[j] < s_ptMin[j] - dWidth * 0.1 ) {
    //            bOutside = true;
    //        }
    //        if ( pt[j] > s_ptMax[j] + dWidth * 0.1 ) {
    //            bOutside = true;
    //        }
    //    }
    //    if ( bOutside == false ) {
    //        s_afaceSurface.addItem( R3Pt_i( in_i1, in_i2, in_i3 ) );
    //    }

    p_ImplicitSurfacer->s_afaceSurface.addItem( R3Pt_i( in_i3, in_i2, in_i1 ) );
    return 1;
}

static void VertProc(VERTICES vs) {
    p_ImplicitSurfacer->s_aptSurface.need( vs.count );
    p_ImplicitSurfacer->s_avecSurface.need( vs.count );
    for ( int i = 0; i < vs.count; i++ ) {
        p_ImplicitSurfacer->s_aptSurface[i] = vs.ptr[i].position;
        p_ImplicitSurfacer->s_avecSurface[i] = vs.ptr[i].normal;
    }
}



void Surfacer::CalSurfacingPara(vector<double>&Vs, vector<int>&labels){

    vector<double>leftcorner(3, DBL_MAX);
    vector<double>rightcorner(3, DBL_MIN);
    vector<double>midpoint(3,0);

    int nv = Vs.size()/3;
    for(int i=0;i<nv;++i){
        auto p_v = Vs.data()+i*3;
        for(int j=0;j<3;++j){
            leftcorner[j] = min(leftcorner[j],p_v[j]);
            rightcorner[j] = max(rightcorner[j],p_v[j]);
            midpoint[j] += p_v[j];
        }
    }
    for(int j=0;j<3;++j){
        midpoint[j] /= nv;
    }

    //double width = MyUtility::_VerticesDistance(leftcorner.data(),rightcorner.data());
    double width = -1;
    for(int j=0;j<3;++j)width = max(width, fabs(rightcorner[j]-leftcorner[j]));
    //dSize = width * 0.02;
    dSize = width * 0.005;

    if(0){
        int bInd = -1;
        double mindist = DBL_MAX;
        for(int i=0;i<nv;++i)if(labels[i]==0){
            auto p_v = Vs.data()+i*3;
            double d;
            if(d = MyUtility::vecSquareDist(p_v, midpoint.data())< mindist){
                bInd = i;mindist = d;
            }
        }
        for(int j=0;j<3;++j){
            st[j] = Vs[bInd*3+j];
        }
        double maxdist = DBL_MIN;
        maxdist = max(maxdist, MyUtility::_VerticesDistance(leftcorner.data(),Vs.data()+bInd*3));
        maxdist = max(maxdist, MyUtility::_VerticesDistance(rightcorner.data(),Vs.data()+bInd*3));
        iBound = (int) (maxdist / (dSize * 0.75 + 1e-6));
    }else{
        for(int j=0;j<3;++j){
            st[j] = midpoint[j];
        }
        iBound = (int) (width / dSize / 2. * 1.75);

    }


}




int GetOffSurfacePoint(vector<double>&offPts, vector<double>&surPts, vector<uint>&surfv, vector<double>&testPts, double offthres){


    int n_offPts = 0;
    int n_surPts = surPts.size()/3;
    int n_face = surfv.size()/3;
    int n_testPts = testPts.size()/3;

    auto p_v = surPts.data();

    auto p_testv = testPts.data();

    double intesetp[3];

    offPts.clear();
    for(int i=0;i<n_testPts;++i){

        vector<double>dist2tri(n_face);
        bool ison = false;
        for(int j=0;j<n_face;++j){
            auto p_fv = surfv.data() + j*3;
            dist2tri[j] = MyUtility::point2TriDistance(p_testv+i*3, p_v+p_fv[0]*3, p_v+p_fv[1]*3, p_v+p_fv[2]*3, intesetp);
            if(dist2tri[j]<offthres){
                ison = true;break;
            }
        }
        if(!ison){
            n_offPts++;
            for(int j=0;j<3;++j)offPts.push_back(p_testv[i*3+j]);
        }

    }

    return n_offPts;


}





double Surfacer::Surfacing_Implicit(vector<double>&Vs, vector<int>&labels, bool ischeckall,
                                    double (*function)(const R3Pt &in_pt)){

    p_ImplicitSurfacer = this;
    ClearBuffer();

    CalSurfacingPara(Vs, labels);


    //polygonize(function, size, bounds, st, triproc, vertproc);
    double thresDist = 1e-3;

    double re_time;
    cout<<"Implicit Surfacing: "<<endl;

    auto t1 = Clock::now();




    if(!ischeckall){
        polygonize(function, dSize, iBound, st, TriProc, VertProc);
        GetCurSurface(all_v,all_fv);
    }else{

        vector<double>offPts;
        vector<double>surPts;
        vector<uint>surfv;
        vector<double>testPts;
        for(int i=0;i<labels.size();++i)if(labels[i]==0){
            for(int j=0;j<3;++j)testPts.push_back(Vs[i*3+j]);
        }

        int ncomp = 0;
        while(true){
            ClearSingleComponentBuffer();
            if(polygonize(function, dSize, iBound, st, TriProc, VertProc))break;

            GetCurSurface(surPts,surfv);
            InsertToCurSurface(surPts,surfv);
            GetOffSurfacePoint(offPts,surPts,surfv,testPts,5e-2);
            cout<<"ncomp found: "<<++ncomp<<" offpts: "<<offPts.size()/3<<endl;
            if(offPts.size()<=10*3)break;
            break;

            for(int j=0;j<3;++j)st[j] = offPts[j];
            testPts = offPts;
        }


    }

    cout<<"Implicit Surfacing Done."<<endl;
    auto t2 = Clock::now();
    cout << "Total Surfacing time: " <<  (re_time = std::chrono::nanoseconds(t2 - t1).count()/1e9) <<endl;

    return re_time;

}

int cc=0;
FT sphere_function (Point_3 p) {
    ++cc;
    const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
    return x2+y2+z2-1;
}

double Surfacer::Surfacing_CGALImplicit(FT (*function)(const Point_3)){

    double re_time;
    cout<<"CGAL Implicit Mesher"<<endl;
    Tr tr;            // 3D-Delaunay triangulation
    C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation
    // defining the surface
    Surface_3 surface(function,             // pointer to function
                      Sphere_3(CGAL::ORIGIN, 2.0)); // bounding sphere
    // Note that "2." above is the *squared* radius of the bounding sphere!
    // defining meshing criteria
    CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,  // angular bound
                                                       0.03,  // radius bound
                                                       0.03); // distance bound
    // meshing surface
    auto t1 = Clock::now();
    CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
    cout<<"Surfacing time: "<<(re_time = std::chrono::nanoseconds(Clock::now() - t1).count()/1e9)<<endl;
    //cout<<"cc: "<<cc<<endl;
    std::cout << "Final number of points: " << tr.number_of_vertices() << "\n";
    std::ofstream out("/Users/Research/Geometry/RBF/data/out.off");
    CGAL::output_surface_facets_to_off (out, c2t3);

    using CGAL::Surface_mesher::number_of_facets_on_surface;
    typedef typename Tr::Vertex_handle Vertex_handle;
    typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
    typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
    std::unordered_map<Vertex_handle, int> V;
    int inum = 0;
    all_v.clear();
    for(Finite_vertices_iterator vit = tr.finite_vertices_begin();vit != tr.finite_vertices_end();++vit){
        const FT* p_v = &(vit->point().x());
        for(int i=0;i<3;++i)all_v.push_back(p_v[i]);
        V[vit] = inum++;
    }
    all_fv.clear();
    for( Finite_facets_iterator fit = tr.finite_facets_begin();fit != tr.finite_facets_end(); ++fit){
        const typename Tr::Cell_handle cell = fit->first;
        const int& index = fit->second;
        if (cell->is_facet_on_surface(index)==true)
        {
            for(int i=0;i<3;++i)all_fv.push_back(V[cell->vertex(tr.vertex_triple_index(index, i))]);
        }
    }

    return re_time;

}

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Timer.h>
#include <CGAL/trace.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Implicit_reconstruction_function.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Polyhedron_VRML_1_ostream.h>

#include <list>
#include <cstdlib>
#include <fstream>
#include <math.h>

#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include "boost/program_options.hpp"

#include <CGAL/disable_warnings.h>
//double Surfacer::Surfacing_Spectral(vector<double>&v, vector<double>&vn){

//    //    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
//    //    typedef Kernel::Point_3 Point;
//    //    typedef Kernel::Vector_3 Vector;
//    //    typedef std::pair<Point, Vector> Point_with_normal;

//    //    typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
//    //    typedef CGAL::First_of_pair_property_map<Point_with_normal> Point_map;
//    //    typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;
//    //    typedef std::list<Point_with_normal> PointList;
//    //    typedef CGAL::Implicit_reconstruction_function<Kernel, PointList, Normal_map> Implicit_reconstruction_function;

//    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
//    // Simple geometric types
//    typedef Kernel::FT FT;
//    typedef Kernel::Point_3 Point;
//    typedef Kernel::Vector_3 Vector;
//    //typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal;
//    typedef std::pair<Point, Vector> Point_with_normal;
//    typedef Kernel::Sphere_3 Sphere;
//    typedef std::vector<Point_with_normal> PointList;
//    typedef CGAL::First_of_pair_property_map<Point_with_normal> Point_map;
//    typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;

//    // polyhedron
//    typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

//    // Spectral implicit function
//    typedef CGAL::Implicit_reconstruction_function<Kernel, PointList, Normal_map> Implicit_reconstruction_function;

//    //typedef CGAL::Implicit_reconstruction_function<Kernel, PointList, CGAL::Identity_property_map<Point_with_normal> > Implicit_reconstruction_function;

//    // Surface mesher
//    typedef CGAL::Surface_mesh_default_triangulation_3 STr;
//    typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
//    typedef CGAL::Implicit_surface_3<Kernel, Implicit_reconstruction_function> Surface_3;

//    double laplacian = 0.1;
//    double bilaplacian = 1.;
//    double ratio = 10.;
//    double fitting = 0.1;
//    //    if(0){
//    //        string filename("/Users/Research/Geometry/RBF/data/doghead.xyz");
//    //        ofstream outer(filename.data(), ofstream::out);
//    //        if (!outer.good()) {
//    //            cout << "Can not create output file " << filename << endl;
//    //            return -1;
//    //        }

//    //        for(int i=0;i<v.size()/3;++i){
//    //            auto p_v = v.data()+i*3;
//    //            auto p_vn = vn.data()+i*3;
//    //            outer<<p_v[0]<<' '<<p_v[1]<<' '<<p_v[2]<<' '<<p_vn[0]<<' '<<p_vn[1]<<' '<<p_vn[2]<<endl;
//    //        }
//    //        outer.close();
//    //    }


//    //    //std::vector<Point_with_normal> points;//kitten.xyz sphere926.pwn sphere_new.xyz
//    PointList points;
//    if(0){
//        std::ifstream stream("/Users/Research/Geometry/RBF/external/cgal-public/Implicit_surface_reconstruction_3/examples/Implicit_surface_reconstruction_3/data/sphere_new.xyz");
//        if (!stream ||
//                !CGAL::read_xyz_points(
//                    stream,
//                    std::back_inserter(points),
//                    CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Point_with_normal>()).
//                    normal_map(CGAL::Second_of_pair_property_map<Point_with_normal>())))
//        {
//            std::cerr << "Error: cannot read file data/kitten.xyz" << std::endl;
//            return EXIT_FAILURE;
//        }
//        stream.close();
//    }else{
//        int  nv = v.size()/3;
//        for(int i=0;i<nv;++i){
//            auto p_v = v.data()+i*3;
//            auto p_vn = vn.data()+i*3;
//            points.emplace_back(make_pair(Point(p_v[0],p_v[1],p_v[2]), Vector(p_vn[0],p_vn[1],p_vn[2])));
//        }
//    }

//    //    Polyhedron output_mesh;

//    //    double average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>
//    //            (points, 6, CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Point_with_normal>()));

//    //    double re_time;
//    //    auto t1 = Clock::now();
//    //    CGAL::spectral_surface_reconstruction_delaunay(points,
//    //                                                   CGAL::First_of_pair_property_map<Point_with_normal>(),
//    //                                                   CGAL::Second_of_pair_property_map<Point_with_normal>(),
//    //                                                   100., 25., output_mesh, 1, 0., average_spacing);
//    //    cout<<"Surfacing time: "<<(re_time = std::chrono::nanoseconds(Clock::now() - t1).count()/1e9)<<endl;

//    Implicit_reconstruction_function function(points, Point_map(), Normal_map());
//    //compute_poisson_implicit_function() compute_spectral_implicit_function(fitting, ratio, bilaplacian, laplacian)
//    if (! function.compute_poisson_implicit_function() )
//    {
//        std::cerr << "Error: cannot compute implicit function" << std::endl;
//        //accumulated_fatal_err = EXIT_FAILURE;
//        //continue;
//    }

//    string curr_outfile("/Users/Research/Geometry/RBF/external/spectra-0.6.2/out.off");
//    function.marching_tetrahedra(0, curr_outfile);

//    readOffFile(curr_outfile,all_v,all_fv);

//    //    //std::cout << "Final number of points: " << output_mesh.number_of_vertices() << "\n";
//    //    std::ofstream out("/Users/Research/Geometry/RBF/data/sphere_spectral-20-30-0.375.off");
//    //    out << output_mesh;
//    //    std::unordered_map<Polyhedron::Vertex_const_iterator, int> V;
//    //    int inum = 0;

//    //    typedef typename Polyhedron::Vertex_const_iterator                  VCI;
//    //    typedef typename Polyhedron::Facet_const_iterator                   FCI;
//    //    typedef typename Polyhedron::Halfedge_around_facet_const_circulator HFCC;

//    //    all_v.clear();
//    //    for ( Polyhedron::Vertex_const_iterator vit = output_mesh.vertices_begin(); vit != output_mesh.vertices_end(); ++vit){
//    //        const FT* p_v = &(vit->point().x());
//    //        for(int i=0;i<3;++i)all_v.push_back(p_v[i]);
//    //        V[vit] = inum++;
//    //    }
//    //    all_fv.clear();
//    //    for( FCI fi = output_mesh.facets_begin(); fi != output_mesh.facets_end(); ++fi) {
//    //        HFCC hc = fi->facet_begin();
//    //        //std::size_t n = circulator_size( hc);
//    //        CGAL_assertion( circulator_size( hc)== 3);
//    //        for(int i=0; i<3;++i){
//    //            all_fv.push_back(V[ VCI(hc->vertex())]);
//    //            ++hc;
//    //        }

//    //    }





//    //    return EXIT_SUCCESS;




//}

double Surfacer::Surfacing_Spectral(vector<double>&v, vector<double>&vn){


    cout<<"Surfacing_Spectral"<<endl;
    string exepath("/Users/Research/Geometry/RBF/external/cgal-public-dev-Implicit_surface_reconstruction_3-octree_based_refinement_update/Implicit_surface_reconstruction_3/test/Implicit_surface_reconstruction_3/");
    string exeinputfile("input.pwn");
    string scriptname("./run.sh");

    string filename = exepath + exeinputfile;
    writeXYZnormal(filename,v,vn);

    string systemcomd = string("cd ")+exepath+string("; ") +scriptname;

    int result = system(systemcomd.data());


    readOffFile(exepath+"iso_facet_1_out.off", all_v,all_fv);
}




double Surfacer::Surfacing_Poisson(vector<double>&v,vector<double>&vn){


//    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
//    typedef Kernel::Point_3 Point;
//    typedef Kernel::Vector_3 Vector;
//    typedef std::pair<Point, Vector> Pwn;
//    typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

//    if(0){
//        string filename("/Users/Research/Geometry/RBF/data/doghead.xyz");
//        ofstream outer(filename.data(), ofstream::out);
//        if (!outer.good()) {
//            cout << "Can not create output file " << filename << endl;
//            return -1;
//        }

//        for(int i=0;i<v.size()/3;++i){
//            auto p_v = v.data()+i*3;
//            auto p_vn = vn.data()+i*3;
//            outer<<p_v[0]<<' '<<p_v[1]<<' '<<p_v[2]<<' '<<p_vn[0]<<' '<<p_vn[1]<<' '<<p_vn[2]<<endl;
//        }
//        outer.close();
//    }


//    std::vector<Pwn> points;//kitten.xyz sphere926.pwn sphere_new.xyz
//    if(0){
//        std::ifstream stream("/Users/Research/Geometry/RBF/external/cgal-public/Implicit_surface_reconstruction_3/examples/Implicit_surface_reconstruction_3/data/sphere_new.xyz");
//        if (!stream ||
//                !CGAL::read_xyz_points(
//                    stream,
//                    std::back_inserter(points),
//                    CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>()).
//                    normal_map(CGAL::Second_of_pair_property_map<Pwn>())))
//        {
//            std::cerr << "Error: cannot read file data/kitten.xyz" << std::endl;
//            return EXIT_FAILURE;
//        }
//        stream.close();
//    }else{
//        int  nv = v.size()/3;
//        for(int i=0;i<nv;++i){
//            auto p_v = v.data()+i*3;
//            auto p_vn = vn.data()+i*3;
//            points.emplace_back(make_pair(Point(p_v[0],p_v[1],p_v[2]), Vector(p_vn[0],p_vn[1],p_vn[2])));
//        }
//    }

//    Polyhedron output_mesh;

//    double average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>
//            (points, 6, CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>()));

//    double re_time;
//    auto t1 = Clock::now();
//    CGAL::poisson_surface_reconstruction_delaunay
//            (points,
//             CGAL::First_of_pair_property_map<Pwn>(),
//             CGAL::Second_of_pair_property_map<Pwn>(),
//             output_mesh, average_spacing/4, 20, 5, 0.1);
//    cout<<"Surfacing time: "<<(re_time = std::chrono::nanoseconds(Clock::now() - t1).count()/1e9)<<endl;

//    //std::cout << "Final number of points: " << output_mesh.number_of_vertices() << "\n";
//    std::ofstream out("/Users/Research/Geometry/RBF/data/sphere_spectral-20-30-0.375.off");
//    out << output_mesh;
//    std::unordered_map<Polyhedron::Vertex_const_iterator, int> V;
//    int inum = 0;

//    typedef typename Polyhedron::Vertex_const_iterator                  VCI;
//    typedef typename Polyhedron::Facet_const_iterator                   FCI;
//    typedef typename Polyhedron::Halfedge_around_facet_const_circulator HFCC;

//    all_v.clear();
//    for ( Polyhedron::Vertex_const_iterator vit = output_mesh.vertices_begin(); vit != output_mesh.vertices_end(); ++vit){
//        const FT* p_v = &(vit->point().x());
//        for(int i=0;i<3;++i)all_v.push_back(p_v[i]);
//        V[vit] = inum++;
//    }
//    all_fv.clear();
//    for( FCI fi = output_mesh.facets_begin(); fi != output_mesh.facets_end(); ++fi) {
//        HFCC hc = fi->facet_begin();
//        //std::size_t n = circulator_size( hc);
//        CGAL_assertion( circulator_size( hc)== 3);
//        for(int i=0; i<3;++i){
//            all_fv.push_back(V[ VCI(hc->vertex())]);
//            ++hc;
//        }

//    }





//    return EXIT_SUCCESS;

    return -1;

}




double Surfacer::Surfacing_PowerCrust(vector<double>&v){

    cout<<"Surfacing_PowerCrust"<<endl;
    string exepath("/Users/Research/Geometry/RBF/external/powercrust/");
    string exeinputfile("input.pts");
    string scriptname("./run.sh");

    string filename = exepath + exeinputfile;
    writeXYZ(filename,v);

    string systemcomd = string("cd ")+exepath+string("; ") +scriptname;

    int result = system(systemcomd.data());


    readOffFile("/Users/Research/Geometry/RBF/external/powercrust/pc.off", all_v,all_fv);

}

void readply(std::wstring filename, vector<double>& vertices, vector<uint>& triangles)
{
    libply::File file(filename);
    const auto& definitions = file.definitions();

    const auto vertexDefinition = definitions.at(0);
    const size_t vertexCount = vertexDefinition.size;
    vertices.reserve(vertexCount);
    libply::ElementReadCallback vertexCallback = [&vertices](libply::ElementBuffer& e)
    {
        for(int i=0;i<3;++i)vertices.push_back(e[i]);
    };

    const auto triangleDefinition = definitions.at(1);
    const size_t triangleCount = triangleDefinition.size;
    triangles.reserve(triangleCount);
    libply::ElementReadCallback triangleCallback = [&triangles](libply::ElementBuffer& e)
    {
        for(int i=0;i<3;++i)triangles.push_back(e[i]);
    };

    file.setElementReadCallback("vertex", vertexCallback);
    file.setElementReadCallback("face", triangleCallback);
    file.read();
}

int readply_vt(std::string filename, vector<double>& v, vector<uint>& t){

    yaply::PlyFile plyFile(filename.c_str());

    std::cout << "loaded file " << filename << ":" << std::endl;
    for (const auto& el : plyFile.elements_) {
        std::cout << "element " << el.name << " " << el.nrElements << std::endl;
        for (const auto& p : el.properties)
            std::cout << "property " << p->name << std::endl;
    }

    auto& vertices = plyFile["vertex"];
    if (vertices.nrElements <= 0) {
        std::cerr << "could not find vertex attribute in the plyfile" << std::endl;
        return EXIT_FAILURE;
    }
    yaply::PLY_ELEMENT& faces = plyFile["face"];
    if (faces.nrElements <= 0) {
        std::cerr << "could not find face attribute in the plyfile" << std::endl;
        return EXIT_FAILURE;
    }
    v.clear();
    std::vector<double> x, y, z;
    if (vertices.getScalarProperty("x", x) && vertices.getScalarProperty("y", y)
            && vertices.getScalarProperty("z", z)) {
        for(int i=0;i<x.size();++i){
            v.push_back(x[i]);v.push_back(y[i]);v.push_back(z[i]);
        }
    }else {
        cout<<"no point xyx!"<<endl;
    }
    t.clear();
    vector<vector<int>> ts;
    {
        //ts = faces.getProperty("vertex_indices");
        int nf;
        faces.getListProperty("vertex_indices",ts,nf);
        for(auto &a:ts)for(auto b:a)t.push_back(b);

    }
}

double Surfacer::Surfacing_ScreenedPoisson(vector<double>&v,vector<double>&vn){

    cout<<"Surfacing_ScreenedPoisson"<<endl;
    string exepath("/Users/Research/Geometry/RBF/external/PoissonRecon/");
    string exeinputfile("input.npts");
    string scriptname("./run.sh");

    string filename = exepath + exeinputfile;
    //cout<<"Surfacing_ScreenedPoisson"<<endl;
    writeXYZnormal(filename,v,vn);

    //cout<<"Surfacing_ScreenedPoisson"<<endl;
    string systemcomd = string("cd ")+exepath+string("; ") +scriptname;

    int result = system(systemcomd.data());

    readply(L"/Users/Research/Geometry/RBF/external/PoissonRecon/Result/output.ply",all_v,all_fv);

    //readOffFile("/Users/Research/Geometry/RBF/external/nasr-1.0/out.off", all_v,all_fv);

}

double Surfacer::Surfacing_NASR(vector<double>&v){

    cout<<"Surfacing_NASR"<<endl;
//    string exepath("/Users/Research/Geometry/RBF/external/nasr-1.0/");
//    string exeinputfile("input.xyz");
//    string scriptname("./run.sh");

//    string filename = exepath + exeinputfile;
//    writeXYZ(filename,v);

//    string systemcomd = string("cd ")+exepath+string("; ") +scriptname;

//    int result = system(systemcomd.data());


//    readOffFile("/Users/Research/Geometry/RBF/external/nasr-1.0/out.off", all_v,all_fv);
    vector<double>vn;
    readObjFile("/Users/Research/Geometry/RBF/data/container/sketch/wang_walsus/surf.obj",all_v,all_fv,vn);


}


void Surfacer::WriteSurface(string fname){

    writeObjFile(fname,all_v,all_fv);

}

void Surfacer::WriteSurface(vector<double> &v, vector<uint>&fv){

    v = all_v;
    fv = all_fv;
}

void Surfacer::WriteSurface(vector<double> **v, vector<uint> **fv){

    *v = &all_v;
    *fv = &all_fv;
}

void Surfacer::ClearBuffer(){

    ClearSingleComponentBuffer();
    all_v.clear();
    all_fv.clear();
}

void Surfacer::ClearSingleComponentBuffer(){
    s_aptSurface.clearcompletely();
    s_avecSurface.clearcompletely();
    s_afaceSurface.clearcompletely();
}

void Surfacer::GetCurSurface(vector<double>&v,vector<uint>&fv){

    int beInd = v.size()/3;
    for(int i=0;i<s_aptSurface.num();++i){
        for(int j=0;j<3;++j)v.push_back( s_aptSurface[i][j] );
    }

    for(int i=0;i<s_afaceSurface.num();++i){
        for(int j=0;j<3;++j)fv.push_back( beInd + s_afaceSurface[i][2-j]);
    }

}


void Surfacer::InsertToCurSurface(vector<double>&v,vector<uint>&fv){

    int beInd = all_v.size()/3;
    all_v.insert(all_v.end(),v.begin(),v.end());

    for(auto a:fv)all_fv.push_back(beInd+a);

}
