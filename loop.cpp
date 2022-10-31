#include "adjacency.h"
#include "element.h"
#include "grid.h"
#include "halfedges.h"
#include "loop.h"
#include "mesh.h"
#include "vec.h"
#include "vertices.h"


namespace flux{

int
getOppositePoint( int p , int q , int t, std::vector<int> triangles ) {
    int k1 = triangles[3*t];
    int k2 = triangles[3*t+1];
    int k3 = triangles[3*t+2];
    if ((k1 == p && k2 == q) || (k1 == q && k2 == p)){
        return k3;
    }else if ((k2 == p && k3 == q) || (k2 == q && k3 == p)){
        return k1;
    }else if ((k3 == p && k1 == q) || (k3 == q && k1 == p)){
        return k2;
    }else{
        return -1;
    }
  
}

int** 
edge_to_triangle( const Mesh<Triangle>& mesh ){
    std::vector<int> edges;
    mesh.get_edges(edges);

    std::vector<int> triangles;
    mesh.get_triangles(triangles);
    
    // this data structure stores edge indeces and their two associated triangle indeces 
    // (i.e. edge 7 is adjacent to triangle 2 and triangle 88)

    int** e2t = 0;
    e2t = new int*[edges.size()/2];
    for(int i = 0; i < edges.size()/2; i++){
        e2t[i] = new int[2];
        int temp = 0;
        int p = edges[i*2];
        int q = edges[i*2+1];
        for(int j = 0; j < triangles.size()/3; j++){
            for(int k = 0; k < 3; k ++){
                if(triangles[3*j+k]== p){
                    for(int l = 0; l < 3; l++){
                        if(triangles[3*j+l]== q){
                            e2t[i][temp] = j;
                            temp++;
                        }
                    }
                }
            }
            if(temp == 2){
                break;
            }
        }
        if(temp == 1){
            e2t[i][1] = -1;
        }
    }
    return e2t;
}

int** 
triangle_to_edge( const Mesh<Triangle>& mesh ){
    std::vector<int> edges;
    mesh.get_edges(edges);

    std::vector<int> triangles;
    mesh.get_triangles(triangles);

    // this data structure stores triangle indeces and their three associated edge indeces 
    // (i.e. triangle 4 has edges indexed at 7, 30, and 99)
    int** t2e = 0;
    t2e = new int*[triangles.size()/3];
    for(int i = 0; i < triangles.size()/3; i++) {
        t2e[i] = new int[3];
        int t1 = triangles[i*3];
        int t2 = triangles[i*3+1];
        int t3 = triangles[i*3+2];
        for(int j = 0; j < edges.size()/2; j ++){
            if(edges[2*j] == t1 && edges[2*j+1] == t2){
                t2e[i][0] = j;
                break;
            }
            else if(edges[2*j] == t2 && edges[2*j+1] == t1){
                t2e[i][0] = j;
                break;
            }
        }
        for(int j = 0; j < edges.size()/2; j ++){
            if(edges[2*j] == t2 && edges[2*j+1] == t3){
                t2e[i][1] = j;
                break;
            }
            else if(edges[2*j] == t3 && edges[2*j+1] == t2){
                t2e[i][1] = j;
                break;
            }
        }
        for(int j = 0; j < edges.size()/2; j ++){
            if(edges[2*j] == t1 && edges[2*j+1] == t3){
                t2e[i][2] = j;
                break;
            }
            else if(edges[2*j] == t3 && edges[2*j+1] == t1){
                t2e[i][2] = j;
                break;
            }
        }
    }
    return t2e;

}

void
loop_subdivision( Mesh<Triangle>& mesh, int levels, int dim ){
    if(dim == 2){
        for(int run = 0; run < levels; run++){
            std::vector<int> edges2d;
            mesh.get_edges(edges2d);

            std::vector<int> triangles2d;
            mesh.get_triangles(triangles2d);

            int** e2t2d = edge_to_triangle(mesh);
            int** t2e2d = triangle_to_edge(mesh);

            
            Vertices vertices2d(2);
            for(int v = 0; v < mesh.vertices().nb(); v++){
                vertices2d.add(mesh.vertices()[v]);
            }

            Vertices replace2d(2);
            for(int v = 0; v < mesh.vertices().nb(); v++){
                double x[2];
                x[0] = mesh.vertices()[v][0];
                x[1] = mesh.vertices()[v][1];

                std::vector<int> points;
                std::vector<int> b_points;
                int flag = 0;
                for(int e = 0; e < edges2d.size()/2; e++){
                    if(edges2d[2*e] == v || edges2d[2*e+1] == v){
                        if( e2t2d[e][0] == -1 || e2t2d[e][1] == -1){
                            if(edges2d[2*e] == v){
                                b_points.push_back(edges2d[2*e+1]);
                            }else{
                                b_points.push_back(edges2d[2*e]);
                            }
                            flag++;
                        }

                        if(edges2d[2*e] == v){
                            points.push_back(edges2d[2*e+1]);
                        }else{
                            points.push_back(edges2d[2*e]);
                        }
                    }
                }
                
                double k2d = 0.0;
                for(int t = 0; t < triangles2d.size()/3; t++){
                    if(triangles2d[3*t] == v || triangles2d[3*t+1] == v || triangles2d[3*t+2] == v){
                        k2d++;
                    }
                }
                
                double beta2d;
                if(k2d==3){
                    beta2d = 3.0/16.0;
                }else{
                    beta2d = 3.0/(8.0*k2d);
                }
                
                double y[2];
                if(flag > 0){
                    double p1 = mesh.vertices()[b_points[0]][0];
                    double p2 = mesh.vertices()[b_points[0]][1];
                    double q1 = mesh.vertices()[b_points[1]][0];
                    double q2 = mesh.vertices()[b_points[1]][1];
                    
                    y[0] = 0.75*x[0] + 0.125*p1 + 0.125*q1;
                    y[1] = 0.75*x[1] + 0.125*p2 + 0.125*q2;
                }
                else{
                    y[0] = (1.0-k2d*beta2d)*x[0];
                    y[1] = (1.0-k2d*beta2d)*x[1];
                    for(int ki = 0; ki < points.size(); ki++){
                        y[0] = y[0] + beta2d*mesh.vertices()[points[ki]][0];
                        y[1] = y[1] + beta2d*mesh.vertices()[points[ki]][1];
                    }
                }
                replace2d.add(y);
            }

            int nb_points2d = mesh.vertices().nb();
            for( int i = 0; i < edges2d.size()/2; i ++){
                double x[2]; 
                int p = edges2d[i*2];
                int q = edges2d[i*2+1];
                
                if( e2t2d[i][0] == -1 || e2t2d[i][1] == -1){
                    x[0] = 0.5*mesh.vertices()[p][0] + 0.5*mesh.vertices()[q][0];
                    x[1] = 0.5*mesh.vertices()[p][1] + 0.5*mesh.vertices()[q][1];
                }else{
                    int v0 = getOppositePoint(p,q,e2t2d[i][0], triangles2d);
                    int v1 = getOppositePoint(p,q,e2t2d[i][1], triangles2d);
                    for (int d = 0; d < 2; d++){
                        x[d] = 0.375*mesh.vertices()[p][d] + 0.375*mesh.vertices()[q][d] +
                        0.125*mesh.vertices()[v0][d] + 0.125*mesh.vertices()[v1][d];
                    }
                }
                vertices2d.add(x);
                replace2d.add(x);
            }
            
            int tri2d[4*triangles2d.size()/3][3];

            for(int t = 0; t < triangles2d.size()/3; t++){
                
                int p0 = triangles2d[3*t];
                int p1 = triangles2d[3*t+1];
                int p2 = triangles2d[3*t+2];
                int q0 = nb_points2d + t2e2d[t][0];
                int q1 = nb_points2d + t2e2d[t][1];
                int q2 = nb_points2d + t2e2d[t][2];

                int t1[3] = {p0,q0,q2};
                int t2[3] = {q0,p1,q1};
                int t3[3] = {q0,q1,q2};
                int t4[3] = {q1,p2,q2};
                
                tri2d[4*t][0] = t1[0];
                tri2d[4*t][1] = t1[1];
                tri2d[4*t][2] = t1[2];
                tri2d[4*t+1][0] = t2[0];
                tri2d[4*t+1][1] = t2[1];
                tri2d[4*t+1][2] = t2[2];
                tri2d[4*t+2][0] = t3[0];
                tri2d[4*t+2][1] = t3[1];
                tri2d[4*t+2][2] = t3[2];
                tri2d[4*t+3][0] = t4[0];
                tri2d[4*t+3][1] = t4[1];
                tri2d[4*t+3][2] = t4[2];
                
            }

            mesh.vertices().clear();
            for(int i = mesh.nb()-1; i >= 0; i--){
                mesh.remove(i);
            }
            
            for(int i = 0; i < replace2d.nb(); i++){
                mesh.vertices().add( replace2d[i] );
            }
            for(int t = 0; t < 4*triangles2d.size()/3; t++){
                mesh.add(tri2d[t]);
            }
        }   
    }else{
        for(int run = 0; run < levels; run++){
    
            std::vector<int> edges;
            mesh.get_edges(edges);

            std::vector<int> triangles;
            mesh.get_triangles(triangles);

            int** e2t = edge_to_triangle(mesh);
            int** t2e = triangle_to_edge(mesh);

            
            Vertices vertices(3);
            for(int v = 0; v < mesh.vertices().nb(); v++){
                vertices.add(mesh.vertices()[v]);
            }

            Vertices replace(3);
            for(int v = 0; v < mesh.vertices().nb(); v++){
                double x[3];
                x[0] = mesh.vertices()[v][0];
                x[1] = mesh.vertices()[v][1];
                x[2] = mesh.vertices()[v][2];

                std::vector<int> points;
                std::vector<int> b_points;
                int flag = 0;
                for(int e = 0; e < edges.size()/2; e++){
                    if(edges[2*e] == v || edges[2*e+1] == v){
                        if( e2t[e][0] == -1 || e2t[e][1] == -1){
                            if(edges[2*e] == v){
                                b_points.push_back(edges[2*e+1]);
                            }else{
                                b_points.push_back(edges[2*e]);
                            }
                            flag++;
                        }

                        if(edges[2*e] == v){
                            points.push_back(edges[2*e+1]);
                        }else{
                            points.push_back(edges[2*e]);
                        }
                    }
                }
                
                double k = 0;
                for(int t = 0; t < triangles.size()/3; t++){
                    if(triangles[3*t] == v || triangles[3*t+1] == v || triangles[3*t+2] == v){
                        k = k + 1.0;
                    }
                }
                
                double beta;
                if(k==3){
                    beta = 3.0/16.0;
                }else{
                    beta = 3.0/(8.0*k);
                }
                
                double y[3];
                if(flag > 0){
                    double p1 = mesh.vertices()[b_points[0]][0];
                    double p2 = mesh.vertices()[b_points[0]][1];
                    double p3 = mesh.vertices()[b_points[0]][2];
                    double q1 = mesh.vertices()[b_points[1]][0];
                    double q2 = mesh.vertices()[b_points[1]][1];
                    double q3 = mesh.vertices()[b_points[1]][2];
                    
                    y[0] = 0.75*x[0] + 0.125*p1 + 0.125*q1;
                    y[1] = 0.75*x[1] + 0.125*p2 + 0.125*q2;
                    y[2] = 0.75*x[2] + 0.125*p3 + 0.125*q3;
                }
                else{
                    y[0] = (1.0-k*beta)*x[0];
                    y[1] = (1.0-k*beta)*x[1];
                    y[2] = (1.0-k*beta)*x[2];
                    for(int ki = 0; ki < points.size(); ki++){
                        y[0] = y[0] + beta*mesh.vertices()[points[ki]][0];
                        y[1] = y[1] + beta*mesh.vertices()[points[ki]][1];
                        y[2] = y[2] + beta*mesh.vertices()[points[ki]][2];
                    }
                }
                replace.add(y);
            }

            int nb_points = mesh.vertices().nb();
            for( int i = 0; i < edges.size()/2; i ++){
                double x[3]; 
                int p = edges[i*2];
                int q = edges[i*2+1];
                
                if( e2t[i][0] == -1 || e2t[i][1] == -1){
                    x[0] = 0.5*mesh.vertices()[p][0] + 0.5*mesh.vertices()[q][0];
                    x[1] = 0.5*mesh.vertices()[p][1] + 0.5*mesh.vertices()[q][1];
                    x[2] = 0.5*mesh.vertices()[p][2] + 0.5*mesh.vertices()[q][2];
                }else{
                    int v0 = getOppositePoint(p,q,e2t[i][0], triangles);
                    int v1 = getOppositePoint(p,q,e2t[i][1], triangles);
                    for (int d = 0; d < 3; d++){
                        x[d] = 0.375*mesh.vertices()[p][d] + 0.375*mesh.vertices()[q][d] +
                        0.125*mesh.vertices()[v0][d] + 0.125*mesh.vertices()[v1][d];
                    }
                }
                vertices.add(x);
                replace.add(x);
            }
            
            int tri[4*triangles.size()/3][3];

            for(int t = 0; t < triangles.size()/3; t++){
                
                int p0 = triangles[3*t];
                int p1 = triangles[3*t+1];
                int p2 = triangles[3*t+2];
                int q0 = nb_points + t2e[t][0];
                int q1 = nb_points + t2e[t][1];
                int q2 = nb_points + t2e[t][2];

                int t1[3] = {p0,q0,q2};
                int t2[3] = {q0,p1,q1};
                int t3[3] = {q0,q1,q2};
                int t4[3] = {q1,p2,q2};
                
                tri[4*t][0] = t1[0];
                tri[4*t][1] = t1[1];
                tri[4*t][2] = t1[2];
                tri[4*t+1][0] = t2[0];
                tri[4*t+1][1] = t2[1];
                tri[4*t+1][2] = t2[2];
                tri[4*t+2][0] = t3[0];
                tri[4*t+2][1] = t3[1];
                tri[4*t+2][2] = t3[2];
                tri[4*t+3][0] = t4[0];
                tri[4*t+3][1] = t4[1];
                tri[4*t+3][2] = t4[2];
            }

            mesh.vertices().clear();
            for(int i = mesh.nb()-1; i >= 0; i--){
                mesh.remove(i);
            }
            
            for(int i = 0; i < replace.nb(); i++){
                mesh.vertices().add( replace[i] );
            }
            for(int t = 0; t < 4*triangles.size()/3; t++){
                mesh.add(tri[t]);
            }
        }
    }
}



} //flux









