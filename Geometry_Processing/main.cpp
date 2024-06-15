#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <cmath>
#include <limits>

//#include "lbfgs.h"
#include "vector-polygon-svg.h"




Vector intersect(const Vector& A, const Vector& B, const Vector& edgeStart, const Vector& edgeEnd) {
    Vector normal = (Vector(edgeEnd[1] - edgeStart[1], edgeStart[0] - edgeEnd[0])).normalized();
    double t = dot(edgeStart - A, normal) / dot(B - A, normal);
    if (0 <= t && t <= 1) return A + t*(B - A);
    return Vector();
}

bool isInside(const Vector& P, const Vector& edgeStart, const Vector& edgeEnd) {
    Vector u = edgeStart;
    Vector v = edgeEnd;
    Vector N = Vector(v[1] - u[1], u[0] - v[0]);
    return dot(P - u, N) <= 0;
}


Polygon sutherlandHodgmanClip(const Polygon& subjectPolygon, const Polygon& clipPolygon) {
    Polygon outPolygon = subjectPolygon;

    for (size_t j = 0; j < clipPolygon.vertices.size(); ++j) {
        Polygon inputPolygon = outPolygon;
        outPolygon.vertices.clear();

        Vector edgeStart = clipPolygon.vertices[j];
        Vector edgeEnd = clipPolygon.vertices[(j + 1) % clipPolygon.vertices.size()];

        for (size_t i = 0; i < inputPolygon.vertices.size(); ++i) {
            Vector curVertex = inputPolygon.vertices[i];
            Vector prevVertex = inputPolygon.vertices[(i > 0) ? (i - 1) : (inputPolygon.vertices.size() - 1)];

            Vector intersection = intersect(prevVertex, curVertex, edgeStart, edgeEnd);

            if (isInside(curVertex, edgeStart, edgeEnd)) {
                if (!isInside(prevVertex, edgeStart, edgeEnd)) {
                    outPolygon.vertices.push_back(intersection);
                }
                outPolygon.vertices.push_back(curVertex);
            } else if (isInside(prevVertex, edgeStart, edgeEnd)) {
                outPolygon.vertices.push_back(intersection);
            }
        }
    }

    return outPolygon;
}


// Function to find intersection point of a line segment with a bisector
Vector computeIntersection(
    const Vector &A,
    const Vector &B, 
    const Vector &Pi, 
    const Vector &Pj, 
    const std::vector<double> &weights) {

    Vector avg((Pi[0] + Pj[0]) / 2, (Pi[1] + Pj[1]) / 2, 0.0);
    Vector weighted = avg + (weights[0] - weights[1])/(2*std::pow((Pi - Pj).norm(), 2))*(Pj - Pi);

    double x = dot(weighted - A, Pi - Pj)/dot(B-A, Pi - Pj);
    //double t = ((M[0] - A[0]) * D[0] + (M[1] - A[1]) * D[1]) / ((B[0] - A[0]) * D[0] + (B[1] - A[1]) * D[1]);
    
    if ((x >= 0) && (x <= 1)) {
        Vector res (A[0] + x * (B[0] - A[0]), A[1] + x * (B[1] - A[1]), 0.0);
        return res;
    }
    else {
    return Vector();
    }
}


bool is_inside(const Vector &A,
        const Vector& Pi,
        const Vector& Pj,
        const std::vector<double> &weights) {

    Vector avg((Pi[0] + Pj[0]) / 2, (Pi[1] + Pj[1]) / 2, 0.0);
    Vector weighted = avg + (weights[0] - weights[1])/(2*std::pow((Pi - Pj).norm(), 2))*(Pj - Pi);

    if (dot(A - weighted, Pj-Pi) < 0) {
        return true;
    }
    else {return false;}
}



std::vector<Polygon> Voronoi(
    const std::vector<Vector> &sites, 
    const Polygon &edges,
    const std::vector<double> &weights = {}
) {
    std::vector<Polygon> output(sites.size());
    #pragma omp parallel for
    for(int i = 0; i < sites.size(); i++) {
        Vector Vi = sites[i];
        Polygon C = edges;
        for(int j = 0; j < sites.size(); j++) {
            if (i == j) {
            continue; // skip that iteration
            }
            Vector Vj = sites[j];
            std::vector<double> w(2);
            w[0] = 1.0;
            w[1] = 1.0;
            std::vector<Vector> border(2);
            border[0] = Vi;
            border[1] = Vj;
            if (weights.size()) {
                w = {weights[i], weights[j]};
            } 
            Polygon out;
            Vector past;
            for(int k = 0; k < C.vertices.size(); k++) {
                
                Vector current = C.vertices[k];
                if (k > 0) {
                    past = C.vertices[k-1];
                }
                else {
                    past = C.vertices[C.vertices.size()-1];
                }
                
                Vector inter_spot = computeIntersection(past, current, border[0], border[1], w);
                if (is_inside(current, border[0], border[1], w)) {
                    if (is_inside(past, border[0], border[1], w) == false) {
                        out.add(inter_spot);
                    }
                    out.add(current);
                } else if (is_inside(past, border[0], border[1], w)) {
                    out.add(inter_spot);
                }
            }
            C = out;
        }
        output[i] = C;
    }
    return output;
}

int main() {
    Polygon firstPolygon({
    Vector(0.3, 0.2), Vector(0.5, 0.3), Vector(0.6, 0.5),
    Vector(0.55, 0.7), Vector(0.4, 0.75), Vector(0.25, 0.6),
    Vector(0.2, 0.4)
});

    Polygon clippingPolygon({
    Vector(0.25, 0.25), Vector(0.75, 0.25),
    Vector(0.75, 0.75), Vector(0.25, 0.75)
    });
    save_svg({firstPolygon, clippingPolygon}, "before.svg");
    Polygon resultPolygon = sutherlandHodgmanClip(firstPolygon, clippingPolygon);
    save_svg({resultPolygon}, "after.svg");

   
    Polygon borders({
        Vector(0, 0), Vector(0, 1),
        Vector(1, 1), Vector(1, 0)
    });

    int number_sites = 1000;


    std::vector<Vector> sites(number_sites);
    srand(time(NULL));
    for (int i = 0; i < number_sites; i++) {
        double a = (double) rand() / RAND_MAX;
        double b = (double) rand() / RAND_MAX;
        sites[i] = Vector(a, b, 0);
    }
    save_svg(Voronoi(sites, borders), "voronoi.svg");



    std::vector<double> weights(number_sites);

    for (int i = 0; i < number_sites; i++) {
        if (sites[i][0] <= 0.1 || sites[i][0] >= 0.9 || sites[i][1] <= 0.1 || sites[i][1] >= 0.9 ) {
            weights[i] = 0.99;
        } else {
            weights[i] = 1;
        }
    }
    std::cout << "saving power diagram..." << std::endl;
    save_svg(Voronoi(sites, borders, weights), "power.svg");

    // std::vector<Vector> sites = {{0, 0}, {100, 0}, {0, 100}, {100, 100}, {50, 50}};
    // std::vector<Polygon> cells;

    // voronoiParallelEnumeration(sites, cells);

    // for (size_t i = 0; i < cells.size(); ++i) {
    //     std::cout << "Voronoi cell for site (" << sites[i].x << ", " << sites[i].y << "): ";
    //     for (const auto &point : cells[i].vertices) {
    //         std::cout << "(" << point.x << ", " << point.y << ") ";
    //     }
    //     std::cout << std::endl;
    // }

    return 0;
}