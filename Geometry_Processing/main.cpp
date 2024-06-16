#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <cmath>
#include <limits>

#include "lbfgs.h"
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




// LBGFS Part -----------------------------------------------------------



objective_function::objective_function(
        const std::vector<double> &lambdas,
        std::function<std::vector<Polygon>(
            const std::vector<Vector>&,
            const Polygon&,
            const std::vector<double>&
        )> Voronoi,
        const std::vector<Vector> &sites,
        const Polygon &bounds
    ) {
    this->m_x = NULL;
    this->lambdas = lambdas;
    this->voronoi_function = Voronoi;
    this->sites = sites;
    this->bounds = bounds;
    this->iterations = 0;
}

objective_function::~objective_function() {
    if (m_x != NULL) {
        lbfgs_free(m_x);
        m_x = NULL;
    }
}

std::vector<Polygon> objective_function::run(int N)
{
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *m_x = lbfgs_malloc(N);

    if (m_x == NULL) {
        printf("ERROR: Failed to allocate a memory block for variables.\n");
        exit(1);
    }
    for (int i = 0; i < N; i ++) {
        m_x[i] = -1;
    }

    int output = lbfgs(N, m_x, &fx, _evaluate, _progress, this, NULL);

    printf("L-BFGS optimization terminated with status code = %d\n", output);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, m_x[0], m_x[1]);

    return this->polygons;
}

lbfgsfloatval_t objective_function::_evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
    ) {
    return reinterpret_cast<objective_function*>(instance)->evaluate(x, g, n, step);
}

lbfgsfloatval_t objective_function::evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
    ) {
    lbfgsfloatval_t fx = 0.0;
    std::vector<double> weights(x, x + n);
    this->polygons = this->Voronoi(sites, bounds, weights);
    if (iterations % 100 == 0) {
        std::string filename = "optimized_" + std::to_string(iterations) + ".svg";
        save_svg(polygons, filename);
    }
    for (uint i = 0; i < n; i++) {
        std::vector<Vector> vertices = this->polygons[i].vertices;
        size_t n = vertices.size();
        Vector point = this->sites[i];
        double area = this->polygons[i].area();
        double tmp = 0;
        if (n > 0) {
            Vector c1 = vertices[0];
            for (uint i = 0; i < n - 2; i++) {
                Vector c2 = vertices[i + 1];
                Vector c3 = vertices[i + 2];
                double T = Polygon({c1,c2,c3}).area();
                tmp += (T/6.) * (
                    dot(c1 - point, c1 - point) +
                    dot(c1 - point, c2 - point) +
                    dot(c1 - point, c3 - point) +
                    dot(c2 - point, c2 - point) +
                    dot(c2 - point, c3 - point) +
                    dot(c3 - point, c3 - point)
                );
            }
        }
        fx += tmp - x[i]*area + this->lambdas[i]*x[i];
        g[i] = area - this->lambdas[i];
    }
    iterations += 1;
    return -fx;
}

int objective_function::_progress(
        void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
    ) {
    return reinterpret_cast<objective_function*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
}

int objective_function::progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
    ) {
    printf("Iteration %d:\n", k);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    return 0;
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