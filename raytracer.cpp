#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <cmath>
#include <limits>
#include <random>
#include <cfloat>
#include <tuple>
#include <omp.h>
#include <list>
#include <chrono>
#include <iostream>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
static std::default_random_engine engine(10);
static std::uniform_real_distribution<double> uniform(0.0, 1.0);





class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    double norm2() const {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }
    double norm() const {
        return sqrt(norm2());
    }
    Vector normalize() {
        double n = norm();
        Vector temp;
        temp[0] = data[0]/n;
        temp[1] = data[1]/n;
        temp[2] = data[2]/n;
        return temp;
    }
    double operator[](int i) const { return data[i]; };
    double& operator[](int i) { return data[i]; };
    double data[3];
};

double dot(const Vector& a, const Vector& b) {
    return a[0]*b[0] + a[1]*b[1]+ a[2]*b[2];
}

Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]);
}

Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0]+b[0], a[1]+b[1], a[2]+b[2]);
}

Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0]-b[0], a[1]-b[1], a[2]-b[2]);
}

Vector operator*(const double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}

Vector operator*(const Vector& a, const double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}

Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}

Vector operator/(const Vector& a, const double b) {
    return Vector(a[0]/b, a[1]/b, a[2]/b);
}

Vector operator/(const double a, const Vector& b) {
    return Vector(b[0]/a, b[1]/a, b[2]/a);
}

Vector operator/(const Vector& a, const Vector& b) {
    return Vector(a[0]/b[0], a[1]/b[1], a[2]/b[2]);
}



class Ray {
    public:
        Vector O;
        Vector u;
        explicit Ray(Vector Origin, Vector dir) {
            O = Origin;
            u = dir;
        }
};


struct Intersection { //defines beahviors of ray-object intersections (to get colors of the pixels at the end)
    Vector P;
    Vector N;
    Vector albedo;
    double t = 0.;
    double refractive_index;
    bool reflective;
    bool intersects = false;

    Intersection(Vector albedo = Vector(0., 0., 0.), double refractive_index = 1., bool reflective = false){
        this->albedo = albedo;
        this->refractive_index = refractive_index;
        this->reflective = reflective;
    }
};




class Geometry {
public:
    virtual Intersection intersect(const Ray& ray) = 0;
};





class BoundingBox {
public:
    Vector B_min;
    Vector B_max;

    explicit BoundingBox(Vector min = Vector(), Vector max = Vector()) {
        B_min = min;
        B_max = max;
    }
};


class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
    };
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;  // indices within the uv coordinates array
    int ni, nj, nk;  // indices within the normals array
    int group;       // face group
};

struct Node {
    BoundingBox bbox;
    int starting_triangle;
    int ending_triangle;
    Node* child_left;
    Node* child_right;
};

class TriangleMesh : public Geometry {
    double scaling_factor;
    Vector translation;
    Vector albedo;
    double refractive_index;
    bool reflective;

public:

    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    BoundingBox bbox;
    Node* root;

    ~TriangleMesh() {}
    TriangleMesh(double scaling_factor,
                 Vector translation,
                 Vector albedo = Vector(0., 0., 0.),
                 double refractive_index = 1.,
                 bool reflective = false) {
        this->scaling_factor = scaling_factor;
        this->translation = translation;
        this->albedo = albedo;
        this->refractive_index = refractive_index;
        this->reflective = reflective;
        this->root = new Node;
    };



// Obj file reader

    void readOBJ(const char* obj) {

        char matfile[255];
        char grp[255];

        FILE* f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f)) {
            char line[255];
            if (!fgets(line, 255, f)) break;

            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());

            if (line[0] == 'u' && line[1] == 's') {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }

            if (line[0] == 'v' && line[1] == ' ') {
                Vector vec;

                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));

                    vertices.push_back(vec);
                    vertexcolors.push_back(col);

                } else {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n') {
                Vector vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't') {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f') {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;

                char* consumedline = line + 1;
                int offset;

                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                    if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                    if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                    if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                    indices.push_back(t);
                } else {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                        if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                        if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3) {
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            indices.push_back(t);
                        } else {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                            if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                            if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }

                consumedline = consumedline + offset;

                while (true) {
                    if (consumedline[0] == '\n') break;
                    if (consumedline[0] == '\0') break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                        if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                        if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                        if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    } else {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                            if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                            if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        } else {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                                if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                                if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            } else {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1) {
                                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                } else {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }

            }

        }
        fclose(f);
        this->BVH(root, 0, indices.size());
    }

    /*  ------------------------ Bounding Box Auxiliary ------------------------   */

    BoundingBox bounding_box(int starting_triangle, int ending_triangle) {
        double minX{DBL_MAX}, minY{DBL_MAX}, minZ{DBL_MAX};
        double maxX{DBL_MIN}, maxY{DBL_MIN}, maxZ{DBL_MIN};
        for (int i=starting_triangle; i<ending_triangle; i++) {
            auto triangle_vertices = {
                this->vertices[this->indices[i].vtxi],
                this->vertices[this->indices[i].vtxj],
                this->vertices[this->indices[i].vtxk]
            };
            for (auto const& v : triangle_vertices) {
                Vector Vertex = scaling_factor * v + translation;
                minX = std::min(minX, Vertex[0]);
                maxX = std::max(maxX, Vertex[0]);
                minY = std::min(minY, Vertex[1]);
                maxY = std::max(maxY, Vertex[1]);
                minZ = std::min(minZ, Vertex[2]);
                maxZ = std::max(maxZ, Vertex[2]);
            }
        }

        return BoundingBox(Vector(minX, minY, minZ), Vector(maxX, maxY, maxZ));
    }

    std::tuple<double, double> get_min(Vector N, Vector Bmin, Vector Bmax, const Ray& ray) {
        double t_B_min = dot(Bmin - ray.O, N) / dot(ray.u, N);
        double t_B_max = dot(Bmax - ray.O, N) / dot(ray.u, N);
        double tx_0 = std::min(t_B_min, t_B_max);
        double tx_1 = std::max(t_B_min, t_B_max);
        return std::make_tuple(tx_0, tx_1);
    }

    bool bounding_box_intersects(const Ray& ray, BoundingBox bounding_box, double* t) {
        // Check for ray-plane intersection
        auto [tx_0, tx_1] = get_min(Vector(1, 0, 0), bounding_box.B_min, bounding_box.B_max, ray);
        auto [ty_0, ty_1] = get_min(Vector(0, 1, 0), bounding_box.B_min, bounding_box.B_max, ray);
        auto [tz_0, tz_1] = get_min(Vector(0, 0, 1), bounding_box.B_min, bounding_box.B_max, ray);

        double min = std::min(std::min(tx_1, ty_1), tz_1);
        double max = std::max(std::max(tx_0, ty_0), tz_0);
        *t = (min > max) ? max : *t;
        return (min > max) ? true : false;
    }

    Vector BaryCenter(int i) {
        Vector vertex1 = (scaling_factor * this->vertices[this->indices[i].vtxi]) + translation;
        Vector vertex2 = (scaling_factor * this->vertices[this->indices[i].vtxj]) + translation;
        Vector vertex3 = (scaling_factor * this->vertices[this->indices[i].vtxk]) + translation;
        return (vertex1 + vertex2 + vertex3)/3.;
    }

    void BVH(Node *node, int starting_triangle, int ending_triangle) {
        node->bbox = this->bounding_box(starting_triangle, ending_triangle);
        node->starting_triangle = starting_triangle;
        node->ending_triangle = ending_triangle;

        Vector diag = node->bbox.B_max - node->bbox.B_min;
        Vector middle_diag = node->bbox.B_min + diag*0.5;
        int longest_axis = (abs(diag[0]) > abs(diag[1]) &&
                            abs(diag[0]) > abs(diag[2])) ? 0 :
                           (abs(diag[1]) > abs(diag[2])) ? 1 : 2;

        int pivot_index = starting_triangle;

        for (int i = starting_triangle; i < ending_triangle; i++) {
            Vector barycenter = this->BaryCenter(i);
            if (barycenter[longest_axis] < middle_diag[longest_axis]) {
                std::swap(indices[i], indices[pivot_index]);
                pivot_index++;
            }
        }

        if (pivot_index <= starting_triangle
            || pivot_index >= ending_triangle - 1
            || ending_triangle-starting_triangle < 5
            )
            return;

        node->child_left = new Node();
        node->child_right = new Node();
        this->BVH(node->child_left, starting_triangle, pivot_index);
        this->BVH(node->child_right, pivot_index, ending_triangle);
    }


Intersection intersect(const Ray &ray) override {

        Intersection I(this->albedo, this->refractive_index, this->reflective);
        double t;
        double t_min{DBL_MAX};
        if (!bounding_box_intersects(ray, this->root->bbox, &t))
            return Intersection();

        std::list<Node*> nodes_to_visit;
        nodes_to_visit.push_front(root);
        while (!nodes_to_visit.empty()) {

            Node *curNode = nodes_to_visit.back();
            nodes_to_visit.pop_back();

            if (curNode->child_left) {
                if (bounding_box_intersects(ray, curNode->child_left->bbox, &t))
                    if (t < t_min)
                        nodes_to_visit.push_back(curNode->child_left);
                if (bounding_box_intersects(ray, curNode->child_right->bbox, &t))
                    if (t < t_min)
                        nodes_to_visit.push_back(curNode->child_right);
            } else {
                Vector A, B, C, e1, e2, N;
                for (int i = curNode->starting_triangle; i < curNode->ending_triangle; i++){
                    TriangleIndices triangle_indices = this->indices[i];
                    A = scaling_factor * vertices[triangle_indices.vtxi] + translation;
                    B = scaling_factor * vertices[triangle_indices.vtxj] + translation;
                    C = scaling_factor * vertices[triangle_indices.vtxk] + translation;
                    e1 = B - A;
                    e2 = C - A;
                    N = cross(e1, e2);

                    double beta = dot(cross(A - ray.O, ray.u), e2) / dot(ray.u, N);
                    double gamma = - dot(cross(A - ray.O, ray.u), e1) / dot(ray.u, N);
                    double alpha = 1. - beta - gamma;
                    double t = dot(A - ray.O, N) / dot(ray.u, N);

                    if (alpha >= 0 && beta >= 0 && gamma >= 0 && t > 0 && t < t_min) {
                        t_min = t;
                        I.intersects = true;
                        I.t = t;
                        I.P = A + beta * e1 + gamma * e2;
                        I.N = N;
                    }
                }
            }
        }
        return I;
    }

};


class Sphere : public Geometry {
    private:
        Vector C;
        Vector albedo;
        double radius;
        double refractive_index;
        bool reflective;
        bool is_hollow;

    public:
        explicit Sphere(Vector C,
                        double radius,
                        Vector albedo,
                        bool reflective = false,
                        double refractive_index = 1.,
                        bool is_hollow = false) {
            this->C = C;
            this->radius = radius;
            this->albedo = albedo;
            this->reflective = reflective;
            this->refractive_index = refractive_index;
            this->is_hollow = is_hollow;
        }

        Intersection intersect(const Ray& ray) override {
            Intersection I(this->albedo, this->refractive_index, this->reflective);

            Vector origin_center = ray.O - this->C;
            double delta = pow(dot(ray.u, origin_center), 2)
                            - (dot(origin_center, origin_center)
                            - pow(radius, 2));

            if (delta >= 0.) {
                double t1 = dot(ray.u, -1.*origin_center) - sqrt(delta);
                double t2 = dot(ray.u, -1.*origin_center) + sqrt(delta);
                I.t = (t1>0) ? t1 : ((t2>0) ? t2 : 0.0);
                I.intersects = (t2 < 0.) ? false : true;
            }

            I.P = ray.O + (I.t * ray.u);
            I.N = (I.P - this->C).normalize();
            I.N = (this->is_hollow) ? -1.*I.N : I.N;
            return I;
        }
};


class Scene {
    private:
        std::vector<Geometry*> geometries;
        Vector S;
        double intensity = 1e5;

    public:
        explicit Scene(Vector light) { S = light; }
        void addGeometry(Geometry* geometry) { geometries.push_back(geometry); }

        Vector random_cos(const Vector &N) {
            double r1 = uniform(engine);
            double r2 = uniform(engine);
            double x = sqrt(1-r2) * cos(2.*M_PI*r1);
            double y = sqrt(1-r2) * sin(2.*M_PI*r1);
            double z = sqrt(r2);

            double min = std::numeric_limits<double>::max();
            int axis = 0;
            for (int i = 0; i < 3; i++)
                if (abs(N[i]) < min) {
                    min = abs(N[i]);
                    axis = i;
                }

            Vector T1 = (axis==0) ? Vector(0., N[2], -N[1]).normalize():
                        (axis==1) ? Vector(N[2], 0., -N[0]).normalize():
                                    Vector(N[1], -N[0], 0.).normalize();
            Vector T2 = cross(N, T1);
            return T1*x + T2*y + N*z;
        }

        Intersection intersect(const Ray& ray) {
            Intersection I_main, I_temp;
            double t = std::numeric_limits<double>::max();
            for (auto& sphere : geometries) {
                I_temp = sphere->intersect(ray);
                if (I_temp.intersects && I_temp.t < t) {
                    t = I_temp.t;
                    I_main = I_temp;
                }
            }
            return I_main;
        }

        Vector getColor(const Ray& ray, int depth) {
            if (depth < 0) {return Vector(0., 0., 0.);}

            Intersection I = intersect(ray);
            Vector Lo(0., 0., 0.);

            if (I.intersects == true) {
                double eps = 1e-10;
                Vector localP = I.P + (eps * I.N);
                Vector localN = I.N;

                if (I.reflective) {
                    Ray reflected_ray = Ray(localP, ray.u - (2*dot(ray.u, localN) * localN));
                    return getColor(reflected_ray, depth - 1);
                }

                if (I.refractive_index != 1.) {
                    double Nu = dot(ray.u, localN);
                    double n1 = (Nu > 0.) ? I.refractive_index : 1.;
                    double n2 = (Nu > 0.) ? 1. : I.refractive_index;
                    localN = (Nu > 0.) ? -1.*localN : localN;

                    localP = I.P - (eps * localN);
                    Nu = dot(ray.u, localN);
                    if (1. - pow(n1/n2, 2) * (1. - pow(Nu, 2)) > 0.) {
                        Vector wt_T = (n1/n2) * (ray.u - Nu * localN);
                        Vector wt_N = -1. * localN * sqrt(1. - pow(n1/n2, 2) * (1 - pow(Nu, 2)));
                        Vector w_t = wt_T + wt_N;
                        double k0 = pow((n1-n2)/(n1+n2), 2);
                        double R = k0 + (1-k0) * pow(1 - abs(dot(localN, w_t)), 5);
                        if (uniform(engine) < R) {
                            Ray reflected_ray = Ray(localP, ray.u - (2*dot(ray.u, I.N) * I.N));
                            return getColor(reflected_ray, depth - 1);
                        } else {
                            Ray refracted_ray = Ray(localP, w_t);
                            return getColor(refracted_ray, depth - 1);
                        }
                    } else {
                        /* total internal reflection */
                        Ray reflected_ray = Ray(localP, ray.u - (2*dot(I.N, ray.u) * I.N));
                        return getColor(reflected_ray, depth - 1);
                    }
                }

                // add direct lighting if diffuse
                double d = (S - localP).norm();
                Vector w_i = (S - localP).normalize();
                Intersection I_light = intersect(Ray(S, w_i*(-1.)));
                double visibility = (!I_light.intersects || I_light.t > d) ? 1. : 0.;
                Lo = intensity / (4*M_PI*d*d) * I.albedo / M_PI * visibility * std::max(0., dot(w_i, localN));

                // indirect lighting part
                Ray random_ray = Ray(localP, random_cos(localN));
                Lo = Lo + I.albedo*getColor(random_ray, depth - 1);
            }

            return Lo;
        }
};

void boxMuller (double stdev, double &x, double &y) {
double r1 = uniform(engine);
double r2 = uniform(engine) ;
x = sqrt(-2 * log(r1)) * cos(2 * M_PI * r2) * stdev;
y = sqrt(-2 * log(r1)) * sin(2 * M_PI * r2) * stdev;
}

int main() {
    auto start = std::chrono::high_resolution_clock::now();

    Scene scene = Scene(Vector(-10, 20, 40));

    // create the large spheres for the edges
    Sphere ceiling = Sphere(Vector(0, 1000, 0), 940, Vector(1, 0, 0)); // red
    Sphere floor = Sphere(Vector(0, -1000, 0), 990, Vector(0, 0, 1)); // blue
    Sphere front = Sphere(Vector(0, 0, -1000), 940, Vector(0, 1, 0)); // green
    Sphere back = Sphere(Vector(0, 0, 1000), 940, Vector(1, 0, 1)); // pink
    Sphere left = Sphere(Vector(1000, 0, 0), 940, Vector(0, 1, 1)); // turquoise
    Sphere right = Sphere(Vector(-1000, 0, 0), 940, Vector(1, 1, 0)); // yellow

    // create the three spheres
    Sphere mirror = Sphere(Vector(20, 0, 0), 10, Vector(1., 1., 1.), true, 1.5);
    Sphere refracted = Sphere(Vector(0, 0, 0), 10, Vector(1., 1., 1.), false, 1.5);
    Sphere hollow_out_part = Sphere(Vector(-20, 0, 0), 10, Vector(1., 1., 1.), false, 1.5);
    Sphere hollow_in_part = Sphere(Vector(-20, 0, 0), 9.8, Vector(1., 1., 1.), false, 1.5, true);

    // add spheres to the scene
    scene.addGeometry(&mirror);
    scene.addGeometry(&refracted);
    scene.addGeometry(&hollow_out_part);
    scene.addGeometry(&hollow_in_part);
    scene.addGeometry(&ceiling);
    scene.addGeometry(&floor);
    scene.addGeometry(&front);
    scene.addGeometry(&back);
    scene.addGeometry(&left);
    scene.addGeometry(&right);

    // add cat to the scene
    // TriangleMesh cat = TriangleMesh(0.6, Vector(0, -10, 0), Vector(1., 1., 1.));
    // cat.readOBJ("cat.obj");
    // scene.addGeometry(&cat);

    int W = 512;
    int H = 512;
    std::vector<unsigned char> image(W*H*3, 0);

    Vector camera_pos = Vector(0, 0, 55);
    double angle = 1.0472; // 60 degrees in radians
    double gamma = 2.2;
    int max_depth = 2;
    int rays_per_pixel = 10;

    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < H; i++)
        for (int j = 0; j < W; j++) {
            Vector color = Vector(0., 0., 0.);
            double x, y;
            for (int k = 0; k < rays_per_pixel; k++) {
                boxMuller(0.5, x, y);
                Vector pixel = Vector(camera_pos[0]+(j+x)+0.5-W/2,
                                      camera_pos[1]-(i+y)-0.5+H/2,
                                      camera_pos[2]-W/(2*tan(angle/2)));
                Ray ray = Ray(camera_pos, (pixel-camera_pos).normalize());
                color = color + scene.getColor(ray, max_depth);
            }

            image[(i * W + j) * 3 + 0] = std::min(255., pow(color[0]/rays_per_pixel, 1./gamma)*255);
            image[(i * W + j) * 3 + 1] = std::min(255., pow(color[1]/rays_per_pixel, 1./gamma)*255);
            image[(i * W + j) * 3 + 2] = std::min(255., pow(color[2]/rays_per_pixel, 1./gamma)*255);
        }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
    std::cout << "Rendering took " << duration.count()/1000 << " seconds." << std::endl;

    stbi_write_png("output.png", W, H, 3, &image[0], 0);
    return 0;
}





