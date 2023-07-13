//
// Created by angel_rarmbt1 on 12/07/2023.
//

#ifndef E5_CURVAS_DTRIANGULATION_H
#define E5_CURVAS_DTRIANGULATION_H

#include <glm/glm.hpp>
#include <vector>
#include <algorithm>

using namespace std;
using namespace glm;

struct Triangle3D {
    vec3 p1, p2, p3;
    Triangle3D(vec3 p1, vec3 p2, vec3 p3) : p1(p1), p2(p2), p3(p3) {}
    bool operator==(const Triangle3D& triangle) const{
        return triangle.p1 == p1 && triangle.p2 == p2 && triangle.p3 == p3;
    }
};

struct Edge {
    vec3 p1, p2;
    Edge(vec3 p1, vec3 p2) : p1(p1), p2(p2) {}
};

//Delaunay triangulation
class DTriangulation{
private:
    // Check if a given point lies inside the circumcircle of a triangle
    bool inCircle(vec3& a, vec3& b, vec3& c, vec3& d) {
        double ax = a.x - d.x;
        double ay = a.y - d.y;
        double az = a.z - d.z;
        double bx = b.x - d.x;
        double by = b.y - d.y;
        double bz = b.z - d.z;
        double cx = c.x - d.x;
        double cy = c.y - d.y;
        double cz = c.z - d.z;

        // Calculate the determinant of the matrix
        double det = (ax * (by * cz - bz * cy)) -
                     (ay * (bx * cz - bz * cx)) +
                     (az * (bx * cy - by * cx));

        return (det > 0);
    }
public:
    DTriangulation(){};
    vector<Triangle3D> performTriangulation(vector<vec3>& points){
        vector<Triangle3D> triangles;
        // Create a super triangle that encompasses all the points
        double minX = points[0].x;
        double minY = points[0].y;
        double minZ = points[0].z;
        double maxX = minX;
        double maxY = minY;
        double maxZ = minZ;

        // Calculate the min and max coords of the points to define the boundaries of the super triangle
        for (vec3 point : points) {
            if (point.x < minX) minX = point.x;
            if (point.y < minY) minY = point.y;
            if (point.z < minZ) minZ = point.z;
            if (point.x > maxX) maxX = point.x;
            if (point.y > maxY) maxY = point.y;
            if (point.z > maxZ) maxZ = point.z;
        }

        double dx = maxX - minX;
        double dy = maxY - minY;
        double dz = maxZ - minZ;
        double deltaMax = std::max(std::max(dx, dy), dz);
        double midX = (minX + maxX) / 2.0;
        double midY = (minY + maxY) / 2.0;
        double midZ = (minZ + maxZ) / 2.0;

        // Vertices of the super triangle
        vec3 p1(midX - 20 * deltaMax, midY - deltaMax, midZ);
        vec3 p2(midX, midY + 20 * deltaMax, midZ);
        vec3 p3(midX + 20 * deltaMax, midY - deltaMax, midZ);

        triangles.emplace_back(p1, p2, p3);

        for (vec3 point : points) {
            vector<Triangle3D> badTriangles;
            // Find bad triangles = if current point lies inside the circumcircle of the triangle
            for (Triangle3D triangle : triangles) {
                vec3 p1 = triangle.p1;
                vec3 p2 = triangle.p2;
                vec3 p3 = triangle.p3;
                if (inCircle(p1, p2, p3, point))
                    badTriangles.push_back(triangle);
            }

            vector<Edge> polygonEdges;
            // Get the edges of the bad triangles
            for (Triangle3D triangle : badTriangles) {
                vec3 p1 = triangle.p1;
                vec3 p2 = triangle.p2;
                vec3 p3 = triangle.p3;

                polygonEdges.emplace_back(p1, p2);
                polygonEdges.emplace_back(p2, p3);
                polygonEdges.emplace_back(p3, p1);
            }

            // Remove bad triangles
            triangles.erase(std::remove_if(triangles.begin(), triangles.end(),
                                           [&](const Triangle3D& triangle) {
                                               return std::find_if(badTriangles.begin(), badTriangles.end(),
                                                                   [&](const Triangle3D& badTriangle) {
                                                                       return triangle == badTriangle;
                                                                   }) != badTriangles.end();
                                           }), triangles.end());

            // Create triangles by connecting the vertices of the polygon edges with the current point
            for (const Edge& edge : polygonEdges) {
                if (std::find_if(polygonEdges.begin(), polygonEdges.end(),
                                 [&](const Edge& e) {
                                     return (edge.p1 == e.p2 && edge.p2 == e.p1);
                                 }) == polygonEdges.end()) {
                    triangles.emplace_back(edge.p1, edge.p2, point);
                }
            }
        }
        // Remove triangles that contain indices beyond the range of the original
        triangles.erase(std::remove_if(triangles.begin(), triangles.end(),
                                       [&](const Triangle3D& triangle) {
                                           return (triangle.p1 == p1 || triangle.p1 == p2 || triangle.p1 == p3 ||
                                                   triangle.p2 == p1 || triangle.p2 == p2 || triangle.p2 == p3 ||
                                                   triangle.p3 == p1 || triangle.p3 == p2 || triangle.p3 == p3);
                                       }), triangles.end());
        return triangles;
    }

};


#endif //E5_CURVAS_DTRIANGULATION_H
