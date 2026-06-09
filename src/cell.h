#pragma once

#include "util.h"
#include "grid.h"

class Vertex;
class Edge;

enum CDirection {NORTHEAST, SOUTHEAST, SOUTH, SOUTHWEST, NORTHWEST, NORTH};
static const CDirection AllDirs[] = {NORTHEAST, SOUTHEAST, SOUTH, SOUTHWEST, NORTHWEST, NORTH};
static const std::map<CDirection, string> CDirStrings = {{NORTHEAST, "NorthEast"}, {SOUTHEAST, "SouthEast"}, {SOUTH, "South"}, {SOUTHWEST, "SouthWest"}, {NORTHWEST, "NorthWest"}, {NORTH, "North"}};

class Cell {

    public:
        Cell();
        Cell(Vector3d centerCart, Vector2d centerUV);
        ~Cell() = default;
        void addVertices(vector<Vertex*> &vertices);
        void addEdges(vector<Edge*> &edges);
        void addEdge(Edge* edge);
        void setRestLengths();
        VectorXd getAgentInput();
        VectorXi getAgentOrdering();
        vector<Edge*> &getEdges();
        vector<Edge*> &getEdges(bool);
        vector<Vertex*> &getVertices();
        num getVolume();
        std::vector<Cell*> getNeighborsVec();
        template<typename T> struct N {
            using value_type = T;
            T *north;
            T *northEast;
            T *southEast;
            T *south;
            T *southWest;
            T *northWest;
            std::map<CDirection, T**> m {{NORTHEAST, &northEast}, {SOUTHEAST, &southEast}, {SOUTH, &south}, {SOUTHWEST, &southWest}, {NORTHWEST, &northWest}, {NORTH, &north}};
        };
        template <typename T> using Neighbors = struct N<T>;
        typedef struct {
            Neighbors<Edge> apical;
            Neighbors<Edge> basal;
        } EdgeMap;
        Neighbors<Cell> neighbors;
        EdgeMap edgeMap;
        Vector3d initialCenter, initialGridCenter, origInitialCenter, initialCenterNormal;
        Vector2d initialCenterUV, origInitialCenterUV;
        bool initialized;
        vector<MatrixXd*> splinePoints, realPoints, tangents;
        vector<num> curvatures, distances;
        int level, index, vSplineIndex, signal;
        vector<int> edgesSetting;
        static inline num maxShrinkage;

    private:
        bool orderingInitialized;
        void setEdgeMap();
        void setAgentOrdering();
        vector<Vertex*> vertices;
        vector<Edge*> edges;
        std::array<Edge*, 60> agentOrdering;
        static inline int currIndex = 0;
        static inline MatrixXd tets = MatrixXd();
};

template<typename T>
struct A 
{
     using value_type = T;
     T someVar;
};

std::ostream& operator<<(std::ostream &s, Cell &c);
