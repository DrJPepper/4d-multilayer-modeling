#pragma once

#include "util.h"
#include "cell.h"

class Cell;
class Edge;
class Grid;
class Widget;

class Container {
    public:
        std::set<Cell*> &getCells();
        std::set<Cell*> &getCells(Cell *excludeCell);
        bool vtkInitialized();
        void setActor(vtkSmartPointer<vtkActor> actor);
        vtkSmartPointer<vtkActor> getActor();
        vtkSmartPointer<vtkActor> actor;
        // Putting this here maybe temporarily
        void removeCell(Cell *cell);
        void setVtkPoints(vtkSmartPointer<vtkPoints> vtkPts);
        std::set<Cell*> cells;

    protected:
        explicit Container(Cell *cell);
        ~Container() = default;
        void addCell(Cell *cell);
        vtkSmartPointer<vtkPoints> vtkPts;
        bool initialized;
        friend Grid;
};

class Vertex : public Container {
    public:
        explicit Vertex(Vector3d position, Cell *cell);
        ~Vertex() = default;
        Vector3d &getPosition();
        void setSphereSource(vtkSmartPointer<vtkSphereSource> source);
        void updateVtkEntities();
        int getIndex();
        json makeEntity();
        json makeEntity(float red, float green, float blue, float opacity);
        void setPtId(vtkIdType id);
        void reset();
        int cartCount;
        Vector3d cartEstimate;
        vector<Vector3d> cartEstVec;
        vector<Edge*> edges;
        Vector3d originalPosition;

    protected:
        void updatePosition(Vector3d &position);
        void updateIndex(int index);
        friend Grid;
        friend Edge;

    private:
        Vector3d position;
        vtkSmartPointer<vtkSphereSource> sphereSource;
        vtkIdType ptId;
        int index;
        static inline int currIndex = -1;
};

// TODO: Try to make these << operations const
std::ostream& operator<<(std::ostream &s, Vertex &v);

class Edge : public Container {
    public:
        num getRestLength();
        num getCurrentLength();
        void setRestLength(num restLength);
        void setRestLength(num restLength, bool avgWithNeighbor);
        Vertex *getV1();
        Vertex *getV2();
        int getIndex();
        void setLineSource(vtkSmartPointer<vtkLineSource> source);
        void updateVtkEntities();
        void colorVtkEntities(bool);
        void colorVtkEntitiesTD();
        json makeEntity();
        json makeEntity(float red, float green, float blue, float opacity);
        void setPtIds(vtkIdType id1, vtkIdType id2);
        void setLineColorsArray(vtkIdType cid, vtkSmartPointer<vtkFloatArray> vtkLineColors);
        num setAction(num, num);
        num L0, origRestLength, curvature, tdFinalLength;
        MatrixXd splinePoints, splinePointsFull;
        bool isVertical, setViaAverage, isSet, isInitialized;
        int type, correctAction;
        std::unordered_map<Cell*, int> agentActions;
        std::deque<int> tdActions;

    protected:
        explicit Edge(Vertex *v1, Vertex *v2, Cell *cell);
        explicit Edge(Vertex *v1, Vertex *v2, Cell *cell, num restLength);
        ~Edge() = default;
        void addCell(Cell *cell);
        void removeCell(Cell *cell);
        void updateIndex(int index);
        friend Grid;
        friend std::ostream& operator<<(std::ostream &s, Edge &e);

    private:
        Vertex *v1;
        Vertex *v2;
        int index;
        num restLength;
        vtkSmartPointer<vtkLineSource> lineSource;
        vtkIdType ptId1;
        vtkIdType ptId2;
        vtkIdType cid;
        vtkSmartPointer<vtkFloatArray> lineColors;
        bool restHasBeenSet;
        static inline int currIndex = 0;
};

std::ostream& operator<<(std::ostream &s, Edge &e);

class Grid {
    public:
        Grid();
        ~Grid() = default;

        SO::Matrix3X &getVertexMat();
        SO::Matrix3X &getVertexMatV();
        vector<std::shared_ptr<SO::Constraint>> &getSOEdges();
        vector<std::shared_ptr<SO::Constraint>> &getSOEdgesV();
        vector<Edge*> &getEdges();
        vector<Edge*> &getEdges(bool);
        std::set<Vertex*> getVertices();
        std::vector<Vertex*> getVerticesV();
        void setVertexMat(SO::Matrix3X &newVertMat);

        vector<Vertex*> &importJSON(json save);

        vector<Vertex*> &addVertices(MatrixXd &vertices, vector<int>&);
        vector<Vertex*> &addVertices(Cell *cell, MatrixXd &vertices);
        vector<Vertex*> &addVertices(Cell *cell, MatrixXd &vertices, vector<int>&);
        vector<Vertex*> &addVertices(MatrixXd &vertices);
        void addCellToVertices(Cell *cell, vector<Vertex*> &vertices);

        vector<Edge*> &addEdges(Cell *cell, vector<vector<Vertex*>> &vertPairs);
        vector<Edge*> &addEdges(Cell *cell, vector<vector<Vertex*>> &vertPairs, num restLength);
        vector<Edge*> &addEdges(vector<vector<Vertex*>> &vertPairs);
        vector<Edge*> &addEdges(vector<vector<Vertex*>> &vertPairs, num restLength);
        Edge *addEdge(Cell *cell, vector<Vertex*> &vertPair);
        Edge *addEdge(Cell *cell, vector<Vertex*> &vertPair, bool setMap);
        Edge *addEdge(vector<Vertex*> &vertPair);
        Edge *addEdge(vector<Vertex*> &vertPair, bool setMap);
        vector<num> getEdgeLengths();
        vector<num> getEdgeLengths(bool);
        num calcEdgeLengthDiff(vector<num>& prevLengths);
        void addCellToEdges(Cell *cell, vector<Edge*> &edges);
        void addCellToEdge(Cell *cell, Edge *edge);
        void cullCellFromGrid(Cell *cell);
        void setWidget(Widget *widget);
        void reset();
        std::vector<Vertex*> verticesV;
        std::set<Vertex*> vertices;
        std::unordered_map<int, Edge*> edgeMap;
        std::unordered_map<int, Vertex*> vertMap;
        std::set<Vector3d*> deletedCCs;
        num kEdge, kAng, kAng90, maxShrinkage;
        uint epochsRun;
        bool freshReset;

    private:
        void flagGridChange();
        SO::Matrix3X vertMat;
        vector<std::shared_ptr<SO::Constraint>> sOEdges;
        vector<std::shared_ptr<SO::Constraint>> sOAngles;
        bool gridChangedM, gridChangedE, tempVThing;
        Widget *widget;
        vector<Edge*> edges;
        vector<Edge*> edgesHoriz;
};
