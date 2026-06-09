#include "grid.h"

Grid::Grid() {
    flagGridChange();
    kEdge = K_EDGE;
    kAng = K_ANG_120;
    kAng90 = K_ANG_90;
    epochsRun = 0;
    freshReset = true;
}

SO::Matrix3X &Grid::getVertexMatV() {
    if (gridChangedM) {
        // ShapeOp wants column major, whereas I've otherwise used row major
        vertMat = SO::Matrix3X::Zero(3, verticesV.size());
        int i = 0;
        for (auto vertex : verticesV) {
            vertMat.col(i) = vertex->getPosition();
            vertex->updateIndex(i);
            i++;
        }
        gridChangedM = false;
    }
    return vertMat;
}

SO::Matrix3X &Grid::getVertexMat() {
    if (gridChangedM) {
        // ShapeOp wants column major, whereas I've otherwise used row major
        vertMat = SO::Matrix3X::Zero(3, vertices.size());
        int i = 0;
        for (auto vertex : vertices) {
            vertMat.col(i) = vertex->getPosition();
            vertex->updateIndex(i);
            i++;
        }
        gridChangedM = false;
    }
    return vertMat;
}

void Grid::flagGridChange() {
    gridChangedM = true;
    gridChangedE = true;
}

vector<Vertex*> &Grid::importJSON(json save) {
    vector<int> indices;
    MatrixXd verts(save["vertices"].size(), 3);
    for (unsigned int i=0; i<save["vertices"].size(); i++) {
        auto v = save["vertices"][i]["position"];
        verts.row(i) << v[0], v[1], v[2];
        indices.push_back(save["vertices"][i]["index"]);
    }
    auto vertObjs = new vector<Vertex*>;
    *vertObjs = addVertices(verts, indices);
    indices.clear();
    for (unsigned int i=0; i<save["edges"].size(); i++) {
        auto e = save["edges"][i];
        vector<Vertex*> pair;
        pair.push_back(vertMap[e["vertices"][0]]);
        pair.push_back(vertMap[e["vertices"][1]]);
        auto eObj = addEdge(pair, false);
        eObj->setRestLength(e["rest_length"]);
        eObj->isVertical = e["type"] == 2;
        eObj->index = e["index"];
        edgeMap[eObj->index] = eObj;
    }
    return (*vertObjs);
}

vector<std::shared_ptr<SO::Constraint>> &Grid::getSOEdges() {
    sOEdges = vector<std::shared_ptr<SO::Constraint>>();
    // Just to make extra sure it's updated
    auto vertMatLocal = getVertexMat();
    for (auto edge : edges) {
        vector<int> pts;
        pts.push_back(edge->getV1()->getIndex());
        pts.push_back(edge->getV2()->getIndex());
        num k = edge->setViaAverage ? kEdge * 0.3 : kEdge;
        k = edge->isVertical ? K_VERT : kEdge;
        auto constraint = std::make_shared<SO::EdgeStrainConstraint>(
            pts, k, vertMatLocal, 1.0, 1.0);
        constraint->setEdgeLength(edge->getRestLength());
        sOEdges.push_back(constraint);
    }
    return sOEdges;
}

vector<std::shared_ptr<SO::Constraint>> &Grid::getSOEdgesV() {
    sOEdges = vector<std::shared_ptr<SO::Constraint>>();
    // Just to make extra sure it's updated
    auto vertMatLocal = getVertexMatV();
    for (auto edge : edges) {
        vector<int> pts;
        pts.push_back(edge->getV1()->getIndex());
        pts.push_back(edge->getV2()->getIndex());
        num k = edge->setViaAverage ? kEdge * 0.3 : kEdge;
        k = edge->isVertical ? K_VERT : kEdge;
        auto constraint = std::make_shared<SO::EdgeStrainConstraint>(
            pts, k, vertMatLocal, 1.0, 1.0);
        constraint->setEdgeLength(edge->getRestLength());
        sOEdges.push_back(constraint);
    }
    return sOEdges;
}

vector<Edge*> &Grid::getEdges(bool removeVertical) {
    if (!removeVertical)
        return edges;
    if (!edgesHoriz.size()) {
        for (auto e : edges)
            if (!e->isVertical)
                edgesHoriz.push_back(e);
    }
    return edgesHoriz;
}

vector<Edge*> &Grid::getEdges() {
    return getEdges(false);
}

std::set<Vertex*> Grid::getVertices() {
    return vertices;
}

std::vector<Vertex*> Grid::getVerticesV() {
    return verticesV;
}

vector<Vertex*> &Grid::addVertices(MatrixXd &vertices, vector<int> &indices) {
    return addVertices(nullptr, vertices, indices);
}

vector<Vertex*> &Grid::addVertices(Cell *cell, MatrixXd &vertices) {
    vector<int> none;
    return addVertices(cell, vertices, none);
}

vector<Vertex*> &Grid::addVertices(MatrixXd &vertices) {
    vector<int> none;
    return addVertices(nullptr, vertices, none);
}

vector<Vertex*> &Grid::addVertices(Cell *cell, MatrixXd &vertices, vector<int> &indices) {
    flagGridChange();
    auto newVerts = new vector<Vertex*>();
    for (unsigned int i=0; i<vertices.rows(); i++) {
        Vector3d vertexVec = vertices.row(i);
        auto currVertex = new Vertex(vertexVec, cell);
        if (indices.size())
            currVertex->index = indices[i];
        vertMap[currVertex->index] = currVertex;
        newVerts->push_back(currVertex);
        this->verticesV.push_back(currVertex);
        this->vertices.insert(currVertex);
    }
    if (cell != nullptr)
        cell->addVertices(*newVerts);
    return *newVerts;
}

vector<Edge*> &Grid::addEdges(vector<vector<Vertex*>> &vertPairs, num restLength) {
    return addEdges(nullptr, vertPairs, restLength);
}

// TODO: Un-duplicate this code
vector<Edge*> &Grid::addEdges(Cell *cell, vector<vector<Vertex*>> &vertPairs, num restLength) {
    flagGridChange();
    auto newEdges = new vector<Edge*>();
    for (auto vertPair : vertPairs) {
        auto currEdge = new Edge(vertPair[0], vertPair[1], cell, restLength);
        newEdges->push_back(currEdge);
        edgeMap[currEdge->index] = currEdge;
        edges.push_back(currEdge);
    }
    if (cell != nullptr)
        cell->addEdges(*newEdges);
    return *newEdges;
}

vector<Edge*> &Grid::addEdges(vector<vector<Vertex*>> &vertPairs) {
    return addEdges(nullptr, vertPairs);
}

vector<Edge*> &Grid::addEdges(Cell *cell, vector<vector<Vertex*>> &vertPairs) {
    flagGridChange();
    auto newEdges = new vector<Edge*>();
    for (auto vertPair : vertPairs) {
        auto currEdge = new Edge(vertPair[0], vertPair[1], cell);
        newEdges->push_back(currEdge);
        edgeMap[currEdge->index] = currEdge;
        edges.push_back(currEdge);
    }
    if (cell != nullptr)
        cell->addEdges(*newEdges);
    return *newEdges;
}

Edge *Grid::addEdge(vector<Vertex*> &vertPair) {
    return addEdge(nullptr, vertPair, true);
}

Edge *Grid::addEdge(vector<Vertex*> &vertPair, bool setMap) {
    return addEdge(nullptr, vertPair, setMap);
}

Edge *Grid::addEdge(Cell *cell, vector<Vertex*> &vertPair) {
    return addEdge(cell, vertPair, true);
}

Edge *Grid::addEdge(Cell *cell, vector<Vertex*> &vertPair, bool setMap) {
    flagGridChange();
    auto edge = new Edge(vertPair[0], vertPair[1], cell);
    if (setMap)
        edgeMap[edge->index] = edge;
    edges.push_back(edge);
    if (cell != nullptr)
        cell->addEdge(edge);
    return edge;
}

void Grid::addCellToVertices(Cell *cell, vector<Vertex*> &vertices) {
    flagGridChange();
    for (auto vertex : vertices) {
        vertex->addCell(cell);
    }
    cell->addVertices(vertices);
}

void Grid::addCellToEdges(Cell *cell, vector<Edge*> &edges) {
    flagGridChange();
    for (auto edge : edges) {
        edge->addCell(cell);
    }
    cell->addEdges(edges);
}

void Grid::addCellToEdge(Cell *cell, Edge *edge) {
    flagGridChange();
    edge->addCell(cell);
    cell->addEdge(edge);
}

void Grid::setWidget(Widget *widget) {
    this->widget = widget;
}

void Grid::setVertexMat(SO::Matrix3X &newVertMat) {
    vertMat = newVertMat;
    for (auto vertex : vertices) {
        Vector3d p = vertMat.col(vertex->getIndex());
        vertex->updatePosition(p);
    }
}

void Grid::reset() {
    gridChangedM = true;
    gridChangedE = true;
    freshReset = true;
    epochsRun = 0;
    for (auto v : vertices)
        v->reset();
    for (auto e : edges)
        e->setRestLength(e->origRestLength);
}

void Grid::cullCellFromGrid(Cell *cell) {
    auto v = new Vector3d();
    *v << cell->initialCenter;
    deletedCCs.insert(v);
    if (cell->neighbors.north)
        cell->neighbors.north->neighbors.south = nullptr;
    if (cell->neighbors.south)
        cell->neighbors.south->neighbors.north = nullptr;
    if (cell->neighbors.northEast)
        cell->neighbors.northEast->neighbors.southWest = nullptr;
    if (cell->neighbors.northWest)
        cell->neighbors.northWest->neighbors.southEast = nullptr;
    if (cell->neighbors.southEast)
        cell->neighbors.southEast->neighbors.northWest = nullptr;
    if (cell->neighbors.southWest)
        cell->neighbors.southWest->neighbors.northEast = nullptr;
    for (auto edge : cell->getEdges()) {
        edge->removeCell(cell);
        if (!edge->cells.size()) {
            edges.erase(std::remove(edges.begin(), edges.end(), edge), edges.end());
        }
    }
    for (auto vertex : cell->getVertices()) {
        vertex->removeCell(cell);
        if (!vertex->cells.size()) {
            vertices.erase(vertex);
            for (unsigned int i=0; i<verticesV.size(); i++) {
                if (verticesV[i] == vertex)
                    verticesV.erase(verticesV.begin()+i);
            }
        }
    }
}

vector<num> Grid::getEdgeLengths() {
    return getEdgeLengths(false);
}

vector<num> Grid::getEdgeLengths(bool setTDFinal) {
    auto lengths = vector<num>();
    for (auto e : getEdges(true)) {
        lengths.push_back(e->getCurrentLength());
        if (setTDFinal)
            e->tdFinalLength = e->getCurrentLength();
    }
    return lengths;
}

num Grid::calcEdgeLengthDiff(vector<num>& prevLengths) {
    num avg = 0.0;
    for (unsigned int i=0; i<getEdges(true).size(); i++) {
        avg += abs(edgesHoriz[i]->getCurrentLength() - prevLengths[i]);
    }
    avg /= edgesHoriz.size();
    return avg;
}
