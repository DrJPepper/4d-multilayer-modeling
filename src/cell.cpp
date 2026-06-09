#include "cell.h"

Cell::Cell() {
    this->orderingInitialized = false;
    this->signal = 5;
    if (!tets.size()) {
        tets = MatrixXd(12, 4);
        tets << 0, 1, 2, 6,
             1, 2, 7, 6,
             2, 8, 7, 6,
             0, 2, 3, 6,
             2, 3, 8, 6,
             3, 9, 8, 6,
             0, 3, 4, 6,
             3, 4, 9, 6,
             4, 10, 9, 6,
             0, 4, 5, 6,
             4, 5, 10, 6,
             5, 11, 10, 6;
    }
}

Cell::Cell(Vector3d centerCart, Vector2d centerUV) : Cell() {
    this->initialCenter = centerCart;
    this->origInitialCenter = centerCart;
    this->initialCenterUV = centerUV;
    this->origInitialCenterUV = centerUV;
    neighbors.north = nullptr;
    neighbors.south = nullptr;
    neighbors.northWest = nullptr;
    neighbors.northEast = nullptr;
    neighbors.southWest = nullptr;
    neighbors.southEast = nullptr;
    initialized = false;
    this->index = ++Cell::currIndex;
}

void Cell::addVertices(vector<Vertex*> &vertices) {
    for (auto vertex : vertices) {
        if (!std::count(this->vertices.begin(), this->vertices.end(), vertex)) {
            this->vertices.push_back(vertex);
        }
    }
}

void Cell::addEdges(vector<Edge*> &edges) {
    vector<Vertex*> verts;
    int i=0, abv, ind;
    string dir;
    for (auto edge : edges) {
        if (!std::count(this->edges.begin(), this->edges.end(), edge)) {
            this->edges.push_back(edge);
            verts.push_back(edge->getV1());
            verts.push_back(edge->getV2());
        }
    }
    addVertices(verts);
    if (this->edges.size() == 18) {
        for (auto edge : this->edges) {
            abv = i % 3;
            ind = (i++ / 3 + 1) % 6;
            switch (abv) {
                case 0:
                    *edgeMap.basal.m[AllDirs[ind]] = edge;
                    dir = "basal";
                    break;
                case 1:
                    *edgeMap.apical.m[AllDirs[ind]] = edge;
                    dir = "apical";
                    break;
            }
        }
    }
}

void Cell::addEdge(Edge* edge) {
    if (!std::count(this->edges.begin(), this->edges.end(), edge))
        this->edges.push_back(edge);
}

std::vector<Cell*> Cell::getNeighborsVec() {
    auto n = neighbors;
    vector<Cell*> nv = {n.northEast, n.southEast, n.south, n.southWest, n.northWest, n.north};
    return nv;
}

vector<Edge*> &Cell::getEdges() {
    return edges;
}

vector<Edge*> &Cell::getEdges(bool removeVertical) {
    if (!removeVertical)
        return edges;
    auto edgesNV = new vector<Edge*>();
    for (auto e : edges)
        if (!e->isVertical)
            edgesNV->push_back(e);
    return *edgesNV;
}

vector<Vertex*> &Cell::getVertices() {
    return vertices;
}

std::ostream& operator<<(std::ostream &s, Cell &c) {
    s << "Cell " << &c << " with:\n";
    auto edges = c.getEdges();
    auto vertices = c.getVertices();
    if (c.curvatures.size())
        s << "- Curvatures: " << c.curvatures[0] << ", " << c.curvatures[1] << ", " << c.curvatures[2] << endl;
    s << "- Center: " << c.initialCenter(0) << ", " << c.initialCenter(1) << ", " << c.initialCenter(2) << endl;
    s << "- Level: " << c.level << endl;
    s << "- Index: " << c.index << endl;
    s << "- V Spline Index: " << c.vSplineIndex << endl;
    for (auto it = edges.begin(); it != edges.end(); it++) {
        s << "- " << **it;
        if (std::next(it) != edges.end()) 
            s << endl;
    }
    return s;
}

num Cell::getVolume() {
    num vol = 0.0;
    for (VectorXd r : tets.rowwise()) {
        vol += tetrahedronVolume(vertices[r(0)]->getPosition(),
                vertices[r(1)]->getPosition(), vertices[r(2)]->getPosition(),
                vertices[r(3)]->getPosition());
    }
    return vol;
}

void Cell::setRestLengths() {
    int i, index;
    num L00, L01, maxCurvShrinkage = 0.0;
    num Ln0, Ln1, R, Rmin, zHeight;
    // Calculating this again for every cell is kind of a waste,
    // although not really enough to bother to cache it somewhere
    zHeight = distance(vertices[0]->getPosition(), vertices[6]->getPosition());

    num angleTemp, angle0, angle1, kappa;
    for (i = 0; i < 3; i++) {
        kappa = curvatures[i];
        Rmin = zHeight / maxShrinkage;
        if (kappa == 0.0) {
            R = -1.0;
            angleTemp = 0.0;
            angle0 = 0.0;
            angle1 = 0.0;
            index = i * 3;
            L00 = edges[index]->L0;
            L01 = edges[index + 9]->L0;
        } else {
            R = 1 / abs(kappa);
            R = R < Rmin ? Rmin : R;
            angleTemp = ((kappa > 0) - (kappa < 0)) *
                        asin(edges[i*3]->L0 / 2.0 / R) * 2.0;
            index = i * 3 + (angleTemp > 0.0);
            L00 = edges[index]->L0;
            L01 = edges[index + 9]->L0;
            angle0 = ((kappa > 0) - (kappa < 0)) *
                        asin(L00 / 2.0 / R) * 2.0;
            angle1 = ((kappa > 0) - (kappa < 0)) *
                        asin(L01 / 2.0 / R) * 2.0;
        }
        if (R < 0.0) {
            Ln0 = L00;
            Ln1 = L01;
        } else {
            Ln0 = L00 - 2 * zHeight * sin(abs(angle0) / 2);
            Ln0 = std::max(Ln0, L00 * static_cast<num>(.5));
            Ln1 = L01 - 2 * zHeight * sin(abs(angle1) / 2);
            Ln1 = std::max(Ln1, L01 * static_cast<num>(.5));
            num shrink = std::max((L00 - Ln0) / L00, (L01 - Ln1) / L01);
            if (shrink > maxCurvShrinkage)
                maxCurvShrinkage = shrink;
        }
        int one = i*3, two = i*3+1;
        edges[one]->setRestLength(L00);
        edges[one + 9]->setRestLength(L00);
        edges[two]->setRestLength(L01);
        edges[two + 9]->setRestLength(L01);
        edges[index]->setRestLength(Ln0);
        edges[index + 9]->setRestLength(Ln1);
    }
}

VectorXd Cell::getAgentInput() {
    if (!orderingInitialized) {
        setAgentOrdering();
        setEdgeMap();
    }
    int s = agentOrdering.size();
    VectorXd out(s);
    for (int i=0; i<s; i++) {
        auto c = agentOrdering[i];
        out(i) = -1.0;
        if (c != nullptr)
            out(i) = c->getCurrentLength();
    }
    return out;
}

VectorXi Cell::getAgentOrdering() {
    if (!orderingInitialized) {
        setAgentOrdering();
        setEdgeMap();
    }
    int s = agentOrdering.size();
    VectorXi out(s);
    for (int i=0; i<s; i++) {
        auto c = agentOrdering[i];
        out(i) = -1;
        if (c != nullptr)
            out(i) = c->getIndex();
    }
    return out;
}

void Cell::setEdgeMap() {
    string dir;
    bool print = false;
    for (auto edge : edges) {
        if (!edge->isVertical) {
            for (auto d : AllDirs) {
                auto n = *neighbors.m[d];
                if (n != nullptr) {
                    auto v = n->getEdges();
                    if(std::find(v.begin(), v.end(), edge) != v.end()) {
                        if (edge->getV1()->getPosition()[2] > 0.0) {
                            if (*edgeMap.apical.m[d] != edge) {
                                dir = "apical";
                                print = true;
                            }
                        } else {
                            if (*edgeMap.basal.m[d] != edge) {
                                dir = "basal";
                                print = true;
                            }
                        }
                        if (print)
                            cout << fmt::format("Error: {} edge in {} direction {} does not match in cell {}\n",
                                    dir, fmt::ptr(edge), CDirStrings.at(d), fmt::ptr(this));
                        print = false;
                    }
                }
            }
        }
    }
}

void Cell::setAgentOrdering() {
    int ind = 0;
    for (auto e : edges)
        if (!e->isVertical)
            agentOrdering[ind++] = e;
    int addStart = 0;
    auto nv = getNeighborsVec();
    for (unsigned int i=0; i<nv.size(); i++) {
        auto curr = nv[i];
        auto next = nv[(i+1)%6];
        if (curr != nullptr) {
            for (int j=0; j<12; j++) {
                agentOrdering[ind++] = curr->edges[(addStart+j++)%18];
                agentOrdering[ind++] = curr->edges[(addStart+j++)%18];
            }
        } else {
            int origInd = ind;
            for (int j=0; j<8; j++) {
                agentOrdering[ind++] = nullptr;
            }
            if (next != nullptr) {
                agentOrdering[origInd+6] = next->edges[(addStart+9)%18];
                agentOrdering[origInd+7] = next->edges[(addStart+10)%18];
            }
        }
        addStart = (addStart + 3) % 18;
    }
    for (auto e : agentOrdering)
        if (e != nullptr)
            if (e->isVertical)
                cout << "Vertical\n";
    orderingInitialized = true;
}
