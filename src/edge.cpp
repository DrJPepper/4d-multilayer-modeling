#include "grid.h"

Edge::Edge(Vertex *v1, Vertex *v2, Cell *cell, num restLength) : Container(cell) {
    L0 = 0.0;
    this->isSet = false;
    this->isInitialized = false;
    this->isVertical = false;
    this->setViaAverage = false;
    this->v1 = v1;
    this->v2 = v2;
    v1->addCell(cell);
    v2->addCell(cell);
    this->restLength = restLength;
    this->origRestLength = restLength;
    restHasBeenSet = false;
    this->index = ++Edge::currIndex;
    correctAction = -1;
}

Edge::Edge(Vertex *v1, Vertex *v2, Cell *cell) :
    Edge(v1, v2, cell, distance(v1->getPosition(), v2->getPosition())) {}

num Edge::getRestLength() {
    return restLength;
}

num Edge::getCurrentLength() {
    return distance(v1->getPosition(), v2->getPosition());
}

void Edge::setRestLength(num restLength) {
    restHasBeenSet = true;
    this->restLength = restLength;
}

void Edge::setRestLength(num restLength, bool avgWithNeighbor) {
    if (restHasBeenSet) {
        this->restLength = (this->restLength + restLength) * 0.5;
    } else {
        setRestLength(restLength);
    }
}

int Edge::getIndex() {
    return this->index;
}

void Edge::updateIndex(int index) {
    this->index = index;
}

void Edge::addCell(Cell *cell) {
    Container::addCell(cell);
    v1->addCell(cell);
    v2->addCell(cell);
}

void Edge::removeCell(Cell *cell) {
    Container::removeCell(cell);
    v1->removeCell(cell);
    v2->removeCell(cell);
}


Vertex *Edge::getV1() {
    return v1;
}

Vertex *Edge::getV2() {
    return v2;
}

void Edge::setLineSource(vtkSmartPointer<vtkLineSource> source) {
    lineSource = source;
}

void Edge::setLineColorsArray(vtkIdType cid, vtkSmartPointer<vtkFloatArray> vtkLineColors) {
    this->cid = cid;
    lineColors = vtkLineColors;
}

void Edge::colorVtkEntitiesTD() {
    if (!isVertical) {
        Vector3d color = {0.0, 0.0, 0.0};
        switch (tdActions.front()) {
            case 0:
                color(1) = 1.0;
                break;
            case 1:
                color(2) = 1.0;
                break;
            case 2:
                color(0) = 1.0;
                break;
            default: 
                color << 1.0, 1.0, 1.0;
                cout << fmt::format("ERROR: Invalid action set for edge {} ({})\n", index, tdActions.front());
        }
        tdActions.pop_front();
        lineColors->SetTuple3(cid, color(0), color(1), color(2));
    }
    vtkPts->Modified();
}

num Edge::setAction(num length, num cutoff) {
    num diff = length - tdFinalLength;
    int corrVal;
    if (abs(diff) < cutoff)
        corrVal = 0;
    else if (diff < 0)
        corrVal = 1;
    else
        corrVal = 2;
    tdActions.push_back(corrVal);
    return abs(diff);
}

void Edge::colorVtkEntities(bool showError) {
    Vector3d color = {0.0, 0.0, 0.0};
    Vector2i actions;
    actions << -1, -1;
    int i = 0;
    for (auto kv : agentActions) {
        actions[i++] = kv.second;
    }
    if (!isVertical) {
        if (showError) {
            for (int i=0; i<2; i++) {
                if (actions[i] == correctAction) {
                    color(1) = 1.0;
                } else if (actions[i] > -1) {
                    color(0) = 1.0;
                }
            }
        } else {
            for (int i=0; i<2; i++) {
                switch (actions[i]) {
                    case 0:
                        color(1) += 1.0;
                        break;
                    case 1:
                        color(2) += 1.0;
                        break;
                    case 2:
                        color(0) += 1.0;
                        break;
                    default: 
                        if (actions[i] == -1 && i == 1)
                            break;
                        color << 1.0, 1.0, 1.0;
                        cout << fmt::format("ERROR: Invalid action set for edge {} ({})\n", index, actions[0]);
                }
            }
        }
        lineColors->SetTuple3(cid, color(0), color(1), color(2));
    }
    vtkPts->Modified();
}

void Edge::updateVtkEntities() {
    auto pt1 = getV1()->getPosition();
    auto pt2 = getV2()->getPosition();
    vtkPts->SetPoint(ptId1, pt1(0), pt1(1), pt1(2));
    vtkPts->SetPoint(ptId2, pt2(0), pt2(1), pt2(2));
    Vector3d color = ratioToRGB(distance(pt1, pt2) / restLength);
    lineColors->SetTuple3(cid, color(0), color(1), color(2));
    vtkPts->Modified();
}

void Edge::setPtIds(vtkIdType id1, vtkIdType id2) {
    ptId1 = id1;
    ptId2 = id2;
}


json Edge::makeEntity() {
    return makeEntity(1.0, 1.0, 1.0, 1.0);
}

json Edge::makeEntity(float red, float green, float blue, float opacity) {
    std::stringstream ss;
    ss << *this;
    MatrixXd m(2, 3);
    m.row(0) = v1->getPosition();
    m.row(1) = v2->getPosition();
    string n = "vector", d = ss.str();
    auto e = makeEntityGeneric(n, red, green, blue, opacity, d, m);
    e["rest_length"] = restLength;
    return e;
}

std::ostream& operator<<(std::ostream &s, Edge &e) {
    return s << "Edge with rest length " << e.getRestLength() << " with:\n\t" << *e.getV1() << "\n\t" << *e.getV2();
}
