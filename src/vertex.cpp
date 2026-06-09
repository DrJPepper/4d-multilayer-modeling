#include "grid.h"

Vertex::Vertex(Vector3d position, Cell *cell) :
        Container(cell) {
    this->position << position;
    originalPosition << position;
    this->index = ++Vertex::currIndex;
    cartCount = 0;
    cartEstimate = Vector3d::Zero();
}

Vector3d &Vertex::getPosition() {
    return position;
}

void Vertex::updatePosition(Vector3d &position) {
    this->position = position;
}

int Vertex::getIndex() {
    return index;
}

void Vertex::updateIndex(int index) {
    this->index = index;
}

void Vertex::setSphereSource(vtkSmartPointer<vtkSphereSource> source) {
    sphereSource = source;
}

void Vertex::updateVtkEntities() {
    auto row = getPosition();
    vtkPts->SetPoint(ptId, row(0), row(1), row(2));
    vtkPts->Modified();
}

json Vertex::makeEntity() {
    return makeEntity(1.0, 1.0, 1.0, 1.0);
}

json Vertex::makeEntity(float red, float green, float blue, float opacity) {
    std::stringstream ss;
    ss << *this;
    string n = "point", d = ss.str();
    MatrixXd m(1, 3);
    m.row(0) = position;
    return makeEntityGeneric(n, red, green, blue, opacity, d, m);
}

void Vertex::setPtId(vtkIdType id) {
    ptId = id;
    initialized = true;
}

void Vertex::reset() {
    position = originalPosition;
}

std::ostream& operator<<(std::ostream &s, Vertex &v) {
    return s << "Vertex at: " << v.getPosition();
}
