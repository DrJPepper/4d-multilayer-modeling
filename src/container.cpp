#include "grid.h"

Container::Container(Cell *cell) {
    if (cell != nullptr)
        cells.insert(cell);
    initialized = false;
}

void Container::addCell(Cell *cell) {
    cells.insert(cell);
}

void Container::setVtkPoints(vtkSmartPointer<vtkPoints> vtkPts) {
    this->vtkPts = vtkPts;
    initialized = true;
}

void Container::removeCell(Cell *cell) {
    cells.erase(cell);
}

std::set<Cell*> &Container::getCells() {
    auto tempCells = new std::set<Cell*>();
    *tempCells = cells;
    return *tempCells;
}

std::set<Cell*> &Container::getCells(Cell *excludeCell) {
    auto tempCells = new std::set<Cell*>();
    *tempCells = getCells();
    tempCells->erase(excludeCell);
    return *tempCells;
}

bool Container::vtkInitialized() {
    return initialized;
}

void Container::setActor(vtkSmartPointer<vtkActor> actor) {
    this->actor = actor;
}

vtkSmartPointer<vtkActor> Container::getActor() {
    return actor;
}
