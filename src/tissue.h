#pragma once

#include "util.h"
#include "cell.h"
#include "grid.h"
#include "model.h"

typedef struct {
    int level;
    int index;
    Cell *cell;
} CellWrapper;

class CellGrid {
    public:
        CellGrid(Model &model, num spacing, bool);
        CellGrid(Model &model, num spacing);
        ~CellGrid() = default;
        Cell &getCell(int x, int y);
        vector<Cell*> &getAllCells();
        vector<vector<CellWrapper*>> grid;
        void adjustCellPlacement();

    private:
        Model model;
        num spacing;
};

class Tissue {
    public:
        Tissue(Model &model, bool);
        Tissue(Model &model);
        ~Tissue() = default;
        void addInitialCell();
        void updateGrid(SO::Matrix3X &newPoints);
        SO::Matrix3X &getVertexMat();
        vector<std::shared_ptr<SO::Constraint>> &getSOEdges();
        vector<std::shared_ptr<SO::Constraint>> &getSOAngles();
        vector<std::shared_ptr<SO::Constraint>> &getSOAnglesCubic();
        vector<std::shared_ptr<SO::Constraint>> &getSOAnglesTemp();
        vector<Cell*> &getCells();
        void importJSON(json save);
        void importJSONFull();
        void importJSONFull(string filename);
        void exportJSON();
        void exportJSON(string filename);
        void buildCubic();
        void addLayersCubic();
        Grid &getGrid();
        Grid grid;
        CellGrid cellGrid;
        num scaleFactor;

    private:
        void altStart(int rows, int cols);
        vector<Cell*> cells;
        vector<std::shared_ptr<SO::Constraint>> sOAnglesT;
        Model model;
        bool cubic;
};
