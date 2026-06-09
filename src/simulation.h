#pragma once

#include "util.h"
#include "tissue.h"
#include "model.h"

struct EvalResults {
    vtkSmartPointer<vtkActor> actor, modelActor;
    MatrixXd result, dst2;
    num avgChamf, error, maxD;
};

class Simulation : public QObject {
    Q_OBJECT
    public:
        Simulation();
        Simulation(bool);
        ~Simulation() = default;
        void runSO();
        void runSO(int);
        void runSO(int, int);
        void runHeadless();
        Tissue &getTissue();
        Model &getModel();
        vtkSmartPointer<vtkActor> triangulateGrid();
        vtkSmartPointer<vtkActor> triangulateGrid(Eigen::Vector4i&);
        vtkSmartPointer<vtkActor> triangulateGridCubic(Eigen::Vector4i&);
        EvalResults &evaluateChamfer();
        Eigen::Array2d getAvgSpringDisp();
        std::ofstream videoTxtFile;
        num secPerEpoch, stepSize, cutoff, divAvg, tdEpochSteps;
        std::shared_ptr<SO::Solver> s;

    signals:
        void renderSignal();

    protected:
        Model model;
        bool firstRun;

    private:
        void runSOCubic(int);
        void *agentPointer;
        Tissue tissue;
        vector<num> finalEdges;
        bool cubic;
};
