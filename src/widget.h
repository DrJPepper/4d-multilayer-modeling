#pragma once

#include "util.h"
#include "simulation.h"
#include <CImg.h>
using namespace cimg_library;

#ifndef WIDGET_H
#define WIDGET_H

namespace Ui {
    class Widget;
}

class ThemeMap {
    public:
        vtkSmartPointer<vtkActor> primaryPoints;
        double primaryPointsR;
        vtkSmartPointer<vtkActor> primarySpline;
        vtkSmartPointer<vtkActor> cellPointsAll;
        double cellPointsAllR;
        vtkSmartPointer<vtkActor> cellPointsUsing;
        double cellPointsUsingR;
        vtkSmartPointer<vtkActor> cellSplines;
        vtkSmartPointer<vtkActor> hexVertices;
        double hexVerticesR;
        vtkSmartPointer<vtkActor> vSplines;

        vtkSmartPointer<vtkActor> paramModelPoints;
        double paramModelPointsR;
        vtkSmartPointer<vtkActor> paramModelFaces;
        vtkSmartPointer<vtkActor> paramPointsPrimary;
        double paramPointsPrimaryR;
        vtkSmartPointer<vtkActor> paramControlPoints;
        double paramControlPointsR;
        vtkSmartPointer<vtkActor> paramSplines;
        vtkSmartPointer<vtkActor> paramPrimarySpline;

        vtkSmartPointer<vtkActor> tris;
};

class Theme {
    public:
        string background;
        string primaryPoints;
        string primarySpline;
        string cellPointsAll;
        string cellPointsUsing;
        string cellSplines;
        string hexVertices;
        string vSplines;

        string paramBackground;
        string paramModelPoints;
        string paramModelFaces;
        string paramPointsPrimary;
        string paramControlPoints;
        string paramSplines;
        string paramPrimarySpline;

        string triBackground;
        string tris;

        string scaleBarFont;
};

class Worker : public QObject {
    Q_OBJECT
    public:
        Worker() = default;
        ~Worker() = default;
    public slots:
        void cellFinder(Widget *widget, double *pickPos);
        void runSO(Widget *, Simulation *sim, int epochs, int genAgentTD);
    signals:
        void finished(const QString&);
        void finishedSO();
};

class Widget : public QWidget {
        Q_OBJECT
    public:
        explicit Widget(QWidget *parent = 0, Simulation *simulation = 0, string inputFileName = "");
        ~Widget();
        void removeActor(vtkSmartPointer<vtkActor> actor);
        void evaluate(bool render);
        num evaluateChamfer(bool render);
        Simulation *getSimulation();
        Ui::Widget *ui;
        vtkSmartPointer<vtkActor> centerActor;
        void writeImageWithText(string fName);

    signals:
        // Fires up the worker that's used in the click callback function.
        // Necessary to use a signal because of how QThreads operate.
        void operate(Widget*, double*);
        void operateSO(Widget*, Simulation*, int, int);
        void evaluateDone();

    private Q_SLOTS:
        void on_lightModeBox_stateChanged();
        void on_continueButton_clicked();
        void on_flipWindowsButton_clicked();
        void on_evaluateButton_clicked();
        void on_centerButton_clicked();
        void on_triToggleButton_clicked();
        void on_resetGridButton_clicked();
        void on_printButton_clicked();
        void on_pngPrintButton_clicked();
        void on_jsonPrintButton_clicked();
        void on_visPrintButton_clicked();
        void evaluateSlotBatch();
        void evaluateSlotSingle();
        void playBell();

    private:
        void changeTheme(Theme &theme);
        void updatePBar(int, int);
        void renderGrid();
        void renderGridAndSave();
        void renderAgency();
        void renderTD();
        void triangulateGrid();
        vtkSmartPointer<vtkActor> triangulateGrid(bool render);
        void runSO();
        void runBatch();
        void renderModel();
        void renderCellMap();
        void wheelEvent(QWheelEvent *event);
        void finishInitCubic();
        void duplicateLayerCubic();
        vtkSmartPointer<vtkActor> renderSpheres(vtkSmartPointer<vtkPoints> points, string color, double radius);
        vtkSmartPointer<vtkPoints> renderSpheres(MatrixXd &points, string color, double radius, vtkSmartPointer<vtkActor>&);
        vtkSmartPointer<vtkActor> renderSpheres(vtkSmartPointer<vtkPoints> points, string color, double radius, double opacity);
        vtkSmartPointer<vtkPoints> renderSpheres(MatrixXd &points, string color, double radius, double opacity, vtkSmartPointer<vtkActor>&);
        vtkSmartPointer<vtkActor> renderPolyline(MatrixXd &points, string color, double radius);
        vtkSmartPointer<vtkActor> renderPolylines(vector<MatrixXd*> lines, string color, double radius, VectorXd &offset);
        vtkSmartPointer<vtkActor> renderPolylines(vector<MatrixXd*> lines, string color, double radius, VectorXd &offset, int modVal);
        vtkSmartPointer<vtkActor> renderLines(vtkSmartPointer<vtkPoints> points, MatrixXi &lines, string color, double radius);
        vtkSmartPointer<vtkActor> renderLines(vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkCellArray> lines, string color, double radius);
        vtkSmartPointer<vtkActor> renderLines(vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkCellArray> lines, string color, double radius, double opacity);
        vtkSmartPointer<vtkActor> renderLines(vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkCellArray> lines, vtkSmartPointer<vtkFloatArray> colors, double radius, double opacity);
        vtkSmartPointer<vtkRenderer> renderer, rendererTri, paramRenderer;
        vtkSmartPointer<vtkGenericOpenGLRenderWindow> vtkWindow, vtkWindowTri, paramWindow;
        vtkSmartPointer<vtkJSONSceneExporter> exporter;
        vtkSmartPointer<vtkOrientationMarkerWidget> axesWidget;
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
        vtkSmartPointer<vtkNamedColors> colors;
        Simulation *sim;
        ThemeMap themeMap;
        bool runAgents, setKs, runningSO, firstBatchRun, pausing, firstWrite, triangulated, windowsFlipped, firstRenderGrid, pastRenderModel;
        QTimer *timer;
        enum ScrollBehavior { rotX, rotY, rotZ, transX, transY, transZ, zoom };
        ScrollBehavior scrollBehavior;
        double defaultScrollValue, scaleFactor, kLinear, kA120, kA90, currErr;
        int counter, startTime, batchIndex, batchTotal;
        vtkSmartPointer<vtkActor> tempPrism;
        vtkSmartPointer<vtkCubeSource> tempCubeSource;
        vtkSmartPointer<vtkScalarBarActor> scalarBar;
        string inputFileName, outStringBatch;
        std::deque<void (Widget::*)()> funDeq;
        Grid *grid;
        Theme theme, light, dark;
        Model &model;
        std::ofstream csvFile;
        json sceneVis;
        bool gridHasBeenRendered, cubic;
};

void clickCallbackFunction(vtkObject* caller, long unsigned int eventId,
        void* clientData, void* callData);

#endif  // WIDGET_H
