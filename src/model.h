#pragma once

#include "util.h"
#include "spline.h"
#include "bezier.h"
#include "cell.h"
#include "BarycentricOCL/mesh.h"

class Model {
    public:
        Model();
        Model(bool);
        ~Model() = default;
	Spline primarySpline;
        MatrixXd primaryUVPoints, primaryRealPoints, primarySplinePoints, primaryTangents, primaryCellCentersCart, primaryCellCentersUV;
        Mesh *mesh;
        vtkSmartPointer<vtkPoints> points;
        vtkSmartPointer<vtkCellArray> lines;
        vtkSmartPointer<vtkActor> modelActor;
        vector<MatrixXd*> vSplines, vTangents, vRealPoints, vCellCentersCart, vCellCentersUV, vCellCentersNormal, offsetSplinesPos, offsetSplinesNeg, offsetTansPos, offsetTansNeg, uSplines;
        vector<Spline*> splines, usplines;
        // Spline "stumps" for configuring the cells
        void populateSpline(bool isUSpline, bool reverse, num uValue, num vValue, num offset, num length, MatrixXd &uvPoints, MatrixXd &realPoints, MatrixXd &splinePoints, MatrixXd &tangents);
        void calcIdealPosition(Vector2d &curCenterUV, Vector3d &curCenterCart, Vector2d &compUV, num targetDist);
        void populateSpline(Spline &spline);
        void populateSpline(Spline &spline, Spline &otherSpline);
        void genSplinePointsAndCCs(Spline &spline, bool reverse);
        void bumpPoint(Vector2d &point, Vector3d &realPoint, num &horizShift, bool isUSpline, Spline &otherSpline, int &stepIn);
        void bumpPoint(Vector2d &point, Vector3d &realPoint, num &horizShift, bool isUSpline, Spline &otherSpline, int &outTri, int &inTri, int &startPoint, int &stepIn);
        num distance3dFrom2d(Vector2d p1, Vector2d p2);
        num surfaceDistance(Vector2d p1, Vector2d p2);
        num surfaceDistance(Vector3d p1, Vector3d p2);
        num spacing, zHeight;
        int startingLineNodeCount, nodeCountU, nodeCountV, countToBump;
        std::string modelFile;
        vector<std::deque<int>*> levels, indexList;

        bool batch, flipPrimary, hasDuplicatedCubic;
        num kAngIn120, kAngIn120Min, kAngIn120Max, kAng120Step;
        num kAngIn90, kAngIn90Min, kAngIn90Max, kAng90Step;
        num kLinearIn, kLinearInMin, kLinearInMax, kLinearStep;
        num surfaceArea, scaleFactor, scale;
        uint epochsIn, epochsInStep, layersCubic;
        Vector4i paramCorners;

    private:
        Vector3d findNextCC(MatrixXd *spline, Vector3d prevCellCenter, int &index, num stop, bool reverse, bool &splineEnded);
        Vector3d findNextCC(MatrixXd *spline, int &index, num stop, bool reverse, bool &splineEnded);

        int centerCC(int splineInd, MatrixXi &indicies, int direction, vector<bool> &doneVec, vector<std::deque<Vector3d>*> &ccs);
        void processK(toml::node_view<toml::node> tbl, num &kIn, num &kMin, num &kMax, num &kStep, num kDefault);
        bool spacingDiv, cubic;
        num startingU, uOffset, vOffset;
        int pointCount, bumpTotal, bumpNum;
        MatrixXd pointMap, vertMat;
        bool printCC;
        void finishInit();
        void setScaleFactor();
        void finishInitCubic();
};
