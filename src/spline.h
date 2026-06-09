#pragma once

#include "bezier.h"
#include "util.h"

enum Direction {FORWARD, REVERSE, CENTER_OUT};
enum SplineType {PRIMARY, VERTICAL, CELL, UNSET};

class Spline {
    public:
        Spline();
        Spline(bool isUSpline, num uvValue, num offset);
        Spline(bool isUSpline, bool setSpacing, num uvValue, num offset);
        Spline(bool isUSpline, Direction direction, num uValue, num vValue, num offset, num length);
        Spline(bool isUSpline, bool setSpacing, Direction direction, bool shift, num uValue, num vValue, num offset, num length);
        Spline duplicate();
        void combine(Spline &splineF, Spline &splineR);
        bool isUSpline, setSpacing, shift, realConstructor, bump;
        num uValue, vValue, offset, length;
        int index;
        MatrixXd uvPoints, realPoints, splinePoints, splinePointsUV, tangents, cellCentersCart, cellCentersUV, cellCenterNormals;
        vector<num> curvatures;
        Direction direction;
        SplineType sType;
        bool halfSpline;

    private:
        void defaultValues();
};

// Builds on bezier.h/cpp and extends those functions to the splines generated
// for CAD models
MatrixXd genSpline(MatrixXd &points, MatrixXd &tangents, int n);
MatrixXd genSpline(MatrixXd &points, MatrixXd &tangents);
MatrixXd initialize(MatrixXd &, MatrixXd*, float, int &shrinkage);
MatrixXd Q(MatrixXd*, float u);
MatrixXd tangent(MatrixXd*, float);
MatrixXd genCurveMat(MatrixXd*, int);
MatrixXd hermiteToBezier(MatrixXd*, MatrixXd*, int);
void smoothPoints(MatrixXd &points, int index);
