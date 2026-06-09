#pragma once

#include "util.h"
#include "BarycentricOCL/bezier.h"

// Basic calculations needed to get positional, distance and curvature
// measurements

MatrixXd QLoop(MatrixXd &points, int count);
int factorial(int);
int kChoosei(int n, int k);
MatrixXd genPatchMat(MatrixXd*, int, int);
MatrixXd genNormMat(MatrixXd*, int, int);
Vector3d qPrime(MatrixXd*, float);
Vector3d qDoublePrime(MatrixXd*, float);
num curvature(Vector3d, Vector3d);
num curvatureFromU(MatrixXd*, float);
num curvatureFromUSigned(MatrixXd*, float, Vector3d normal);
num curvatureInU(MatrixXd *points, float u, float v);
num curvatureInV(MatrixXd *points, float u, float v);
num distanceInU(MatrixXd *points, float ustart, float uend, float v);
num distanceInV(MatrixXd *points, float u, float vstart, float vend);
num curvatureSign(MatrixXd *points, float u);
num fullSplineSignedCurvature(MatrixXd *points, MatrixXd *tangents, Vector3d normal);
vector<num> fullSplineSignedCurvatureList(MatrixXd *points, MatrixXd *tangents, vector<Vector3d> normals, int count);
MatrixXd hermiteToBezier(MatrixXd *points, MatrixXd *tangents);
