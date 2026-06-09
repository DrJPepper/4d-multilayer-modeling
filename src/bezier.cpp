#include "bezier.h"

// Runs Q equation for a whole Bezier curve
MatrixXd QLoop(MatrixXd &points, int count) {
    MatrixXd result = MatrixXd::Zero(count, 3);
    num step = 1.0 / static_cast<num>(count), u;
    for (int w = 0; w < count; w++) {
        MatrixXd acc = MatrixXd::Zero(1, 3);
        u = step * static_cast<num>(w);
        for (int i = 0; i < 4; i++) {
            acc += points.row(i) * kChoosei(3, i) * pow(1 - u, 3 - i) * pow(u, i);
        }
        result.row(w) = acc;
    }
    return result;
}

TEST_CASE("testing the Q function") {
    MatrixXd points(4, 3);
    MatrixXd result(1, 3);
    result << 6, 1.5, 0;
    points << 0, 1.5, 0, 2, 1.5, 0, 4, 1.5, 0, 6, 1.5, 0;
    CHECK(Q(&points, 1.0) == result);
}

TEST_CASE("testing the factorial function") {
    CHECK(factorial(0) == 1);
    CHECK(factorial(1) == 1);
    CHECK(factorial(2) == 2);
    CHECK(factorial(3) == 6);
    CHECK(factorial(10) == 3628800);
}

// Converts a Bezier patch into a u x v point matrix
MatrixXd genPatchMat(MatrixXd *points, int u, int v) {
    MatrixXd patch = MatrixXd::Zero(u * v, 3);
    float du, dv;
    int count = 0;
    for (int i = 0; i < v; i++) {
        for (int j = 0; j < u; j++) {
            du = j / static_cast<float>(u - 1);
            dv = i / static_cast<float>(v - 1);
            patch.row(count) = p(points, du, dv);
            count++;
        }
    }
    return patch;
}

// Populates a matrix with the normal values at each point
MatrixXd genNormMat(MatrixXd *points, int u, int v) {
    MatrixXd normals = MatrixXd::Zero(u * v, 3);
    float du, dv;
    int count = 0;
    for (int i = 0; i < v; i++) {
        for (int j = 0; j < u; j++) {
            du = 1.0 / static_cast<float>(u - 1);
            dv = 1.0 / static_cast<float>(v - 1);
            normals.row(count) = normal(points, j * du, i * dv).normalized();
            count++;
        }
    }
    return normals;
}

TEST_CASE("testing the dQdu function") {
    MatrixXd points(16, 3);
    Vector3d result;
    result << 6, 0, 0;
    points << 0, 0, 0, 2, 0, 0, 4, 0, 0, 6, 0, 0, 0, 2, 0, 2, 2, 0, 4, 2, 0, 6,
        2, 0, 0, 4, 0, 2, 4, 0, 4, 4, 0, 6, 4, 0, 0, 6, 0, 2, 6, 0, 4, 6, 0, 6,
        6, 0;
    CHECK(dQdu(&points, 0.75, 0.25) == result);
}

// Runs Q' equation
Vector3d qPrime(MatrixXd *points, float u) {
    Vector3d acc;
    acc << 0, 0, 0;
    for (int i = 0; i < 3; i++)
        acc += kChoosei(2, i) * pow(u, i) * pow(1 - u, 2 - i) *
               (3 * (Vector3d)(points->row(i + 1) - points->row(i)));
    return acc;
}

TEST_CASE("testing the qPrime function") {
    MatrixXd points(4, 3);
    Vector3d result;
    result << 6, 0, 0;
    points << 0, 1.5, 0, 2, 1.5, 0, 4, 1.5, 0, 6, 1.5, 0;
    // I guess it makes sense this is the same as the dQdu test case
    CHECK(qPrime(&points, 0.75) == result);
}

// Runs Q'' equation
Vector3d qDoublePrime(MatrixXd *points, float u) {
    Vector3d acc;
    acc << 0, 0, 0;
    for (int i = 0; i < 2; i++)
        acc += pow(u, i) * pow(1 - u, 1 - i) * 6 *
               (Vector3d)(points->row(i + 2) - 2 * points->row(i + 1) +
                          points->row(i));
    return acc;
}

TEST_CASE("testing the qDoublePrime function") {
    MatrixXd points(4, 3);
    Vector3d result;
    result << 0, 0, 0;
    points << 0, 1.5, 0, 2, 1.5, 0, 4, 1.5, 0, 6, 1.5, 0;
    CHECK(qDoublePrime(&points, 0.75) == result);
}

// Actual direct curvature calculation
num curvature(Vector3d prime, Vector3d numPrime) {
    return (prime.cross(numPrime)).norm() / pow(prime.norm(), 3);
}

// Curvature from points including the direction of curvature as a sign
num curvatureFromUSigned(MatrixXd *points, float u, Vector3d normal) {
    Vector3d prime = qPrime(points, u);
    Vector3d numPrime = qDoublePrime(points, u);
    num result = normal.dot(numPrime);
    return sgn(result) * curvature(prime, numPrime);
}

// Curvature from points without directionality
num curvatureFromU(MatrixXd *points, float u) {
    Vector3d prime = qPrime(points, u);
    Vector3d numPrime = qDoublePrime(points, u);
    return curvature(prime, numPrime);
}

MatrixXd hermiteToBezier(MatrixXd *points, MatrixXd *tangents) {
    MatrixXd result = MatrixXd::Zero(4, 3);
    result.row(0) = points->row(0);
    result.row(1) = points->row(0) + (tangents->row(0) / 3.0);
    result.row(2) = points->row(1) - (tangents->row(1) / 3.0);
    result.row(3) = points->row(1);
    return result;
}

num distanceInV(MatrixXd *points, float u, float vstart, float vend) {
    MatrixXd one = p(points, u, vstart);
    MatrixXd two = p(points, u, vend);
    return distance(one, two);
}

num distanceInU(MatrixXd *points, float ustart, float uend, float v) {
    MatrixXd one = p(points, ustart, v);
    MatrixXd two = p(points, uend, v);
    return distance(one, two);
}

TEST_CASE("testing the distanceInU function") {
    MatrixXd points(16, 3);
    points << 0, 0, 0, 2, 0, 0, 4, 0, 0, 6, 0, 0, 0, 2, 0, 2, 2, 0, 4, 2, 0, 6,
        2, 0, 0, 4, 0, 2, 4, 0, 4, 4, 0, 6, 4, 0, 0, 6, 0, 2, 6, 0, 4, 6, 0, 6,
        6, 0;
    CHECK(distanceInU(&points, 0.0, 0.5, 0.5) == 3.0);
}

std::vector<num> fullSplineSignedCurvatureList(MatrixXd *points, MatrixXd
        *tangents, vector<Vector3d> normals, int count) {
    auto curvatures = new vector<num>();
    int steps = count / (points->rows() - 1);
    for (int i = 0; i < points->rows()-1; i++) {
        Vector3d normal = normals[i];
        MatrixXd hermitePoints = points->block(i, 0, 2, 3);
        MatrixXd hermiteTans = tangents->block(i, 0, 2, 3);
        MatrixXd bezierPoints = hermiteToBezier(&hermitePoints, &hermiteTans);
        for (int j = 0; j < steps; j++) {
            num currCurv = curvatureFromUSigned(&bezierPoints, 0.5, normal);
            curvatures->push_back(currCurv);
        }
    }
    return *curvatures;
}

num fullSplineSignedCurvature(MatrixXd *points, MatrixXd *tangents, Vector3d normal) {
    num avgCurv = 0.0;
    for (int i = 0; i < points->rows()-1; i++) {
        MatrixXd hermitePoints = points->block(i, 0, 2, 3);
        MatrixXd hermiteTans = tangents->block(i, 0, 2, 3);
        MatrixXd bezierPoints = hermiteToBezier(&hermitePoints, &hermiteTans);
        num curCurv = curvatureFromUSigned(&bezierPoints, 0.5, normal);
        avgCurv += curCurv;
    }
    avgCurv /= points->rows() - 1;
    return avgCurv;
}

/*
 * points: 16 Bezier patch control points
 * u, v: u and v values
 * inU: a flag to signal which direction we want it run in
 */
num curvatureSign(MatrixXd *points, float u, float v, bool inU) {
    num result;
    Vector3d N, qDP;
    MatrixXd bezPtsU, bezPtsV;
    auto ps = p(points, u, v);
    // Generate 4 Bezier curve points for V
    MatrixXd hermitePts = MatrixXd::Zero(2, 3);
    MatrixXd hermiteTans = MatrixXd::Zero(2, 3);
    hermitePts.row(0) = p(points, u, 0.0);
    hermitePts.row(1) = p(points, u, 1.0);
    hermiteTans.row(0) = dQdv(points, u, 0.0);
    hermiteTans.row(1) = dQdv(points, u, 1.0);
    bezPtsV = hermiteToBezier(&hermitePts, &hermiteTans);
    // Generate 4 Bezier curve points for U
    hermitePts = MatrixXd::Zero(2, 3);
    hermiteTans = MatrixXd::Zero(2, 3);
    hermitePts.row(0) = p(points, 0.0, v); 
    hermitePts.row(1) = p(points, 1.0, v); 
    hermiteTans.row(0) = dQdu(points, 0.0, v);
    hermiteTans.row(1) = dQdu(points, 1.0, v);
    bezPtsU = hermiteToBezier(&hermitePts, &hermiteTans);
    // Calculate N and Q (which should be equivalent to P)
    // num prime in the requested direction
    N = qPrime(&bezPtsU, u).cross(qPrime(&bezPtsV, v));
    if (inU) {
        qDP = qDoublePrime(&bezPtsU, u);
    } else {
        qDP = qDoublePrime(&bezPtsV, v);
    }
    // Get the final dot product result
    result = N.dot(qDP);
    // Return the sign of the result as a floating point
    return static_cast<num>(-1.0 * ((result > 0) - (result < 0)));
}

num curvatureInV(MatrixXd *points, float u, float v) {
    num curvature;
    MatrixXd bezPts;
    MatrixXd hermitePts = MatrixXd::Zero(2, 3);
    MatrixXd hermiteTans = MatrixXd::Zero(2, 3);
    hermitePts.row(0) = p(points, u, 0.0);
    hermitePts.row(1) = p(points, u, 1.0);
    hermiteTans.row(0) = dQdv(points, u, 0.0);
    hermiteTans.row(1) = dQdv(points, u, 1.0);
    bezPts = hermiteToBezier(&hermitePts, &hermiteTans);
    curvature = curvatureFromU(&bezPts, v);
    return curvatureSign(points, u, v, false) * curvature;
}

num curvatureInU(MatrixXd *points, float u, float v) {
    num curvature;
    MatrixXd bezPts;
    MatrixXd hermitePts = MatrixXd::Zero(2, 3);
    MatrixXd hermiteTans = MatrixXd::Zero(2, 3);
    hermitePts.row(0) = p(points, 0.0, v);
    hermitePts.row(1) = p(points, 1.0, v);
    hermiteTans.row(0) = dQdu(points, 0.0, v);
    hermiteTans.row(1) = dQdu(points, 1.0, v);
    bezPts = hermiteToBezier(&hermitePts, &hermiteTans);
    curvature = curvatureFromU(&bezPts, u);
    return curvatureSign(points, u, v, true) * curvature;
}

TEST_CASE("testing the curvatureInU function") {
    MatrixXd points(16, 3);
    points << 0, 0, 0, 2, 0, 0, 4, 0, 0, 6, 0, 0, 0, 2, 0, 2, 2, 0, 4, 2, 0, 6,
        2, 0, 0, 4, 0, 2, 4, 0, 4, 4, 0, 6, 4, 0, 0, 6, 0, 2, 6, 0, 4, 6, 0, 6,
        6, 0;
    CHECK(curvatureInU(&points, 0.5, 0.5) == 0.0);
}
