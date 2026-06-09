#include "spline.h"

void Spline::defaultValues() {
    this->isUSpline = false;
    this->setSpacing = false;
    this->direction = FORWARD;
    this->shift = false;
    this->uValue = -1.0;
    this->vValue = -1.0;
    this->offset = 0.0;
    this->length = -1.0;
    this->realConstructor = true;
    this->bump = false;
    this->sType = UNSET;
    this->index = -1;
    this->halfSpline = false;
}

Spline::Spline() {
    defaultValues();
    this->realConstructor = false;
}

Spline::Spline(bool isUSpline, num uvValue, num offset) {
    defaultValues();
    this->isUSpline = isUSpline;
    this->uValue = uvValue;
    this->vValue = uvValue;
    this->offset = offset;
}

Spline::Spline(bool isUSpline, bool setSpacing, num uvValue, num offset) {
    defaultValues();
    this->isUSpline = isUSpline;
    this->setSpacing = setSpacing;
    this->uValue = uvValue;
    this->vValue = uvValue;
    this->offset = offset;
}

Spline::Spline(bool isUSpline, Direction direction, num uValue, num vValue, num offset, num length) {
    defaultValues();
    this->isUSpline = isUSpline;
    this->direction = direction;
    this->uValue = uValue;
    this->vValue = vValue;
    this->offset = offset;
    this->length = length;
}

Spline::Spline(bool isUSpline, bool setSpacing, Direction direction, bool shift, num uValue, num vValue, num offset, num length) {
    defaultValues();
    this->isUSpline = isUSpline;
    this->setSpacing = setSpacing;
    this->direction = direction;
    this->shift = shift;
    this->uValue = uValue;
    this->vValue = vValue;
    this->offset = offset;
    this->length = length;
}

Spline Spline::duplicate() {
    auto s = Spline(isUSpline, setSpacing, direction, shift, uValue, vValue, offset, length);
    s.bump = bump;
    s.sType = sType;
    s.index = index;
    return s;
}

void Spline::combine(Spline &splineF, Spline &splineR) {
    realPoints = MatrixXd::Zero(splineF.realPoints.rows()+splineR.realPoints.rows()-1, 3);
    realPoints << splineR.realPoints.colwise().reverse(), splineF.realPoints(Eigen::seq(1,last),all);
    uvPoints = MatrixXd::Zero(splineF.uvPoints.rows()+splineR.uvPoints.rows()-1, 3);
    uvPoints << splineR.uvPoints.colwise().reverse(), splineF.uvPoints(Eigen::seq(1,last),all);
    tangents = MatrixXd::Zero(splineR.tangents.rows()+splineF.tangents.rows()-1, 3);
    tangents << splineR.tangents.colwise().reverse(), splineF.tangents(Eigen::seq(1,last),all);
    splinePoints = MatrixXd::Zero(splineR.splinePoints.rows()+splineF.splinePoints.rows()-1, 3);
    splinePoints << splineR.splinePoints.colwise().reverse(), splineF.splinePoints(Eigen::seq(1,last),all);
    splinePointsUV = MatrixXd::Zero(splineR.splinePointsUV.rows()+splineF.splinePointsUV.rows()-1, 2);
    splinePointsUV << splineR.splinePointsUV.colwise().reverse(), splineF.splinePointsUV(Eigen::seq(1,last),all);
    // NOTE: This is currently being overridden by different code within the
    // populateSpline method.
    if (!shift) {
        Vector2d pp;
        pp << uValue, vValue;
        cellCentersCart = MatrixXd::Zero(splineR.cellCentersCart.rows()+splineF.cellCentersCart.rows()+1, 3);
        cellCentersCart << splineR.cellCentersCart.colwise().reverse(), splinePoints.row(0), splineF.cellCentersCart;
        cellCentersUV = MatrixXd::Zero(splineR.cellCentersUV.rows()+splineF.cellCentersUV.rows()+1, 2);
        cellCentersUV << splineR.cellCentersUV.colwise().reverse(), pp, splineF.cellCentersUV;
    } else {
        cellCentersCart = MatrixXd::Zero(splineR.cellCentersCart.rows()+splineF.cellCentersCart.rows(), 3);
        cellCentersCart << splineR.cellCentersCart.colwise().reverse(), splineF.cellCentersCart;
        cellCentersUV = MatrixXd::Zero(splineR.cellCentersUV.rows()+splineF.cellCentersUV.rows(), 2);
        cellCentersUV << splineR.cellCentersUV.colwise().reverse(), splineF.cellCentersUV;
    }
}


MatrixXd genSpline(MatrixXd &points, MatrixXd &tangents) {
    return genSpline(points, tangents, 50);
}

// Takes the points that need intersected and returns the stitched spline
MatrixXd genSpline(MatrixXd &points, MatrixXd &tangents, int n) {
    float tension = 0.0;
    int shrinkage;

    // Run primary functions
    points = initialize(points, &tangents, tension, shrinkage);
    MatrixXd curve = MatrixXd::Zero((points.rows() - 1) * (n + 1), 3);
    for (int i = 0; i < points.rows() - 1; i++) {
        MatrixXd temp = hermiteToBezier(&points, &tangents, i);
        curve.block(i * (n + 1), 0, 1 + n, 3) = genCurveMat(&temp, n);
    }
    return curve;
}

// Blends the controls points together from one curve to the next
void smoothPoints(MatrixXd &points, int index) {
    Vector3d zero(static_cast<num>(0), static_cast<num>(0), static_cast<num>(0));
    Vector3d p0, p1 = zero, p2 = zero;
    num t, d;
    int counter = index - 1;
    if (index > 1 && index < points.rows()) {
        p0 = points.row(index);
        if (p0 != zero) {
            while (p1 == zero && counter >= 0 && (index - counter) < 7) {
                p1 = points.row(counter);
                counter--;
            }
            while (p2 == zero && counter >= 0 && (index - counter) < 7) {
                p2 = points.row(counter);
                counter--;
            }
            if (p1 != zero && p2 != zero) {
                t = -(p1 - p0).dot(p2 - p1) / pow2((p2 - p1).norm());
                d = sqrt(pow2((p1(0) - p0(0)) + (p2(0) - p1(0)) * t) +
                        pow2((p1(1) - p0(1)) + (p2(1) - p1(1)) * t) +
                        pow2((p1(2) - p0(2)) + (p2(2) - p1(2)) * t));
                if (d >= 1.0) {
                    points.row(index) = zero;
                }
            }
        }
    }
}

// Populates matrix of points
MatrixXd initialize(MatrixXd &points, MatrixXd *tangents, float tension,
        int &shrinkage) {
    // Store the last tangent in a temp file and first into our real matrix
    MatrixXd zero = MatrixXd::Zero(1, 3);
    shrinkage = 0;
    Vector3d vec, startRow, endRow;
    MatrixXd pointsNew = points;
    int counter = 0, start, end, div;
    bool tempB;
    while (false) {
        if (pointsNew.rows() > 4)
            smoothPoints(pointsNew, counter);
        tempB = pointsNew.row(counter) == zero;
        if (tempB) {
            start = counter;
            while (tempB && counter < pointsNew.rows() - 1) {
                counter++;
                if (pointsNew.rows() > 4)
                    smoothPoints(pointsNew, counter);
                tempB = pointsNew.row(counter) == zero;
            }
            if (counter != pointsNew.rows() - 1) {
                counter--;
            }
            end = counter;
            if (start == 0 && end == pointsNew.rows() - 1) {
                shrinkage = points.rows();
                return zero;
            } else if (start == 0) {
                shrinkage += end + 1;
                counter = -1;
                MatrixXd tempTemp =
                    pointsNew(Eigen::seq(end + 1, pointsNew.rows() - 1), all);
                pointsNew = tempTemp;
            } else if (end == (pointsNew.rows() - 1)) {
                shrinkage += end - start + 1;
                MatrixXd tempTemp2 = pointsNew(Eigen::seq(0, start - 1), all);
                pointsNew = tempTemp2;
                counter = pointsNew.rows();
            } else {
                startRow = pointsNew.row(start - 1);
                endRow = pointsNew.row(end + 1);
                vec = endRow - startRow;
                div = end - start + 2;
                for (int q = 1; q <= end - start + 1; q++) {
                    pointsNew.row(start + q - 1) = startRow + vec / div * q;
                }
                counter = end;
            }
        }
        counter++;
    }
    points = pointsNew;
    MatrixXd lastTangent = MatrixXd::Zero(1, 3);
    tangents->conservativeResize(tangents->rows() + 1, tangents->cols());
    tangents->row(tangents->rows() - 1) = points.row(1) - points.row(0);
    tangents->row(tangents->rows() - 1) *= (1 - tension);
    lastTangent = points.row(points.rows() - 1) - points.row(points.rows() - 2);
    lastTangent *= (1 - tension);
    points = pointsNew;

    // TODO: I may be generating an extra point and tan here, but it's working
    // correctly with graddesc so I'm not going to mess with it further for now
    MatrixXd temp = MatrixXd::Zero(2, 3);
    // Calculate tangents
    for (int i = 1; i < points.rows() - 1; i++) {
        tangents->conservativeResize(tangents->rows() + 1, tangents->cols());
        temp.row(0) = points.row(i - 1);
        temp.row(1) = points.row(i + 1);
        tangents->row(tangents->rows() - 1) = tangent(&temp, tension);
    }

    // Move over the final tangent
    tangents->conservativeResize(tangents->rows() + 1, tangents->cols());
    tangents->row(tangents->rows() - 1) = lastTangent;
    return pointsNew;
}

MatrixXd tangent(MatrixXd *inTans, float tension) {
    return (inTans->row(1) - inTans->row(0)) / 2.0 * (1 - tension);
}

MatrixXd hermiteToBezier(MatrixXd *points, MatrixXd *tangents, int row1) {
    MatrixXd result = MatrixXd::Zero(4, 3);
    result.row(0) = points->row(row1);
    result.row(1) = points->row(row1) + (tangents->row(row1) / 3.0);
    result.row(2) = points->row(row1 + 1) - (tangents->row(row1 + 1) / 3.0);
    result.row(3) = points->row(row1 + 1);
    return result;
}

num tempCounter = 0;

// Populates a matrix with values from the Q calculations
MatrixXd genCurveMat(MatrixXd *points, int n) {
    MatrixXd curve = MatrixXd::Zero(n + 1, 3);
    float du = 1.0 / (float)n;
    for (int i = 0; i <= n; i++) {
        curve.row(i) = Q(points, i * du);
    }
    return curve;
}
