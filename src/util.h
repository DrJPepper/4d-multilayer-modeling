#pragma once

#include "const.h"
#include "BarycentricOCL/eigen.h"
#include "BarycentricOCL/util.h"

#include "Constraint.h"
#include "Solver.h"
#include "doctest.h"
#include "tqdm.h"

#include <QtCore/QtGlobal>
#include <QApplication>
#include <QWidget>
#include <QWheelEvent>
#include <QSlider>
#include <QVTKOpenGLNativeWidget.h>
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
#include <QSound>
#else
#include <QSoundEffect>
#endif

#include <qsurfaceformat.h>
#include <vtkActor.h>
#include <vtkTriangle.h>
#include <vtkTransform.h>
#include <vtkActor2D.h>
#include <vtkActorCollection.h>
#include <vtkBoxWidget.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkColor.h>
#include <vtkConeSource.h>
#include <vtkCylinderSource.h>
#include <vtkDataSetMapper.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkLight.h>
#include <vtkLineSource.h>
#include <vtkNamedColors.h>
#include <vtkOBJExporter.h>
#include <vtkX3DExporter.h>
#include <vtkOBJImporter.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkProperty2D.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkCubeSource.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkTubeFilter.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vtkGlyph3DMapper.h>
#include <vtkGlyph3D.h>
#include <vtkLine.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkInteractorStyle.h>
#include <vtkBMPWriter.h>
#include <vtkImageWriter.h>
#include <vtkJPEGWriter.h>
#include <vtkPNGWriter.h>
#include <vtkPNMWriter.h>
#include <vtkPostScriptWriter.h>
#include <vtkTIFFWriter.h>
#include <vtkWindowToImageFilter.h>
#include <vtkScalarBarActor.h>
#include <vtkLookupTable.h>
#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkJSONSceneExporter.h>
#include <vtkVRMLExporter.h>

#include <vtkCallbackCommand.h>
#include <vtkCommand.h>
#include <vtkPropPicker.h>
#include <vtkRendererCollection.h>
#include <vtkPointPicker.h>
#include <vtkRenderWindowInteractor.h>

#include <filesystem>
#include <deque>
#include <chrono>
#include <csignal>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/range/combine.hpp>
#include <toml++/toml.hpp>
#include <curl/curl.h>

#include <lsqcpp/lsqcpp.hpp>
#include <nlohmann/json.hpp>

using json = nlohmann::json;
namespace SO = ShapeOp;
namespace po = boost::program_options;

using Pointcloud = Eigen::Matrix<num, 3, Eigen::Dynamic>;
using Vector6 = Eigen::Matrix<num, 6, 1>;
using JacobiMatrix = Eigen::Matrix<num, Eigen::Dynamic, 6>;
using ColVectorX = Eigen::Matrix<num, Eigen::Dynamic, 1>;

extern toml::table setts;
// Temporary
extern ColVector3d transFinal;
extern Matrix3d rotFinal;

Matrix3d rotationMatrix(num x, num y, num z);
Vector3d ratioToRGB(num ratio);
num vectAngle(Vector3d one, Vector3d two);
num addNoise(num original, num noiseMax);
num sgn(num);
json makeEntityGeneric(string type, float red, float green, float blue, MatrixXd position, num radius);
json makeEntityGeneric(string type, float red, float green, float blue, string description, MatrixXd position);
json makeEntityGeneric(string type, float red, float green, float blue, float opacity, string description, MatrixXd position);
json makeEntityGeneric(string type, float red, float green, float blue, float opacity, string description, MatrixXd position, num radius);
vtkSmartPointer<vtkPoints> accessEachVertex(const vtkSmartPointer<vtkPolyData>& mesh, MatrixXd &vertMat);
vtkSmartPointer<vtkPoints> accessEachVertex(const vtkSmartPointer<vtkPolyData>& mesh, MatrixXd &vertMat, Vector3d offset);
vtkSmartPointer<vtkCellArray> accessEachFace(const vtkSmartPointer<vtkPolyData>& mesh);
num angle3Points(Vector3d a, Vector3d b, Vector3d c);
num triangleCheck(Vector3d a, Vector3d b, Vector3d c);
num stdev(vector<num> &arr);
num mean(vector<num> &arr);
num roundDouble(num x, int place);
void writeImage(std::string const& fileName, vtkRenderWindow* renWin, bool rgba = true);
long startTimer();
long stopTimer();
long stopTimer(bool);
long stopTimer(long);
long stopTimer(bool, long);
bool btwn0and1(num x);

template <typename T> bool btwn0and1(T x) {
    int c;
    if (std::is_same<T, Vector2d>::value)
        c = 2;
    else if (std::is_same<T, Vector3d>::value)
        c = 3;
    bool out = true;
    for (int i=0; i<c; i++)
        out &= btwn0and1(x(i));
    return out;
}

struct Callback
{
    Callback(Pointcloud &pointcloudA, Pointcloud &pointcloudB)
        : pointcloudA(&pointcloudA), pointcloudB(&pointcloudB)
    { }

    Pointcloud *pointcloudA = nullptr;
    Pointcloud *pointcloudB = nullptr;

    bool operator()(const lsqcpp::Index iteration,
            const Vector6& xval,
            const ColVectorX&,
            const JacobiMatrix&,
            const Vector6&,
            const Vector6&)
    {
        ColVector3d trans = xval.segment(0, 3);
        Matrix3d rot = lsqcpp::parameter::decodeRotation(xval.segment(3, 3));
        transFinal = trans;
        rotFinal = rot;

        Pointcloud cloud(pointcloudB->rows(), pointcloudB->cols());
        for(lsqcpp::Index i = 0; i < pointcloudB->cols(); ++i)
        {
            cloud.col(i) = rot * pointcloudB->col(i) + trans;
        }

        std::stringstream ss;
        ss << std::setw(3) << std::setfill('0') << iteration << "_pointcloud.b.csv";

        ss.str("");
        ss << std::setw(3) << std::setfill('0') << iteration << "_pointcloud.a.csv";

        return true;
    }
};

struct Objective1
{
    constexpr static bool ComputesJacobian = false;

    Objective1() = default;

    Objective1(Pointcloud &pointcloudA, Pointcloud &pointcloudB)
        : pointcloudA(&pointcloudA), pointcloudB(&pointcloudB)
    { }

    Pointcloud *pointcloudA = nullptr;
    Pointcloud *pointcloudB = nullptr;

    template<typename I, typename O>
        void operator()(const Eigen::MatrixBase<I> &xval,
                Eigen::MatrixBase<O> &fval) const
        {
            ColVector3d translation = xval.segment(0, 3);
            Matrix3d rotation = lsqcpp::parameter::decodeRotation(xval.segment(3, 3));

            fval.derived().resize(pointcloudA->cols());
            for(lsqcpp::Index i = 0; i < pointcloudA->cols(); ++i)
            {
                fval(i) = (pointcloudA->col(i) - (rotation * pointcloudB->col(i) + translation)).norm();
            }
        }
};

struct Objective
{
    constexpr static bool ComputesJacobian = false;

    Objective() = default;

    Objective(Pointcloud &pointcloudA, Pointcloud &pointcloudB)
        : pointcloudA(&pointcloudA), pointcloudB(&pointcloudB)
    { }

    Pointcloud *pointcloudA = nullptr;
    Pointcloud *pointcloudB = nullptr;

    template<typename I, typename O>
        void operator()(const Eigen::MatrixBase<I> &xval,
                Eigen::MatrixBase<O> &fval) const
        {
            ColVector3d translation = xval.segment(0, 3);
            Matrix3d rotation = lsqcpp::parameter::decodeRotation(xval.segment(3, 3));

            fval.derived().resize(pointcloudA->cols());
            for(lsqcpp::Index i = 0; i < pointcloudA->cols(); ++i)
            {
                fval(i) = std::numeric_limits<num>::max();
                for (int j=0; j<pointcloudB->cols(); ++j) {
                    fval(i) = std::min(fval(i), (pointcloudA->col(i) - (rotation * pointcloudB->col(j) + translation)).norm());
                }
            }
        }
};
