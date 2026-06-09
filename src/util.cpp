#include "util.h"

long _startTime = 0;

Vector3d ratioToRGB(num ratio) {
    num cap = 0.2, green, other;
    num minCap = 1.0 - cap;
    num maxCap = 1.0 + cap;
    Vector3d result = Vector3d::Zero(3);

    if (ratio > maxCap)
        result << 255.0, 0.0, 0.0;
    else if (ratio < minCap)
        result << 0.0, 0.0, 255.0;
    else if (ratio > 1.0) {
        other = (ratio - 1.0) / cap * 255.0;
        green = 255 - other;
        result << other, green, 0.0;
    } else if (ratio <= 1.0) {
        other = (1.0 - ratio) / cap * 255.0;
        green = 255 - other;
        result << 0.0, green, other;
    } else {
        cerr << "WARNING: NaN in ratioToRGB" << endl;
    }

    return result / 255.0;
}

num distance2d(Vector2d one, Vector2d two) {
    MatrixXd oneMat(1, 3);
    oneMat << one(0), one(1), 0.0;
    MatrixXd twoMat(1, 3);
    twoMat << two(0), two(1), 0.0;
    return distance(oneMat, twoMat);
}

num vectAngle(Vector3d one, Vector3d two) {
    return acos((one.dot(two)) / (one.norm() * two.norm()));
}

num addNoise(num original, num noiseMax) {
    num lower = original * (1.0 - noiseMax);
    num upper = original * (1.0 + noiseMax);
    std::random_device r;
    std::uniform_real_distribution<num> unif(lower, upper);
    std::default_random_engine re(r());
    return unif(re);
}

json makeEntityGeneric(string type, float red, float green, float blue, MatrixXd position, num radius) {
    return makeEntityGeneric(type, red, green, blue, 1.0, "", position, radius);
}

json makeEntityGeneric(string type, float red, float green, float blue, string description, MatrixXd position) {
    return makeEntityGeneric(type, red, green, blue, 1.0, description, position, -1.0);
}

json makeEntityGeneric(string type, float red, float green, float blue, float opacity, string description, MatrixXd position) {
    return makeEntityGeneric(type, red, green, blue, opacity, description, position, -1.0);
}

json makeEntityGeneric(string type, float red, float green, float blue, float opacity, string description, MatrixXd position, num radius) {
    num a = 4, b = 2;
    bool pline = !type.compare("polyline") || !type.compare("y");
    if (pline)
        a = 7;
    MatrixXd pos(position.rows(), position.cols());
    for (int i=0; i<pos.rows(); i++) {
        for (int j=0; j<pos.cols(); j++) {
            pos(i,j) = roundDouble(position(i,j), a);
        }
    }
    json entity;
    entity["t"] = type;
    if (red != 1.0 || green != 1.0 || blue != 1.0)
        entity["c"] = {roundDouble(red, b), roundDouble(green, b), roundDouble(blue, b)};
    if (opacity != 1.0)
        entity["o"] = roundDouble(opacity, b);
    if (description.compare(""))
        entity["d"] = description;
    if (radius > 0.0)
        entity["r"] = roundDouble(radius, b);
    if (!type.compare("p") || !type.compare("position")) {
        entity["p"] = {pos(0), pos(1), pos(2)};
    } else if (!type.compare("vector") || !type.compare("line") || !type.compare("v")) {
        entity["p"] = {pos(0, 0), pos(0, 1), pos(0, 2), pos(1, 0), pos(1, 1), pos(1, 2)};
    } else if (pline) {
        for (auto row : pos.rowwise()) {
            entity["p"].push_back({row(0), row(1), row(2)});
        }
    } else {
        entity = {"Errored", fmt::format("Invalid type {}", type)};
    }

    return entity;
}

vtkSmartPointer<vtkPoints> accessEachVertex(const vtkSmartPointer<vtkPolyData>& mesh, MatrixXd &vertMat) {
    Vector3d offset;
    offset << 0, 0, 0;
    return accessEachVertex(mesh, vertMat, offset);
}

vtkSmartPointer<vtkPoints> accessEachVertex(const vtkSmartPointer<vtkPolyData>& mesh, MatrixXd &vertMat, Vector3d offset) {
    vtkSmartPointer<vtkPoints> vertices = mesh->GetPoints();
    vtkSmartPointer<vtkDataArray> verticesArray = vertices->GetData();

    int numberOfVertices = vertices->GetNumberOfPoints();
    vertMat.resize(numberOfVertices, 3);
    // Access 3D coordinate [x, y, z] of each vertex 
    vtkNew<vtkPoints> points;
    for( int i = 0; i < numberOfVertices; i++ )
    {
        float x = verticesArray->GetComponent(i, 0) + offset(0);
        float y = verticesArray->GetComponent(i, 1) + offset(1);
        float z = verticesArray->GetComponent(i, 2) + offset(2);
        points->InsertNextPoint(x, y, z);

        vertMat.row(i) << x, y, z;
    }

    return points;
}

vtkSmartPointer<vtkCellArray> accessEachFace(const vtkSmartPointer<vtkPolyData>& mesh) {
    int numberOfFaces = mesh->GetNumberOfCells();
    vtkNew<vtkCellArray> lines;

    // Access each mesh face. A face is defined by the indices 
    // of the participating vertices.
    // this mesh has triangle faces - therefore three vertices.
    for( int i = 0; i < numberOfFaces; i++)
    {
        vtkSmartPointer<vtkIdList> face = vtkSmartPointer<vtkIdList>::New();
        mesh->GetCellPoints(i,face);

        for (int j=0; j<3; j++) {
            vtkNew<vtkLine> line;
            line->GetPointIds()->SetId(0, face->GetId(j));
            line->GetPointIds()->SetId(1, face->GetId((j+1)%3));
            lines->InsertNextCell(line);
        }
    }
    return lines;
}

num angle3Points(Vector3d a, Vector3d b, Vector3d c) {
    Vector3d v1 = a - b;
    Vector3d v2 = c - b;
    return acos(v1.dot(v2) / (v1.norm() * v2.norm()));
}

num triangleCheck(Vector3d a, Vector3d b, Vector3d c) {
    num ang1 = angle3Points(b, a, c);
    num ang2 = angle3Points(b, c, a);
    return abs(ang1 - ang2);
}

num mean(vector<num> &arr) {
    num out = 0.0;
    for (num d : arr)
        out += d;
    return out / arr.size();
}

num roundDouble(num x, int place) {
    // avoiding the pow function
    num mult = 1.0;
    for (int i=0; i<place; i++)
        mult *= 10.0;
    return static_cast<num>(round(x * mult)) / mult;
}

num stdev(vector<num> &arr) {
    num sum = 0.0, std = 0.0, mean;

    for (num a : arr)
        sum += a;

    mean = sum / arr.size();

    for (num a : arr)
        std += pow2(a - mean);

    return sqrt(std / arr.size());
}

num sgn(num x) {
    return static_cast<num>(-1.0 * ((x > 0) - (x < 0)));
}

// From VTK Demos
void writeImage(std::string const& fileName, vtkRenderWindow* renWin, bool rgba) {
    if (!fileName.empty()) {
        std::string fn = fileName;
        std::string ext;
        auto found = fn.find_last_of(".");
        if (found == std::string::npos) {
            ext = ".png";
            fn += ext;
        }
        else {
            ext = fileName.substr(found, fileName.size());
        }
        std::locale loc;
        std::transform(ext.begin(), ext.end(), ext.begin(),
                [=](char const& c) { return std::tolower(c, loc); });
        auto writer = vtkSmartPointer<vtkImageWriter>::New();
        if (ext == ".bmp") {
            writer = vtkSmartPointer<vtkBMPWriter>::New();
        } else if (ext == ".jpg") {
            writer = vtkSmartPointer<vtkJPEGWriter>::New();
        } else if (ext == ".pnm") {
            writer = vtkSmartPointer<vtkPNMWriter>::New();
        } else if (ext == ".ps") {
            if (rgba) {
                rgba = false;
            }
            writer = vtkSmartPointer<vtkPostScriptWriter>::New();
        } else if (ext == ".tiff") {
            writer = vtkSmartPointer<vtkTIFFWriter>::New();
        } else {
            writer = vtkSmartPointer<vtkPNGWriter>::New();
        }
        vtkNew<vtkWindowToImageFilter> window_to_image_filter;
        window_to_image_filter->SetInput(renWin);
        window_to_image_filter->SetScale(1); // image quality
        if (rgba) {
            window_to_image_filter->SetInputBufferTypeToRGBA();
        } else {
            window_to_image_filter->SetInputBufferTypeToRGB();
        }
        // Read from the front buffer.
        window_to_image_filter->ReadFrontBufferOff();
        window_to_image_filter->Update();

        writer->SetFileName(fn.c_str());
        writer->SetInputConnection(window_to_image_filter->GetOutputPort());
        writer->Write();
    } else {
        cerr << "No filename provided." << std::endl;
    }
}

bool btwn0and1(num x) {
    return x >= 0 && x <= 1;
}

long startTimer() {
    _startTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch()).count();
    return _startTime;
}

long stopTimer() {
    return stopTimer(true, -1);
}

long stopTimer(bool printTime) {
    return stopTimer(printTime, -1);
}

long stopTimer(long start) {
    return stopTimer(true, start);
}

long stopTimer(bool printTime, long start) {
    long endTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch()).count();
    if (start < 0)
        start = _startTime;
    long diff = endTime - start;
    if (printTime)
        cout << "Timer Duration: " << diff << endl;
    return diff;
}

Matrix3d rotationMatrix(num x, num y, num z) {
    Matrix3d *m = new Matrix3d();
    *m << cos(z)*cos(y), cos(z)*sin(y)*sin(x)-sin(z)*cos(x), cos(z)*sin(y)*cos(x)+sin(z)*sin(x),
         sin(z)*cos(y), sin(z)*sin(y)*sin(x)+cos(z)*cos(x), sin(z)*sin(y)*cos(x)-cos(z)*sin(x),
         -sin(y), cos(y)*sin(x), cos(y)*cos(x);
    return *m;
}
