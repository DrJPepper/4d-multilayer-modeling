#include "widget.h"
#include "ui_widget.h"

// Being lazy just going to use a global variable instead
// of making a mess of the functions
vtkSmartPointer<vtkRenderer> currentRenderer;

Widget::Widget(QWidget *parent, Simulation *simulation, string inputFileName) : QWidget(parent), ui(new Ui::Widget), model(simulation->getModel()) {
    scrollBehavior = zoom;
    ui->setupUi(this);
    pastRenderModel = false;

    dark.background = "Black";
    dark.primaryPoints = "White";
    dark.primarySpline = "Aqua";
    dark.cellPointsAll = "GreenYellow";
    dark.cellPointsUsing = "LightSalmon";
    dark.cellSplines = "Orange";
    dark.hexVertices = "Red";
    dark.vSplines = "Fuchsia";

    dark.paramBackground = "Black";
    dark.paramModelPoints = "red";
    dark.paramModelFaces = "Green";
    dark.paramPointsPrimary = "White";
    dark.paramPrimarySpline = "Aqua";
    dark.paramControlPoints = "Orange";
    dark.paramSplines = "Fuchsia";

    dark.scaleBarFont = "white";

    dark.triBackground = "Black";
    dark.tris = "Cornsilk";


    light.background = "White";
    light.primaryPoints = "silver";
    light.primarySpline = "slateblue";
    light.cellPointsAll = "tomato";
    light.cellPointsUsing = "red";
    light.cellSplines = "fuchsia";
    light.hexVertices = "greenyellow";
    light.vSplines = "darkviolet";

    light.paramBackground = light.background;
    light.paramModelPoints = "red";
    light.paramModelFaces = "green";
    light.paramPointsPrimary = light.primaryPoints;
    light.paramPrimarySpline = "aqua";
    light.paramControlPoints = light.primaryPoints;
    light.paramSplines = "red";

    light.scaleBarFont = "black";

    light.triBackground = light.background;
    light.tris = "gray";
    
    //theme = ui->lightModeBox->isChecked() ? light : dark;
    theme = light;
    gridHasBeenRendered = false;

    colors = vtkSmartPointer<vtkNamedColors>::New();
    scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
    renderer = vtkSmartPointer<vtkRenderer>::New();
    rendererTri = vtkSmartPointer<vtkRenderer>::New();
    currentRenderer = renderer;

    renderer->SetBackground(colors->GetColor3d(theme.background).GetData());
    rendererTri->SetBackground(colors->GetColor3d(theme.triBackground).GetData());

    vtkWindow = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    vtkWindowTri = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    vtkWindow->AddRenderer(renderer);
    vtkWindowTri->AddRenderer(rendererTri);
    ui->qvtkWidget->setRenderWindow(vtkWindow);

    renderWindowInteractor = vtkWindow->GetInteractor();
    defaultScrollValue = dynamic_cast<vtkInteractorStyle*>(
      renderWindowInteractor->GetInteractorStyle())
      ->GetMouseWheelMotionFactor();

    exporter = vtkSmartPointer<vtkJSONSceneExporter>::New();
    exporter->SetActiveRenderer(renderer);
    exporter->SetRenderWindow(vtkWindow);

    if (!inputFileName.compare("")) 
        this->inputFileName = "../models/result.obj";
    else
        this->inputFileName = inputFileName;

    sim = simulation;
    grid = &sim->getTissue().getGrid();
    grid->setWidget(this);

    startTime = time(0);
    auto name = fmt::format("../exports/export_{}", startTime);
    int epochs = setts["simulation"]["epochs"].value_or(EPOCHS_PER_CYCLE);
    if (model.epochsInStep > 0)
        epochs = model.epochsInStep;
    ui->epochTextEdit->setText(QString::fromStdString(fmt::format("{}", epochs)));
    ui->jsonTextEdit->setText(QString::fromStdString(name));
    ui->objTextEdit->setText(QString::fromStdString(name));
    ui->pngTextEdit->setText(QString::fromStdString(name));
    ui->linearKEdit->setText(QString::fromStdString(fmt::format("{}", sim->getTissue().getGrid().kEdge)));
    ui->angularKEdit->setText(QString::fromStdString(fmt::format("{}", sim->getTissue().getGrid().kAng)));
    ui->angular90KEdit->setText(QString::fromStdString(fmt::format("{}", sim->getTissue().getGrid().kAng90)));
    ui->cutoffTextEdit->setText(QString::fromStdString(fmt::format("{}", sim->cutoff)));
    ui->stepSizeTextEdit->setText(QString::fromStdString(fmt::format("{}", sim->divAvg)));

    vtkNew<vtkAxesActor> axes;
    axesWidget = vtkSmartPointer<vtkOrientationMarkerWidget>::New();;
    num rgba[4]{0.0, 0.0, 0.0, 0.0};
    colors->GetColor("Black", rgba);
    axesWidget->SetOutlineColor(rgba[0], rgba[1], rgba[2]);
    axesWidget->SetOrientationMarker(axes);
    axesWidget->SetInteractor(renderWindowInteractor);
    axesWidget->SetViewport(0.0, model.batch ? 0.0333 : 0.0, 0.0784, model.batch ? 0.13: 0.097);
    axesWidget->EnabledOn();
    axesWidget->InteractiveOn();

    this->connect(sim->s.get(), &SO::Solver::progress, this, &Widget::updatePBar);
    this->connect(sim, &Simulation::renderSignal, this, &Widget::renderGridAndSave, Qt::BlockingQueuedConnection);

    setKs = true;
    runningSO = false;
    triangulated = false;
    firstRenderGrid = true;
    scaleFactor = model.scaleFactor;
    cubic = !static_cast<string>(setts["grid"]["type"].value_or("epi")).compare("cubic");
    runAgents = true;
    if (cubic) {
        runAgents = false;
        finishInitCubic();
        funDeq.push_back(&Widget::renderGrid);
        funDeq.push_back(&Widget::runSO);
        funDeq.push_back(&Widget::duplicateLayerCubic);
    }

    bool imported = setts["global"]["import_from_json"].value_or(false);
    if (!setts["global"]["scene_saving"].value_or(false))
        ui->visPrintButton->setVisible(false);
    if (imported) {
        ui->qvtkWidgetParam->setVisible(false);
        scaleFactor = sim->getTissue().scaleFactor;
    }
    if (model.batch) {
        std::filesystem::create_directory(fmt::format("../exports/batch_{}", startTime));
        csvFile.open(fmt::format("../exports/batch_{}/results.csv", startTime));
        if (!csvFile.good()) {
            cout << "ERROR: CSV output file cannot be created\n";
            exit(1);
        }
        csvFile << "kLinear,k90,k120,Error\n";
        csvFile.flush();
        renderModel();
        auto actors = renderer->GetActors();
        actors->InitTraversal();
        while (actors->GetNumberOfItems())
            renderer->RemoveActor(actors->GetNextItem());
        firstBatchRun = true;
        timer = new QTimer(this);
        timer->singleShot(1000, this, &Widget::runBatch);
    } else if (!cubic) {
        if (setts["global"]["scene_saving"].value_or(false)) {
            sceneVis["list"].push_back({});
            sceneVis["reset"] = true;
            sceneVis["glyph"] = true;
        }
        if (imported) {
            renderGrid();
        } else {
            if (setts["global"]["scene_saving"].value_or(false)) {
                sceneVis["list"].push_back({});
                sceneVis["reset"] = true;
                sceneVis["glyph"] = true;
            }
            vtkNew<vtkCallbackCommand> clickCallback;
            clickCallback->SetCallback(clickCallbackFunction);
            clickCallback->SetClientData(this);
            renderWindowInteractor->AddObserver(vtkCommand::LeftButtonPressEvent,
                    clickCallback);
            funDeq.push_back(&Widget::renderCellMap);
            funDeq.push_back(&Widget::renderGrid);
            renderModel();
        }
    }
}

void Widget::runBatch() {
    int waitTime = 0;
    if (firstBatchRun) {
        batchIndex = 0;
        batchTotal = static_cast<int>((model.kLinearInMax - model.kLinearInMin) / model.kLinearStep);
        batchTotal += static_cast<int>((model.kAngIn120Max - model.kAngIn120Min) / model.kAng120Step);
        batchTotal += static_cast<int>((model.kAngIn90Max - model.kAngIn90Min) / model.kAng90Step);
        ui->overallPBar->setMinimum(0);
        ui->overallPBar->setMaximum(batchTotal);
        ui->overallPBar->setValue(0);
        firstBatchRun = false;
        firstWrite = true;
        kLinear = model.kLinearInMin;
        kA120 = model.kAngIn120Min;
        kA90 = model.kAngIn90Min;
        counter = 0;
        grid->epochsRun = 0;
        outStringBatch = fmt::format("../exports/batch_{}/{:04d}_k90-{:.2f}_k120-{:.2f}_kL-{:.2f}",
                startTime, counter, kA90, kA120, kLinear);
        pausing = false;
        ui->linearKEdit->setReadOnly(true);
        ui->angularKEdit->setReadOnly(true);
        ui->angular90KEdit->setReadOnly(true);
        ui->epochTextEdit->setReadOnly(true);

        if (cubic) {
            renderGrid();
            sim->runSO(boost::lexical_cast<int>(ui->epochTextEdit->text().toStdString()), false);
            duplicateLayerCubic();
            grid->epochsRun = 0;
            for (auto v : grid->getVertices())
                v->originalPosition << v->getPosition();
        }
        renderGrid();
        currentRenderer->ResetCamera();
        currentRenderer->GetRenderWindow()->Render();
        std::filesystem::create_directory(outStringBatch);
        waitTime = 1000;
    } else {
        ui->overallPBar->setValue(batchIndex);
    }

    ui->infoBox->setPlainText(QString::fromStdString(fmt::format("Epochs done: {}", grid->epochsRun)));
    if (grid->epochsRun >= model.epochsIn && !pausing) {
        pausing = true;
        waitTime = round(boost::lexical_cast<num>(ui->epochTimer->text().toStdString())
            * model.epochsInStep * 1000);
        waitTime = waitTime > 3000 ? 3000 : waitTime;
        csvFile << fmt::format("{},{},{},{}\n", kLinear, kA90, kA120, currErr);
        csvFile.flush();
        batchIndex++;
        timer->singleShot(waitTime, this, &Widget::runBatch);
    } else {
        if (pausing) {
            if (model.kLinearInMax < 0)
                return;
            pausing = false;
            kA90 += model.kAng90Step;
            if (kA90 >= model.kAngIn90Max) {
                kA90 = model.kAngIn90;
                kA120 += model.kAng120Step;
                if (kA120 > model.kAngIn120Max * 1.001) {
                    kA120 = model.kAngIn120;
                    kLinear += model.kLinearStep;
                }
            }
            if (model.kLinearInMax > 0 && kLinear > model.kLinearInMax * 1.001) {
                csvFile.close();
                return;
            }
            grid->epochsRun = 0;
            ui->angular90KEdit->setText(QString::fromStdString(fmt::format("{}", kA90)));
            ui->angularKEdit->setText(QString::fromStdString(fmt::format("{}", kA120)));
            ui->linearKEdit->setText(QString::fromStdString(fmt::format("{}", kLinear)));
            ui->infoBox->setPlainText(QString::fromStdString(fmt::format("Epochs done: {}", grid->epochsRun)));
            outStringBatch = fmt::format("../exports/batch_{}/{:04d}_k90-{:.2f}_k120-{:.2f}_kL-{:.2f}",
                    startTime, ++counter, kA90, kA120, kLinear);
            std::filesystem::create_directory(outStringBatch);
            grid->reset();
            renderGrid();
        }
        cout << "epochsRun: " << grid->epochsRun << endl;
        grid->kAng = kA120;
        grid->kAng90 = kA90;
        grid->kEdge = kLinear;
        timer->singleShot(waitTime, this, &Widget::runSO);
    }
}

void Worker::runSO(Widget *widget, Simulation *sim, int epochs, int genAgentTD) {
    try {
        sim->runSO(epochs, genAgentTD);
    } catch (std::exception& e) {
        cout << "Error during runSO:\n" << e.what() << endl;
    }
    emit finishedSO();
}

void Widget::evaluateSlotSingle() {
    evaluateChamfer(false);
}

void Widget::evaluateSlotBatch() {
    evaluateChamfer(false);
    runBatch();
}

void Widget::playBell() {
    // https://assets.mixkit.co/active_storage/sfx/586/586.wav
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    QSound::play("../lib/586.wav");
#else
    auto qs = QSoundEffect();
    qs.setSource(QUrl("../lib/586.wav"));
    qs.play();
#endif
}

void Widget::runSO() {
    if (!runningSO) {
        ui->shapeOpRunning->setText(QString::fromStdString("ShapeOp Running"));
        ui->shapeOpRunning->setText(QString::fromStdString("ShapeOp Running"));
        ui->shapeOpRunning->setStyleSheet("QLabel { background-color : red; color : white; }");
        sim->stepSize = boost::lexical_cast<num>(ui->stepSizeTextEdit->text().toStdString());
        sim->cutoff = boost::lexical_cast<num>(ui->cutoffTextEdit->text().toStdString());
        runningSO = true;
        Worker *worker = new Worker;
        QThread *workerThread = new QThread;
        worker->moveToThread(workerThread);
        auto func = &Widget::renderGrid;
        int runAgent = static_cast<int>(ui->trainingDataRadio->isChecked()) +
            static_cast<int>(ui->runAgentRadio->isChecked())*2;
        if (runAgent == 1)
            func = &Widget::renderTD;
        else if (runAgent == 2)
            func = &Widget::renderAgency;
        this->connect(worker, &Worker::finishedSO, this, func);
        this->connect(worker, &Worker::finishedSO, workerThread, &QThread::quit);
        if (model.batch) {
            this->connect(this, &Widget::operateSO, worker, &Worker::runSO);
            this->connect(worker, &Worker::finishedSO, this, &Widget::evaluateSlotBatch);
        } else {
            this->connect(this, &Widget::operateSO, worker, &Worker::runSO);
        }
        if (ui->chimeCheckBox->isChecked()) {
            this->connect(worker, &Worker::finishedSO, this, &Widget::playBell);
        }
        workerThread->start();
        workerThread->setPriority(QThread::TimeCriticalPriority);
        emit this->operateSO(this, sim, boost::lexical_cast<int>(ui->epochTextEdit->text().toStdString()), runAgent);
    }
}

vtkSmartPointer<vtkActor> Widget::renderSpheres(vtkSmartPointer<vtkPoints> points, string color, num radius) {
    return renderSpheres(points, color, radius, 1.0);
}

vtkSmartPointer<vtkActor> Widget::renderSpheres(vtkSmartPointer<vtkPoints> points, string color, num radius, num opacity) {
    auto sphereSource = vtkSmartPointer<vtkSphereSource>::New();
    sphereSource->SetRadius(radius);

    auto pd = vtkSmartPointer<vtkPolyData>::New();
    pd->SetPoints(points);

    auto mapper = vtkSmartPointer<vtkGlyph3DMapper>::New();
    mapper->SetInputData(pd);
    mapper->SetSourceConnection(sphereSource->GetOutputPort());
    mapper->ScalarVisibilityOff();
    mapper->ScalingOff();

    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    auto gc = colors->GetColor3d(color);
    auto colorVals = gc.GetData();
    actor->GetProperty()->SetOpacity(opacity);
    currentRenderer->AddActor(actor);
    actor->GetProperty()->SetColor(colorVals);
    if (setts["global"]["scene_saving"].value_or(false)) {
        num x[3];
        for (int i=0; i<points->GetNumberOfPoints(); i++) {
            points->GetPoint(i, x);
            MatrixXd loc(1, 3);
            loc << x[0], x[1], x[2];
            if (currentRenderer == renderer) {
                json pt = makeEntityGeneric("p", colorVals[0], colorVals[1], colorVals[2], loc, radius);
                sceneVis["list"][sceneVis["list"].size()-1]["e"].push_back(pt);
            }
        }
    }
    return actor;
}

vtkSmartPointer<vtkPoints> Widget::renderSpheres(MatrixXd &points, string color, num radius, vtkSmartPointer<vtkActor> &dataOut) {
    return renderSpheres(points, color, radius, 1.0, dataOut);
}

vtkSmartPointer<vtkPoints> Widget::renderSpheres(MatrixXd &points, string color, num radius, num opacity, vtkSmartPointer<vtkActor> &dataOut) {
    vtkNew<vtkPoints> vtkPoints;
    for (auto point : points.rowwise())
        vtkPoints->InsertNextPoint(point(0), point(1), point(2));
    dataOut = renderSpheres(vtkPoints, color, radius, opacity);
    return vtkPoints;
}

vtkSmartPointer<vtkActor> Widget::renderLines(vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkCellArray> lines, vtkSmartPointer<vtkFloatArray> colors, num radius, num opacity) {
    vtkNew<vtkPolyData> linesPolyData;
    linesPolyData->SetPoints(points);
    linesPolyData->SetLines(lines);
    linesPolyData->GetCellData()->SetScalars(colors);
    vtkNew<vtkTubeFilter> tubeFilter;
    tubeFilter->SetInputData(linesPolyData);
    tubeFilter->SetNumberOfSides(8);
    tubeFilter->SetRadius(radius);
    tubeFilter->Update();
    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputConnection(tubeFilter->GetOutputPort());
    mapper->SetColorMode(2);

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetLineWidth(4);

    renderer->AddActor(actor);
    if (setts["global"]["scene_saving"].value_or(false)) {
        num x1[3], x2[3], colorVals[3];
        for (int i=0; i<points->GetNumberOfPoints(); i+=2) {
            points->GetPoint(i, x1);
            points->GetPoint(i+1, x2);
            MatrixXd loc(2, 3);
            loc << x1[0], x1[1], x1[2], x2[0], x2[1], x2[2];
            colors->GetTuple(i/2, colorVals);
            if (currentRenderer == renderer) {
                json pt = makeEntityGeneric("v", colorVals[0], colorVals[1], colorVals[2], loc, radius);
                sceneVis["list"][sceneVis["list"].size()-1]["e"].push_back(pt);
            }
        }
    }
    return actor;
}

vtkSmartPointer<vtkActor> Widget::renderLines(vtkSmartPointer<vtkPoints> points, MatrixXi &lines, string color, num radius) {
    vtkNew<vtkCellArray> vtkLines;

    for( int i = 0; i < lines.rows(); i++)
    {
        auto row = lines.row(i);
        for (int j=0; j<3; j++) {
            vtkNew<vtkLine> line;
            line->GetPointIds()->SetId(0, row(j));
            line->GetPointIds()->SetId(1, row((j+1)%3));
            vtkLines->InsertNextCell(line);
        }
    }
    return renderLines(points, vtkLines, color, radius);
}

vtkSmartPointer<vtkActor> Widget::renderLines(vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkCellArray> lines, string color, num radius) {
    return renderLines(points, lines, color, radius, 1.0);
}

vtkSmartPointer<vtkActor> Widget::renderLines(vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkCellArray> lines, string color, num radius, num opacity) {
    vtkNew<vtkPolyData> linesPolyData;
    linesPolyData->SetPoints(points);
    linesPolyData->SetLines(lines);
    vtkNew<vtkTubeFilter> tubeFilter;
    tubeFilter->SetInputData(linesPolyData);
    tubeFilter->SetNumberOfSides(8);
    tubeFilter->SetRadius(radius);
    tubeFilter->Update();
    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputConnection(tubeFilter->GetOutputPort());

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetLineWidth(4);
    auto gc = colors->GetColor3d(color);
    auto cd = gc.GetData();
    actor->GetProperty()->SetColor(cd);
    actor->GetProperty()->SetOpacity(opacity);

    currentRenderer->AddActor(actor);
    if (setts["global"]["scene_saving"].value_or(false)) {
        num x1[3], x2[3];
        for (int i=0; i<points->GetNumberOfPoints(); i+=2) {
            points->GetPoint(i, x1);
            points->GetPoint(i+1, x2);
            MatrixXd loc(2, 3);
            loc << x1[0], x1[1], x1[2], x2[0], x2[1], x2[2];
            if (currentRenderer == renderer) {
                json pt = makeEntityGeneric("v", cd[0], cd[1], cd[2], opacity, "", loc, radius);
                sceneVis["list"][sceneVis["list"].size()-1]["e"].push_back(pt);
            }
        }
    }
    return actor;
}

vtkSmartPointer<vtkActor> Widget::renderPolyline(MatrixXd &points, string color, num radius) {
    vector<MatrixXd*> linesVec;
    linesVec.push_back(&points);
    VectorXd offset(3);
    offset << 0, 0, 0;
    return renderPolylines(linesVec, color, radius, offset);
}

vtkSmartPointer<vtkActor> Widget::renderPolylines(vector<MatrixXd*> lines, string color, num radius, VectorXd &offset) {
    return renderPolylines(lines, color, radius, offset, 1);
}

vtkSmartPointer<vtkActor> Widget::renderPolylines(vector<MatrixXd*> lines, string color, num radius, VectorXd &offset, int modVal) {
    vtkNew<vtkPoints> vtkPoints;
    vtkNew<vtkCellArray> cells;
    int ii = 0, k, ptCount;
    for (auto line : lines) {
        if (line->rows() > 1) {
            k = 0, ptCount=0;
            for (auto point : line->rowwise()) {
                if (!(k % modVal)) {
                    if (line->cols() == 2)
                        vtkPoints->InsertNextPoint(point(0) + offset(0), point(1) + offset(1), offset(2));
                    else
                        vtkPoints->InsertNextPoint(point(0) + offset(0), point(1) + offset(1), point(2) + offset(2));
                    ptCount++;
                }
                k++;
            }
            vtkNew<vtkPolyLine> polyLine;
            polyLine->GetPointIds()->SetNumberOfIds(ptCount);
            for (int i = 0; i < ptCount; i++) {
                polyLine->GetPointIds()->SetId(i, ii);
                ii++;
            }
            cells->InsertNextCell(polyLine);
        }
    }
    // Create a polydata to store everything in
    vtkNew<vtkPolyData> polyDataLine;

    // Add the points to the dataset
    polyDataLine->SetPoints(vtkPoints);

    // Add the lines to the dataset
    polyDataLine->SetLines(cells);

    vtkNew<vtkTubeFilter> tubes;
    tubes->SetInputData(polyDataLine);
    tubes->SetRadius(radius);
    tubes->SetNumberOfSides(8);

    // Set up actor and mapper
    vtkNew<vtkPolyDataMapper> mapperLine;
    mapperLine->SetInputConnection(tubes->GetOutputPort());

    vtkNew<vtkActor> actorLine;
    actorLine->SetMapper(mapperLine);
    auto gc = colors->GetColor3d(color);
    auto cd = gc.GetData();
    actorLine->GetProperty()->SetColor(cd);
    currentRenderer->AddActor(actorLine);
    if (setts["global"]["scene_saving"].value_or(false)) {
        for (auto line : lines) {
            if (currentRenderer == renderer && line->rows() > 1) {
                json pt = makeEntityGeneric("y", cd[0], cd[1], cd[2], 1.0, "", *line, radius);
                sceneVis["list"][sceneVis["list"].size()-1]["e"].push_back(pt);
            }
        }
    }
    return actorLine;
}

void Widget::renderCellMap() {
    auto cellGrid = sim->getTissue().cellGrid;
    vtkNew<vtkPoints> cellPoints, linePointsN, linePointsS, linePointsNW, linePointsNE, linePointsSW, linePointsSE;
    vtkNew<vtkCellArray> linesN, linesS, linesNW, linesNE, linesSW, linesSE;
    int n = 0, s = 0, ne = 0, nw = 0, se = 0, sw = 0;
    pastRenderModel = true;

    for (auto cell : sim->getTissue().getCells()) {
        auto point = cell->initialCenter;
        cellPoints->InsertNextPoint(point(0), point(1), point(2));
        auto northC = cell->neighbors.north;
        auto southC = cell->neighbors.south;
        auto northWC = cell->neighbors.northWest;
        auto southWC = cell->neighbors.southWest;
        auto northEC = cell->neighbors.northEast;
        auto southEC = cell->neighbors.southEast;
        if (northC != nullptr) {
            auto point2 = northC->initialCenter;
            linePointsN->InsertNextPoint(point(0), point(1), point(2));
            linePointsN->InsertNextPoint(point2(0), point2(1), point2(2));
            vtkNew<vtkLine> line;
            line->GetPointIds()->SetId(0, n);
            line->GetPointIds()->SetId(1, n+1);
            linesN->InsertNextCell(line);
            n += 2;
        }
        if (southC != nullptr) {
            auto point2 = southC->initialCenter;
            linePointsS->InsertNextPoint(point(0), point(1), point(2));
            linePointsS->InsertNextPoint(point2(0), point2(1), point2(2));
            vtkNew<vtkLine> line;
            line->GetPointIds()->SetId(0, s);
            line->GetPointIds()->SetId(1, s+1);
            linesS->InsertNextCell(line);
            s += 2;
        }
        if (northEC != nullptr) {
            auto point2 = northEC->initialCenter;
            linePointsNE->InsertNextPoint(point(0), point(1), point(2));
            linePointsNE->InsertNextPoint(point2(0), point2(1), point2(2));
            vtkNew<vtkLine> line;
            line->GetPointIds()->SetId(0, ne);
            line->GetPointIds()->SetId(1, ne+1);
            linesNE->InsertNextCell(line);
            ne += 2;
        }
        if (southEC != nullptr) {
            auto point2 = southEC->initialCenter;
            linePointsSE->InsertNextPoint(point(0), point(1), point(2));
            linePointsSE->InsertNextPoint(point2(0), point2(1), point2(2));
            vtkNew<vtkLine> line;
            line->GetPointIds()->SetId(0, se);
            line->GetPointIds()->SetId(1, se+1);
            linesSE->InsertNextCell(line);
            se += 2;
        }
        if (northWC != nullptr) {
            auto point2 = northWC->initialCenter;
            linePointsNW->InsertNextPoint(point(0), point(1), point(2));
            linePointsNW->InsertNextPoint(point2(0), point2(1), point2(2));
            vtkNew<vtkLine> line;
            line->GetPointIds()->SetId(0, nw);
            line->GetPointIds()->SetId(1, nw+1);
            linesNW->InsertNextCell(line);
            nw += 2;
        }
        if (southWC != nullptr) {
            auto point2 = southWC->initialCenter;
            linePointsSW->InsertNextPoint(point(0), point(1), point(2));
            linePointsSW->InsertNextPoint(point2(0), point2(1), point2(2));
            vtkNew<vtkLine> line;
            line->GetPointIds()->SetId(0, sw);
            line->GetPointIds()->SetId(1, sw+1);
            linesSW->InsertNextCell(line);
            sw += 2;
        }
    }
    
    num v = 0.02 * scaleFactor;
    renderLines(linePointsN, linesN, "Blue", v, 0.5);
    renderLines(linePointsS, linesS, "Red", v, 0.5);
    renderLines(linePointsNE, linesNE, "Blue", v, 0.5);
    renderLines(linePointsSE, linesSE, "Red", v, 0.5);
    renderLines(linePointsNW, linesNW, "Blue", v, 0.5);
    renderLines(linePointsSW, linesSW, "Red", v, 0.5);
    renderSpheres(cellPoints, "Yellow", 0.04 * scaleFactor);
}

void Widget::duplicateLayerCubic() {
    sim->getTissue().addLayersCubic();
    // TODO: Render the new vertices and edges, might just want to throw
    // everything away and add everything back in again
    firstRenderGrid = true;
    renderGrid();
}

void Widget::finishInitCubic() {
    currentRenderer = renderer;
    vtkSmartPointer<vtkActor> dummy;
    for (auto s : model.splines) {
        MatrixXd sp = s->splinePoints(Eigen::seq(0,last,20),all);
        renderSpheres(sp, "white", 0.03 * scaleFactor, dummy);
    }
    for (auto s : model.usplines) {
        MatrixXd sp = s->splinePoints(Eigen::seq(0,last,20),all);
        renderSpheres(sp, "white", 0.03 * scaleFactor, dummy);
    }
    // TODO: Clean this up/figure out why the initial import doesn't work
    vtkNew<vtkOBJImporter> importer2;
    importer2->SetFileName(model.modelFile.c_str());
    importer2->Update();
    vtkActorCollection *allActors2 = importer2->GetRenderer()->GetActors();
    auto modelActor2 = allActors2->GetLastActor();
    modelActor2->SetScale(model.scale);
    modelActor2->GetProperty()->LightingOff();
    modelActor2->GetProperty()->SetOpacity(0.75);
    modelActor2->GetProperty()->SetColor(
            colors->GetColor3d("darkgreen").GetData());
    renderer->AddActor(modelActor2);
    ////////////////////////
    currentRenderer->ResetCamera();
    currentRenderer->GetRenderWindow()->Render();
}

void Widget::renderModel() {
    vtkNew<vtkOBJImporter> importer;
    importer->SetFileName(model.modelFile.c_str());
    importer->Update();
    vtkActorCollection *allActors = importer->GetRenderer()->GetActors();
    vtkSmartPointer<vtkActor> modelActor;
    vtkSmartPointer<vtkPoints> points;
    vtkSmartPointer<vtkCellArray> lines;
    vector<MatrixXd*> offsetPoints;
    // For shifting the model away from the parameterization
    VectorXd offset(3);
    offset << 0, 0, 0;
    offset = offset.transpose();
    MatrixXd vertMat = MatrixXd::Zero(0, 3);
    modelActor = allActors->GetLastActor();
    modelActor->SetScale(model.scale);
    auto mesh = dynamic_cast<vtkPolyData*>(modelActor->GetMapper()->GetInput());
    points = accessEachVertex(mesh, vertMat);
    vtkNew<vtkTransform> transform;
    transform->Scale(model.scale, model.scale, model.scale);
    transform->TransformPoints(points, points);
    lines = accessEachFace(mesh);
    string pointString = "point";

    // Real points need offset, and cell points need put into a
    // single data structure and offset
    vtkNew<vtkPoints> cellPoints, cellPointsShifted, hexPoints, cartEsts, pSplinePoints;

    for (auto point : grid->deletedCCs)
        cellPoints->InsertNextPoint((*point)(0), (*point)(1), (*point)(2));

    for (auto cell : sim->getTissue().getCells()) {
        auto point = cell->initialCenter;
        cellPointsShifted->InsertNextPoint(point(0), point(1), point(2));
    }

    for (auto col : sim->getTissue().cellGrid.grid) {
        for (auto wrapper : col)
            for (MatrixXd *mat : wrapper->cell->splinePoints)
                offsetPoints.push_back(mat);
    }

    int k;
    vector<MatrixXd*> hexEdges;
    for (auto edge : grid->getEdges()) {
        k = 0;
        vector<Vector3d> currEdge;
        for (auto point : edge->splinePoints.rowwise()) {
            if (!(k % 8) || point == edge->splinePoints.row(edge->splinePoints.rows()-1))
                currEdge.push_back(point);
            k++;
        }
        MatrixXd *currEMat = new MatrixXd();
        *currEMat = MatrixXd::Zero(currEdge.size(), 3);
        for (long unsigned int i=0; i<currEdge.size(); i++)
            currEMat->row(i) = currEdge[i];
        hexEdges.push_back(currEMat);
    }

    k = 0;
    for (auto row : model.primarySplinePoints.rowwise()) {
        if (!(k % 20))
            pSplinePoints->InsertNextPoint(row(0), row(1), row(2));
        k++;
    }

    for (auto vert : grid->getVertices()) {
        if (vert->cartCount) {
            Vector3d p = vert->cartEstimate;
            cartEsts->InsertNextPoint(p(0), p(1), p(2));
        }
    }

    Vector3d minPoint, maxPoint;
    num ma = std::numeric_limits<num>::max();
    num mi = std::numeric_limits<num>::min();
    minPoint << ma, ma, ma;
    maxPoint << mi, mi, mi;
    for (Vector3d row : model.mesh->realPoints.rowwise()) {
        for (int z=0; z<3; z++) {
            if (row(z) < minPoint(z))
                minPoint(z) = row(z);
            if (row(z) > maxPoint(z))
                maxPoint(z) = row(z);
        }
    }
    vector<MatrixXd*> uvPoints;
    for (auto s : model.splines) {
        uvPoints.push_back(&s->splinePointsUV);
    }
    vtkNew<vtkPoints> uvControlPoints;

    for (auto spline : model.splines) {
        for (auto point : spline->uvPoints.rowwise()) {
            uvControlPoints->InsertNextPoint(point(0), point(1), point(2));
        }
    }

    // I came up with the numbers below for the face model
    // where this value is about 6.67, hence the divisor
    //scaleFactor = distance(minPoint, maxPoint) / 6.67 * (15.0 / static_cast<num>(model.startingLineNodeCount));

    paramRenderer = vtkSmartPointer<vtkRenderer>::New();
    paramRenderer->SetBackground(colors->GetColor3d(theme.paramBackground).GetData());

    paramWindow = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    paramWindow->AddRenderer(paramRenderer);
    ui->qvtkWidgetParam->setRenderWindow(paramWindow);

    if (!cubic && !model.batch) {
        MatrixXd paramPointsMoved = model.mesh->paramPoints.rowwise() + offset;
        currentRenderer = paramRenderer;
        //themeMap.paramModelPointsR = 0.005;
        themeMap.paramModelPointsR = 0.0;
        auto vtkParamPoints = renderSpheres(paramPointsMoved, theme.paramModelPoints, themeMap.paramModelPointsR, themeMap.paramModelPoints);
        MatrixXi paramFaces = MatrixXi::Zero(model.mesh->faces.size(), 3);
        for (long unsigned int i=0; i<model.mesh->faces.size(); i++)
            paramFaces.row(i) << model.mesh->faces[i]->vertices;
        themeMap.paramModelFaces = renderLines(vtkParamPoints, paramFaces, theme.paramModelFaces, 0.0013);
        paramPointsMoved = model.primaryUVPoints.rowwise() + offset;
        themeMap.paramPointsPrimaryR = 0.006;
        themeMap.paramControlPointsR = 0.0045;
        themeMap.paramSplines = renderPolylines(uvPoints, theme.paramSplines, 0.003, offset, 20);
        themeMap.paramPrimarySpline = renderPolyline(model.primarySpline.splinePointsUV, theme.paramPrimarySpline, 0.003);
        currentRenderer = renderer;
    }
    themeMap.primaryPointsR = 0.03 * scaleFactor;
    renderSpheres(model.primaryRealPoints, theme.primaryPoints, themeMap.primaryPointsR, themeMap.primaryPoints);
    offset << 0, 0, 0;
    if (!static_cast<string>(setts["model"]["type"].value_or("obj")).compare("obj")) {
        themeMap.cellPointsAllR = 0.045 * scaleFactor;
        themeMap.cellPointsUsingR = 0.065 * scaleFactor;
    } else {
        themeMap.cellPointsAllR = 0.025 * scaleFactor;
        themeMap.cellPointsUsingR = 0.03 * scaleFactor;
    }
    themeMap.cellPointsAll = renderSpheres(cellPoints, theme.cellPointsAll, themeMap.cellPointsAllR);
    centerActor = renderSpheres(cellPointsShifted, theme.cellPointsUsing, themeMap.cellPointsUsingR);
    themeMap.cellPointsUsing = centerActor;
    themeMap.cellSplines = renderPolylines(hexEdges, theme.cellSplines, 0.015 * scaleFactor, offset);
    themeMap.hexVerticesR = 0.03 * scaleFactor;
    themeMap.hexVertices = renderSpheres(cartEsts, theme.hexVertices, themeMap.hexVerticesR);
    themeMap.primarySpline = renderPolyline(model.primarySplinePoints, theme.primarySpline, 0.02 * scaleFactor);
    themeMap.vSplines = renderPolylines(model.vSplines, theme.vSplines, 0.02 * scaleFactor, offset, 20);
    // TODO: Clean this up/figure out why the initial import doesn't work
    vtkNew<vtkOBJImporter> importer2;
    importer2->SetFileName(model.modelFile.c_str());
    importer2->Update();
    vtkActorCollection *allActors2 = importer2->GetRenderer()->GetActors();
    auto modelActor2 = allActors2->GetLastActor();
    modelActor2->SetScale(model.scale);
    modelActor2->GetProperty()->LightingOff();
    modelActor2->GetProperty()->SetOpacity(0.35);
    modelActor2->GetProperty()->SetColor(
            colors->GetColor3d("darkgreen").GetData());
    renderer->AddActor(modelActor2);
    renderer->ResetCamera();
    renderer->GetRenderWindow()->Render();
}

void Widget::renderTD() {
    auto edges = grid->getEdges(true);
    auto checkEdge = edges[0];

    for (auto edge : edges) {
        edge->colorVtkEntitiesTD();
    }
    
    if (checkEdge->tdActions.size())
        funDeq.push_back(&Widget::renderTD);
    else
        funDeq.push_back(&Widget::renderGrid);
    renderer->ResetCamera();
    renderer->GetRenderWindow()->Render();
}

void Widget::renderAgency() {
    auto edges = grid->getEdges();

    bool checked = ui->vizErrorCheckBox->isChecked();
    for (auto edge : edges) {
        edge->colorVtkEntities(checked);
    }

    funDeq.push_back(&Widget::renderGrid);
    renderer->ResetCamera();
    renderer->GetRenderWindow()->Render();
}

void Widget::renderGridAndSave() {
    auto dirName = "../exports/temp";
    auto pngName = fmt::format("{:05d}.png", grid->epochsRun); 
    auto fName = fmt::format("{}/{}", dirName, pngName);
    cout << fName << endl;
    renderGrid();
    if (grid->freshReset) {
        grid->freshReset = false;
        std::filesystem::remove_all(dirName);
        std::filesystem::create_directory(dirName);
    }
    writeImageWithText(fName);
}

void Widget::renderGrid() {
    if (cubic)
        gridHasBeenRendered = true;
    string fName;
    if (model.batch)
        fName = fmt::format("{}/epoch_{:09d}.png", outStringBatch, grid->epochsRun);
    ui->shapeOpRunning->setText(QString::fromStdString("ShapeOp not Running"));
    ui->shapeOpRunning->setStyleSheet("QLabel {}");
    ui->epochTimer->setText(QString::fromStdString(fmt::format("{:.8f}", sim->secPerEpoch)));
    ui->epochsRunLabel->setText(QString::fromStdString(fmt::format("{}", grid->epochsRun)));
    runningSO = false;
    vtkNew<vtkNamedColors> colors;
    int c=0;
    auto vertices = grid->getVertices();
    vtkNew<vtkPoints> points;
    for (auto vertex : vertices) {
        if (vertex->vtkInitialized()) {
            // Calling here instead of automatically with updatePosition
            // to save cycles when not re rendering
            vertex->updateVtkEntities();
        } else {
            auto row = vertex->getPosition();

            vertex->setVtkPoints(points);
            vertex->setPtId(points->InsertNextPoint(row(0), row(1), row(2)));
        }
    }

    auto edges = grid->getEdges();

    if (firstRenderGrid) {
        renderSpheres(points, "White", 0.05 * scaleFactor);
        for (auto e : edges)
            e->origRestLength = e->getRestLength();
    }


    vtkNew<vtkCellArray> lines;
    vtkNew<vtkPoints> pts;
    vtkNew<vtkFloatArray> lineColors;
    lineColors->SetNumberOfComponents(3);
    lineColors->SetName("Colors");

    c = 0;
    for (auto edge : edges) {
        if (edge->vtkInitialized()) {
            edge->updateVtkEntities();
        } else {
            auto pt1 = edge->getV1()->getPosition();
            auto pt2 = edge->getV2()->getPosition();

            auto id1 = pts->InsertNextPoint(pt1(0), pt1(1), pt1(2));
            auto id2 = pts->InsertNextPoint(pt2(0), pt2(1), pt2(2));
            edge->setVtkPoints(pts);
            edge->setPtIds(id1, id2);
            vtkNew<vtkLine> line;
            line->GetPointIds()->SetId(0, c);
            line->GetPointIds()->SetId(1, c+1);
            lines->InsertNextCell(line);
            Vector3d color = ratioToRGB(distance(pt1, pt2) / edge->getRestLength());
            lineColors->InsertNextTuple3(color(0), color(1), color(2));
            edge->setLineColorsArray(c / 2, lineColors);
            c += 2;
        }
    }

    if (firstRenderGrid)
        renderLines(pts, lines, lineColors, 0.025 * scaleFactor, 1.0);

    renderer->ResetCamera();
    renderer->GetRenderWindow()->Render();
    if (model.batch) {
        if (firstWrite) {
            timer = new QTimer(this);
            auto f = std::bind(&Widget::writeImageWithText, this, fName);
            timer->singleShot(500, f);
        } else {
            writeImageWithText(fName);
        }
    }
    if (setts["global"]["scene_saving"].value_or(false) && !firstRenderGrid) {
        for (auto vert : vertices) {
            MatrixXd loc(1, 3);
            loc << vert->getPosition().transpose();
            if (currentRenderer == renderer) {
                json pt = makeEntityGeneric("p", 1.0, 1.0, 1.0, loc, 0.05 * scaleFactor);
                sceneVis["list"][sceneVis["list"].size()-1]["e"].push_back(pt);
            }
        }
        for (auto edge : edges) {
            auto pt1 = edge->getV1()->getPosition();
            auto pt2 = edge->getV2()->getPosition();
            Vector3d color = ratioToRGB(distance(pt1, pt2) / edge->getRestLength());
            MatrixXd loc(2, 3);
            loc << pt1.transpose(), pt2.transpose();
            if (currentRenderer == renderer) {
                json pt = makeEntityGeneric("v", color(0), color(1), color(2), 1.0, "", loc, 0.025 * scaleFactor);
                sceneVis["list"][sceneVis["list"].size()-1]["e"].push_back(pt);
            }
        }
    }
    firstRenderGrid = false;
}

void Widget::triangulateGrid() {
    triangulateGrid(true);
}

vtkSmartPointer<vtkActor> Widget::triangulateGrid(bool render) {
    auto actor = sim->triangulateGrid();
    actor->GetProperty()->SetColor(
            colors->GetColor3d(theme.tris).GetData());

    if (render) {
        auto actors = rendererTri->GetActors();
        actors->InitTraversal();
        while (actors->GetNumberOfItems())
            rendererTri->RemoveActor(actors->GetNextItem());

        rendererTri->AddActor(actor);
        themeMap.tris = actor;
        rendererTri->ResetCamera();
        vtkWindowTri->Render();
    }
    return actor;
}

void Widget::writeImageWithText(string fName) {
    unsigned char black[]  = {0, 0, 0};
    writeImage(fName, vtkWindow, false);
    CImg<unsigned char> img(fName.c_str());
    img.draw_text(10, img.height() - 40, "kA90: %.2f, kA120: %.2f, kLinear: %.2f, Epochs: %d",
            black, 0, 1, 32, sim->getTissue().getGrid().kAng90,
            sim->getTissue().getGrid().kAng, sim->getTissue().getGrid().kEdge,
            firstWrite ? 0 : grid->epochsRun);
    firstWrite = false;
    img.save_png(fName.c_str());
}

void Widget::removeActor(vtkSmartPointer<vtkActor> actor) {
    renderer->RemoveActor(actor);
}

Simulation *Widget::getSimulation() {
    return sim;
}

void Worker::cellFinder(Widget *widget, num* pickPos) {
    auto cells = widget->getSimulation()->getTissue().getCells();
    Vector3d pos;
    pos << pickPos[0], pickPos[1], pickPos[2];
    num bestDist = std::numeric_limits<num>::max();
    Cell *bestCell = nullptr;
    for (auto cell : cells) {
        Vector3d c = cell->initialCenter;
        num dist = distance(pos, c);
        if (dist < bestDist) {
            bestDist = dist;
            bestCell = cell;
        }
    }
    std::stringstream ss;
    ss << *bestCell;
    QString out = QString::fromStdString(ss.str());
    emit finished(out);
}

void clickCallbackFunction(vtkObject* caller,
        long unsigned int eventId,
        void* clientData,
        void* vtkNotUsed(callData))
{
    vtkNew<vtkNamedColors> colors;
    auto interactor = reinterpret_cast<vtkRenderWindowInteractor*>(caller);
    auto widget = reinterpret_cast<Widget*>(clientData);
    int* clickPos = interactor->GetEventPosition();
    // Pick from this location.
    vtkNew<vtkPropPicker> picker;
    picker->Pick(clickPos[0], clickPos[1], 0, interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());

    auto pickedActor = picker->GetActor();
    if (pickedActor == widget->centerActor) {
        // I feel like there should be a simpler way to do this, but this is what I came
        // up with after going through the Qt docs for a while. It does work and even
        // allows multiple background threads in tandem without blocking GUI updates,
        // just seems overly complicated.
        num* pos = picker->GetPickPosition();
        Worker *worker = new Worker;
        QThread *workerThread = new QThread;
        worker->moveToThread(workerThread);
        widget->connect(widget, &Widget::operate, worker, &Worker::cellFinder);
        widget->connect(worker, &Worker::finished, widget->ui->infoBox, &QTextBrowser::setPlainText);
        widget->connect(worker, &Worker::finished, workerThread, &QThread::quit);
        workerThread->start();
        emit widget->operate(widget, pos);
    }
}

void Widget::on_continueButton_clicked() {
    if (model.batch) {
        return;
    }
    if (setts["global"]["scene_saving"].value_or(false)) {
        if (!ui->concatCheckBox->isChecked()) {
            sceneVis.clear();
            sceneVis["reset"] = true;
            sceneVis["glyph"] = true;
        }
        sceneVis["list"].push_back({});
    }
    sim->getTissue().getGrid().kEdge = boost::lexical_cast<num>(ui->linearKEdit->text().toStdString());
    sim->getTissue().getGrid().kAng = boost::lexical_cast<num>(ui->angularKEdit->text().toStdString());
    sim->getTissue().getGrid().kAng90 = boost::lexical_cast<num>(ui->angular90KEdit->text().toStdString());
    if (funDeq.size()) {
        if (!gridHasBeenRendered) {
            auto actors = renderer->GetActors();
            actors->InitTraversal();
            while (actors->GetNumberOfItems())
                renderer->RemoveActor(actors->GetNextItem());
        }
        (this->*funDeq.front())();
        funDeq.pop_front();
        renderer->ResetCamera();
        renderer->GetRenderWindow()->Render();
    } else {
        gridHasBeenRendered = true;
        runSO();
        renderer->ResetCamera();
        renderer->GetRenderWindow()->Render();
    }
}

void Widget::on_flipWindowsButton_clicked() {
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> dummy;
    auto temp1 = ui->qvtkWidget->renderWindow();
    auto temp2 = ui->qvtkWidgetParam->renderWindow();
    exporter->SetActiveRenderer(temp2->GetRenderers()->GetNextItem());
    exporter->SetRenderWindow(temp2);
    ui->qvtkWidgetParam->setRenderWindow(dummy);
    ui->qvtkWidget->setRenderWindow(temp2);
    ui->qvtkWidgetParam->setRenderWindow(temp1);
    temp1->Render();
    temp2->Render();
}

num Widget::evaluateChamfer(bool render) {
    auto ev = sim->evaluateChamfer();
    if (render) {
        auto actors = renderer->GetActors();
        actors->InitTraversal();
        while (actors->GetNumberOfItems())
            renderer->RemoveActor(actors->GetNextItem());

        vtkSmartPointer<vtkActor> dummy;

        // Add the geometry and topology to the polydata
        auto triPD = dynamic_cast<vtkPolyData*>(ev.actor->GetMapper()->GetInput());
        MatrixXd vv = MatrixXd::Zero(0, 3);
        vtkSmartPointer<vtkPoints> triPts1 = accessEachVertex(triPD, vv);
        vtkNew<vtkPoints> triPts;

        vtkNew<vtkTransform> transform;
        transform->Scale(model.scale, model.scale, model.scale);
        transform->TransformPoints(triPts1, triPts);
        for (int i=0; i<triPts->GetNumberOfPoints(); i++) {
            Vector3d x;
            x << ev.result.row(i);
            triPts->SetPoint(i, x(0), x(1), x(2));
        }
        // Create a polydata object
        vtkSmartPointer<vtkPolyData> polyData =
            vtkSmartPointer<vtkPolyData>::New();

        vtkNew<vtkFloatArray> ptColors;
        ptColors->SetNumberOfComponents(3);
        ptColors->SetName("Colors");

        int c = triPts->GetNumberOfPoints();
        for (int i=0; i<c; i++) {
            Vector3d color, ptV;
            auto pt = triPts->GetPoint(i);
            ptV << pt[0], pt[1], pt[2];
            num minD = std::numeric_limits<num>::max(), currD;
            for (Vector3d r2 : ev.dst2.rowwise()) {
                currD = distance(ptV, r2);
                if (currD < minD)
                    minD = currD;
            }
            Vector3d colorV;
            num x = minD / ev.maxD;
            colorV << std::min(1.0, 2.0 * x), std::min(1.0, 2.0 * (1.0 - x)), 0.0;
            ptColors->InsertNextTuple3(colorV(0), colorV(1), colorV(2));
        }

        vtkNew<vtkLookupTable> lut;
        lut->SetNumberOfTableValues(4);
        lut->SetTableRange(0.0, ev.maxD);
        int tableSize = 256;
        lut->SetNumberOfTableValues(tableSize);
        lut->Build();

        for (auto i = 0; i < lut->GetNumberOfColors(); ++i)
        {
            num s = static_cast<num>(i) / lut->GetNumberOfColors();
            std::array<num, 4> rgba{std::min(1.0, 2.0 * s), std::min(1.0, 2.0 * (1.0 - s)), 0.0, 1.0};
            lut->SetTableValue(static_cast<vtkIdType>(i), rgba.data());
        }

        // Add the geometry and topology to the polydata
        polyData->SetPoints(triPts);
        polyData->SetPolys(triPD->GetPolys());
        polyData->GetPointData()->SetScalars(ptColors);
        vtkNew<vtkPolyDataMapper> mapper;
        mapper->SetInputData(polyData);
        mapper->SetColorMode(2);
        mapper->SetScalarModeToUsePointData();
        mapper->SetScalarRange(0.0, ev.maxD);

        vtkNew<vtkActor> actor2;
        actor2->SetMapper(mapper);


        auto scaleSize = 1750;

        scalarBar->SetNumberOfLabels(4);

        mapper->SetLookupTable(lut);
        scalarBar->SetLookupTable(lut);
        scalarBar->SetTitle("Chamfer Distance");
        scalarBar->GetTitleTextProperty()->SetFontSize(36);
        scalarBar->GetLabelTextProperty()->SetFontSize(24);
        scalarBar->GetAnnotationTextProperty()->SetFontSize(24);
        auto gc3d = colors->GetColor3d(theme.scaleBarFont);
        auto cl = gc3d.GetData();
        scalarBar->GetTitleTextProperty()->SetColor(cl);
        scalarBar->GetLabelTextProperty()->SetColor(cl);
        scalarBar->GetAnnotationTextProperty()->SetColor(cl);
        scalarBar->UnconstrainedFontSizeOn();
        scalarBar->SetMaximumWidthInPixels(scaleSize / 8);
        scalarBar->SetMaximumHeightInPixels(scaleSize / 3);
        scalarBar->SetPosition(0.8, 0.25);

        ev.modelActor->GetProperty()->SetColor(
                colors->GetColor3d("darkmagenta").GetData());
        ev.modelActor->GetProperty()->SetOpacity(0.3);
        renderer->AddActor(ev.modelActor);

        renderer->AddActor(actor2);
        renderer->AddViewProp(scalarBar);

        renderer->ResetCamera();
        renderer->GetRenderWindow()->Render();
    }

    currErr = ev.avgChamf;
    return ev.avgChamf;
}

void Widget::evaluate(bool render) {
    if (render) {
        auto actors = renderer->GetActors();
        actors->InitTraversal();
        while (actors->GetNumberOfItems())
            renderer->RemoveActor(actors->GetNextItem());
    }

    auto cells = sim->getTissue().getCells();
    vector<Vector3d> pointVec;
    std::set<Vertex*> vertSet;
    vector<Vector3d> vertVecSrc, vertVecDst;
    vtkNew<vtkPoints> targetPoints; 
    for (auto cell : cells) {
        auto verts = cell->getVertices();
        for (int i=0; i<6; i++) {
            if (!vertSet.count(verts[i]) && !vertSet.count(verts[i+6]) && verts[i]->cartCount) {
                Vector3d v = verts[i]->getPosition() + verts[i+6]->getPosition();
                v /= 2.0;
                vertVecSrc.push_back(v);
                vertSet.insert(verts[i]);
                vertSet.insert(verts[i+6]);
                Vector3d v2 = verts[i]->cartEstimate;
                targetPoints->InsertNextPoint(v2(0), v2(1), v2(2));
                vertVecDst.push_back(v2);
            }
        }
    }
    SO::Matrix3X src(3, vertVecSrc.size()), dst(3, vertVecSrc.size());
    for (long unsigned int i=0; i<vertVecSrc.size(); i++) {
        src.col(i) = vertVecSrc[i];
        dst.col(i) = vertVecDst[i];
    }

    Pointcloud pointcloudA = dst;
    Pointcloud pointcloudB = src;

    auto callback = Callback(pointcloudA, pointcloudB);

    lsqcpp::GaussNewton<num, 6, Eigen::Dynamic, Objective, lsqcpp::ArmijoBacktracking, lsqcpp::DenseCholeskySolver> optimizer;
    optimizer.setMinimumGradientLength(1e-3F);
    optimizer.setMinimumStepLength(1e-3F);
    optimizer.setObjective({pointcloudA, pointcloudB});
    optimizer.setCallback(callback);
    optimizer.setMaximumIterations(1000);
    optimizer.setVerbosity(0);

    Vector6 xval(6);
    xval.setZero();
    xval.segment(3, 3) = lsqcpp::parameter::encodeRotation(Matrix3d::Identity());

    optimizer.minimize(xval);
    MatrixXd result(src.cols(), src.rows());
    for(int i = 0; i < src.cols(); i++)
    {
        result.row(i) = (rotFinal * src.col(i) + transFinal).transpose();
    }
    MatrixXd src2 = src.transpose();
    MatrixXd dst2 = dst.transpose();
    currErr = (result - dst2).rowwise().norm().sum();
    cout << "CurrErr: " << currErr << endl;

    if (render) {
        for (Vector3d r : result.rowwise())
            pointVec.push_back(r);

        MatrixXd gridMat = MatrixXd::Zero(pointVec.size(), 3);
        for (long unsigned int i=0; i<pointVec.size(); i++)
            gridMat.row(i) = pointVec[i];
        vtkSmartPointer<vtkActor> dummy;
        renderSpheres(gridMat, "Red", 0.035 * scaleFactor, dummy);
        renderSpheres(targetPoints, "Yellow", 0.035 * scaleFactor);

        vtkNew<vtkCellArray> lines;
        vtkNew<vtkPoints> pts;
        vtkNew<vtkFloatArray> lineColors;
        lineColors->SetNumberOfComponents(3);
        lineColors->SetName("Colors");

        int c = 0;
        for (int i=0; i<dst2.rows(); i++) {
            vtkNew<vtkLine> line;
            line->GetPointIds()->SetId(0, c);
            line->GetPointIds()->SetId(1, c+1);
            lines->InsertNextCell(line);
            lineColors->InsertNextTuple3(1.0, 1.0, 1.0);
            c += 2;
        }

        renderLines(pts, lines, lineColors, 0.012 * scaleFactor, 1.0);

        renderer->ResetCamera();
        renderer->GetRenderWindow()->Render();
    }
}

void Widget::updatePBar(int iter, int total) {
    ui->cyclePBar->setMinimum(0);
    ui->cyclePBar->setMaximum(total);
    ui->cyclePBar->setValue(iter);
}

void Widget::on_evaluateButton_clicked() {
    evaluateChamfer(true);
}

void Widget::on_resetGridButton_clicked() {
    grid->reset();
    renderGrid();
}

void Widget::on_printButton_clicked() {
    exporter->SetFileName(ui->objTextEdit->text().toLocal8Bit());
    exporter->Update();
}

void Widget::on_pngPrintButton_clicked() {
    string fPrefix = ui->pngTextEdit->text().toStdString();
    int count = 0;
    string fName = fPrefix + ".png";
    if (std::filesystem::exists(fName)) {
        count++;
        fName = fmt::format("{}-{}.png", fPrefix, count);
        while (std::filesystem::exists(fName)) {
            count++;
            fName = fmt::format("{}-{}.png", fPrefix, count);
        }
    }
    cout << fName << endl;
    writeImage(fName, ui->qvtkWidget->renderWindow(), false);
}

void Widget::on_visPrintButton_clicked() {
    std::ofstream o("out_vis.json");
    o << sceneVis;
    o.close();
}

void Widget::on_jsonPrintButton_clicked() {
    json output;
    for (auto v : grid->vertices) {
        json entry;
        auto p = v->getPosition();
        entry["position"] = {p(0), p(1), p(2)};
        output["vertices"].push_back(entry);
    }
    for (auto e : grid->getEdges()) {
        json entry;
        entry["rest_length"] = e->getRestLength();
        entry["vertices"] = {e->getV1()->getIndex(), e->getV2()->getIndex()};
        output["edges"].push_back(entry);
    }
    for (auto c : sim->getTissue().getCells()) {
        json entry;
        for (auto e : c->getEdges())
            entry["edges"].push_back(e->getIndex());
        for (auto v : c->getVertices())
            entry["vertices"].push_back(v->getIndex());
        output["cells"].push_back(entry);
    }
    output["scale_factor"] = scaleFactor;

    std::ofstream o(fmt::format("{}.json", ui->jsonTextEdit->text().toStdString()));
    o << output;
    o.close();
}

void Widget::on_lightModeBox_stateChanged() {
    if (ui->lightModeBox->isChecked())
        changeTheme(light);
    else
        changeTheme(dark);
}

void Widget::on_centerButton_clicked() {
    currentRenderer->ResetCamera();
    currentRenderer->GetRenderWindow()->Render();
}

void Widget::on_triToggleButton_clicked() {
    if (!triangulated) {
        currentRenderer = rendererTri;
        ui->triToggleButton->setText(QString::fromStdString("Show Grid"));
        ui->qvtkWidget->setRenderWindow(vtkWindowTri);
        exporter->SetActiveRenderer(rendererTri);
        exporter->SetRenderWindow(vtkWindowTri);
        triangulateGrid();
        vtkWindowTri->Render();
    } else {
        currentRenderer = renderer;
        ui->triToggleButton->setText(QString::fromStdString("Triangulate"));
        ui->qvtkWidget->setRenderWindow(vtkWindow);
        exporter->SetActiveRenderer(renderer);
        exporter->SetRenderWindow(vtkWindow);
        vtkWindow->Render();
    }
    triangulated = !triangulated;
}

void Widget::wheelEvent(QWheelEvent *event) {
    renderer->GetRenderWindow()->Render();
}

Widget::~Widget() {
    delete ui;
}

void Widget::changeTheme(Theme &theme) {
    this->theme = theme;
    paramRenderer->SetBackground(colors->GetColor3d(theme.paramBackground).GetData());
    renderer->SetBackground(colors->GetColor3d(theme.background).GetData());
    rendererTri->SetBackground(colors->GetColor3d(theme.triBackground).GetData());

    scalarBar->GetTitleTextProperty()->SetColor(colors->GetColor3d(theme.scaleBarFont).GetData());
    scalarBar->GetLabelTextProperty()->SetColor(colors->GetColor3d(theme.scaleBarFont).GetData());
    scalarBar->GetAnnotationTextProperty()->SetColor(colors->GetColor3d(theme.scaleBarFont).GetData());

    vector<vtkSmartPointer<vtkActor>> actors;
    vector<vtkSmartPointer<vtkRenderer>> rends;
    vector<num> radii;
    vector<string> themeVals;
    if (!pastRenderModel) {
        actors.push_back(themeMap.primaryPoints);
        radii.push_back(themeMap.primaryPointsR);
        themeVals.push_back(theme.primaryPoints);
        rends.push_back(renderer);
        actors.push_back(themeMap.cellPointsAll);
        radii.push_back(themeMap.cellPointsAllR);
        themeVals.push_back(theme.cellPointsAll);
        rends.push_back(renderer);
        actors.push_back(themeMap.cellPointsUsing);
        radii.push_back(themeMap.cellPointsUsingR);
        themeVals.push_back(theme.cellPointsUsing);
        rends.push_back(renderer);
        actors.push_back(themeMap.hexVertices);
        radii.push_back(themeMap.hexVerticesR);
        themeVals.push_back(theme.hexVertices);
        rends.push_back(renderer);
    }

    for (long unsigned int i=0; i<actors.size(); i++) {
        auto actor = actors[i];
        num radius = radii[i];
        auto ren = rends[i];
        vtkNew<vtkSphereSource> sphereNew;
        sphereNew->SetRadius(radius);
        vtkNew<vtkGlyph3DMapper> mapperNew;
        mapperNew->SetInputData(actor->GetMapper()->GetInput());
        mapperNew->SetSourceConnection(sphereNew->GetOutputPort());
        vtkNew<vtkActor> NewA;
        NewA->SetMapper(mapperNew);
        ren->RemoveActor(actor);
        ren->AddActor(NewA);
        actor = NewA;
        actor->GetProperty()->SetColor(colors->GetColor3d(themeVals[i]).GetData());
    }

    themeMap.primarySpline->GetProperty()->SetColor(colors->GetColor3d(theme.primarySpline).GetData());
    themeMap.cellSplines->GetProperty()->SetColor(colors->GetColor3d(theme.cellSplines).GetData());
    themeMap.vSplines->GetProperty()->SetColor(colors->GetColor3d(theme.vSplines).GetData());
    themeMap.paramModelFaces->GetProperty()->SetColor(colors->GetColor3d(theme.paramModelFaces).GetData());
    themeMap.paramSplines->GetProperty()->SetColor(colors->GetColor3d(theme.paramSplines).GetData());
    themeMap.paramPrimarySpline->GetProperty()->SetColor(colors->GetColor3d(theme.paramPrimarySpline).GetData());

    if (themeMap.tris != nullptr)
        themeMap.tris->GetProperty()->SetColor(colors->GetColor3d(theme.tris).GetData());

    ui->qvtkWidget->renderWindow()->Render();
    ui->qvtkWidgetParam->renderWindow()->Render();
}
