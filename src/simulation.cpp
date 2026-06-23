#include "simulation.h"
#include "agent.h"

Simulation::Simulation() :
        model( Model(false) ), tissue( Tissue(model, false) ) {
    Simulation(false);
}

Simulation::Simulation(bool cubic) :
        model( Model(cubic) ), tissue( Tissue(model, cubic) ) {
    s = std::make_shared<SO::Solver>();
    firstRun = true;
    this->cubic = cubic;
    tdEpochSteps = 14;

    if (!cubic) {
        if (setts["global"]["import_from_json"].value_or(false)) {
            tissue.importJSONFull();
        } else {
            tissue.exportJSON();
        }
        Agent *agent = new Agent(tissue.grid);
        agentPointer = agent;
    }
}

Tissue &Simulation::getTissue() {
    return tissue;
}

Model &Simulation::getModel() {
    return model;
}

void Simulation::runSO() {
    runSO(EPOCHS_PER_CYCLE);
}

void Simulation::runSO(int epochs) {
    runSO(epochs, 0);
}

void Simulation::runSOCubic(int epochs) {
    auto edges = tissue.grid.getEdges();
    auto constraints = tissue.grid.getSOEdgesV();
    auto angles = tissue.getSOAnglesCubic();
    for (auto c : angles)
        constraints.push_back(c);
    SO::Matrix3X points = tissue.grid.getVertexMatV();
    s->setPoints(points);
    s->clearConstraints();
    for (auto c : constraints)
        s->addConstraint(c);

    s->initialize();
    s->solve(epochs);
    tissue.grid.epochsRun += epochs;
    SO::Matrix3X newPoints = s->getPoints();
    tissue.updateGrid(newPoints);
}

void Simulation::runSO(int epochs, int agencyFlag) {
    epochs += 1;
    bool genAgentTD, runAgent;
    num currDiff = 0.0, divAvgMod = 0.0, initialDiff = 0.0;
    genAgentTD = agencyFlag == 1;
    runAgent = agencyFlag == 2;
    if (!genAgentTD && cubic) {
        runSOCubic(epochs);
        return;
    }
    auto cells = tissue.getCells();
    int steps;
    if (genAgentTD) {
        steps = epochs;
    } else {
        steps = 1;
    }
    int epochSteps = epochs / steps;
    auto start = std::chrono::high_resolution_clock::now();
    json td;
    std::unordered_map<int, vector<num>> edgeMap;
    std::unordered_map<int, vector<num>> cellMap;
    string layer;
    Cell::Neighbors<Edge> layerMap;
    string outString = fmt::format("k90-{:.2f}_k120-{:.2f}_kEdge-{:.2f}_{}",
            tissue.grid.kAng90, tissue.grid.kAng, tissue.grid.kEdge, time(0));
    vector<int> frames;
    vector<string> durations;
    if (genAgentTD && cubic) {
        tissue.grid.reset();
        runSO(epochs, 0);
        finalEdges = tissue.grid.getEdgeLengths(true);
        tissue.grid.reset();
        videoTxtFile.open(fmt::format("../exports/video/input.txt"));
        emit renderSignal();
        initialDiff = tissue.grid.calcEdgeLengthDiff(finalEdges);
        divAvg = initialDiff / static_cast<num>(tdEpochSteps);
        cout << initialDiff << endl;
        stepSize = divAvg;
        cutoff = divAvg * 0.6;
        int i;
        for (i=0; i<10; i++) {
            frames.push_back(i);
            durations.push_back("1.0");
            cout << i << ", 1.0\n";
        }
        int currStep = static_cast<int>((epochs - 10) * 0.05 / 100.0);
        cout << currStep << endl;
        for (i=i+currStep; i<=static_cast<int>(epochs * 0.05); i+=currStep) {
            frames.push_back(i);
            durations.push_back("0.05");
            cout << i << ", 0.05\n";
        }
        currStep = static_cast<int>((epochs - i) / 15.0);
        for (i+=currStep;i<epochs-1; i+=currStep) {
            frames.push_back(i);
            durations.push_back("1.0");
            cout << i << ", 1.0\n";
        }
        // For some reason ffmpeg ignores the last frame
        // so need to add a duplicate
        frames.push_back(epochs - 1);
        durations.push_back("10.0");
        frames.push_back(epochs - 1);
        durations.push_back("0.01");
        cout << "Save Number,Epoch,Target Change,Actual Change,Percent Difference\n";
    } else if (genAgentTD) {
        for (auto edge : tissue.grid.getEdges(true))
            edgeMap[edge->getIndex()] = vector<num>();
        for (auto cell : cells) {
            json cellJS;
            cellJS["index"] = cell->index;
            cellJS["center"] = cell->initialGridCenter;
            for (auto edge : cell->getEdges(true))
                cellJS["edges"].push_back(edge->getIndex());
            for (int i=0; i<2; i++) {
                if (i) {
                    layer = "apical";
                    layerMap = cell->edgeMap.apical;
                } else {
                    layer = "basal";
                    layerMap = cell->edgeMap.basal;
                }
                for (auto d : AllDirs) {
                    cellJS["edge_map"][layer][CDirStrings.at(d)] = (*layerMap.m[d])->getIndex();
                }
            }
            auto cs = cell->getNeighborsVec();
            for (auto neighbor : cs) {
                if (neighbor != nullptr) {
                    cellJS["neighbors"].push_back(neighbor->index);
                } else {
                    cellJS["neighbors"].push_back(-1);
                }
            }
            for (auto d : AllDirs) {
                auto n = *cell->neighbors.m[d];
                if (n != nullptr)
                    cellJS["neighbor_map"][CDirStrings.at(d)] = n->index;
                else
                    cellJS["neighbor_map"][CDirStrings.at(d)] = -1;
            }
            for (int i : cell->getAgentOrdering())
                cellJS["process_order"].push_back(i);
            td["cells"].push_back(cellJS);
        }
        tissue.grid.reset();
        runSO(epochs, 0);
        finalEdges = tissue.grid.getEdgeLengths(true);
        tissue.grid.reset();
        videoTxtFile.open(fmt::format("../exports/video/input.txt"));
        emit renderSignal();
        initialDiff = tissue.grid.calcEdgeLengthDiff(finalEdges);
        divAvg = initialDiff / static_cast<num>(tdEpochSteps);
        cout << initialDiff << endl;
        stepSize = divAvg;
        cutoff = divAvg * 0.6;
        int i;
        for (i=0; i<10; i++) {
            frames.push_back(i);
            durations.push_back("1.0");
            cout << i << ", 1.0\n";
        }
        int currStep = static_cast<int>((epochs - 10) * 0.05 / 100.0);
        cout << currStep << endl;
        for (i=i+currStep; i<=static_cast<int>(epochs * 0.05); i+=currStep) {
            frames.push_back(i);
            durations.push_back("0.05");
            cout << i << ", 0.05\n";
        }
        currStep = static_cast<int>((epochs - i) / 15.0);
        for (i+=currStep;i<epochs-1; i+=currStep) {
            frames.push_back(i);
            durations.push_back("1.0");
            cout << i << ", 1.0\n";
        }
        // For some reason ffmpeg ignores the last frame
        // so need to add a duplicate
        frames.push_back(epochs - 1);
        durations.push_back("10.0");
        frames.push_back(epochs - 1);
        durations.push_back("0.01");
        cout << "Save Number,Epoch,Target Change,Actual Change,Percent Difference\n";
    }
    unsigned int j, k;
    vector<num> prevEdges = tissue.grid.getEdgeLengths(),
        currEdges;
    td["rest_length_step"] = divAvg;
    td["cutoff"] = cutoff;
    int m=0;
    auto constraints = tissue.grid.getSOEdges();
    auto angles = tissue.getSOAngles();
    for (auto c : angles)
        constraints.push_back(c);
    SO::Matrix3X points = tissue.getVertexMat();
    s->setPoints(points);
    s->clearConstraints();
    for (auto c : constraints) {
        s->addConstraint(c);
    }

    s->initialize();
    int c=0;
    for (int i=0; i<steps; i++) {
        if (i || !genAgentTD) {
            if (runAgent) {
                for (j=0; j<cells.size(); j++) {
                    auto cell = cells[j];
                    json inp;
                    inp["cell"] = cell->index;
                    inp["epoch"] = tissue.grid.epochsRun;

                    auto cs = cell->getNeighborsVec();
                    cs.push_back(cell);
                    std::set<int> done;
                    for (auto neighbor : cs) {
                        if (neighbor != nullptr) {
                            inp["signal"].push_back({neighbor->index, neighbor->signal});
                            auto edges = neighbor->getEdges(true);
                            for (k=0; k<edges.size(); k++) {
                                auto edge = edges[k];
                                if (!done.contains(edge->getIndex())) {
                                    done.insert(edge->getIndex());
                                    json ee;
                                    ee["index"] = edge->getIndex();
                                    ee["length"] = edge->getCurrentLength();
                                    inp["edges"].push_back(ee);
                                }
                            }
                        }
                    }
                    json agentStats = ((Agent*)agentPointer)->statAgent();
                    stepSize = agentStats["step_size"];
                    cutoff = agentStats["cutoff"];
                    json agentOut = ((Agent*)agentPointer)->runAgent(inp);
                    json vals = agentOut["edges"];
                    if (agentOut["signal"] == 1)
                        cell->signal += 1;
                    else if (agentOut["signal"] == 2)
                        cell->signal -= 1;
                    cell->signal = min(max(cell->signal, 0), 10);
                    for (auto v : vals) {
                        int ind, action;
                        v["index"].get_to(ind);
                        v["action"].get_to(action);
                        tissue.grid.edgeMap[ind]->agentActions[cell] = action;
                    }
                }
                int total = 0, correct = 0, c;
                int totalC = 0, correctC = 0;
                int totalD = 0, correctD = 0;
                num diffs = 0.0;
                for (auto e : tissue.grid.getEdges(true)) {
                    num diff = e->getCurrentLength() - e->origRestLength;
                    diffs += diff;
                    int corrVal;
                    if (abs(diff) < cutoff) {
                        corrVal = 0;
                    } else if (diff < 0) {
                        corrVal = 1;
                    } else {
                        corrVal = 2;
                    }
                    e->correctAction = corrVal;
                    num sSMod = stepSize, rlMod = 0.0;
                    if (e->agentActions.size() == 2)
                        sSMod /= 2.0;
                    for (auto kv : e->agentActions) {
                        total++;
                        c = kv.second == corrVal;
                        correct += c;
                        if (corrVal) {
                            totalC++;
                            correctC += c;
                        } else {
                            totalD++;
                            correctD += c;
                        }
                        switch (kv.second) {
                            case 0:
                                rlMod = 0.0;
                                break;
                            case 1:
                                rlMod = sSMod;
                                break;
                            case 2:
                                rlMod = -sSMod;
                                break;
                            default:
                                cout << "Error: invalid action when updating rest lengths (" << kv.second << endl;
                                exit(1);
                        }
                        e->setRestLength(e->getCurrentLength()+rlMod);
                    }
                }
                cout << fmt::format("Don't change correct: {} / {} ({})\n", correctD, totalD, static_cast<num>(correctD) / totalD);
                cout << fmt::format("Change correct: {} / {} ({})\n", correctC, totalC, static_cast<num>(correctC) / totalC);
                cout << fmt::format("Overall: {} / {} ({})\n", correct, total, static_cast<num>(correct) / total);
                cout << fmt::format("Average Diff: {}\n", diffs / static_cast<num>(total));
                constraints = tissue.grid.getSOEdges();
                angles = tissue.getSOAngles();
                for (auto c : angles)
                    constraints.push_back(c);
                points = tissue.getVertexMat();
                s->setPoints(points);
                s->clearConstraints();
                for (auto c : constraints)
                    s->addConstraint(c);

                s->initialize();
            }
            s->solve(runAgent ? 1000 : epochSteps);
            tissue.grid.epochsRun += epochSteps;
            SO::Matrix3X newPoints = s->getPoints();

            tissue.updateGrid(newPoints);
        }
        if (genAgentTD) {
            if (tissue.grid.epochsRun == static_cast<uint>(frames[c])) {
                string durationStr = durations[c++];
                videoTxtFile << fmt::format(
                        "file ../temp/{:05d}.png\nduration {}\n",
                        tissue.grid.epochsRun, durationStr); 
                emit renderSignal();
            }
            currDiff = initialDiff - tissue.grid.calcEdgeLengthDiff(finalEdges);
            if (!i || currDiff > divAvgMod || i+1 == steps) {
                td["epochs"].push_back(tissue.grid.epochsRun);
                for (auto e : tissue.grid.getEdges()) {
                    if (!e->isVertical) {
                        edgeMap[e->getIndex()].push_back(e->getCurrentLength());
                    }
                }
                for (auto cell : tissue.getCells()) {
                    cellMap[cell->index].push_back(cell->getVolume());
                }
                prevEdges = tissue.grid.getEdgeLengths();
                cout << fmt::format("{},{},{:.6f},{:.6f},{:.2f}%\n", m+1,
                        tissue.grid.epochsRun, divAvgMod, currDiff,
                        (currDiff-divAvgMod)/divAvgMod*100.0);
                initialDiff = tissue.grid.calcEdgeLengthDiff(finalEdges);
                divAvgMod = tissue.grid.calcEdgeLengthDiff(finalEdges) /
                    (static_cast<num>(tdEpochSteps - m++));
            }
        }
    }
    if (genAgentTD && cubic) {
        videoTxtFile.close();
        cout << endl;
        string commandS = fmt::format("cd ../exports/video; ffmpeg -f concat -safe 0 -i input.txt -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' -c:v libvpx-vp9 -c:a libopus -y -an {}.mp4; echo 'Video {} rendered\n'", outString, outString);
        cout << "Avg distance from rest length: " << tissue.grid.calcEdgeLengthDiff(finalEdges) << endl;
        system(commandS.c_str());
    } else if (genAgentTD) {
        for (auto kv : edgeMap) {
            json ee;
            auto edge = tissue.grid.edgeMap[kv.first];
            ee["index"] = kv.first;
            ee["target"] = edge->origRestLength;
            Vector3d p1, p2;
            p1 = edge->getV1()->originalPosition;
            p2 = edge->getV2()->originalPosition;
            std::vector<num> v1(p1.data(), p1.data() + p1.rows() * p1.cols());
            std::vector<num> v2(p2.data(), p2.data() + p2.rows() * p2.cols());
            ee["v1"] = v1;
            ee["v2"] = v2;
            for (auto l : kv.second) {
                ee["lengths"].push_back(l);
                edge->setAction(l, cutoff);
            }
            td["edges"].push_back(ee);
        }
        for (auto kv : cellMap) {
            // TODO: could be made a lot more efficient
            for (unsigned int i=0; i<td["cells"].size(); i++) {
                if (td["cells"][i]["index"] == kv.first) {
                    for (auto v : kv.second) {
                        td["cells"][i]["volumes"].push_back(v);
                    }
                    break;
                }
            }
        }
        std::ofstream outFile("training_data.json");
        outFile << td;
        outFile.close();
        videoTxtFile.close();
        cout << endl;
        string commandS = fmt::format("cd ../exports/video; ffmpeg -f concat -safe 0 -i input.txt -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' -c:v libvpx-vp9 -c:a libopus -y -an {}.mp4; echo 'Video {} rendered\n'", outString, outString);
        cout << "Avg distance from rest length: " << tissue.grid.calcEdgeLengthDiff(finalEdges) << endl;
        system(commandS.c_str());
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    auto dur_num = static_cast<num>(duration.count()) / 1e6;
    secPerEpoch = dur_num / epochs;
}

void Simulation::runHeadless() {
    if (!setts["simulation"]["batch"]) {
        model.kLinearInMin = model.kLinearIn;
        model.kLinearInMax = model.kLinearIn;
        model.kLinearStep = model.kLinearIn;
        model.kAngIn120Max = model.kAngIn120;
        model.kAngIn120Min = model.kAngIn120;
        model.kAng120Step = model.kAngIn120;
        model.kAngIn90Max = model.kAngIn90;
        model.kAngIn90Min = model.kAngIn90;
        model.kAng90Step = model.kAngIn90;
    }
    int totalLin = static_cast<int>((model.kLinearInMax - model.kLinearInMin) / model.kLinearStep) + 1;
    int total120 = static_cast<int>((model.kAngIn120Max - model.kAngIn120Min) / model.kAng120Step) + 1;
    int total90 = static_cast<int>((model.kAngIn90Max - model.kAngIn90Min) / model.kAng90Step) + 1;
    int startTime = time(0);
    int steps = setts["simulation"]["epoch_steps"].value_or(500);

    string batchDir = fmt::format("../exports/batch_{}", startTime);
    std::ofstream csvFile;

    json output;
    output["scale_factor"] = model.scaleFactor;

    for (auto e : tissue.grid.getEdges()) {
        json entry;
        entry["rest_length"] = e->getRestLength();
        entry["vertices"] = {e->getV1()->getIndex(), e->getV2()->getIndex()};
        output["edges"].push_back(entry);
    }

    if (cubic) {
        runSO();
        tissue.addLayersCubic();
        for (auto v : tissue.grid.getVertices())
            v->originalPosition << v->getPosition();
    }

    for (auto e : tissue.grid.getEdges())
        e->origRestLength = e->getRestLength();

    for (int i=0; i<totalLin; i++) {
        for (int j=0; j<total120; j++) {
            for (int k=0; k<total90; k++) {
                tissue.grid.reset();
                tissue.grid.kEdge = model.kLinearInMin + i * model.kLinearStep;
                tissue.grid.kAng = model.kAngIn120Min + j * model.kAng120Step;
                tissue.grid.kAng90 = model.kAngIn90Min + k * model.kAng90Step;

                for (int w=0; w<setts["simulation"]["epochs"].value_or(10000); w+=steps) {
                    cout << "\tEpoch: " << w << endl;
                    runSO(steps, false);
                }
                if (!i && !j && !k) {
                    std::filesystem::create_directories(batchDir);
                    csvFile.open(fmt::format("../exports/batch_{}/results.csv", startTime));
                    if (!csvFile.good()) {
                        cout << "ERROR: CSV output file cannot be created\n";
                        exit(1);
                    }
                    csvFile << "kLinear,k90,k120,Error\n";
                    csvFile.flush();
                }
                csvFile << fmt::format("{:.3f},{:.3f},{:.3f},{:.6f}\n", tissue.grid.kEdge, tissue.grid.kAng90, tissue.grid.kAng, evaluateChamfer().avgChamf);
                csvFile.flush();
            }
        }
    }
}

EvalResults &Simulation::evaluateChamfer() {
    getAvgSpringDisp();
    auto ev = new EvalResults();
    auto cells = getTissue().getCells();
    Eigen::Vector4i corners;
    ev->actor = triangulateGrid(corners);

    vtkSmartPointer<vtkPoints> points, points2;
    MatrixXd vertMat1 = MatrixXd::Zero(0, 3);
    MatrixXd vertMat2 = MatrixXd::Zero(0, 3);

    auto mesh = dynamic_cast<vtkPolyData*>(ev->actor->GetMapper()->GetInput());
    points = accessEachVertex(mesh, vertMat1);

    vtkNew<vtkTransform> transform;
    transform->Scale(model.scale, model.scale, model.scale);
    transform->TransformPoints(points, points);

    vtkNew<vtkOBJImporter> importer;
    string inputFileName;
    if (!static_cast<string>(setts["model"]["type"].value_or("obj")).compare("obj"))
        inputFileName = setts["model"]["model"].value_or("../obj_models/model.obj");
    else
        inputFileName = "./temp.obj";
    importer->SetFileName(inputFileName.c_str());
    importer->Update();
    vtkActorCollection *allActors = importer->GetRenderer()->GetActors();
    ev->modelActor = allActors->GetLastActor();
    ev->modelActor->SetScale(model.scale);
    ev->modelActor->GetProperty()->LightingOff();
    auto mesh2 = dynamic_cast<vtkPolyData*>(ev->modelActor->GetMapper()->GetInput());
    points2 = accessEachVertex(mesh2, vertMat2);

    vector<Vector3d> vertVecSrc, vertVecDst;
    Vector3d midDst, midSrc, translation, elemsSum;
    elemsSum << 0.0, 0.0, 0.0;
    for (Vector3d r : vertMat1.rowwise()) {
        elemsSum += r;
        vertVecSrc.push_back(r);
    }
    auto s = static_cast<num>(vertVecSrc.size());
    midSrc << elemsSum(0)/s, elemsSum(1)/s, elemsSum(2)/s;

    elemsSum << 0.0, 0.0, 0.0;
    int a = 50;
    MatrixXd p2ds = MatrixXd::Zero(a*a, 2), p3ds;
    vector<string> colors;
    for (int i=0; i<a; i++) {
        for (int j=0; j<a; j++) {
            int ind = i*a+j;
            p2ds.row(ind) << bound(i / static_cast<num>(a)), bound(j / static_cast<num>(a));
            if (!ind || ind == a-1 || ind == a*(a-1) || ind == a*a-1) {
                colors.push_back(" Green");
            } else {
                colors.push_back(" Blue");
            }
        }
    }
    p3ds = model.mesh->barycentricLoopingF(p2ds);
    elemsSum = p3ds.colwise().sum();
    s = static_cast<num>(a*a);
    midDst << elemsSum(0)/s, elemsSum(1)/s, elemsSum(2)/s;
    Eigen::Matrix3Xd srcC(3,4);
    if (cubic) {
        srcC <<
            vertVecSrc[corners(0)].transpose(),
            vertVecSrc[corners(1)].transpose(),
            vertVecSrc[corners(2)].transpose(),
            vertVecSrc[corners(3)].transpose();
    } else {
        for (int i=0; i<4; i++) {
            auto cVerts = cells[corners(i)]->getVertices();
            srcC.col(i) << 0.0, 0.0, 0.0;
            for (long unsigned int j=0; j<cVerts.size(); j++)
                srcC.col(i) += cVerts[j]->getPosition();
            srcC.col(i) /= cVerts.size();
        }
    }
    
    Eigen::Matrix3Xd dstC(3,4);
    dstC << p3ds.row(0).transpose(),
        p3ds.row(a-1).transpose(),
        p3ds.row(a*(a-1)).transpose(),
        p3ds.row(a*a-1).transpose();

    Eigen::Matrix4d T = Eigen::umeyama(srcC, dstC, false);

    Eigen::Matrix3d R = T.block<3, 3>(0, 0);
    Eigen::Vector3d t = T.block<3, 1>(0, 3);

    SO::Matrix3X src(3, vertVecSrc.size()), dst = p3ds.transpose();

    Matrix3d rotTemp = R;

    translation = t.transpose();

    for (long unsigned int i=0; i<vertVecSrc.size(); i++) {
        src.col(i) = rotTemp * vertVecSrc[i].transpose() + translation.transpose();
    }

    MatrixXd srct, dstt;
    srct = src.transpose();
    dstt = dst.transpose();

    Pointcloud pointcloudA = dst;
    Pointcloud pointcloudB = src;

    auto callback = Callback(pointcloudA, pointcloudB);

    lsqcpp::GaussNewton<num, 6, Eigen::Dynamic, Objective, lsqcpp::ArmijoBacktracking, lsqcpp::DenseCholeskySolver> optimizer;
    optimizer.setMinimumStepLength(0.015);
    optimizer.setObjective({pointcloudA, pointcloudB});
    optimizer.setCallback(callback);
    optimizer.setMaximumIterations(25);
    optimizer.setVerbosity(1);

    Vector6 xval(6);
    xval.setZero();
    xval.segment(3, 3) = lsqcpp::parameter::encodeRotation(Matrix3d::Identity());

    optimizer.minimize(xval);
    ev->result = MatrixXd(src.cols(), src.rows());
    for(int i = 0; i < src.cols(); i++)
    {
        ev->result.row(i) = (rotFinal * src.col(i) + transFinal).transpose();
    }
    ev->dst2 = dst.transpose();

    vector<num> errVals, linErr, angErr90, angErr120;
    ev->error = 0.0;
    ev->maxD = 0.0;
#pragma omp parallel for
    for (Vector3d r1: ev->result.rowwise()) {
        num minD = std::numeric_limits<num>::max(), currD;
        for (Vector3d r2 : ev->dst2.rowwise()) {
            currD = distance(r1, r2);
            if (currD < minD)
                minD = currD;
        }
#pragma omp critical
        {
            errVals.push_back(minD);
            ev->error += minD;
            if (minD > ev->maxD)
                ev->maxD = minD;
        }
    }


    /*for (auto e : tissue.grid.getEdges()) {
        if (e->isVertical) {
            for (auto eh : tissue.grid.getEdges(true)) {
                Vertex *v1, *v2, *h1, *h2, *a1, *a2, *a3;
                bool set = false;
                v1 = e->getV1();
                v2 = e->getV2();
                h1 = eh->getV1();
                h2 = eh->getV2();

                if (v1->getIndex() == h1->getIndex()) {
                    a1 = v2;
                    a2 = v1;
                    a3 = h2;
                    set = true;
                }
                if (v1->getIndex() == h2->getIndex()) {
                    a1 = v2;
                    a2 = v1;
                    a3 = h1;
                    set = true;
                }
                if (v2->getIndex() == h1->getIndex()) {
                    a1 = v1;
                    a2 = v2;
                    a3 = h2;
                    set = true;
                }
                if (v2->getIndex() == h2->getIndex()) {
                    a1 = v1;
                    a2 = v2;
                    a3 = h1;
                    set = true;
                }

                if (set) {
                    Vector3d one = (a1->getPosition() - a2->getPosition()).normalized();
                    Vector3d two = (a3->getPosition() - a2->getPosition()).normalized();
                    //num  += abs(acos(zudir.dot(ldir)) - H_PI) / H_PI;
                }
            }
        }
    }*/

    Vector3d mmin, mmax, bbox;
    mmin = ev->result.colwise().minCoeff();
    mmax = ev->result.colwise().maxCoeff();
    bbox = (mmax - mmin).cwiseAbs();
    ev->avgChamf = ev->error / ev->result.rows();
    cout << "Bounding Box Size: " << bbox << " = " << bbox[0]*bbox[1]*bbox[2] << endl;
    cout << "Surface area: " << model.surfaceArea << endl;
    cout << "Average Chamfer Distance: " << ev->avgChamf << endl;
    cout << "Max Chamfer Distance: " << ev->maxD << endl;
    cout << "Chamfer Distance StDev: " << stdev(errVals) << endl;

    return *ev;
}


vtkSmartPointer<vtkActor> Simulation::triangulateGridCubic(Eigen::Vector4i &corners) {
    vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::New();
    vtkNew<vtkCellArray> triangles;
    auto verts = tissue.grid.getVerticesV();
    int lVerts = verts.size() / model.layersCubic;
    std::ofstream o("cubic_triangulation.obj");
    int rows = model.nodeCountU, cols = model.nodeCountV;
    corners << 0, cols-3, lVerts-cols-3, lVerts-1-3;
    for (int i=0; i<lVerts; i++) {
        Vector3d p = verts[i]->getPosition();
        for (uint j=1; j<model.layersCubic; j++) {
            p += verts[i+j*lVerts]->getPosition();
        }
        p /= model.layersCubic;
        points->InsertNextPoint(p(0), p(1), p(2));
        o << "v " << p << endl;
    }
    for (int i = 0; i < rows-1; i++) {
        for (int j = 0; j < cols-1; j++) {
            vtkSmartPointer<vtkTriangle> triangle =
                vtkSmartPointer<vtkTriangle>::New();
            Vector3i f = {i * cols + j,
                          (i + 1) * cols + j,
                          i * cols + j + 1};
            triangle->GetPointIds()->SetId(0, f(0));
            triangle->GetPointIds()->SetId(1, f(1));
            triangle->GetPointIds()->SetId(2, f(2));
            o << fmt::format("f {0}//{0} {1}//{1} {2}//{2}\n", f(0)+1, f(1)+1, f(2)+1);
            triangles->InsertNextCell(triangle);
            Vector3i f2 = {i * cols + j + 1,
                          (i + 1) * cols + j,
                          (i + 1) * cols + j + 1};
            vtkSmartPointer<vtkTriangle> triangle2 =
                vtkSmartPointer<vtkTriangle>::New();
            // Center
            triangle2->GetPointIds()->SetId(0, f2(0));
            // Edge points
            triangle2->GetPointIds()->SetId(1, f2(1));
            triangle2->GetPointIds()->SetId(2, f2(2));
            triangles->InsertNextCell(triangle2);
            o << fmt::format("f {0}//{0} {1}//{1} {2}//{2}\n", f2(0)+1, f2(1)+1, f2(2)+1);
        }
    }
    o.close();
    // Create a polydata object
    vtkSmartPointer<vtkPolyData> polyData =
        vtkSmartPointer<vtkPolyData>::New();

    // Add the geometry and topology to the polydata
    polyData->SetPoints(points);
    polyData->SetPolys(triangles);
    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputData(polyData);

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    return actor;
}

vtkSmartPointer<vtkActor> Simulation::triangulateGrid() {
    Eigen::Vector4i corners;
    return triangulateGrid(corners);
}

vtkSmartPointer<vtkActor> Simulation::triangulateGrid(Eigen::Vector4i &corners) {
    if (cubic)
        return triangulateGridCubic(corners);
    vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::New();
    vtkNew<vtkCellArray> triangles;
    int startInd = 0;
    int cRows = model.vCellCentersCart.size(), l, cInd;
    corners << std::numeric_limits<int>::min(),
            std::numeric_limits<int>::max(), std::numeric_limits<int>::max(),
            std::numeric_limits<int>::min();
    auto cells = tissue.getCells();
    for (unsigned int i=0; i<cells.size(); i++) {
    //for (auto cell : tissue.getCells()) {
        auto cell = cells[i];
        l = cell->level;
        cInd = cell->vSplineIndex;
        if (cInd == 1) {
            if (l > corners(0))
                corners(0) = i;
            if (l < corners(2))
                corners(2) = i;
        } else if (cInd == cRows - 2) {
            if (l < corners(1))
                corners(1) = i;
            if (l > corners(3))
                corners(3) = i;
        }
        MatrixXd positions = MatrixXd::Zero(6, 3);
        auto verts = cell->getVertices();
        Vector3d center = Vector3d::Zero();
        for (int i=0; i<6; i++) {
            positions.row(i) = verts[i]->getPosition() + verts[i+6]->getPosition();
            positions.row(i) /= 2.0;
            center += positions.row(i);
            points->InsertNextPoint(positions(i, 0), positions(i, 1), positions(i, 2));
        }
        center /= 6.0;
        points->InsertNextPoint(center(0), center(1), center(2));
        for (int i=0; i<6; i++) {
            vtkSmartPointer<vtkTriangle> triangle =
                vtkSmartPointer<vtkTriangle>::New();
            // Center
            triangle->GetPointIds()->SetId(0, startInd+6);
            // Edge points
            triangle->GetPointIds()->SetId(1, startInd+i);
            triangle->GetPointIds()->SetId(2, startInd+((i+1)%6));
            triangles->InsertNextCell(triangle);
        }
        startInd += 7;
    }
    for (int i=0; i<4; i++)
        cout << cells[corners(i)]->initialCenterUV << endl;
    cout << endl;

    // Create a polydata object
    vtkSmartPointer<vtkPolyData> polyData =
        vtkSmartPointer<vtkPolyData>::New();

    // Add the geometry and topology to the polydata
    polyData->SetPoints(points);
    polyData->SetPolys(triangles);
    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputData(polyData);

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    return actor;
}

// Returns how far linear and angular springs are from their rest lengths
Eigen::Array2d Simulation::getAvgSpringDisp() {
    Eigen::Array2d disp;
    disp << 0.0, 0.0;
    Eigen::Array2i count;
    count << 0, 0;
    vector<num> linErr;
    num err = 0;
    for (auto e : tissue.grid.getEdges(true)) {
        linErr.push_back(abs((e->getRestLength() - e->getCurrentLength()) / e->getRestLength()));
    }
    for (num e : linErr)
        err += e;
    cout << "Edge length error percent: " << err / linErr.size() * 100 << "%\n";
    cout << "Edge length stdev: " << stdev(linErr) << endl;
    disp[1] /= count[1];
    tissue.getSOAngles();
    return disp;
}
