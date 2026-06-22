#include "model.h"

Model::Model() {
    Model(false);
}

Model::Model(bool cubic) {
    // TODO: Deal with processing the config values differently (ks, grid sizing)
    startingLineNodeCount = setts["isolines"]["initial_node_count"].value_or(50);
    countToBump = setts["isolines"]["spline_bump_count"].value_or(0);
    startingU = setts["isolines"]["primary_u_value"].value_or(0.5);
    uOffset = setts["isolines"]["u_spline_offset"].value_or(0.0);
    vOffset = setts["isolines"]["v_spline_offset"].value_or(0.0);
    layersCubic = setts["grid"]["layers"].value_or(3);
    batch = setts["simulation"]["batch"].value_or(false);
    hasDuplicatedCubic = false;
    num zh = setts["grid"]["z_height"].value_or(0.0);
    num zhm = setts["grid"]["z_height_mult"].value_or(0.0);
    num zhd = setts["grid"]["z_height_div"].value_or(0.0);
    spacingDiv = true;
    if (zh) {
        zHeight = zh;
        spacingDiv = false;
    } else if (zhm) {
        zHeight = 1.0 / zhm;
    } else if (zhd) {
        zHeight = zhd;
    } else {
        zHeight = 2.0;
    }
    processK(setts["grid"]["k_ang_120"], kAngIn120, kAngIn120Min, kAngIn120Max, kAng120Step, K_ANG_120);
    processK(setts["grid"]["k_ang_90"], kAngIn90, kAngIn90Min, kAngIn90Max, kAng90Step, K_ANG_90);
    processK(setts["grid"]["k_linear"], kLinearIn, kLinearInMin, kLinearInMax, kLinearStep, K_EDGE);
    epochsIn = setts["simulation"]["epochs"].value_or(EPOCHS_PER_CYCLE);
    epochsInStep = setts["simulation"]["epoch_step_size"].value_or(epochsIn);

    scale = setts["model"]["scale"].value_or(1.0);

    modelFile = setts["model"]["model"].value_or("../obj_models/model.obj");
    string inputType = setts["model"]["type"].value_or("obj");

    mesh = new Mesh(modelFile, inputType);
    mesh->scale = setts["model"]["scale"].value_or(1.0);
    mesh->flipPrimary = setts["model"]["flip_primary"].value_or(0);

    if (!inputType.compare("raytrace") || !inputType.compare("bezier")) {
        mesh->triRows = setts["triangulation"]["x"].value_or(20);
        mesh->triCols = setts["triangulation"]["y"].value_or(20);
        modelFile = "./temp.obj";
        if (!inputType.compare("raytrace")) {
            auto rb1 = setts["ray_trace"]["ray_start"];
            auto rb2 = setts["ray_trace"]["ray_box_size"];
            mesh->rayBox << rb1[0].value_or(0.0), rb1[1].value_or(0.0), rb2[0].value_or(10.0), rb2[1].value_or(10.0);
            mesh->rtsettings->scale = setts["ray_trace"]["scale"].value_or(1.0);
            auto rot = setts["ray_trace"]["rotate"];
            mesh->rtsettings->rotate << rot[0].value_or(0), rot[1].value_or(0), rot[2].value_or(0);
            auto Zv = setts["ray_trace"]["Zv"];
            mesh->rtsettings->Zv << Zv[0].value_or(1), Zv[1].value_or(0), Zv[2].value_or(0);
            auto Vup = setts["ray_trace"]["Vup"];
            mesh->rtsettings->Vup << Vup[0].value_or(0), Vup[1].value_or(1), Vup[2].value_or(0);
        }
    } else {
        auto pc = setts["model"]["param_corners"];
        mesh->paramCorners << pc[0].value_or(0), pc[1].value_or(0),
                     pc[2].value_or(0), pc[3].value_or(0);
    }

    mesh->initializeMesh();

    surfaceArea = 0.0;
    for (auto f : mesh->faces) {
        Vector3d p1 = mesh->realPoints.row(f->vertices(0));
        Vector3d p2 = mesh->realPoints.row(f->vertices(1));
        Vector3d p3 = mesh->realPoints.row(f->vertices(2));
        surfaceArea += scaleneArea(p1, p2, p3);
    }

    if (setts["global"]["import_from_json"].value_or(false)) {
        return;
    }

    this->cubic = cubic;
    if (cubic)
        finishInitCubic();
    else
        finishInit();
}

void Model::finishInitCubic() {
    num maxU = 0.0, maxV = 0.0, zHeightBackup = 0.0, zHeightDivisor = zHeight, spacingBackup = 0.0;
    splines = vector<Spline*>(startingLineNodeCount);
    usplines = vector<Spline*>(startingLineNodeCount);
    for (int i=0; i<startingLineNodeCount; i++) {
        num uv = static_cast<num>(i) / (startingLineNodeCount - 1);
        auto splineU = new Spline(true, true, uv, 0.0);
        auto splineV = new Spline(false, true, uv, 0.0);
        splineU->sType = PRIMARY;
        splineV->sType = PRIMARY;
        populateSpline(*splineU);
        if (spacing > spacingBackup) {
            spacingBackup = spacing;
            zHeightBackup = zHeight;
        }
        zHeight = zHeightDivisor;
        populateSpline(*splineV);
        if (spacing > spacingBackup) {
            spacingBackup = spacing;
            zHeightBackup = zHeight;
        }
        zHeight = zHeightDivisor;
        splines[i] = splineV;
        usplines[i] = splineU;
        num distU = 0.0, distV = 0.0;
        for (int i=0; i<splineU->splinePoints.rows()-1; i++) {
            distU += distance(splineU->splinePoints.row(i), splineU->splinePoints.row(i+1));
            distV += distance(splineV->splinePoints.row(i), splineV->splinePoints.row(i+1));
        }
        // NOTE: I don't know why but the V and U need flipped for this to work right,
        // still trying to figure that one out
        maxU = std::max(distV, maxU);
        maxV = std::max(distU, maxV);
    }
    zHeight = zHeightBackup;
    num aspectRatio = maxU / maxV;
    if (aspectRatio > 1.0) {
        nodeCountU = startingLineNodeCount;
        nodeCountV = ceil(nodeCountU / aspectRatio);
        splines = vector<Spline*>(nodeCountV);
        for (int i=0; i<nodeCountV; i++) {
            num uv = static_cast<num>(i) / (nodeCountV - 1);
            auto splineV = new Spline(false, false, uv, 0.0);
            splineV->sType = PRIMARY;
            populateSpline(*splineV);
            vector<Vector3d> normals;
            int steps = splineV->splinePoints.rows() / splineV->realPoints.rows();
            for (int k=0; k<splineV->realPoints.rows()-1; k++) {
                Vector3d normal, tCopy = splineV->splinePoints.row(steps / 2 + steps * k);
                mesh->barycentricLoopingR(tCopy, normal);
                normals.push_back(normal);
            }
            splineV->curvatures = fullSplineSignedCurvatureList(
                    &splineV->realPoints, &splineV->tangents,
                    normals,
                    static_cast<int>(splineV->splinePoints.rows()));
            splines[i] = splineV;
        }
        for (int i=0; i<nodeCountU; i++) {
            auto splineU = usplines[i];
            vector<Vector3d> normals;
            int steps = splineU->splinePoints.rows() / splineU->realPoints.rows();
            for (int k=0; k<splineU->realPoints.rows()-1; k++) {
                Vector3d normal, tCopy = splineU->splinePoints.row(steps / 2 + steps * k);
                mesh->barycentricLoopingR(tCopy, normal);
                normals.push_back(normal);
            }
            splineU->curvatures = fullSplineSignedCurvatureList(
                    &splineU->realPoints, &splineU->tangents,
                    normals,
                    static_cast<int>(splineU->splinePoints.rows()));
        }
    } else {
        nodeCountV = startingLineNodeCount;
        nodeCountU = ceil(nodeCountV * aspectRatio);
        usplines = vector<Spline*>(nodeCountU);
        for (int i=0; i<nodeCountU; i++) {
            num uv = static_cast<num>(i) / (nodeCountU - 1);
            auto splineU = new Spline(true, false, uv, 0.0);
            splineU->sType = PRIMARY;
            populateSpline(*splineU);
            vector<Vector3d> normals;
            int steps = splineU->splinePoints.rows() / splineU->realPoints.rows();
            for (int k=0; k<splineU->realPoints.rows()-1; k++) {
                Vector3d normal, tCopy = splineU->splinePoints.row(steps / 2 + steps * k);
                mesh->barycentricLoopingR(tCopy, normal);
                normals.push_back(normal);
            }
            splineU->curvatures = fullSplineSignedCurvatureList(
                    &splineU->realPoints, &splineU->tangents,
                    normals,
                    static_cast<int>(splineU->splinePoints.rows()));
            usplines[i] = splineU;
        }
        for (int i=0; i<nodeCountV; i++) {
            auto splineV = splines[i];
            vector<Vector3d> normals;
            int steps = splineV->splinePoints.rows() / splineV->realPoints.rows();
            for (int k=0; k<splineV->realPoints.rows()-1; k++) {
                Vector3d normal, tCopy = splineV->splinePoints.row(steps / 2 + steps * k);
                mesh->barycentricLoopingR(tCopy, normal);
                normals.push_back(normal);
            }
            splineV->curvatures = fullSplineSignedCurvatureList(
                    &splineV->realPoints, &splineV->tangents,
                    normals,
                    static_cast<int>(splineV->splinePoints.rows()));
        }
    }
    cout << "U length: " << maxU << endl;
    cout << "V length: " << maxV << endl;
    cout << "Aspect ratio: " << std::min(maxU, maxV) / std::max(maxU, maxV) << endl;
    cout << fmt::format("Node count (u, v): ({}, {})\n", nodeCountU-1, nodeCountV-1);
    setScaleFactor();
}

void Model::finishInit() {
    primarySpline = Spline(true, true, startingU, uOffset);
    primarySpline.sType = PRIMARY;
    populateSpline(primarySpline);
    primaryUVPoints = primarySpline.uvPoints;
    primaryRealPoints = primarySpline.realPoints;
    primarySplinePoints = primarySpline.splinePoints;
    primaryTangents = primarySpline.tangents;
    primaryCellCentersCart = primarySpline.cellCentersCart;
    primaryCellCentersUV = primarySpline.cellCentersUV;

    MatrixXd uvPointsV, tempRealPointsVF, tempRealPointsVR, tempSPointsF, tempSPointsR, tempTanF, tempTanR, tempCCCF, tempCCCR, tempCCVF, tempCCVR;

    int midpt = startingLineNodeCount / 2;
    vector<vector<int>> toDo;
    vector<int> midptv, left, right;
    midptv.push_back(midpt);
    for (int i=midpt-1; i>=0; i--)
        left.push_back(i);
    for (int i=midpt+1; i<startingLineNodeCount; i++)
        right.push_back(i);
    toDo.push_back(midptv);
    toDo.push_back(left);
    toDo.push_back(right);
    tqdm bar;
    bar.set_label("Initializing Splines");

    int c=0;
    bar.progress(c, startingLineNodeCount);
    Spline midSpline, otherSpline;
    Spline *spline;

    vRealPoints = vector<MatrixXd*>(startingLineNodeCount);
    vSplines = vector<MatrixXd*>(startingLineNodeCount);
    vTangents = vector<MatrixXd*>(startingLineNodeCount);
    vCellCentersCart = vector<MatrixXd*>(startingLineNodeCount);
    vCellCentersUV = vector<MatrixXd*>(startingLineNodeCount);
    vCellCentersNormal = vector<MatrixXd*>(startingLineNodeCount);
    splines = vector<Spline*>(startingLineNodeCount);

    bumpTotal = countToBump * 2;
    bumpNum = 1;
    for (auto stack : toDo) {
        if (stack.size() != 1) {
            otherSpline = midSpline;
        }
        for (int i : stack) {
            bool shift = (i + 1) % 2;
            auto realPointsV = new MatrixXd;
            auto tangentsV = new MatrixXd;
            auto splinePointsV = new MatrixXd;
            auto cellCentersVCart = new MatrixXd;
            auto cellCentersVUV = new MatrixXd;
            spline = new Spline(false, false, CENTER_OUT, shift, startingU, primaryCellCentersUV(i,1), vOffset, 1.0);
            spline->bump = i < countToBump || i >= startingLineNodeCount - countToBump;
            spline->sType = VERTICAL;
            spline->index = i;
            populateSpline(*spline, otherSpline);
            if (spline->bump && !otherSpline.bump) {
                auto realPointsVi = new MatrixXd;
                auto tangentsVi = new MatrixXd;
                auto splinePointsVi = new MatrixXd;
                auto cellCentersVCarti = new MatrixXd;
                auto cellCentersVUVi = new MatrixXd;
                (*realPointsVi) = otherSpline.realPoints;
                (*tangentsVi) = otherSpline.tangents;
                (*splinePointsVi) = otherSpline.splinePoints;
                (*cellCentersVCarti) = otherSpline.cellCentersCart;
                (*cellCentersVUVi) = otherSpline.cellCentersUV;
                vRealPoints[otherSpline.index] = realPointsVi;
                vSplines[otherSpline.index] = splinePointsVi;
                vTangents[otherSpline.index] = tangentsVi;
                vCellCentersCart[otherSpline.index] = cellCentersVCarti;
                vCellCentersUV[otherSpline.index] = cellCentersVUVi;
            }
            (*realPointsV) = spline->realPoints;
            (*tangentsV) = spline->tangents;
            (*splinePointsV) = spline->splinePoints;
            (*cellCentersVCart) = spline->cellCentersCart;
            (*cellCentersVUV) = spline->cellCentersUV;
            vRealPoints[i] = realPointsV;
            vSplines[i] = splinePointsV;
            splines[i] = spline;

            if (stack.size() == 1) {
                midSpline = *spline;
            } else {
                otherSpline = *spline;
            }
            vTangents[i] = tangentsV;
            vCellCentersCart[i] = cellCentersVCart;
            vCellCentersUV[i] = cellCentersVUV;
            c++;
            bar.progress(c, startingLineNodeCount);
        }
    }
    bar.reset();

    vector<std::deque<Vector3d>*> ccs;
    vector<bool> doneUp, doneDown;
    MatrixXi indicies = MatrixXi::Zero(startingLineNodeCount, 2);
    for (int z = 0; z < startingLineNodeCount; z++) {
        ccs.push_back(new std::deque<Vector3d>);
        levels.push_back(new std::deque<int>);
        indexList.push_back(new std::deque<int>);
        doneUp.push_back(false);
        doneDown.push_back(false);
    }
    int index;
    bool splineEnded;
    for (int z = 0; z < startingLineNodeCount; z+=2) {
        Vector3d pc = primaryCellCentersCart.row(z);
        ccs[z]->push_front(findNextCC(vSplines[z], pc, index, spacing / 2, false, splineEnded));
        levels[z]->push_front(-1);
        indexList[z]->push_front(index);
        indicies(z, 0) = index;
        ccs[z]->push_back(findNextCC(vSplines[z], index, spacing, true, splineEnded));
        levels[z]->push_back(1);
        indexList[z]->push_back(index);
        indicies(z, 1) = index;
        if (z + 1 < startingLineNodeCount) {
            ccs[z+1]->push_back(primaryCellCentersCart.row(z+1));
            levels[z+1]->push_front(0);
            if (splines[z+1]->bump) {
                // TODO: Make this a function and reuse here and in bumpPoint
                num distMin = std::numeric_limits<num>::max();
                Vector3d center = primaryCellCentersCart.row(z+1);
                for (int i=0; i<vSplines[z+1]->rows(); i++) {
                    num currDist = distance(center, vSplines[z+1]->row(i));
                    if (currDist < distMin) {
                        index = i;
                        distMin = currDist;
                    }
                }
            } else {
                index = (indicies(z, 0) + indicies(z, 1)) / 2.0;
            }
            indexList[z+1]->push_front(index);
            indicies.row(z+1) << index, index;
        }
    }
    bool done = false;
    int shift, i = 1; 
    Vector3d cc;
    MatrixXi oldIndicies = MatrixXi::Zero(indicies.rows(), indicies.cols());
    oldIndicies << indicies;
    while (!done) {
        done = true;
        shift = (i + 1) % 2;
        // Add the new cells for every other spline
        num maxDiff = 0;
        for (int z = !shift; z < startingLineNodeCount; z+=2) {
            for (int w = 0; w < 2; w++) {
                for (int k = 0; k < 2; k++) {
                    if ((!k && !doneUp[z]) || (k && !doneDown[z])) {
                        cc = findNextCC(vSplines[z], indicies(z, k), spacing, k, splineEnded);
                        if (!splineEnded) {
                            if (k) {
                                ccs[z]->push_back(cc);
                                levels[z]->push_back((*levels[z])[levels[z]->size()-1]+2);
                                indexList[z]->push_back(indicies(z, k));
                            } else {
                                ccs[z]->push_front(cc);
                                levels[z]->push_front((*levels[z])[0]-2);
                                indexList[z]->push_front(indicies(z, k));
                            }
                            num dif = abs(surfaceDistance((*ccs[z])[ccs[z]->size()-1], (*ccs[z])[ccs[z]->size()-2]) - spacing);
                            if (dif > maxDiff) {
                                maxDiff = dif;
                            }
                        } else {
                            if (k) {
                                doneDown[z] = true;
                            } else {
                                doneUp[z] = true;
                            }
                        }
                    }
                }
            }
        }
        for (int z = !shift; z < startingLineNodeCount - 1; z+=2) {
            // The k stuff handles populating both column 0 and 2 when we're on z = 1
            //int k = - (z == 1);
            int k = 0;
            while (k < 1) {
                int splineInd = k ? z - 1 : z + 1;
                auto spline = vSplines[splineInd];
                if (!doneUp[splineInd] &&
                        ((splineInd-1>=0 && !doneUp[splineInd-1]) ||
                         (splineInd+1 < startingLineNodeCount && !doneUp[splineInd+1]))) {
                    index = centerCC(splineInd, indicies, 0, doneUp, ccs);
                    if (index >= 0 && distance(spline->row(indicies(splineInd, 0)), spline->row(index))
                            > 0.8 * spacing) {
                        indicies(splineInd, 0) = index;
                        ccs[splineInd]->push_front(spline->row(index));
                        levels[splineInd]->push_front((*levels[splineInd])[0]-2);
                        indexList[splineInd]->push_front(index);
                    } else {
                        cc = findNextCC(vSplines[splineInd], indicies(splineInd, 0), spacing, false, splineEnded);
                        if (!splineEnded) {
                            ccs[splineInd]->push_front(cc);
                            levels[splineInd]->push_front((*levels[splineInd])[0]-2);
                            indexList[splineInd]->push_front(indicies(splineInd, 0));
                        } else {
                            doneUp[splineInd] = true;
                        }
                    }
                } else {
                    cc = findNextCC(vSplines[splineInd], indicies(splineInd, 0), spacing, false, splineEnded);
                    if (!splineEnded) {
                        ccs[splineInd]->push_front(cc);
                        levels[splineInd]->push_front((*levels[splineInd])[0]-2);
                        indexList[splineInd]->push_front(indicies(splineInd, 0));
                    } else {
                        doneUp[splineInd] = true;
                    }
                }
                int s = -1;
                if (splineInd-1 >= 0)
                    s = ccs[splineInd-1]->size()-1;
                if (!doneDown[splineInd] &&
                        ((splineInd-1>=0 && !doneDown[splineInd-1]) ||
                         (splineInd+1<startingLineNodeCount && !doneDown[splineInd+1])) &&
                        static_cast<int>(ccs.size()) > splineInd+1 && static_cast<int>(ccs[splineInd+1]->size()) > s) {
                    index = centerCC(splineInd, indicies, 1, doneDown, ccs);
                    if (index >= 0 && distance(spline->row(indicies(splineInd, 1)), spline->row(index))
                            > 0.8 * spacing) {
                        indicies(splineInd, 1) = index;
                        ccs[splineInd]->push_back(spline->row(index));
                        levels[splineInd]->push_back((*levels[splineInd])[levels[splineInd]->size()-1]+2);
                        indexList[splineInd]->push_back(index);
                    } else {
                        cc = findNextCC(vSplines[splineInd], indicies(splineInd, 1), spacing, true, splineEnded);
                        if (!splineEnded) {
                            ccs[splineInd]->push_back(cc);
                            levels[splineInd]->push_back((*levels[splineInd])[levels[splineInd]->size()-1]+2);
                            indexList[splineInd]->push_back(indicies(splineInd, 1));
                        } else {
                            doneDown[splineInd] = true;
                        }
                    }
                } else {
                    cc = findNextCC(vSplines[splineInd], indicies(splineInd, 1), spacing, true, splineEnded);
                    if (!splineEnded) {
                        ccs[splineInd]->push_back(cc);
                        levels[splineInd]->push_back((*levels[splineInd])[levels[splineInd]->size()-1]+2);
                        indexList[splineInd]->push_back(indicies(splineInd, 1));
                    } else {
                        doneDown[splineInd] = true;
                    }
                }
                k++;
            }
        }
        i++;
        done = doneUp[0] && doneDown[0];
        for (int z = 1; z < startingLineNodeCount - 1; z++) {
            done = done && doneUp[z] && doneDown[z];
        }
    }
    cout << "\33[2K\r";
    bar.set_label("Running Reverse Barycentric");
    int s = static_cast<int>(ccs.size());
    vector<vector<Cell*>> tempCCGrid;
    for (int z = 0; z<s; z++) {
        auto col = *ccs[z];
        MatrixXd t(col.size(), 3), c(col.size(), 2), n(col.size(), 3);
        vector<Cell*> currVec;
        tempCCGrid.push_back(currVec);
        for (int k = 0; k<static_cast<int>(col.size()); k++) {
            t.row(k) = vSplines[z]->row((*indexList[z])[k]);
            Vector3d tCopy = t.row(k), normal;
            c.row(k) = mesh->barycentricLoopingR(tCopy, normal);
            n.row(k) << normal;
        }
        *vCellCentersCart[z] = t;
        *vCellCentersUV[z] = c;
        vCellCentersNormal[z] = new MatrixXd;
        *vCellCentersNormal[z] = n;
        MatrixXd f = splines[z]->splinePointsUV;
    }
    num maxU = 0, maxV = 0;
    int maxVCells = 0;
    auto sp = primarySpline.splinePoints;
    for (int i=0; i<sp.rows()-1; i++) {
        maxU += distance(sp.row(i), sp.row(i+1));
    }
    for (auto s : splines) {
        num ddist = 0;
        sp = s->splinePoints;
        for (int i=0; i<sp.rows()-1; i++) {
            ddist += distance(sp.row(i), sp.row(i+1));
        }
        maxV = std::max(maxV, ddist);
        maxVCells = std::max(maxVCells, static_cast<int>(s->cellCentersUV.rows()));
    }
    cout << "Primary spline length: " << maxU << endl;
    cout << "Max vertical spline length: " << maxV << endl;
    cout << "Aspect ratio: " << std::min(maxU, maxV) / std::max(maxU, maxV) << endl;
    cout << fmt::format("Node count (u, v (max)): ({}, {})\n", primarySpline.cellCentersUV.rows(), maxVCells);
    bar.reset();
    setScaleFactor();
}

void Model::setScaleFactor() {
    Vector3d minPoint, maxPoint;
    double ma = std::numeric_limits<double>::max();
    double mi = std::numeric_limits<double>::min();
    minPoint << ma, ma, ma;
    maxPoint << mi, mi, mi;
    for (Vector3d row : mesh->realPoints.rowwise()) {
        for (int z=0; z<3; z++) {
            if (row(z) < minPoint(z))
                minPoint(z) = row(z);
            if (row(z) > maxPoint(z))
                maxPoint(z) = row(z);
        }
    }

    // I came up with the numbers below for the face model
    // where this value is about 6.67, hence the divisor
    //scaleFactor = distance(minPoint, maxPoint) / 6.67 * (15.0 / static_cast<double>(startingLineNodeCount));

    scaleFactor = .8*(maxPoint - minPoint).cwiseAbs().minCoeff() / log(static_cast<double>(startingLineNodeCount));
}

void Model::processK(toml::node_view<toml::node> tbl, num &kIn, num &kMin, num &kMax, num &kStep, num kDefault) {
    auto tblArr = tbl.as_array();
    if (tblArr) {
        if (tblArr->size() != 3) {
            cerr << "ERROR: k values must be either singleton or length 3 array\n";
            exit(1);
        }
        setts["simulation"].as_table()->insert("batch", true);
        kMin = tbl[0].value_or(kDefault);
        kIn = kMin;
        kMax = tbl[1].value_or(kDefault);
        kStep = tbl[2].value_or(kMax - kMin);
    } else {
        kIn = tbl.value_or(kDefault);
        kMin = kIn;
        kMax = kIn;
        kStep = kIn;
    }
}

int Model::centerCC(int splineInd, MatrixXi &indicies, int direction, vector<bool> &doneVec, vector<std::deque<Vector3d>*> &ccs) {
    bool check1, check2;
    auto spline = vSplines[splineInd];
    num minA = std::numeric_limits<num>::max(), currA;
    int s1, s2;
    if (direction) {
        s1 = -1;
        if (splineInd-1 >= 0)
            s1 = ccs[splineInd-1]->size()-1;
        s2 = s1 - 1;
    } else {
        s1 = 0;
        s2 = 1;
    }
    int index = -1;
    check1 = splineInd-1 >= 0 && !doneVec[splineInd-1];
    check2 = splineInd + 1 < startingLineNodeCount && !doneVec[splineInd+1];
    for (int w = indicies(splineInd, direction); direction ? w > 0 : w < spline->rows(); direction ? w-- : w++) {
        vector<num> dists;
        currA = 0.0;
        if (check1) {
            dists.push_back(distance((*ccs[splineInd-1])[s1], spline->row(w)));
            dists.push_back(distance(spline->row(w), (*ccs[splineInd-1])[s2]));
        }
        if (check2) {
            dists.push_back(distance((*ccs[splineInd+1])[s1], spline->row(w)));
            dists.push_back(distance(spline->row(w), (*ccs[splineInd+1])[s2]));
        }
        num angErr = 0.0;
        if (check1 && check2) {
            currA = stdev(dists);
            num err1 = abs(angle3Points((*ccs[splineInd-1])[s1], spline->row(w), (*ccs[splineInd+1])[s2]) - PI);
            num err2 = abs(angle3Points((*ccs[splineInd-1])[s2], spline->row(w), (*ccs[splineInd+1])[s1]) - PI);
            angErr = err1 + err2;
        } else {
            currA = dists[0] + dists[1];
        }
        if ((angErr < 2.0) && (currA < minA)) {
            minA = currA;
            index = w;
        }
    }
    return index;
}

Vector3d Model::findNextCC(MatrixXd *spline, Vector3d prevCellCenter, int &index, num stop, bool reverse, bool &splineEnded) {
    num minD = std::numeric_limits<num>::max(), currD;
    index = -1;
    for (int i = 0; i < spline->rows(); i++) {
        currD = distance(spline->row(i), prevCellCenter);
        if (currD < minD) {
            minD = currD;
            index = i;
        }
    }
    return findNextCC(spline, index, stop, reverse, splineEnded);
}

Vector3d Model::findNextCC(MatrixXd *spline, int &index, num stop, bool reverse, bool &splineEnded) {
    int step = reverse ? -1 : 1;
    num dist = 0;

    while (dist < stop && (index <= spline->rows() - 5) && (index >= 5)) {
        dist += distance(spline->row(index), spline->row(index+step));
        index += step;
    }
    splineEnded = dist < 0.8 * stop;
    return spline->row(std::min(static_cast<int>(spline->rows())-1, index + step));
}

num Model::distance3dFrom2d(Vector2d p1, Vector2d p2) {
    Vector3d p12, p22;
    p12 = mesh->barycentricLoopingF(p1);
    p22 = mesh->barycentricLoopingF(p2);
    return distance(p12, p22);
}

num Model::surfaceDistance(Vector3d p1, Vector3d p2) {
    Vector2d p12, p22;
    p12 = mesh->barycentricLoopingR(p1);
    p22 = mesh->barycentricLoopingR(p2);
    return surfaceDistance(p12, p22);
}

num Model::surfaceDistance(Vector2d p1, Vector2d p2) {
    num uDiff = p2(0) - p1(0);
    num vDiff = p2(1) - p1(1);
    num steps = 30;
    num uStep = uDiff / steps;
    num vStep = vDiff / steps;
    num dist = 0.0;
    MatrixXd pts(static_cast<int>(steps)+1, 3);
    pts.row(0) << p1, 0.0;
    Vector3d point;
    point << pts.row(0);
    for (int i=1; i<static_cast<int>(steps)+1; i++) {
        point(0) += uStep;
        point(1) += vStep;
        pts.row(i) << point;
    }
    MatrixXd pps = mesh->barycentricLoopingF(pts);
    for (int i=1; i<static_cast<int>(steps)+1; i++) {
        dist += distance(pps.row(i-1), pps.row(i));
    }
    return dist;
}

void Model::bumpPoint(Vector2d &point, Vector3d &realPoint, num &horizShift, bool isUSpline, Spline &otherSpline, int &stepIn) {
    int outTri = -1, inTri = -1, startPoint = -1;
    bumpPoint(point, realPoint, horizShift, isUSpline, otherSpline, outTri, inTri, startPoint, stepIn);
}

void Model::bumpPoint(Vector2d &point, Vector3d &realPoint, num &horizShift, bool isUSpline, Spline &otherSpline, int &outTri, int &inTri, int &startPoint, int &stepIn) {
    num mdist = std::numeric_limits<num>::max(), minDist = std::numeric_limits<num>::max(), currDist;
    int ind = -1, j, currInd, steps = 100, rows =
        otherSpline.splinePointsUV.rows(), incr = rows / steps, minInd = -1;
    realPoint = mesh->barycentricLoopingF(point);
    for (j=0; j<steps; j++) {
        currInd = j*incr;
        currDist = surfaceDistance(point, otherSpline.splinePointsUV.row(currInd));
        if (currDist < minDist) {
            minInd = currInd;
            minDist = currDist;
        }
    }
    for (j=std::max(0, minInd-incr); j<min(rows-1, minInd+incr); j++) {
        currDist = surfaceDistance(point, otherSpline.splinePointsUV.row(j));
        if (currDist < minDist) {
            minInd = j;
            minDist = currDist;
        }
    }
    ind = minInd;
    mdist = minDist;
    Vector2d ppUV = otherSpline.splinePointsUV.row(ind);
    num dd;
    startPoint = ind;
    bool done = false, flipped = false;
    num ss = 0.001, currShift = 0.0, sign = 1.0, tt = abs(mdist - spacing);
    Vector2d np, prevNP;
    while (!done) {
        currShift += sign * ss;
        if (isUSpline)
            np << point(0) + currShift, point(1);
        else
            np << point(0), point(1) + currShift;
        done = np(0) > 1.0 || np(0) < 0.0 || np(1) > 1.0 || np(1) < 0.0;
        np << bound(np(0)), bound(np(1));
        Vector3d newRealPoint = mesh->barycentricLoopingF(np);
        dd = abs(surfaceDistance(np, ppUV) - spacing);
        if (dd > tt) {
            if (flipped)
                done = true;
            flipped = true;
            sign *= -1.0;
        } else {
            tt = dd;
            realPoint << newRealPoint;
        }
    }
    point << np;
    horizShift += currShift;
}

void Model::populateSpline(Spline &spline) {
    auto s = Spline();
    populateSpline(spline, s);
}

void Model::populateSpline(Spline &spline, Spline &otherSpline) {
    Vector2d point;
    bool fullLen = spline.length <= 0.0;
    num horizShift = 0.0;
    int cPtCount = cubic ? 20 : MAX_CONTROL_PT_COUNT ;
    int maxPts = static_cast<int>(round(cPtCount /
                (spline.halfSpline ? 2.0 : 1.0)) *
                (fullLen ? 1.0 : spline.length));
    if (spline.direction == CENTER_OUT) {
        Vector2d uvPoint;
        uvPoint << spline.uValue, spline.vValue;
        Spline splineF = spline.duplicate(), splineR = spline.duplicate();
        splineF.halfSpline = true;
        splineR.halfSpline = true;
        splineF.direction = FORWARD;
        splineR.direction = REVERSE;
        populateSpline(splineF, otherSpline);
        populateSpline(splineR, otherSpline);
        spline.combine(splineF, splineR);
        if (spline.bump && !otherSpline.bump) {
            int diff = otherSpline.index - spline.index;
            auto thirdSpline = *splines[spline.index + diff * 2];
            for (int i=0; i<std::min(std::min(otherSpline.uvPoints.rows(),
                            thirdSpline.uvPoints.rows()),
                            spline.uvPoints.rows()); i++) {
                otherSpline.uvPoints.row(i) = (spline.uvPoints.row(i) + thirdSpline.uvPoints.row(i)) / 2.0;
            }
            otherSpline.realPoints = mesh->barycentricLoopingF(otherSpline.uvPoints);
            genSplinePointsAndCCs(otherSpline, false);
        }
    } else {
        bool reverse = spline.direction == REVERSE;
        num revMod = reverse ? -1 : 1, revMod2 = reverse ? 1 : 0, stepSize;
        string pointString = "point";
        string vecString = "vector";
        vector<json> entries;
        spline.length = fullLen ? 1.0 : spline.length;
        int pointCount = std::max(10, maxPts);
        spline.uvPoints = MatrixXd::Zero(pointCount, 3);
        spline.realPoints = MatrixXd::Zero(pointCount, 3);
        bool done = false;
        int i = 0;
        if (fullLen) {
            num incr = 0.000001, step;
            while (true) {
                step = revMod * incr * i;
                if (spline.isUSpline)
                    point << spline.uValue + step * sin(spline.offset), revMod2 + step * cos(spline.offset);
                else
                    point << revMod2 + step * cos(spline.offset), spline.vValue + step * sin(spline.offset);
                if (point(0) > 1.0 || point(0) < 0.0 || point(1) > 1.0 || point(1) < 0.0) {
                    break;
                }
                i++;
            }
            stepSize = step / static_cast<num>(pointCount - 1);
        } else {
            stepSize = spline.length / static_cast<num>(pointCount - 1);
        }

        int startPoint = -1;
        tqdm bar;
        if (spline.bump) {
            num toPrint = 1.0 + ((static_cast<num>(bumpNum) - 1.0) / 2.0);
            bar.set_label(fmt::format("Bumping Spline {}/{}", toPrint, bumpTotal));
            bar.progress(0, pointCount);
        }
        int stepIn = 0;
        int rs = static_cast<int>(otherSpline.splinePointsUV.rows());
        for (i = 0; i < pointCount && !done; i++) {
            num step = revMod * stepSize * i;
            if (!fullLen)
                if (spline.isUSpline)
                    point << spline.uValue + horizShift + step * sin(spline.offset), spline.vValue + step * cos(spline.offset);
                else
                    point << spline.uValue + step * cos(spline.offset), spline.vValue + horizShift + step * sin(spline.offset);
            else
                if (spline.isUSpline)
                    point << spline.uValue + horizShift + step * sin(spline.offset), revMod2 + step * cos(spline.offset);
                else
                    point << revMod2 + step * cos(spline.offset), spline.vValue + horizShift + step * sin(spline.offset);
            if ((point(0) > 1.0 || point(0) < 0.0 || point(1) > 1.0 || point(1) < 0.0)) {
                done = i;
                if (i > 1 && i < maxPts - 1) {
                    break;
                }
            }
            point << bound(point(0)), bound(point(1));
            Vector2d pointBackup;
            pointBackup << point;
            Vector3d realPoint;
            // TODO: Second part of this check didn't do what I wanted but leaving for now
            if (otherSpline.realConstructor && spline.bump && (startPoint < 0 || !(startPoint < 10 || startPoint > rs - 10))) {
                int inTri=-1, outTri=-1;
                bumpPoint(point, realPoint, horizShift, spline.isUSpline, otherSpline, outTri, inTri, startPoint, stepIn);
                bar.progress(i, pointCount);
            }

            if ((!cubic && point(0) > 1e-5 && point(1) > 1e-5) || point(0) >= 0) {
                spline.uvPoints.row(i) << point, 0.0;
            } else {
                done = i;
                if (!done) {
                    spline.uvPoints.row(i) << point, 0.0;
                } else {
                    spline.uvPoints.row(i) << spline.uvPoints.row(i-1);
                }
            }
        }
        spline.realPoints = mesh->barycentricLoopingF(spline.uvPoints);
        bumpNum += otherSpline.realConstructor && spline.bump;
        spline.realPoints.conservativeResize(i, spline.realPoints.cols());
        spline.uvPoints.conservativeResize(i, spline.uvPoints.cols());
        genSplinePointsAndCCs(spline, reverse);
    }
}

void Model::genSplinePointsAndCCs(Spline &spline, bool reverse) {
    int i;
    spline.tangents = MatrixXd::Zero(0, 3);
    if (spline.realPoints.rows() < 2) {
        cerr << "Error: spline with fewer than 2 control points passed to genSplinePoints\n";
        exit(1);
    }
    if (reverse) {
        // Was having weird issues if I set any of these directly as their
        // reverse so using temp buffers for all this
        MatrixXd t = spline.realPoints.colwise().reverse();
        spline.splinePoints = genSpline(t, spline.tangents, 250);
    } else {
        MatrixXd rp = MatrixXd::Zero(spline.uvPoints.rows(), 3);
        spline.splinePoints = genSpline(spline.realPoints, spline.tangents, 250);
    }
    if (reverse) {
        MatrixXd t = spline.tangents.colwise().reverse();
        spline.tangents << t;
        t = spline.splinePoints.colwise().reverse();
        spline.splinePoints << t;
    }
    num splineLen = 0;
    for (i=0; i<spline.splinePoints.rows() - 1; i++) {
        splineLen += distance(spline.splinePoints.row(i), spline.splinePoints.row(i+1));
    }
    if (spline.setSpacing) {
        // I think this is right now but still might be worth triple checking (1/2)
        spacing = splineLen / (static_cast<num>(startingLineNodeCount));
        if (spacingDiv)
            zHeight = splineLen / FIXED_ZHEIGHT_DIV / zHeight;
    }
    num stop = spacing;
    num dist = 0, currDist;
    if (spline.length == 1.0 || spline.shift) {
        stop = spacing / 2.0;
    }
    // Pairs with this (2/2)
    int centerCount = round(std::max(static_cast<num>(0.0), splineLen / spacing));
    spline.cellCentersCart = MatrixXd::Zero(centerCount, 3);
    spline.cellCentersUV = MatrixXd::Zero(centerCount, 2);
    int c = 0;
    i = 0;
    vector<int> iVals;
    while (c < centerCount) {
        currDist = distance(spline.splinePoints.row(i), spline.splinePoints.row(i+1));
        dist += currDist;
        i++;
        if (dist >= stop || (i > spline.splinePoints.rows() - 5)) {
            if (dist - stop < currDist / 2.0)
                i--;
            dist = 0;
            stop = spacing;
            spline.cellCentersCart.row(c) = spline.splinePoints.row(i);
            iVals.push_back(i);
            c++;
        }
    }
    MatrixXd zz(2, 3);
    zz << 0.5, 0.0, 0.0,
          0.5, 1.0, 0.0;
    MatrixXd yy = mesh->barycentricLoopingF(zz);
    spline.splinePointsUV = mesh->barycentricLoopingR(spline.splinePoints).block(0,0,spline.splinePoints.rows(),2);
    if (centerCount) {
        c=0;
        for (int i : iVals) {
            spline.cellCentersUV.row(c) = spline.splinePointsUV.row(i);
            c++;
        }
    }
}

void Model::calcIdealPosition(Vector2d &curCenterUV, Vector3d &curCenterCart, Vector2d &compUV, num targetDist) {
    Vector2d diff = curCenterUV - compUV;
    num stepSize = 0.0002;
    Vector2d step = diff.normalized() * stepSize;
    Vector3d compCart = mesh->barycentricLoopingF(compUV);
    int i = 20;
    num currDist = 0;
    Vector2d currUV(0, 0);
    Vector3d currCart(0, 0, 0);
    while (currDist < targetDist) {
        Vector2d currUV = compUV + step * i;
        const num MAX = 0.99999, MIN = 0.00001;
        currUV = currUV.cwiseMin(MAX);
        currUV = currUV.cwiseMax(MIN);
        Vector3d currCart = mesh->barycentricLoopingF(currUV);
        currDist = distance(compCart, currCart);
        i++;
        if (currUV(0) == MIN || currUV(1) == MIN || currUV(0) == MAX || currUV(1) == MAX) {
            break;
        }
    }
    curCenterUV = currUV;
    curCenterCart = currCart;
}
