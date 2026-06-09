#include "tissue.h"

template <typename T>
string vecToString(vector<T*> vec) {
    std::stringstream ss;
    for (auto t : vec)
        ss << (*t) << endl << endl;
    return ss.str();
}
template <typename T>
string vecToStringOneLine(vector<T> vec) {
    std::stringstream ss;
    bool b = true;
    for (auto t : vec) {
        if (b)
            b = false;
        else
            ss << ", ";
        ss << t;
    }
    return ss.str();
}

Tissue::Tissue(Model &model) : cellGrid(CellGrid(model, model.spacing, false)), model(model) {
    Tissue(model, false);
}

void Tissue::addLayersCubic() {
    uint i, j;
    num maxShrinkage = grid.maxShrinkage, kappa, Rmin, R, angle, L0, LN;
    num maxKappa = 0.0, distanceMax = 0.0;
    cout << "Z height: " << model.zHeight << endl;
    model.zHeight /= (model.layersCubic - 1);
    cout << "Interlayer spacing: " << model.zHeight << endl;
    for (auto e : grid.getEdges()) {
        if (e->curvature > maxKappa)
            maxKappa = e->curvature;
        if (e->getRestLength() > distanceMax)
            distanceMax = e->getRestLength();
    }
    cout << "Max curvature: " << maxKappa << endl;
    R = 1.0 / maxKappa;
    num T = R * maxShrinkage;
    num maxZHeight = T / (model.layersCubic - 1);
    if (model.zHeight > maxZHeight) {
        cout << "WARN: User specified zHeight greater than maximum allowable for model, defaulting to max of " << maxZHeight << endl;
        model.zHeight = maxZHeight;
    }
    MatrixXd vMatOrig = grid.getVertexMatV().transpose();
    MatrixXd vMatNew(vMatOrig.rows()*(model.layersCubic-1), 3);
    vector<vector<int>> vertPairsVertInt;
    for (i=1; i<model.layersCubic; i++) {
        for (j=0; j<vMatOrig.rows(); j++) {
            int ind = (i-1)*vMatOrig.rows()+j;
            vMatNew.row(ind) << vMatOrig.row(j);
            vMatNew.row(ind)(2) = i * model.zHeight;
            vector<int> pair = {ind, static_cast<int>(i*vMatOrig.rows()+j)};
            vertPairsVertInt.push_back(pair);
        }
    }
    grid.addVertices(vMatNew);
    auto verts = grid.getVerticesV();
    for (i=0; i<verts.size()/model.layersCubic; i++) {
        for (j=1; j<model.layersCubic; j++) {
            auto v = verts[i+j*vMatOrig.rows()];
            auto oP = verts[i]->originalPosition;
            v->originalPosition << oP(0), oP(1), v->originalPosition(2);
        }
    }
    vector<vector<Vertex*>> vertPairs;
    vector<num> restLengths;
    num maxShrink = 0.0;
    for (auto e : grid.getEdges()) {
        int i1 = e->getV1()->getIndex();
        int i2 = e->getV2()->getIndex();

        VectorXd LNs(model.layersCubic);
        kappa = e->curvature;
        Rmin = model.zHeight / maxShrinkage;
        L0 = e->getRestLength();
        if (kappa == 0.0) {
            R = -1.0;
            angle = 0.0;
        } else {
            R = 1 / abs(kappa);
            R = std::max(R, Rmin);
            angle = ((kappa > 0) - (kappa < 0)) *
                        asin(L0 / 2.0 / R) * 2.0;
        }
        for (i=0; i<model.layersCubic; i++) {
            if (R < 0.0) {
                LNs(i) = L0;
            } else {
                LN = L0 - i * 2.0 * model.zHeight * sin(abs(angle) / 2.0);
                LN = std::max(LN, L0 * (1.0 - maxShrinkage));
                LNs(i) = LN;
                num shrink = 1.0 - LN / L0;
                if (shrink > maxShrink) {
                    maxShrink = shrink;
                }
            }
        }

        if (angle < 0.0)
            LNs.reverseInPlace();

        e->setRestLength(LNs(0));
        for (i=1; i<model.layersCubic; i++) {
            vector<Vertex*> pair = {verts[i1+i*vMatOrig.rows()], verts[i2+i*vMatOrig.rows()]};
            vertPairs.push_back(pair);
            restLengths.push_back(LNs(i));
        }
    }
    cout << "Max shrinkage L0->Ln: " << maxShrink << endl;
    auto newEdges = grid.addEdges(vertPairs);
    for (i=0; i<newEdges.size(); i++) {
        newEdges[i]->setRestLength(restLengths[i]);
    }
    for (auto p : vertPairsVertInt) {
            vector<Vertex*> pair = {verts[p[0]], verts[p[1]]};
            auto e = grid.addEdge(pair);
            e->isVertical = true;
    }
    model.hasDuplicatedCubic = true;
}

void Tissue::buildCubic() {
    int i, j, k, w, nU = model.nodeCountU, nV = model.nodeCountV;
    int rows = model.splines[0]->splinePoints.rows();
    num sp = static_cast<num>(model.splines[0]->splinePoints.rows());
    num L0Max = 0.0, stepsU = sp / (nV - 1), stepsV = sp / (nU - 1);
    MatrixXd vertMat = MatrixXd::Zero(nU*nV, 3);
    num x = 0.0, y = 0.0, z = 0.0;
    vector<num> distsU, distsV, curvsU, curvsV;
    for (w=0; w<2; w++) {
        for (i=0; i<nU; i++) {
            for (j=0; j<nV; j++) {
                if (!w) {
                    if (j < nV - 1) {
                        int start = j * stepsU, count = 0;
                        num distU = 0.0, curvU = 0.0;
                        for (k=start; k<std::min(static_cast<int>(start+stepsU), rows-1); k++) {
                            distU += distance(
                                    model.usplines[i]->splinePoints.row(k),
                                    model.usplines[i]->splinePoints.row(k+1));
                            curvU += model.usplines[i]->curvatures[k];
                            count++;
                        }
                        curvsU.push_back(curvU / count);
                        distsU.push_back(distU);
                        L0Max = std::max(L0Max, distU);
                    }
                    if (i < nU - 1) {
                        int start = i * stepsV, count = 0;
                        num distV = 0.0, curvV = 0.0;
                        for (k=start; k<std::min(static_cast<int>(start+stepsV), rows-1); k++) {
                            distV += distance(
                                    model.splines[j]->splinePoints.row(k),
                                    model.splines[j]->splinePoints.row(k+1));
                            curvV += model.splines[j]->curvatures[k];
                            count++;
                        }
                        curvsV.push_back(curvV / count);
                        distsV.push_back(distV);
                        L0Max = std::max(L0Max, distV);
                    }
                } else {
                    x = i * L0Max;
                    y = j * L0Max;
                    vertMat.row(i*nV+j) << x, y, z;
                }
            }
        }
    }
    auto verts = grid.addVertices(vertMat);
    num maxShrink = 0.0;
    for (i=0; i<nU; i++) {
        for (j=0; j<nV-1; j++) {
            vector<Vertex*> pair = {verts[i*nV+j], verts[i*nV+j+1]};
            auto e = grid.addEdge(pair);
            e->setRestLength(distsU[i*(nV-1)+j]);
            num shrink = 1.0 - (e->getRestLength() / L0Max);
            if (shrink > maxShrink)
                maxShrink = shrink;
            e->curvature = curvsU[i*(nV-1)+j];
        }
    }
    for (j=0; j<nV; j++) {
        for (i=0; i<nU-1; i++) {
            vector<Vertex*> pair = {verts[j+i*nV], verts[j+(i+1)*nV]};
            auto e = grid.addEdge(pair);
            e->setRestLength(distsV[j+i*nV]);
            num shrink = 1.0 - (e->getRestLength() / L0Max);
            if (shrink > maxShrink)
                maxShrink = shrink;
            e->curvature = curvsV[j+i*nV];
        }
    }
    cout << "Ls: " << L0Max << endl;
    cout << "Max shrinkage Ls->L0: " << maxShrink << endl;
}

Tissue::Tissue(Model &model, bool cubic) : cellGrid(CellGrid(model, model.spacing, cubic)), model(model) {
    grid = Grid();
    grid.maxShrinkage = setts["grid"]["max_shrinkage"].value_or(0.8);
    Cell::maxShrinkage = grid.maxShrinkage;
    this->cubic = cubic;
    if (model.kAngIn120 >= 0)
        grid.kAng = model.kAngIn120;
    if (model.kAngIn90 >= 0)
        grid.kAng90 = model.kAngIn90;
    if (model.kLinearIn >= 0)
        grid.kEdge = model.kLinearIn;
    if (setts["global"]["import_from_json"].value_or(false))
        return;
    if (cubic) {
        buildCubic();
        return;
    }

    auto vertices = new std::set<Vertex*>();
    tqdm bar;
    bar.set_label("Setting cell neighbors");
    int c = 0;
    for (auto col : cellGrid.grid) {
        bar.progress(c, cellGrid.grid.size());
        c++;
        for (auto cellWrapper : col) {
            auto cell = cellWrapper->cell;
            cell->level = cellWrapper->level;
            cells.push_back(cell);
            MatrixXd initPts = MatrixXd::Zero(12, 3);
            for (int i=0; i<6; i++) {
                num ang = i * PI / 3, s = model.spacing / sqrt(3.0);
                Vector3d cc = cell->initialGridCenter;
                initPts.row(i) << cc(0) + cos(ang) * s, cc(1) + sin(ang) * s, 0;
                initPts.row(i+6) << cc(0) + cos(ang) * s, cc(1) + sin(ang) * s, model.zHeight;
            }
            vector<Vertex*> vertices2(12);
            Cell *neighbor;
            for (int i = 0; i<6; i++) {
                int ind1, ind2;
                switch (i) {
                    case 0:
                        // TODO: West and East are like backwards, this works
                        // but the naming is wrong
                        neighbor = cell->neighbors.southEast;
                        ind1 = 4;
                        ind2 = 3;
                        break;
                    case 1:
                        neighbor = cell->neighbors.south;
                        ind1 = 5;
                        ind2 = 4;
                        break;
                    case 2:
                        neighbor = cell->neighbors.southWest;
                        ind1 = 0;
                        ind2 = 5;
                        break;
                    case 3:
                        neighbor = cell->neighbors.northWest;
                        ind1 = 1;
                        ind2 = 0;
                        break;
                    case 4:
                        neighbor = cell->neighbors.north;
                        ind1 = 2;
                        ind2 = 1;
                        break;
                    default:
                        neighbor = cell->neighbors.northEast;
                        ind1 = 3;
                        ind2 = 2;
                }
                if (neighbor != nullptr && neighbor->initialized) {
                    auto nVerts = neighbor->getVertices();
                    vertices2[i] = nVerts[ind1];
                    vertices2[(i+1)%6] = nVerts[ind2];
                    vertices2[i+6] = nVerts[ind1+6];
                    vertices2[(i+1)%6+6] = nVerts[ind2+6];
                } else if (vertices2[i] == nullptr) {
                    Vector3d temp = initPts.row(i);
                    vertices2[i] = new Vertex(temp, cell);
                    temp = initPts.row(i+6);
                    vertices2[i+6] = new Vertex(temp, cell);
                }
            }
            for (auto v : vertices2) {
                vertices->insert(v);
            }
            grid.addCellToVertices(cell, vertices2);
            cell->initialized = true;
            // This should only need called once but putting it here for now
            grid.vertices = *vertices;
            for (int i=0; i<6; i++) {
                switch (i) {
                    case 0:
                        // TODO: West and East are backwards, this works
                        // but the naming is wrong
                        neighbor = cell->neighbors.southEast;
                        break;
                    case 1:
                        neighbor = cell->neighbors.south;
                        break;
                    case 2:
                        neighbor = cell->neighbors.southWest;
                        break;
                    case 3:
                        neighbor = cell->neighbors.northWest;
                        break;
                    case 4:
                        neighbor = cell->neighbors.north;
                        break;
                    default:
                        neighbor = cell->neighbors.northEast;
                }
                if (neighbor != nullptr && neighbor->initialized) {
                    int ind = ((i + 3) % 6) * 3;
                    vector<Edge*> edges;
                    auto nEdges = neighbor->getEdges();
                    edges.push_back(nEdges[ind]);
                    edges.push_back(nEdges[ind+1]);
                    edges.push_back(nEdges[ind+2]);
                    grid.addCellToEdges(cell, edges);
                } else {
                    vector<vector<Vertex*>> vertPairs;
                    auto pair1 = new vector<Vertex*>();
                    auto pair2 = new vector<Vertex*>();
                    auto pair3 = new vector<Vertex*>();
                    pair1->push_back(vertices2[i]);
                    pair1->push_back(vertices2[(i+1)%6]);
                    pair2->push_back(vertices2[i+6]);
                    pair2->push_back(vertices2[(i+1)%6+6]);
                    pair3->push_back(vertices2[i]);
                    pair3->push_back(vertices2[i+6]);
                    vertPairs.push_back(*pair1);
                    vertPairs.push_back(*pair2);
                    vertPairs.push_back(*pair3);
                    auto newEdges = grid.addEdges(cell, vertPairs);
                    newEdges[0]->type = 0;
                    newEdges[1]->type = 1;
                    newEdges[2]->type = 2;
                    newEdges[2]->isVertical = true;
                }
            }
        }
    }
    bar.reset();
    bar.set_label("Initializing cell splines");
    vector<num> offsets = {PI / 6.0, PI / 2.0, -PI / 6.0};
    c = 0;
    cout << "Cell count: " << cells.size() << endl;
    for (auto cell : cells) {
        {
            bar.progress(c, cells.size());
            c++;
        }
        // TODO Old: Fix speedup stuff
        for (num offset : offsets) {
            Vector2d neighborUV;
            if (cell->neighbors.south != nullptr) {
                neighborUV = cell->neighbors.south->initialCenterUV;
            } else if (cell->neighbors.north != nullptr) {
                neighborUV = cell->neighbors.north->initialCenterUV;
            } else {
                cerr << "Error: cell column of length 1 encountered\n";
                exit(1);
            }
            num dist = 1.0 / model.startingLineNodeCount * 6;
            auto splinePointsO = new MatrixXd;
            auto realPointsO = new MatrixXd;
            auto tangentsO = new MatrixXd;
            MatrixXd uvPoints, tempRealPointsF, tempRealPointsR, tempSPointsF, tempSPointsR, tempTanF, tempTanR;
            auto splineF = Spline(false, FORWARD, cell->initialCenterUV(0), cell->initialCenterUV(1), offset, dist / 2.0);
            auto splineR = Spline(false, REVERSE, cell->initialCenterUV(0), cell->initialCenterUV(1), offset, dist / 2.0);
            model.populateSpline(splineF);
            model.populateSpline(splineR);
            tempRealPointsF = splineF.realPoints;
            tempSPointsF = splineF.splinePoints;
            tempTanF = splineF.tangents;
            tempRealPointsR = splineR.realPoints;
            tempSPointsR = splineR.splinePoints;
            tempTanR = splineR.tangents;
            MatrixXd tempRP, tempTan;
            tempTan = MatrixXd::Zero(tempTanF.rows()+tempTanR.rows(), 3);
            tempTan << tempTanR.colwise().reverse(), tempTanF;
            tempRP = MatrixXd::Zero(tempRealPointsF.rows()+tempRealPointsR.rows(), 3);
            tempRP << tempRealPointsR.colwise().reverse(), tempRealPointsF;
            (*splinePointsO) = genSpline(tempRP, tempTan, 5);

            int rws = splinePointsO->rows();
            int one = round(rws * 0.033);
            int two = round(rws / 3.0),
                three = two * 2,
                four = one ? rws - one : rws - 1;

            tempRP = MatrixXd::Zero(4, 3);
            tempTan = MatrixXd::Zero(0, 3);
            tempRP << splinePointsO->row(one), splinePointsO->row(two), splinePointsO->row(three), splinePointsO->row(four);

            genSpline(tempRP, tempTan);
            (*tangentsO) = MatrixXd::Zero(2, 3);
            (*tangentsO) << tempTan.row(1), tempTan.row(2);
            (*realPointsO) = MatrixXd::Zero(2, 3);
            (*realPointsO) << tempRP.row(1), tempRP.row(2);

            Vector3d normal = cell->initialCenterNormal;

            (*splinePointsO) = MatrixXd::Zero(3, 3);
            (*splinePointsO) << realPointsO->row(0), cell->initialCenter, realPointsO->row(1);
            num disty = 0.0;
            for (int i=0; i<splinePointsO->rows()-1; i++) {
                disty += distance(splinePointsO->row(i), splinePointsO->row(i+1));
            }
            if (disty > 10.0) {
                cout << splineF.uvPoints << endl << splineR.uvPoints << endl;
                cout << *splinePointsO << endl << endl;
            }
            cell->splinePoints.push_back(splinePointsO);
            cell->curvatures.push_back(fullSplineSignedCurvature(realPointsO, tangentsO, normal));
        }
    }
    bar.reset();
    bar.set_label("Analyzing splines");
    c = 0;
    for (auto cell : cells) {
        bar.progress(++c, cells.size());
        int ptCount = 200;
        vector<Cell*> nVec = cell->getNeighborsVec();
        for (int i=0; i<6; i++) {
            num u1, u2, v1, v2, uMult, vMult, u, v;
            auto e1 = cell->getEdges()[i*3], e2 = cell->getEdges()[i*3+1];
            int j = (i + 2) % 6;
            if (!e1->isInitialized && !e2->isInitialized && nVec[i] != nullptr && nVec[j] != nullptr) {
                e1->isInitialized = true;
                e2->isInitialized = true;
                e1->setViaAverage = false;
                e2->setViaAverage = false;
                u1 = nVec[i]->initialCenterUV(0);
                u2 = nVec[j]->initialCenterUV(0);
                v1 = nVec[i]->initialCenterUV(1);
                v2 = nVec[j]->initialCenterUV(1);
                uMult = (u2 - u1) / ptCount;
                vMult = (v2 - v1) / ptCount;
                MatrixXd puv(ptCount, 3), pcart;
                for (int k=0; k<ptCount; k++) {
                    u = uMult * k + u1;
                    v = vMult * k + v1;
                    puv.row(k) << u, v, 0.0;
                }
                e1->splinePointsFull = model.mesh->barycentricLoopingF(puv);
                e2->splinePointsFull = e1->splinePointsFull;
            } else if (!e1->isInitialized && !e2->isInitialized) {
                e1->setViaAverage = true;
                e2->setViaAverage = true;
            }
        }
#pragma omp parallel for num_threads(6)
        for (int i=0; i<6; i++) {
            int minInd;
            auto e1 = cell->getEdges()[i*3], e2 = cell->getEdges()[i*3+1],
            e3 = cell->getEdges()[((i+1)*3)%18], e4 = cell->getEdges()[(((i-1)%6+6)%6)*3];
            Vertex *v11, *v12, *v21, *v22;
            if (e1->getV1() == e3->getV1() || e1->getV1() == e3->getV2()) {
                v11 = e1->getV1();
                v12 = e1->getV2();
                v21 = e2->getV1();
                v22 = e2->getV2();
            } else {
                v11 = e1->getV2();
                v12 = e1->getV1();
                v21 = e2->getV2();
                v22 = e2->getV1();
            }
            if (!e1->isSet && !e1->setViaAverage) {
                minInd = -1;
                num bestDist = std::numeric_limits<num>::max(), currDist;
                if (e3->splinePointsFull.rows()) {
                    for (int k=0; k<e1->splinePointsFull.rows(); k++) {
                        Vector3d pt = e1->splinePointsFull.row(k);
                        for (int w=0; w<e3->splinePointsFull.rows(); w++) {
                            currDist = distance(pt, e3->splinePointsFull.row(w));
                            if (currDist < bestDist) {
                                bestDist = currDist;
                                minInd = k;
                            }
                        }
                    }
                }
                if (minInd < (e1->splinePointsFull.rows() / 2))
                    minInd = 2 * e1->splinePointsFull.rows() / 3;

                v11->cartEstimate += e1->splinePointsFull.row(minInd);
#pragma omp critical
                {
                    v11->cartEstVec.push_back(e1->splinePointsFull.row(minInd));
                }
                v21->cartEstimate += e1->splinePointsFull.row(minInd);
                v11->cartCount++;
                v21->cartCount++;
                e1->splinePoints = e1->splinePointsFull.block(0, 0,
                        minInd+1, 3).eval();
                e2->splinePoints = e1->splinePoints;
                minInd = -1;
                if (e4->splinePointsFull.rows()) {
                    bestDist = std::numeric_limits<num>::max();
                    for (int k=0; k<e1->splinePoints.rows(); k++) {
                        Vector3d pt = e1->splinePoints.row(k);
                        for (int w=0; w<e4->splinePointsFull.rows(); w++) {
                            currDist = distance(pt, e4->splinePointsFull.row(w));
                            if (currDist < bestDist) {
                                bestDist = currDist;
                                minInd = k;
                            }
                        }
                    }
                }
                if (minInd == -1)
                    minInd = e1->splinePoints.rows() / 2;

                v12->cartEstimate += e1->splinePoints.row(minInd);
#pragma omp critical
                {
                    v12->cartEstVec.push_back(e1->splinePoints.row(minInd));
                }
                v22->cartEstimate += e1->splinePoints.row(minInd);
                v12->cartCount++;
                v22->cartCount++;
                e1->splinePoints = e1->splinePoints.block(minInd, 0,
                        e1->splinePoints.rows()-minInd, 3).eval();
                e2->splinePoints = e1->splinePoints;
                e1->isSet = true;
                e2->isSet = true;
            }
        }
    }
    bar.reset();
    for (auto vert : grid.getVertices()) {
        if (vert->cartCount)
            vert->cartEstimate /= static_cast<num>(vert->cartCount);
    }
    for (auto edge : grid.getEdges()) {
        num origDist = distance(edge->getV1()->getPosition(), edge->getV2()->getPosition());
        if (edge->splinePoints.size()) {
            edge->L0 = 0.0;
            for (int i=0; i<edge->splinePoints.rows()-2; i++) {
                edge->L0 += distance(edge->splinePoints.row(i), edge->splinePoints.row(i+1));
            }
            edge->L0 = std::max(origDist * static_cast<num>(0.4), edge->L0);
        } else {
            if (!edge->setViaAverage) {
                edge->isVertical = true;
                edge->type = 2;
                edge->isSet = true;
                edge->L0 = origDist;
            }
        }
    }
    vector<Cell*> cellsToRedo;
    for (auto cell : cells) {
        vector<Edge*> edgesToSet, edges = cell->getEdges();
        num avg = 0.0;
        int count = 0;
        for (int i=0; i<18; i+=3) {
            auto edge = edges[i];
            if (edge->setViaAverage && !edge->isSet) {
                edgesToSet.push_back(edge);
                edgesToSet.push_back(edges[i+1]);
                count++;
            } else {
                avg += edge->L0;
            }
        }
        if (count == 6) {
            cellsToRedo.push_back(cell);
        } else {
            avg /= static_cast<num>(count);
            for (auto edge : edgesToSet) {
                edge->L0 = avg;
                edge->isSet = true;
            }
            cell->setRestLengths();
        }
    }
    for (auto cell : cellsToRedo) {
        vector<Edge*> edgesToSet, edges = cell->getEdges();
        num avg = 0.0;
        int count = 0;
        for (int i=0; i<18; i+=3) {
            auto edge = edges[i];
            if (edge->setViaAverage && !edge->isSet) {
                edgesToSet.push_back(edge);
                edgesToSet.push_back(edges[i+1]);
                count++;
            } else {
                avg += edge->L0;
            }
        }
        avg /= static_cast<num>(count);
        for (auto edge : edgesToSet) {
            edge->L0 = avg;
            edge->isSet = true;
        }
    }
    for (auto edge : grid.getEdges()) {
        num dist = distance(edge->getV1()->getPosition(), edge->getV2()->getPosition());
        edge->L0 = bound(edge->L0, dist * 0.7, dist * 1.5);
    }
    vector<Edge*> nonVertEdgesTop, nonVertEdgesBot, redoTop, redoBot;
    vector<vector<Edge*>*> top, bot;
    top.push_back(&nonVertEdgesTop);
    top.push_back(&redoTop);
    bot.push_back(&nonVertEdgesBot);
    bot.push_back(&redoBot);
    for (auto cell : cells) {
        auto currEdges = cell->getEdges();
        for (unsigned int i=0; i<currEdges.size(); i+=3) {
            if (!currEdges[i]->splinePoints.size() &&
                    !std::count(nonVertEdgesTop.begin(), nonVertEdgesTop.end(), currEdges[i])) {
                nonVertEdgesTop.push_back(currEdges[i]);
                nonVertEdgesBot.push_back(currEdges[i+1]);
            }
        }
    }
    std::set<Cell*> toDelete;
    for (int k=0; k<2; k++) {
        auto arrTop = top[k];
        auto arrBot = bot[k];
        for (unsigned int i=0; i<arrTop->size(); i++) {
            auto edgeTop = (*arrTop)[i];
            auto edgeBot = (*arrBot)[i];
            c = 0;
            num avg = 0.0;
            for (auto cell : edgeTop->cells) {
                for (unsigned int j=0; j<cell->getEdges().size(); j+=3) {
                    auto edge = cell->getEdges()[j];
                    if (!std::count(arrTop->begin(), arrTop->end(), edge)) {
                        avg += edge->L0;
                        c++;
                    }
                }
            }
            if (!c) {
                if (k) {
                    cerr << "WARNING: An edge was not able to be set on the second pass\n";
                    for (auto cell : edgeTop->cells)
                        toDelete.insert(cell);
                } else {
                    redoTop.push_back(edgeTop);
                    redoBot.push_back(edgeBot);
                }
            } else {
                avg /= c;
                edgeTop->L0 = std::min(avg, edgeTop->L0);
                edgeBot->L0 = std::min(avg, edgeBot->L0);
            }
        }
    }
    for (auto cell : cells) {
        c = 0;
        for (auto edge : cell->getEdges()) {
            c += edge->setViaAverage;
        }
        if (c > 6)
            toDelete.insert(cell);
    }
    for (auto cell : cells) {
        cell->setRestLengths();
    }
    for (auto cell : toDelete) {
        grid.cullCellFromGrid(cell);
        cells.erase(std::remove(cells.begin(), cells.end(), cell), cells.end());
    }
    cout << endl << "Ls: " << grid.getEdges()[0]->getCurrentLength() << endl;
    cout << "Interlayer spacing: " << model.zHeight << endl;
}

CellGrid::CellGrid(Model &model, num spacing) : model(model) {
    CellGrid(model, spacing, false);
}

CellGrid::CellGrid(Model &model, num spacing, bool cubic) : model(model) {
    if (cubic)
        return;
    int offset = 0, lastOffset = 0;
    for (int i=0; i<static_cast<int>(model.vSplines.size()); i++) {
        auto spline = *model.vSplines[i];
        auto tans = *model.vTangents[i];
        auto cellCentersUV = *model.vCellCentersUV[i];
        auto cellCentersNormal = *model.vCellCentersNormal[i];
        auto cellCentersCart = *model.vCellCentersCart[i];
        auto realPoints = *model.vRealPoints[i];
        vector<CellWrapper*> colVec;
        num bestInd = -1;
        for (int j=0; j<cellCentersCart.rows(); j++) {
            auto cell = new Cell(cellCentersCart.row(j), cellCentersUV.row(j));
            cell->initialCenterNormal = cellCentersNormal.row(j);
            cell->vSplineIndex = model.splines[i]->index;
            if (i && !j) {
                Vector3d bestP;
                num bestDist = std::numeric_limits<num>::max(), currDist;
                auto topWestCellWrapper = grid[i-1][0];
                auto topWestCell = topWestCellWrapper->cell;
                for (auto p : spline.rowwise()) {
                    currDist = distance(topWestCell->initialCenter, p);
                    if (currDist < bestDist) {
                        bestDist = currDist;
                        bestP = p;
                    }
                }
                Vector2d uv;
                uv << model.mesh->barycentricLoopingR(bestP);
                num diff = uv(0) - cellCentersUV(j, 0);
                bestDist = std::numeric_limits<num>::max();
                if (diff > 0.0) {
                    auto cellCenters = *model.vCellCentersCart[i];
                    for (int z=0; z<cellCenters.rows(); z++) {
                        currDist = distance(cellCenters.row(z), bestP);
                        if (currDist < bestDist) {
                            bestInd = z;
                            bestDist = currDist;
                        }
                    }
                    bestInd -= diff < 0.0;

                    offset = lastOffset - 1 - bestInd * 2;
                } else {
                    auto cellCenters = *model.vCellCentersCart[i-1];
                    for (int z=0; z<cellCenters.rows(); z++) {
                        currDist = distance(cellCenters.row(z), bestP);
                        if (currDist < bestDist) {
                            bestInd = z;
                            bestDist = currDist;
                        }
                    }
                    offset = lastOffset + 1 + bestInd * 2;
                }
            }
            auto wrapper = new CellWrapper;
            wrapper->level = (*model.levels[i])[j];
            wrapper->index = (*model.indexList[i])[j];
            wrapper->cell = cell;
            colVec.push_back(wrapper);
            lastOffset = offset;
        }

        grid.push_back(colVec);
    }

    for (int i=0; i<static_cast<int>(grid.size()); i++) {
        auto col = grid[i];
        for (int j=0; j<static_cast<int>(col.size()); j++) {
            auto wrapper = col[j];
            int level = wrapper->level;
            auto cell = wrapper->cell;
            if (j < static_cast<int>(col.size()) - 1) {
                cell->neighbors.south = col[j+1]->cell;
                col[j+1]->cell->neighbors.north = cell;
            }
            if (j) {
                cell->initialGridCenter = cell->neighbors.north->initialGridCenter;
                cell->initialGridCenter(1) += spacing;
            } else if (i && !j) { 
                cell->initialGridCenter(0) = spacing * sqrt(3.0) / 2.0 * i;
                cell->initialGridCenter(1) = spacing / 2.0 * level;
            } else if (!j && !i) {
                auto temp = new Vector3d;
                (*temp) << 0, spacing / 2.0 * level, 0;
                cell->initialGridCenter = (*temp);
            }
            cell->initialGridCenter(2) = 0.0;

            if (i) {
                auto prevCol = grid[i-1];
                for (int k=0; k < static_cast<int>(grid[i-1].size()); k++) {
                    auto westWrapper = grid[i-1][k];
                    int westLevel = westWrapper->level;
                    auto westCell = westWrapper->cell;
                    if (westLevel == (level - 1)) {
                        cell->neighbors.northWest = westCell;
                        westCell->neighbors.southEast = cell;
                    } else if (westLevel == (level + 1)) {
                        cell->neighbors.southWest = westCell;
                        westCell->neighbors.northEast = cell;
                    }
                }
            }
        }
    }
}

void CellGrid::adjustCellPlacement() {
    int midpt = grid.size() / 2;
    auto col = grid[midpt];
    vector<CellWrapper*> otherCol;
    int colmidpt = col.size() / 2;

    // Spacing middle column relative to itself
    vector<int> toDo;
    for (int i=colmidpt-1; i>=0; i--)
        toDo.push_back(i);
    for (unsigned int i=colmidpt+1; i<col.size(); i++)
        toDo.push_back(i);
    for (int i : toDo) {
        auto comp = col[i-1];
        if (i < colmidpt)
            comp = col[i+1];
        model.calcIdealPosition(col[i]->cell->initialCenterUV,
                col[i]->cell->initialCenter, comp->cell->initialCenterUV, model.spacing);
    }

    tqdm bar;
    bar.set_label("Recentering Cells");
    // Spacing the remaining columns relative to their more central neighbor
    toDo.clear();
    num s = model.spacing * COS_30;
    for (int i=midpt-1; i>=0; i--)
        toDo.push_back(i);
    for (unsigned int i=midpt+1; i<grid.size(); i++)
        toDo.push_back(i);
    int k=0;
    for (int i : toDo) {
        bar.progress(k, toDo.size());
        k++;
        col = grid[i];
        if (i < midpt)
            otherCol = grid[i+1];
        else
            otherCol = grid[i-1];
        unsigned int ind = 0, otherColInd = 0;
        CellWrapper *currCell = col[0];
        if (currCell->level < otherCol[0]->level) {
            ind = (otherCol[0]->level + 1 - currCell->level) / 2;
        } else if (otherCol[0]->level + 1 != currCell->level) {
            otherColInd = (currCell->level - otherCol[0]->level - 1) / 2;
        }
        bool done = false;
        while (!done) {
            Vector2d pt = (otherCol[otherColInd]->cell->initialCenterUV +
                    otherCol[otherColInd+1]->cell->initialCenterUV) / 2.0;
            model.calcIdealPosition(col[ind]->cell->initialCenterUV,
                    col[ind]->cell->initialCenter, pt, s);
            ind++;
            otherColInd++;
            done = otherColInd >= (otherCol.size() - 1) ||
                ind >= col.size();
        }
    }
}

void Tissue::addInitialCell() {
    auto cell = new Cell();
    cells.push_back(cell);
    
    MatrixXd initPts = MatrixXd::Zero(12, 3);
    for (int i=0; i<6; i++) {
        num ang = i * PI / 3;
        initPts.row(i) << cos(ang), sin(ang), 0;
        initPts.row(i+6) << cos(ang), sin(ang), 1;
    }
    auto vertices = grid.addVertices(cell, initPts);
    vector<vector<Vertex*>> vertPairs;
    for (int i=0; i<6; i++) {
        auto pair1 = new vector<Vertex*>();
        auto pair2 = new vector<Vertex*>();
        auto pair3 = new vector<Vertex*>();
        pair1->push_back(vertices[i]);
        pair1->push_back(vertices[(i+1)%6]);
        pair2->push_back(vertices[i+6]);
        pair2->push_back(vertices[(i+1)%6+6]);
        pair3->push_back(vertices[i]);
        pair3->push_back(vertices[i+6]);
        vertPairs.push_back(*pair1);
        vertPairs.push_back(*pair2);
        vertPairs.push_back(*pair3);
    }
    grid.addEdges(cell, vertPairs);
}

vector<Cell*> &Tissue::getCells() {
    return cells;
}

SO::Matrix3X &Tissue::getVertexMat() {
    return grid.getVertexMat();
}

vector<std::shared_ptr<SO::Constraint>> &Tissue::getSOEdges() {
    return grid.getSOEdges();
}

vector<std::shared_ptr<SO::Constraint>> &Tissue::getSOAnglesTemp() {
    return sOAnglesT;
}

vector<std::shared_ptr<SO::Constraint>> &Tissue::getSOAnglesCubic() {
    int i, j, nU = model.nodeCountU, nV = model.nodeCountV;
    uint l;
    auto sOAngles = new vector<std::shared_ptr<SO::Constraint>>();
    auto points = grid.getVertexMatV();
    auto verts = grid.getVerticesV();
    num err;
    std::vector<num> errVert, errHoriz;
    int vs = verts.size();
    const num angTol = 0.2, angle = H_PI;
    int lVerts = vs / model.layersCubic;
    for (l=0; l<(model.hasDuplicatedCubic ? model.layersCubic : 1); l++) {
        int layerOffset = l * lVerts;
        for (i=0; i<nU-1; i++) {
            for (j=0; j<nV-1; j++) {
                vector<int> pts, pts2, pts3, pts4;
                auto one = verts[i*nV+j+layerOffset];
                auto two = verts[i*nV+j+1+layerOffset];
                auto three = verts[(i+1)*nV+j+layerOffset];
                auto four = verts[(i+1)*nV+j+1+layerOffset];
                pts.push_back(one->getIndex());
                pts.push_back(two->getIndex());
                pts.push_back(three->getIndex());
                pts2.push_back(two->getIndex());
                pts2.push_back(four->getIndex());
                pts2.push_back(one->getIndex());
                pts3.push_back(three->getIndex());
                pts3.push_back(one->getIndex());
                pts3.push_back(four->getIndex());
                pts4.push_back(four->getIndex());
                pts4.push_back(three->getIndex());
                pts4.push_back(two->getIndex());
                for (auto a : { pts, pts2, pts3, pts4 }) {
                    errHoriz.push_back(abs(vectAngle(points.col(a[1]) -
                            points.col(a[0]), points.col(a[2]) -
                            points.col(a[0])) - H_PI) / H_PI);
                }
                auto constraint = make_shared<SO::AngleConstraint>(
                    pts, grid.kAng, points, angle - angTol, angle + angTol);
                auto constraint2 = make_shared<SO::AngleConstraint>(
                    pts2, grid.kAng, points, angle - angTol, angle + angTol);
                auto constraint3 = make_shared<SO::AngleConstraint>(
                    pts3, grid.kAng, points, angle-angTol, angle+angTol);
                auto constraint4 = make_shared<SO::AngleConstraint>(
                    pts4, grid.kAng, points, angle-angTol, angle+angTol);
                sOAngles->push_back(constraint2);
            }
        }
    }
    for (l=0; l<(model.hasDuplicatedCubic ? model.layersCubic : 0); l++) {
        for (i=0; i<nU; i++) {
            for (j=0; j<nV; j++) {
                vector<int> indsH, indsV;
                int currI = i*nV+j+lVerts;
                auto one = verts[currI];
                if (j<nV-1)
                    indsH.push_back(currI+1);
                if (j)
                    indsH.push_back(currI-1);
                if (i<nU-1)
                    indsH.push_back((i+1)*nV+j+lVerts);
                if (i)
                    indsH.push_back((i-1)*nV+j+lVerts);
                if (l<model.layersCubic-2)
                    indsV.push_back(currI+lVerts);
                if (l)
                    indsV.push_back(currI-lVerts);
                for (int h : indsH) {
                    for (int v : indsV) {
                        auto two = verts[h];
                        auto three = verts[v];
                        vector<int> pts;
                        pts.push_back(one->getIndex());
                        pts.push_back(two->getIndex());
                        pts.push_back(three->getIndex());
                        errVert.push_back(abs(vectAngle(points.col(pts[1]) -
                                points.col(pts[0]), points.col(pts[2]) -
                                points.col(pts[0])) - H_PI) / H_PI);
                        auto constraint = make_shared<SO::AngleConstraint>(
                            pts, grid.kAng90, points, angle - angTol, angle + angTol);
                        sOAngles->push_back(constraint);
                    }
                }
            }
        }
    }
    err = 0.0;
    for (num e : errVert)
        err += e;
    if (err > 1e-2) {
        cout << "Vert degree error percent: " << err / errVert.size() * 100 << "%\n";
        cout << "Vert degree stdev: " << stdev(errVert) << endl;
    }

    err = 0.0;
    for (num e : errHoriz)
        err += e;
    if (err > 1e-2) {
        cout << "Horiz degree error percent: " << err / errHoriz.size() * 100 << "%\n";
        cout << "Horiz degree stdev: " << stdev(errHoriz) << endl;
    }
    return *sOAngles;
}

vector<std::shared_ptr<SO::Constraint>> &Tissue::getSOAngles() {
    if (cubic)
        return getSOAnglesCubic();
    auto sOAngles = new vector<std::shared_ptr<SO::Constraint>>();
    auto points = grid.getVertexMat();
    const num angTol = 0.01;
    num angle = 2.0944, err;
    std::vector<num> angles = { angle, H_PI, H_PI };
    std::vector<num> err90, err120;
    std::vector<std::vector<num>*> errs = { &err120, &err90, &err90 };
    int k = 0;
    for (auto cell : cells) {
        vector<vector<int>> hexAngs;
        vector<num> percentsCell;
        auto verts = cell->getVertices();
        for (int i=0; i<CELL_POLYGON+1; i+=CELL_POLYGON) {
            for (int j=0; j<CELL_POLYGON; j++) {
                vector<int> pts, pts2, pts3, pts4;
                auto one = verts[i+j];
                auto two = verts[i+(j+1)%CELL_POLYGON];
                auto three = verts[i+(j+2)%CELL_POLYGON];
                pts.push_back(two->getIndex());
                pts.push_back(one->getIndex());
                pts.push_back(three->getIndex());
                hexAngs.push_back(pts);
                int modResult = (j-1)%6;
                modResult = modResult < 0 ? 5 : modResult;
                int otherLayer = !i ? 6 : -6;
                pts3.push_back(verts[i+j]->getIndex());
                pts3.push_back(verts[i+j+otherLayer]->getIndex());
                pts3.push_back(verts[i+modResult]->getIndex());
                pts4.push_back(verts[i+j]->getIndex());
                pts4.push_back(verts[i+j+otherLayer]->getIndex());
                pts4.push_back(verts[i+(j+1)%6]->getIndex());
                k = 0;
                for (auto a : { pts, pts3, pts4 }) {
                    errs[k]->push_back(abs(vectAngle(points.col(a[1]) -
                            points.col(a[0]), points.col(a[2]) -
                            points.col(a[0])) - angles[k]) / angles[k]);
                    k++;
                }
                auto constraint = make_shared<SO::AngleConstraint>(
                    pts, grid.kAng, points, angle - angTol, angle + angTol);
                auto constraint3 = make_shared<SO::AngleConstraint>(
                    pts3, grid.kAng90, points, H_PI-angTol, H_PI+angTol);
                auto constraint4 = make_shared<SO::AngleConstraint>(
                    pts4, grid.kAng90, points, H_PI-angTol, H_PI+angTol);
                sOAngles->push_back(constraint);
                sOAngles->push_back(constraint3);
                sOAngles->push_back(constraint4);
            }
        }
    }
    err = 0.0;
    for (num e : err90)
        err += e;
    if (err > 1e-2) {
        cout << "90 degree error percent: " << err / err90.size() * 100 << "%\n";
        cout << "90 degree stdev: " << stdev(err90) << endl;
    }

    err = 0.0;
    for (num e : err120)
        err += e;
    if (err > 1e-2) {
        cout << "120 degree error percent: " << err / err120.size() * 100 << "%\n";
        cout << "120 degree stdev: " << stdev(err120) << endl;
    }
    return *sOAngles;
}

void Tissue::updateGrid(SO::Matrix3X &newPoints) {
    grid.setVertexMat(newPoints);
}

void Tissue::importJSONFull() {
    importJSONFull(setts["model"]["save_file"].value_or(fmt::format("../exports/saves/model_{}.json", time(0))));
}

void Tissue::importJSONFull(string filename) {
    std::ifstream f(filename);
    json save = json::parse(f);
    f.close();
    bool same = true;
    vector<num> pCorners;
    string inpType = setts["model"]["type"].value_or("obj");
    if (!inpType.compare("obj")) {
        for (long unsigned int i=0; i<setts["model"]["param_corners"].as_array()->size(); i++) {
            same &= setts["model"]["param_corners"][i].value_or(-1.0) == save["param_corners"][i];
        }
    }
    same &= save["initial_node_count"] == setts["isolines"]["initial_node_count"].value_or(-1);
    same &= save["spline_bump_count"] == setts["isolines"]["spline_bump_count"].value_or(-1);
    same &= save["primary_u_value"] == setts["isolines"]["primary_u_value"].value_or(-1);
    same &= save["u_spline_offset"] == setts["isolines"]["u_spline_offset"].value_or(-1);
    same &= save["v_spline_offset"] == setts["isolines"]["v_spline_offset"].value_or(-1);
    same &= save["z_height"] == setts["grid"]["z_height"].value_or(-1);
    same &= save["z_height_mult"] == setts["grid"]["z_height_mult"].value_or(-1);
    same &= save["z_height_div"] == setts["grid"]["z_height_div"].value_or(-1);
    if (same) {
        importJSON(save);
    } else {
        // TODO: Fix this so it just regenerates the tissue
        cerr << "ERROR: Saved model's config does not match current config\n";
        exit(1);
    }
}

void Tissue::importJSON(json save) {
    // For some reason if I try to set model.scaleFactor here it gets deleted. There's some issue
    // with the reference to model being different in this method call that I can't figure out,
    // so this is a work around.
    scaleFactor = save["scale_factor"];
    grid.importJSON(save);
    std::unordered_map<int, Cell*> cellMap;
    cellMap[-1] = nullptr;
    for (auto c : save["cells"]) {
        auto cell = new Cell();
        vector<Edge*> cellEdges;
        vector<Vertex*> cellVerts;
        for (auto i : c["vertices"])
            cellVerts.push_back(grid.vertMap[i]);
        for (auto i : c["edges"])
            cellEdges.push_back(grid.edgeMap[i]);
        cell->addVertices(cellVerts);
        cell->index = c["index"];
        cell->initialGridCenter << c["cc"][0], c["cc"][1], c["cc"][2];
        cellMap[cell->index] = cell;
        grid.addCellToEdges(cell, cellEdges);
        cells.push_back(cell);
    }
    for (uint i=0; i<save["cells"].size(); i++) {
        auto n = save["cells"][i]["neighbors"];
        auto cell = cells[i];
        cell->neighbors.northEast = cellMap[n[0]];
        cell->neighbors.southEast = cellMap[n[1]];
        cell->neighbors.south = cellMap[n[2]];
        cell->neighbors.southWest = cellMap[n[3]];
        cell->neighbors.northWest = cellMap[n[4]];
        cell->neighbors.north = cellMap[n[5]];
    }
}

Grid &Tissue::getGrid() {
    return grid;
}

void Tissue::exportJSON() {
    exportJSON(setts["model"]["save_file"].value_or(fmt::format("../exports/saves/model_{}.json", time(0))));
}

void Tissue::exportJSON(string filename) {
    json output;
    for (auto v : grid.vertices) {
        json entry;
        auto p = v->getPosition();
        entry["position"] = {p(0), p(1), p(2)};
        entry["index"] = v->getIndex();
        output["vertices"].push_back(entry);
    }
    for (auto e : grid.getEdges()) {
        json entry;
        entry["rest_length"] = e->getRestLength();
        entry["index"] = e->getIndex();
        entry["vertices"] = {e->getV1()->getIndex(), e->getV2()->getIndex()};
        entry["type"] = e->type;
        output["edges"].push_back(entry);
    }
    for (auto c : cells) {
        json entry;
        entry["index"] = c->index;
        entry["cc"] = c->initialGridCenter;
        for (auto e : c->getEdges())
            entry["edges"].push_back(e->getIndex());
        for (auto v : c->getVertices())
            entry["vertices"].push_back(v->getIndex());
        for (auto n : c->getNeighborsVec()) {
            int ind = n != nullptr ? n->index : -1;
            entry["neighbors"].push_back(ind);
        }
        output["cells"].push_back(entry);
    }

    vector<num> pCorners;
    for (num p : model.paramCorners) {
        pCorners.push_back(p);
    }
    output["param_corners"] = pCorners;
    output["scale_factor"] = model.scaleFactor;
    output["initial_node_count"] = setts["isolines"]["initial_node_count"].value_or(-1);
    output["spline_bump_count"] = setts["isolines"]["spline_bump_count"].value_or(-1);
    output["primary_u_value"] = setts["isolines"]["primary_u_value"].value_or(-1);
    output["u_spline_offset"] = setts["isolines"]["u_spline_offset"].value_or(-1);
    output["v_spline_offset"] = setts["isolines"]["v_spline_offset"].value_or(-1);
    output["z_height"] = setts["grid"]["z_height"].value_or(-1);
    output["z_height_mult"] = setts["grid"]["z_height_mult"].value_or(-1);
    output["z_height_div"] = setts["grid"]["z_height_div"].value_or(-1);
}
