import json
import os
import re
import signal
import subprocess
import shlex
from math import inf
from pathlib import Path
from time import sleep
from typing import Optional, List, Any, Tuple

from fastapi import BackgroundTasks, FastAPI
from pydantic import BaseModel, model_validator
from uvicorn.server import Server
import logging

from model import *
from dataset import *

logging.basicConfig(level=logging.WARNING)
uvicorn_access = logging.getLogger("uvicorn.access")
uvicorn_access.disabled = True

class Directions(BaseModel):
    North: int
    NorthEast: int
    NorthWest: int
    South: int
    SouthEast: int
    SouthWest: int


class EdgeMap(BaseModel):
    apical: Directions
    basal: Directions


class Cell(BaseModel):
    edge_map: EdgeMap
    edges: List[int] = []
    index: int
    neighbor_map: Directions
    neighbors: List[int] = []
    process_order: Optional[List[int]] = []


class Edge(BaseModel):
    index: int
    length: Optional[float] = None
    lengths: Optional[List[float]] = []

    @model_validator(mode='after')
    @classmethod
    def check_one_length_value(cls, model):
        if (model.length is None) == (len(model.lengths) == 0):
            raise ValueError('Exactly one of "length" or "lengths" must be provided for an edge')
        return model

class EdgeBasic(BaseModel):
    index: int
    length: float

class EdgeOutput(BaseModel):
    index: int
    action: int

class InputData(BaseModel):
    cell: int
    edges: List[EdgeBasic] = []
    signal: List[Tuple[int, int]]
    epoch: int

class Test(BaseModel):
    hi: str

class Edges(BaseModel):
    edges: List[Edge] = []

app = FastAPI()
app.should_exit = False

@app.get("")
@app.get("/")
def test_route():
    return {"message": "Server is up"}

@app.get("/stat_model")
@app.get("/stat_model/")
def stat_model_route():
    return {'step_size': ld.rest_length_step, 'cutoff': ld.cutoff}

@app.post("/run_model")
@app.post("/run_model/")
def run_model_route(inp: InputData):
    if inp.epoch != run_model_route.epochs:
        run_model_route.epochs = inp.epoch
        run_model_route.total = 0
        run_model_route.correct = 0
    data = ld.cell_map[inp.cell].to(DEVICE)
    edge_order = ld.edge_map[inp.cell]
    edge_lengths = [-1] * len(edge_order.keys())
    for e in inp.edges:
        edge_lengths[edge_order[e.index]] = (e.length - ld.data_mean) / ld.data_std
    signals = dict(inp.signal)
    inds = data['cell'].indices.tolist()
    for i in range(len(inds)):
        ind = inds[i]
        sig = signals[ind] / 10.0
        data['cell'].x[i, 0] = sig
    inds = data['edge'].indices.flatten()[data.train_mask].tolist()
    yl = []
    for i in range(len(inds)):
        ind = inds[i]
        c = edge_lengths[edge_order[ind]]
        n = ld.targets[ind]
        diff = c - n
        if abs(diff) < ld.cutoff:
            yl.append(0)
        elif diff < 0.0:
            yl.append(1)
        else:
            yl.append(2)
    data['edge'].x = ft([[i] for i in edge_lengths]).to(DEVICE)
    target = torch.tensor(yl).type(torch.LongTensor).to(DEVICE)
    with torch.no_grad():
        out, signal_out = model(data)

    run_model_route.total += len(out)
    run_model_route.correct += (out == target).count_nonzero()
    indices = data['edge'].indices[data.train_mask].flatten().tolist()
    output = {'signal': signal_out.item(), 'edges': []}
    for i in range(len(out)):
        output['edges'].append(EdgeOutput(index=indices[i], action=out[i]))
    print(output)
    return output

config = CONFIG
ld = LocalDataset('./training_data.json')
dataset, y_list = ld.build_td()
count = 10
model = HeteroGNN(dataset[0].metadata(), config, dataset[0]).to(DEVICE)
model.load_state_dict(torch.load('./graph_model_agency.pth'))
model.eval()
run_model_route.total = 0
run_model_route.correct = 0
run_model_route.epochs = 0

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("__main__:app", host="0.0.0.0", port=8000, reload=True, log_config=None)
