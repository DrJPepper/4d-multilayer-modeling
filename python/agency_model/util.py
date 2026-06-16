import sys
import json
from collections import defaultdict
from copy import deepcopy

import requests
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

import torch
import torch_geometric.nn as pygnn
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.nn import ModuleDict, LazyLinear
from torch_geometric.data import Data, HeteroData, InMemoryDataset
import torch_geometric.transforms as T
from torch_geometric.loader import DataLoader
from torch_geometric.nn import HeteroConv, Linear, to_hetero, NNConv, MessagePassing
from torch_geometric.nn.pool import global_mean_pool
from torch_geometric.nn import AttentionalAggregation
from torch_geometric.utils import to_dense_batch

DEVICE = 'cuda'

DIRS = ["NorthEast", "SouthEast", "South", "SouthWest", "NorthWest", "North"]
LAYERS = ['apical', 'basal']

ft = lambda x: torch.tensor(x, dtype=torch.float32)

def o7(i):
    o = [0] * 7
    o[i] = 1
    return o
