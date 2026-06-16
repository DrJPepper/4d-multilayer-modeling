# Generates an OBJ based on a given function

import sys
import math
import random
import pickle
import os
import argparse

from lib.simpleeval import SimpleEval, simple_eval

import numpy as np
from pathlib import Path
from functools import reduce
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

DEFAULT_TESS_LINES = 5

def iter_dict(d, ks):
    return reduce(lambda a, b: a[b], ks, d)

def parse_args():
    parser = argparse.ArgumentParser(
              prog = 'python3 gen_surface.py',
              description = 'Generate a tessellated based on a function',
              )
    default_config_file = 'surface.toml'
    inputs = [
            ('c', 'config_file', default_config_file, None),
            ('x', 'expression', 'sin(x)', ['function', 'expression']),
            ('s', 'start', 0, ['function', 'start']),
            ('e', 'end', '4*pi', ['function', 'end']),
            ('t', 'thickness', 4, ['surface', 'thickness']),
            ('l', 'tess_lines', DEFAULT_TESS_LINES, ['surface', 'tess_lines']),
            ('a', 'line_spacing', None, ['surface', 'line_spacing']),
            ('p', 'point_count', 10, ['surface', 'point_count'])
        ]
    for i in inputs:
        parser.add_argument(f'-{i[0]}', f'--{i[1]}', required=False)
    args, extras = parser.parse_known_args()
    args = vars(args)
    if not args['config_file']:
        if len(extras):
            args['config_file'] = extras[0]
        elif Path(default_config_file).is_file():
            args['config_file'] = default_config_file
        else:
            args['config_file'] = ''
    if args['config_file']:
        setts = tomllib.loads(Path(args['config_file']).read_text(encoding='utf-8'))
        for i in inputs:
            if i[3] and args[i[1]] is None:
                try:
                    args[i[1]] = iter_dict(setts, i[3])
                except KeyError:
                    pass
    for i in inputs:
        if args[i[1]] is None:
            args[i[1]] = i[2]
    return args

def triangle_normal(p1, p2, p3):
    """
    Calculate the normal vector of a triangle defined by points p1, p2, p3.

    Parameters:
        p1, p2, p3 : array-like of shape (3,)
            The 3D coordinates of the triangle's vertices.

    Returns:
        normal : ndarray of shape (3,)
            A normalized vector perpendicular to the triangle's face.
    """
    p1, p2, p3 = np.array(p1), np.array(p2), np.array(p3)
    
    v1 = p2 - p1
    v2 = p3 - p1
    
    normal = np.cross(v1, v2)
    
    norm = np.linalg.norm(normal)
    if norm == 0:
        raise ValueError("The points are collinear; the normal is undefined.")
    return normal / norm

def flip(x):
    if x >= 2.0:
        return -1
    return 1

def step(x):
    return x%2

def main():
    args = parse_args()
    s = SimpleEval()
    funcs = {"sin": math.sin, "sqrt": math.sqrt, "step": step, "flip": flip}
    names = {"pi": math.pi}
    s.functions = funcs
    expr = s.parse(args['expression'])
    to_parse = ['start', 'end']
    for i in to_parse:
        args[i] = simple_eval(str(args[i]), functions=funcs, names=names)
    if args['line_spacing']:
        # NOTE: If the config file is set to the default this won't print,
        # which I don't think is a big enough deal to account for
        if args['tess_lines'] != DEFAULT_TESS_LINES:
            print('WARNING: both line_spacing and tess_lines given,' +
                  ' line_spacing will be used')
        tess_lines = args['thickness'] / args['line_spacing']
        i_tess_lines = int(tess_lines)
        thickness = i_tess_lines * args['line_spacing']
        if round(tess_lines, 3) != float(i_tess_lines):
            print('WARNING: thickness divided by line_spacing is not an' +
                  f' integer, thickness will be shrunk to {thickness}')
        args['thickness'] = thickness
        args['tess_lines'] = i_tess_lines
    tlr = range(args['tess_lines'])
    lines = [[] for _ in tlr]
    points = []
    for i in range(args['point_count']):
        x = args['start'] + i / args['point_count'] * (args['end'] - args['start'])
        s.names = {'x': x} | names
        y = s.eval(args["expression"], previously_parsed=expr)
        print(f'{x}, {y}')
        for j in tlr:
            z = float(j) / (args['tess_lines'] - 1) * args['thickness']
            lines[j].append([y, x, z])
    for l in lines:
        points += l
    faces = []
    normals = []
    for i in range(args['tess_lines']-1):
        for j in range(args['point_count']-1):
            faces.append([i * args['point_count'] + j,
                          (i + 1) * args['point_count'] + j,
                          i * args['point_count'] + j + 1])
            normals.append(triangle_normal(points[faces[-1][0]], points[faces[-1][1]], points[faces[-1][2]]))
            faces.append([i * args['point_count'] + j + 1,
                          (i + 1) * args['point_count'] + j,
                          (i + 1) * args['point_count'] + j + 1])
            normals.append(triangle_normal(points[faces[-1][0]], points[faces[-1][1]], points[faces[-1][2]]))
    v_normals = [np.array([0.0, 0.0, 0.0]) for _ in points]
    vc = [0] * len(points)
    for i in range(len(faces)):
        f = faces[i]
        n = normals[i]
        for j in f:
            v_normals[j] += n
            vc[j] += 1
    for i in range(len(vc)):
        v_normals[i] /= vc[i]
    v_normals = [(-1.0*i).tolist() for i in v_normals]
    with open('out.obj', 'w') as f:
        for i in points:
            f.write(f'v {i[1]} {i[0]} {i[2]}\n')
        for i in v_normals:
            f.write(f'vn {i[0]} {i[1]} {i[2]}\n')
        for i in faces:
            f.write(f'f {i[0]+1}//{i[0]+1} {i[1]+1}//{i[1]+1} {i[2]+1}//{i[2]+1}\n')

if __name__ == '__main__':
    main()
