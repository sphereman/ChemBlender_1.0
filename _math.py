import numpy as np
from . import mesh
from math import pi

def get_length(vector):
    length = np.linalg.norm(vector)
    return length

def normalize(vector):
    return vector/get_length(vector) if get_length(vector) else vector

# rotate a vector around an axis
def rotate_vec(vec, axis, angle):
    u = normalize(axis)
    x, y, z = u[0], u[1], u[2]
    cos = np.cos(angle)
    sin = np.sin(angle)
    R = np.array([[cos+x*x*(1-cos), x*y*(1-cos)-z*sin, x*z*(1-cos)+y*sin],
                  [x*y*(1-cos)+z*sin, cos+y*y*(1-cos), y*z*(1-cos)-x*sin],
                  [x*z*(1-cos)-y*sin, y*z*(1-cos)+x*sin, cos+z*z*(1-cos)],])
    return np.dot(R, vec)

# calculate a vertical vector of the vec
def vertical_vec(pos1, pos2):
    v = normalize(pos2-pos1)
    if v[0] == 0 and v[1] == 0:
        vertical = np.array((1,0,0))
    else:
        vertical = normalize((-v[1], v[0], 0))
    return vertical

# get heading, pitch, bank vectors of a vert
def vert_hpb(vert, rotate_angle):
    vecs = []
    for edge in vert.link_edges:
        link_vert = mesh.get_link_vert(vert, edge)
        vecs.append(normalize(np.array(vert.co) - np.array(link_vert.co)))
    bank = np.array((0.0, 0.0, 0.0))
    for vec in vecs: bank += vec
    bank = normalize(bank)
    # 此处有待优化
    # -------------
    if len(vecs)>=3:
        if abs(np.dot(np.cross(vecs[0],vecs[1]),vecs[2])) <= 0.0025:   # coplanar vecs
            bank = np.cross(vecs[0],vecs[1])
    # -------------
    if len(vecs) == 0:
        bank = np.array((1, 0, 0))
        pitch = np.array((0, 1, 0))
    elif len(vecs) == 1:
        pitch = vertical_vec(np.array((0,0,0)), bank)
    elif len(vecs) >= 2:
        pitch = normalize(vecs[1]-vecs[0])
    pitch = rotate_vec(pitch, bank, rotate_angle)
    heading = np.cross(bank, pitch)
    return (bank, pitch, heading)


def branches_dir(vert, branches, deflection_angle):
    bank, pitch, heading = vert_hpb(vert, rotate_angle=pi/2)
    banks = []
    pitches = []
    if branches == 1:
        bank = rotate_vec(bank, heading, deflection_angle*pi/180)
        pitch = np.cross(heading, bank)
        banks.append(bank)
        pitches.append(pitch)
    elif branches == 2:
        bank_1 = rotate_vec(bank, heading, deflection_angle*pi/180)
        bank_2 = rotate_vec(bank_1, bank, pi)
        pitch_1 = np.cross(heading, bank_1)
        pitch_2 = np.cross(heading, bank_2)
        banks.append(bank_1)
        banks.append(bank_2)
        pitches.append(pitch_1)
        pitches.append(pitch_2)
    elif branches == 3:
        bank_1 = normalize(bank+2.824*pitch)
        bank_2 = rotate_vec(bank_1, bank, pi*2/3)
        bank_3 = rotate_vec(bank_1, bank, pi*4/3)
        pitch_1 = np.cross(heading, bank_1)
        pitch_2 = rotate_vec(pitch_1, bank, pi*2/3)
        pitch_3 = rotate_vec(pitch_1, bank, pi*4/3)
        banks.append(bank_1)
        banks.append(bank_2)
        banks.append(bank_3)
        pitches.append(pitch_1)
        pitches.append(pitch_2)
        pitches.append(pitch_3)
    elif branches == 4:
        bank_1 = pitch
        bank_2 = rotate_vec(bank_1, bank, pi*2/3)
        bank_3 = rotate_vec(bank_1, bank, pi*4/3)
        pitch_1 = np.cross(bank, bank_1)
        pitch_2 = np.cross(bank, bank_2)
        pitch_3 = np.cross(bank, bank_3)
        banks.append(bank_1)
        banks.append(bank_2)
        banks.append(bank_3)
        banks.append(bank)
        pitches.append(pitch_1)
        pitches.append(pitch_2)
        pitches.append(pitch_3)
        pitches.append(pitch)
    elif branches == 5:
        bank_1 = pitch
        bank_2 = rotate_vec(bank_1, bank, pi/2)
        bank_3 = rotate_vec(bank_1, bank, pi)
        bank_4 = rotate_vec(bank_1, bank, pi*3/2)
        pitch_1 = np.cross(bank, bank_1)
        pitch_2 = np.cross(bank, bank_2)
        pitch_3 = np.cross(bank, bank_3)
        pitch_4 = np.cross(bank, bank_4)
        banks.append(bank_1)
        banks.append(bank_2)
        banks.append(bank_3)
        banks.append(bank_4)
        banks.append(bank)
        pitches.append(pitch_1)
        pitches.append(pitch_2)
        pitches.append(pitch_3)
        pitches.append(pitch_4)
        pitches.append(pitch)
    return (banks, pitches)
# type: ignore