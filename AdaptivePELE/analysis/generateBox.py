from __future__ import print_function, unicode_literals
import os
import argparse
import numpy as np
from AdaptivePELE.atomset import atomset
from AdaptivePELE.utilities import utilities


def parse_arguments():
    """
        Create command-line interface
    """
    desc = "Create a visualization for a cylindrical box"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-t", "--top", type=float, nargs="*", default=None, help="Coordinates for the bottom base")
    parser.add_argument("-b", "--bottom", type=float, nargs="*", default=None, help="Coordinates for the top base")
    parser.add_argument("-c", "--center", type=float, nargs="*", default=None, help="Coordinates for the center")
    parser.add_argument("--type", type=str, default="cylinder", help="Type of box to build (cylinder or sphere)")
    parser.add_argument("--radius", type=float, default=10, help="Radius of the cylinder base")
    parser.add_argument("--n_points", type=int, default=10, help="Number of points to write, the larger it is the finer the discretization")
    parser.add_argument("--title", type=str, default="box.pdb", help="Name to store the output")
    parser.add_argument("--output", type=str, default="", help="Path to store the output")
    parser.add_argument("--from_pdb", type=str, default=None, help="Construct the box from a input pdb")
    args = parser.parse_args()
    return args.top, args.bottom, args.radius, args.title, args.output, args.n_points, args.from_pdb, args.center, args.type


def generate_perpendicular_vector(input_vec):
    dim = len(input_vec)
    # generate two random coordinates
    p = np.random.rand(dim-1)
    # select a non-zero component in the input vector
    ind = np.where(input_vec)[0][0]
    random_indices = list(range(dim))
    random_indices.remove(ind)
    new_vec = np.zeros_like(input_vec)
    new_vec[ind] = -input_vec[random_indices].dot(p)/input_vec[ind]
    new_vec[random_indices] = p
    new_vec /= np.linalg.norm(new_vec)
    return new_vec


def generate_cylinder_box(top, bot, radius, n_points):
    top = np.array(top)
    bot = np.array(bot)
    l = np.linalg.norm(top-bot)
    alpha = top-bot
    alpha /= np.linalg.norm(alpha)
    alpha_rot = generate_perpendicular_vector(alpha)*radius
    angle = 2*np.pi/n_points
    rot_mat = utilities.generateRotationMatrixAroundAxis(alpha, angle)
    rotated_vectors = []
    new_alpha = alpha_rot
    for i in range(n_points):
        new_alpha = np.dot(rot_mat, new_alpha)
        assert np.abs(np.linalg.norm(new_alpha) - radius) < 1e-4, (np.linalg.norm(new_alpha), radius)
        rotated_vectors.append(new_alpha)

    coords = [top, bot, bot+alpha_rot, bot-alpha_rot, top+alpha_rot, top-alpha_rot]
    points = [bot+i*alpha for i in np.linspace(0, l, num=n_points)]
    for point in points:
        coords += [point+alpha_rot, point-alpha_rot]+[point+alpha_2 for alpha_2 in rotated_vectors]
    return coords


def generate_sphere_box(center, radius, n_points):
    center = np.array(center)
    coords = [center]
    # start from x axis
    point = np.zeros(3)
    point[0] = 1
    angle = 2*np.pi/n_points
    angle_ax = 2*np.pi/n_points
    coords.append(center+point*radius)
    ax_rot_mat = utilities.generateRotationMatrixAroundAxis(point, angle_ax)
    # rotate around y axis
    alpha = np.zeros(3)
    alpha[1] = 1
    for k in range(n_points):
        if k:
            # rotate around x axis to obtain new rotation axis
            alpha = np.dot(ax_rot_mat, alpha)
        rot_mat = utilities.generateRotationMatrixAroundAxis(alpha, angle)
        for j in range(n_points):
            point = np.dot(rot_mat,point)
            coords.append(center+radius*point)
    return coords

def main(top, bot, radius, title="box.pdb", output="", n_points=10, input_pdb=None, center=None, type_box="cylinder"):
    if output:
        if not os.path.exists(output):
            os.makedirs(output)
        title = os.path.join(output, title)
    if input_pdb is not None:
        in_PDB = atomset.PDB()
        in_PDB.initialise(input_pdb, resname="DUM")
        for atom in in_PDB.atoms.values():
            if atom.name == "DUMT":
                top = atom.getAtomCoords()
            if atom.name == "DUMB":
                bot = atom.getAtomCoords()
            if atom.name == "DUM":
                center = atom.getAtomCoords()
    if type_box == "cylinder":
        if top is None or bot is None:
            raise utilities.RequiredParameterMissingException("Coordinates for the top and bottom base center not specified!!!")
        coords = generate_cylinder_box(top, bot, radius, n_points)
    elif type_box == "sphere":
        if center is None:
            raise utilities.RequiredParameterMissingException("Coordinates for the center not specified!!!")
        coords = generate_sphere_box(center, radius, n_points)
    coords = np.array(coords)
    utilities.write_PDB_clusters(coords, title)


if __name__ == "__main__":
    t, b, rad, file_name, out, n, pdb_in, c, type_box = parse_arguments()
    main(t, b, rad, file_name, out, n, pdb_in, c, type_box)
