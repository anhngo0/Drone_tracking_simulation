# import json
# import numpy as np
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib.animation import FuncAnimation

# from matplotlib.lines import Line2D

# from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# # ===== load json =====
# filename = "../pf_json/0_pf.json"   # đổi thành file của bạn

# with open(filename, "r") as f:
#     data = json.load(f)

# # ===== extract particles =====
# px = []
# py = []
# pz = []

# for p in data["particles"]:
#     pos = p["pos"]
#     px.append(pos[0])
#     py.append(pos[1])
#     pz.append(pos[2])

# # ===== extract positions =====
# real = data["pos"]
# est = data["est_pos"]

# # ==== poly fov =====
# def draw_poly(poly, color):
#     V = {v["id"]: v["pos"] for v in poly["vertices"]}
#     for f in poly["faces"]:
#         xs, ys, zs = zip(*(V[i] for i in f["v"] + [f["v"][0]]))
#         ax.plot(xs, ys, zs, color=color)

# # ===== plot =====

# fig = plt.figure(figsize=(30, 30))
# ax = fig.add_subplot(111, projection="3d")

# # particles
# ax.scatter(px, py, pz, color="lightblue", s=5, label="particles")

# # real position
# ax.scatter(real[0], real[1], real[2], color="red", s=80, label="pos")

# # estimated position
# ax.scatter(est[0], est[1], est[2], color="black", s=80, label="est_pos")

# # fov polys
# polyhedra = data["polyhedra"]
# for cam_id, poly in polyhedra.items():
#     draw_poly(data["polyhedra"][cam_id], "lightgray")

# # labels
# ax.set_xlabel("X")
# ax.set_ylabel("Y")
# ax.set_zlabel("Z")

# ax.legend()
# plt.show()

import json
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ===== danh sách file =====
files = [f"../pf_json/{i}_pf.json" for i in range(21)]

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(111, projection="3d")

# ===== poly drawing =====
def draw_poly(ax, poly, color):
    V = {v["id"]: v["pos"] for v in poly["vertices"]}
    for f in poly["faces"]:
        xs, ys, zs = zip(*(V[i] for i in f["v"] + [f["v"][0]]))
        ax.plot(xs, ys, zs, color=color)

# ===== update function =====
def update(frame):

    ax.cla()  # clear axis

    filename = files[frame]

    with open(filename, "r") as f:
        data = json.load(f)

    # particles
    px, py, pz = [], [], []
    for p in data["particles"]:
        pos = p["pos"]
        px.append(pos[0])
        py.append(pos[1])
        pz.append(pos[2])

    # positions
    real = data["pos"]
    est = data["est_pos"]

    # draw particles
    ax.scatter(px, py, pz, color="lightblue", s=5, label="particles")

    # real position
    ax.scatter(real[0], real[1], real[2], color="red", s=80, label="pos")

    # estimate
    ax.scatter(est[0], est[1], est[2], color="black", s=80, label="est_pos")

    # draw FOV polyhedra
    polyhedra = data["polyhedra"]
    for cam_id, poly in polyhedra.items():
        draw_poly(ax, poly, "gray")

    ax.set_title(f"Frame {frame}")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    ax.legend()

# ===== animation =====
ani = FuncAnimation(
    fig,
    update,
    frames=len(files),
    interval=2000,  # 5 seconds
    repeat=False
)

plt.show()