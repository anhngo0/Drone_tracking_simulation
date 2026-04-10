import json
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

from matplotlib.lines import Line2D

from mpl_toolkits.mplot3d.art3d import Poly3DCollection

with open("../intersect.json") as f:
    J = json.load(f)

fig = plt.figure(figsize=(16,16))
ax = fig.add_subplot(111, projection="3d")
ax.set_xlabel("X",fontsize=30, labelpad=25)
ax.set_ylabel("Z",fontsize=30, labelpad=25)
ax.set_zlabel("Y",fontsize=30, labelpad=25)
ax.tick_params(axis='x', pad=15, labelsize=25)
ax.tick_params(axis='y', pad=15, labelsize=25)
ax.tick_params(axis='z', pad=15, labelsize=25)

def plot3(ax, x, y, z, *args, **kwargs):
    return ax.plot(x, z, y, *args, **kwargs)

def scatter3(ax, x, y, z, *args, **kwargs):
    return ax.scatter(x, z, y, *args, **kwargs)

# -------------------------------------------------
# Draw polyhedron
# -------------------------------------------------
def draw_poly(poly, color):
    V = {v["id"]: v["pos"] for v in poly["vertices"]}
    for f in poly["faces"]:
        xs, ys, zs = zip(*(V[i] for i in f["v"] + [f["v"][0]]))
        plot3(ax,xs, ys, zs, color=color)
# -------------------------------------------------
# Draw P and Q
# -------------------------------------------------
polyhedra = J["polyhedra"]
# for cam_id, poly in polyhedra.items():
#     draw_poly(J["polyhedra"][cam_id], "blue")

def draw_solid_poly(poly, face_color, edge_color="k", alpha=0.4):
    # Lấy vertex array
    V = np.array([v["pos"] for v in poly["vertices"]])

    # Lấy các mặt
    faces = [
        [V[i] for i in f["v"]]
        for f in poly["faces"]
    ]

    solid = Poly3DCollection(
        faces,
        facecolor=face_color,
        edgecolor=edge_color,
        linewidths=0.5,
        alpha=alpha
    )
    ax.add_collection3d(solid)

# draw_solid_poly(
#     J["intersection"]["intersect"],
#     face_color="lightgray",
#     alpha=0.7
# )

#---------------------------------------------

# -------------------------------------------------
# Intersection point
# -------------------------------------------------
x, y, z = J["intersection"]["O"]
# scatter3(ax,[x], [y], [z], color="black", s=60)

V_int = np.array([
v["pos"] for v in J["intersection"]["intersect"]["vertices"]
])

# scatter3(ax,
#     V_int[:, 0],
#     V_int[:, 1],
#     V_int[:, 2],
#     color="red",
#     s=20,
#     depthshade=True
# )


with open("../trajectory.json", "r") as f:
    data = json.load(f)

# === estimate data === 
with open("../two_traj.json", "r") as f:
    est_data = json.load(f)

points = np.array(data["points"])
segments = data["segments"]

# ===== Load motion =====
with open("../trajectory_motion.json") as f:
    motion = json.load(f)

frames = motion["frames"]

# interpolation points
scatter3(ax,points[:,0], points[:,1], points[:,2],
           color="black", s=60)

# #========== Camera =========
cameras = np.array([
    [0, 2.5, 0],    # Camera 1
    [60, 2.5, 0],   # Camera 2
    [60, 2.5, 60],  # Camera 3
    [0, 2.5, 60]    # Camera 4
])

scatter3(ax,
    cameras[:, 0],
    cameras[:, 1],
    cameras[:, 2],
    color='blue',
    s=100,
    label='Camera'
)

# # ================= Straight axis line =================
P0 = points[0]
P1 = points[-1]

plot3(ax,[P0[0], P1[0]],
        [P0[1], P1[1]],
        [P0[2], P1[2]],
        color="lightgray", linewidth=1, label="Axis")

# # ================= Cubic spline eval =================
def eval_cubic(coeff, tau):
    a, b, c, d = coeff
    return a + b*tau + c*tau**2 + d*tau**3

distance = np.linalg.norm(P1-P0)
n= len(segments)
dx = distance/n

traj = []
for i, seg in enumerate(segments):
    # 0 <= tau <= dx = distance/n
    taus = np.linspace(0, dx, 40)
    for tau in taus:
        x = eval_cubic(seg["x"], tau)
        y = eval_cubic(seg["y"], tau)
        z = eval_cubic(seg["z"], tau)
        traj.append([x, y, z])

traj = np.array(traj)

# spline trajectory
plot3(ax,traj[:,0], traj[:,1], traj[:,2],
        color="black", linewidth=2, label="Interpolated trajectory")

# ================= estimate data ===========================
est_traj = []
for x, y, z in est_data["est_pos"]:
    est_traj.append([x, y, z])

est_traj = np.array(est_traj)

# spline trajectory
# scatter3(ax,est_traj[:,0], est_traj[:,1], est_traj[:,2],
#         color="blue", s = 5, label="estimate trajectory")
# plot3(ax,est_traj[:,0], est_traj[:,1], est_traj[:,2],
#         color="blue", linewidth=2, label="estimate trajectory")
# ================= 10 evenly spaced centers =================

N = 10
centers = np.linspace(P0, P1, N)

# direction vector
w = P1 - P0
w = w / np.linalg.norm(w)

# find perpendicular basis (u, v)
tmp = np.array([0, 0, 1])
if abs(np.dot(tmp, w)) > 0.9:
    tmp = np.array([0, 1, 0])

u = np.cross(w, tmp)
u = u / np.linalg.norm(u)
v = np.cross(w, u)

scatter3(ax,
    centers[:,0],
    centers[:,1],
    centers[:,2],
    color="lightgray",
    s=30,           
    depthshade=False,
    label="Circle centers"
)

# # ================= Draw circles =================
R = 50.0
theta = np.linspace(0, 2*np.pi, 120)

for C in centers:
    circle = np.array([
        C + R*(np.cos(t)*u + np.sin(t)*v)
        for t in theta
    ])
    plot3(ax,circle[:,0], circle[:,1], circle[:,2],
            color="gray", linewidth=0.8)


# moving object
dot, = plot3(ax,[], [], [], "ro", markersize=8)

# ===== Animation =====
# def update(frame_id):
#     if frame_id >= len(frames):
#         return dot,
#     p = frames[frame_id]["pos"]
#     dot.set_data([p[0]], [p[1]])
#     dot.set_3d_properties([p[2]])
#     print(frame_id)
#     return dot,

# ani = FuncAnimation(
#     fig,
#     update,
#     frames=len(frames),
#     interval=1*1000,
#     blit=False,
#     # repeat=False
# )

# ---- Custom legend (black dot & gray line) ----
# legend_elements = [
#     Line2D([0], [0], marker='o', color='gray',
#            linestyle='None', markersize=2, label='Estimate Position'),
#     Line2D([0], [0], color='black', lw=2, label='Trajectory')
# ]

legend_elements = [
    Line2D([0], [0], marker='o', color='blue',
           linestyle='None', markersize=10, label='Camera'),
    Line2D([0], [0], color='black', lw=3, label='Trajectory')
]

ax.legend(handles=legend_elements, fontsize=30)

# ===== View setup =====

# ===== Legend =====
legend_elements = [
    Line2D([0], [0], marker='o', color='blue',
           linestyle='None', markersize=10, label='Camera'),
    Line2D([0], [0], color='black', lw=3, label='Trajectory')
]

ax.legend(handles=legend_elements, fontsize=30)

# ===== Save image =====

ax.set_box_aspect([1,1,1])
# ax.autoscale_view()

# ax.legend(handles=legend_elements)
plt.grid(False)
# ax.grid(False)

# ax.xaxis._axinfo["grid"]["linewidth"] = 0
# ax.yaxis._axinfo["grid"]["linewidth"] = 0
# ax.zaxis._axinfo["grid"]["linewidth"] = 0

ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False

ax.zaxis._axinfo['tick']['inward_factor'] = 0
ax.zaxis._axinfo['tick']['outward_factor'] = 0.4
ax.yaxis._axinfo['tick']['inward_factor'] = 0
ax.yaxis._axinfo['tick']['outward_factor'] = 0.4

from matplotlib.ticker import MaxNLocator
ax.zaxis.set_major_locator(MaxNLocator(4))
ax.yaxis.set_major_locator(MaxNLocator(4))


plt.tight_layout()
# plt.subplots_adjust(left=0, right=1, bottom=0, top=1)

# ax.view_init(elev=15, azim=-90)
# ax.view_init(elev=90, azim=-90)
# plt.savefig("trajectory_plot.pdf", bbox_inches='tight', pad_inches=0)

plt.show()
