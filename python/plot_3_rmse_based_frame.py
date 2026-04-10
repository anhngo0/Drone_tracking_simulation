import json
import numpy as np
import matplotlib.pyplot as plt

# ===== Config =====
files = [
    ("../rmse_based_frame_tracking_0_00873_8.json", "σ = 0.5°", "black"),
    ("../rmse_based_frame_tracking_0_01745_8.json", "σ = 1°", "red"),
    ("../rmse_based_frame_tracking_0_035_8.json", "σ = 2°", "blue"),
]

plt.figure(figsize=(12, 8))

# ===== Loop qua từng file =====
for file_path, label, color in files:
    with open(file_path, "r") as f:
        data = json.load(f)

    rmse_data = data["rmse_per_frame"]

    frames = [entry["frame"] for entry in rmse_data]
    rmse_values = [entry["rmse"] for entry in rmse_data]

    # full frame range
    full_frames = list(range(min(frames), max(frames) + 1))

    rmse_dict = dict(zip(frames, rmse_values))
    rmse_full = [rmse_dict.get(f, np.nan) for f in full_frames]

    plt.plot(full_frames, rmse_full, linewidth=2, color=color, label=label)

# ===== Style =====
plt.xlabel("Frame", fontsize=32)
plt.ylabel("RMSE (m)", fontsize=32)
    
plt.title("RMSE per Frame (Tracking)", fontsize=32)

plt.xticks(fontsize=28)
plt.yticks(fontsize=28)

plt.legend(fontsize=28)
plt.grid(True)

plt.tight_layout()

plt.savefig("tracking_plot.pdf", bbox_inches='tight', pad_inches=0)
# ===== Show =====
plt.show()