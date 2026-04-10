import json
import numpy as np
import matplotlib.pyplot as plt

# ===== 1. Load JSON file =====
with open("../rmse_based_frame_tracking_0_035_6.json", "r") as f:
    data = json.load(f)

rmse_data = data["rmse_per_frame"]

# ===== 2. Extract frame & rmse =====
frames = [entry["frame"] for entry in rmse_data]
rmse_values = [entry["rmse"] for entry in rmse_data]

# ===== 3. Create full frame range =====
full_frames = list(range(min(frames), max(frames) + 1))

# map frame → rmse
rmse_dict = dict(zip(frames, rmse_values))

# chèn NaN nếu thiếu frame
rmse_full = [rmse_dict.get(f, np.nan) for f in full_frames]

# ===== 4. Plot =====
plt.figure(figsize=(8,5))
plt.plot(full_frames, rmse_full, linewidth=1, label="Proposed Algorithm")

plt.xlabel("Frame")
plt.ylabel("RMSE (m)")
plt.title("RMSE per Frame")

plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()