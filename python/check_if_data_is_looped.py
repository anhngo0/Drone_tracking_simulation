import json
from collections import Counter

# ===== Load file =====
with open("../rmse_based_frame_tracking_0_035.json", "r") as f:
    data = json.load(f)

rmse_data = data["rmse_per_frame"]

# ===== Extract frames =====
frames = [entry["frame"] for entry in rmse_data]

# ===== Count occurrences =====
frame_counts = Counter(frames)

# ===== Find duplicates =====
duplicates = {frame: count for frame, count in frame_counts.items() if count > 1}

# ===== Print result =====
if duplicates:
    print("Duplicate frames found:")
    for frame, count in duplicates.items():
        print(f"Frame {frame} appears {count} times")
else:
    print("No duplicate frames found.")