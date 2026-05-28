import os
import json
import numpy as np

data_folder = "Data"

master_dict = {}
decoder = json.JSONDecoder()

# --------------------------------------------------
# LOOP THROUGH ALL FILES
# --------------------------------------------------
for filename in os.listdir(data_folder):

    filepath = os.path.join(data_folder, filename)

    if not os.path.isfile(filepath):
        continue

    # only process known file types
    if not (filename.endswith(".json") or filename.endswith(".txt")):
        continue

    # --------------------------------------------------
    # scrape metadata dynamically from filename
    # --------------------------------------------------
    name_no_ext = os.path.splitext(filename)[0]
    parts = name_no_ext.split("__")

    metadata = {}

    for part in parts:
        if "_" in part:
            key, value = part.split("_", 1)
            metadata[key] = value

    # skip files without required keys
    if not all(k in metadata for k in ["J", "Hz", "I"]):
        continue

    try:
        J = float(metadata["J"])
        Hz = float(metadata["Hz"])
        I = int(metadata["I"])
    except ValueError:
        print(f"Skipping malformed filename: {filename}")
        continue

    # --------------------------------------------------
    # create nested structure
    # master_dict[J][Hz][I]
    # --------------------------------------------------
    if J not in master_dict:
        master_dict[J] = {}

    if Hz not in master_dict[J]:
        master_dict[J][Hz] = {}

    if I not in master_dict[J][Hz]:
        master_dict[J][Hz][I] = {
            "skyrmions": {},
            "info": {}
        }

    entry = master_dict[J][Hz][I]

    # ==================================================
    # SKYRMION JSON FILE
    # ==================================================
    if filename.endswith(".json"):

        print(f"Loading skyrmions: {filename}")

        with open(filepath, "r") as f:
            content = f.read().strip()

        pos = 0

        while pos < len(content):

            while pos < len(content) and content[pos].isspace():
                pos += 1

            if pos >= len(content):
                break

            try:
                obj, idx = decoder.raw_decode(content[pos:])
                pos += idx

                step = obj["step"]

                # store original object
                entry["skyrmions"][step] = obj

            except json.JSONDecodeError as e:
                print(f"Error in {filename} near position {pos}: {e}")
                break

    # ==================================================
    # ENERGY / INFO TXT FILE
    # ==================================================
    elif filename.endswith(".txt"):

        print(f"Loading energy info: {filename}")

        try:
            data = np.genfromtxt(filepath, skip_header=1)

            if data.ndim == 1:
                data = data.reshape(1, -1)

            entry["info"] = {
                "steps": np.array(data[:, 0]),
                "Energy": np.array(data[:, 1]),
                "acceptance": np.array(data[:, 2]),
                "beta": np.array(data[:, 3]),
                "Q": np.array(data[:, 4]),
                "J": np.array(data[:, 5]),
                "D": np.array(data[:, 6]),
                "Hz": np.array(data[:, 7]),
                "I": np.array(data[:, 8]),
            }

        except Exception as e:
            print(f"Failed to load {filename}: {e}")

#%%

