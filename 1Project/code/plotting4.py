import os
import json

data_folder = "Data"

master_dict = {}
decoder = json.JSONDecoder()

# ==========================================================
# loop through all json files
# ==========================================================
filepaths = [
    os.path.join(data_folder, f)
    for f in os.listdir(data_folder)
    if os.path.isfile(os.path.join(data_folder, f))
    and f.endswith(".json")
]

print(f"Found {len(filepaths)} JSON files")

for filepath in filepaths:

    filename = os.path.basename(filepath)

    # ------------------------------------------------------
    # scrape metadata dynamically from filename
    # ------------------------------------------------------
    name_no_ext = os.path.splitext(filename)[0]
    parts = name_no_ext.split("__")

    metadata = {}

    for part in parts:
        if "_" in part:
            key, value = part.split("_", 1)
            metadata[key] = value

    # skip malformed filenames
    if not all(k in metadata for k in ["D", "Hz", "I"]):
        print(f"Skipping {filename}")
        continue

    try:
        J = float(metadata["D"])
        Hz = float(metadata["Hz"])
        I = int(metadata["I"])
    except ValueError:
        print(f"Malformed metadata in {filename}")
        continue

    print(f"Loading {filename}")

    # ------------------------------------------------------
    # initialize nested dictionary
    # master_dict[J][Hz][I]
    # ------------------------------------------------------
    if J not in master_dict:
        master_dict[J] = {}

    if Hz not in master_dict[J]:
        master_dict[J][Hz] = {}

    master_dict[J][Hz][I] = {}

    # ------------------------------------------------------
    # load concatenated json objects
    # ------------------------------------------------------
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
            master_dict[J][Hz][I][step] = obj

        except json.JSONDecodeError as e:
            print(f"Error in {filename} near position {pos}: {e}")
            break

print("Finished loading.")