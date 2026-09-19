"""
Step 1 — Inspect KU-Optofil (read-only). Prints facts we must VERIFY before designing composites.
Usage (from the thesis folder):  python scripts/01_inspect_kuoptofil.py --data data/ku_optofil
Paste the whole console output back into the chat.
"""
import argparse, os, zipfile, collections, random
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--data", default="data/ku_optofil")
a = ap.parse_args()
D = a.data
pd.set_option("display.width", 200); pd.set_option("display.max_columns", 50)

for name in ["metadata.csv", "metadata_with_patient_level_splits.csv"]:
    p = os.path.join(D, name)
    print("\n" + "=" * 80 + f"\n{name}")
    if not os.path.exists(p):
        print("  MISSING"); continue
    df = pd.read_csv(p)
    print(f"  rows={len(df)}  cols={list(df.columns)}")
    print(df.head(5).to_string())
    for c in df.columns:
        nun = df[c].nunique()
        if nun <= 40:
            print(f"\n  value counts of '{c}' ({nun} unique):")
            print(df[c].value_counts(dropna=False).to_string())
        else:
            print(f"\n  '{c}': {nun} unique values, e.g. {df[c].dropna().astype(str).unique()[:3].tolist()}")
    # Look for columns that could encode field / slide / coordinates
    hints = [c for c in df.columns if any(h in c.lower() for h in
             ["patient", "slide", "smear", "field", "image", "x", "y", "bbox", "coord", "split", "fold", "stain", "date"])]
    print(f"\n  columns possibly relevant to patients/fields/coords/splits: {hints}")

# Zip / folder structure
zp = os.path.join(D, "dataset.zip")
root = os.path.join(D, "dataset")
print("\n" + "=" * 80)
files = []
if os.path.isdir(root):
    for r, _, fs in os.walk(root):
        for f in fs: files.append(os.path.relpath(os.path.join(r, f), root))
    print(f"extracted folder: {len(files)} files")
elif os.path.exists(zp):
    with zipfile.ZipFile(zp) as z: files = [n for n in z.namelist() if not n.endswith("/")]
    print(f"dataset.zip (not extracted): {len(files)} files")
ext = collections.Counter(os.path.splitext(f)[1].lower() for f in files)
print("extensions:", dict(ext))
top = collections.Counter(f.replace("\\", "/").split("/")[0] for f in files)
print("top-level entries:", dict(list(top.items())[:30]))
lvl2 = collections.Counter("/".join(f.replace("\\", "/").split("/")[:2]) for f in files)
print("2nd-level entries (first 40):"); [print("  ", k, v) for k, v in list(lvl2.items())[:40]]
print("random file names:", random.Random(0).sample(files, min(8, len(files))))

# Image sizes on a sample
try:
    from PIL import Image
    imgs = [f for f in files if os.path.splitext(f)[1].lower() in (".png", ".jpg", ".jpeg", ".tif", ".tiff", ".bmp")]
    samp = random.Random(1).sample(imgs, min(200, len(imgs)))
    sizes = collections.Counter()
    if os.path.isdir(root):
        for f in samp: sizes[Image.open(os.path.join(root, f)).size] += 1
    else:
        with zipfile.ZipFile(zp) as z:
            for f in samp: sizes[Image.open(z.open(f)).size] += 1
    print("image sizes (sample of 200):", sizes.most_common(10))
except Exception as e:
    print("image size check skipped:", e)
print("\nDONE")
