"""
Crop annotated cells out of RV-PBS fields (RQ4 cross-source analysis).
Each CVAT polygon -> square crop centred on the polygon's bounding box, side = 1.5 x the longer bbox side
(CHOICE: KU-Optofil's own cropping rule is not documented; 1.5x keeps the whole cell plus a margin of background),
resized to 368 x 368 to match KU-Optofil crop size. Output: data/rv_pbs_crops/<id>.jpg + rv_pbs_crops.csv
Labels are RV-PBS's own; 'blast cell' is mapped to is_blast = 1. Image id is kept so analyses can cluster by image.
"""
import glob, os, xml.etree.ElementTree as ET
import numpy as np, pandas as pd
from PIL import Image
OUT = "data/rv_pbs_crops"; os.makedirs(OUT, exist_ok=True); rows = []
for x in sorted(glob.glob("data/rv_pbs/*/annotations.xml")):
    folder = os.path.dirname(x); root = ET.parse(x).getroot()
    for im in root.findall("image"):
        path = os.path.join(folder, im.get("name"))
        if not os.path.exists(path): continue
        img = None
        for j, poly in enumerate(im.findall("polygon")):
            lab = poly.get("label")
            if lab in (None, "bg"): continue
            pts = np.array([[float(v) for v in p.split(",")] for p in poly.get("points").split(";")])
            x0, y0 = pts.min(0); x1, y1 = pts.max(0); cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
            half = 0.75 * max(x1 - x0, y1 - y0)
            if img is None: img = Image.open(path).convert("RGB")
            box = tuple(int(round(v)) for v in (cx - half, cy - half, cx + half, cy + half))
            crop = img.crop(box).resize((368, 368), Image.BILINEAR)
            cid = f"{os.path.basename(folder).replace(' ', '_')}__{os.path.splitext(im.get('name'))[0]}__{j}"
            crop.save(f"{OUT}/{cid}.jpg", quality=95)
            rows.append(dict(crop_id=cid, image=f"{os.path.basename(folder)}/{im.get('name')}", label=lab,
                             is_blast=int(lab == "blast cell"), bbox_side=float(max(x1 - x0, y1 - y0)), folder=os.path.basename(folder)))
pd.DataFrame(rows).to_csv("data/rv_pbs_crops.csv", index=False)
print(len(rows), "crops;", sum(r["is_blast"] for r in rows), "blasts;", len({r['image'] for r in rows}), "images")
