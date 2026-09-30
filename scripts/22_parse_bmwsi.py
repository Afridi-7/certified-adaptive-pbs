"""
Parse BoneMarrowWSI-PediatricLeukemia annotation series (DICOM Microscopy Bulk Simple Annotations) into tables.
Input: data/bmwsi/ann (extracted from data/bmwsi/bmwsi_ann.zip, produced by notebooks/21_bmwsi_annotations_colab.ipynb).
Output: data/bmwsi/cells.csv (one row per annotated box: patient, source series, label, box, centre, roi index)
        data/bmwsi/rois.csv  (one row per monolayer ROI rectangle)
Coordinates are in slide pixel units as stored in the annotation (2D, 0.11 um/px per the dataset description).
"""
import glob, os
import numpy as np, pandas as pd, pydicom, highdicom as hd

D = "data/bmwsi/ann"; cells, rois = [], []
for f in sorted(glob.glob(f"{D}/bonemarrowwsi_pediatricleukemia/**/*.dcm", recursive=True)):
    ds = pydicom.dcmread(f); desc = ds.SeriesDescription; pid = ds.PatientID
    ann = hd.ann.MicroscopyBulkSimpleAnnotations.from_dataset(ds)
    src = "roi" if desc.startswith("Monolayer") else ("consensus" if desc.startswith("Consensus") else "unlabeled")
    for g in ann.get_annotation_groups():
        for c in g.get_graphic_data(coordinate_type=ds.AnnotationCoordinateType):
            c = np.asarray(c, float); x0, y0 = c[:, 0].min(), c[:, 1].min(); x1, y1 = c[:, 0].max(), c[:, 1].max()
            row = dict(patient=pid, series=ds.SeriesInstanceUID, x0=x0, y0=y0, x1=x1, y1=y1)
            if src == "roi": rois.append(row)
            else: cells.append(dict(row, source=src, label=g.label, meaning=g.annotated_property_type.meaning,
                                    cx=(x0 + x1) / 2, cy=(y0 + y1) / 2))
C = pd.DataFrame(cells); R = pd.DataFrame(rois)
R["roi"] = R.groupby("patient").cumcount()
C["roi"] = -1
for pid, g in R.groupby("patient"):
    m = C.patient == pid
    for _, r in g.iterrows():
        inside = m & (C.cx >= r.x0) & (C.cx <= r.x1) & (C.cy >= r.y0) & (C.cy <= r.y1)
        C.loc[inside, "roi"] = r.roi
C.to_csv("data/bmwsi/cells.csv", index=False); R.to_csv("data/bmwsi/rois.csv", index=False)
print("cells", len(C), "| rois", len(R), "| patients with rois", R.patient.nunique())
print(C.groupby("source").agg(boxes=("label", "size"), patients=("patient", "nunique")))
