import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from joblib import dump
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    confusion_matrix, accuracy_score, precision_score,
    recall_score, f1_score
)
import seaborn as sns

# ======================================================
# CONFIG
# ======================================================

DATA_PATH = r"C:\Users\Bellamkonda Sowmya\OneDrive\Desktop\2nd yr notes\ps1\authml\authml\SMILES_Big_Data_Set.csv"

LOGP_MODEL = "logp_model.pkl"
PIC50_MODEL = "pic50_model.pkl"

TEST_SPLIT = 0.2
RANDOM_SEED = 42


# ======================================================
# EXACT DESCRIPTORS USED BY app.py
# ======================================================

def desc_logp(mol):
    return [
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        Descriptors.NumHDonors(mol),
        Descriptors.NumHAcceptors(mol),
        Descriptors.TPSA(mol),
    ]


def desc_pic50(mol):
    return [
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        Descriptors.NumHDonors(mol),
        Descriptors.NumHAcceptors(mol),
        Descriptors.NumRotatableBonds(mol),
        Descriptors.TPSA(mol),
    ]


# ======================================================
# SAFE FEATURE BUILDER
# ======================================================

def build_features(smiles_series):
    clean_smiles = []
    descX_logp = []
    descX_pic = []
    bad_rows = []

    for idx, smi in enumerate(smiles_series):
        mol = Chem.MolFromSmiles(smi)

        if mol is None:
            bad_rows.append(idx)
            continue

        clean_smiles.append(smi)
        descX_logp.append(desc_logp(mol))
        descX_pic.append(desc_pic50(mol))

    return (
        np.array(descX_logp),
        np.array(descX_pic),
        np.array(clean_smiles),
        bad_rows
    )


def train_rf(X, y):
    xtr, xts, ytr, yts = train_test_split(
        X, y,
        test_size=TEST_SPLIT,
        random_state=RANDOM_SEED,
    )

    model = RandomForestRegressor(
        n_estimators=700,
        max_depth=None,
        n_jobs=-1,
        random_state=RANDOM_SEED,
    )
    model.fit(xtr, ytr)
    pred = model.predict(xts)

    rmse = np.sqrt(mean_squared_error(yts, pred))
    r2 = r2_score(yts, pred)

    return model, rmse, r2, yts, pred


def show_plot(true_vals, predicted_vals, title):
    plt.figure()
    plt.scatter(true_vals, predicted_vals)
    plt.plot(
        [true_vals.min(), true_vals.max()],
        [true_vals.min(), true_vals.max()],
    )
    plt.xlabel("True Values")
    plt.ylabel("Predicted Values")
    plt.title(title)
    plt.grid(True)
    plt.show()


# ======================================================
# LOAD CSV AND CLEAN DATA
# ======================================================

df = pd.read_csv(DATA_PATH)

required = ["SMILES", "logP", "pIC50"]
for col in required:
    if col not in df.columns:
        raise Exception(f"Missing column '{col}' in CSV")

df = df[required]

# Drop rows with NaN targets
df = df.dropna(subset=["logP", "pIC50"]).reset_index(drop=True)


# Build descriptors, filter out invalid SMILES
X_logp_raw, X_pic_raw, clean_SMILES, bad_rows = build_features(df["SMILES"])

df_clean = df.drop(index=bad_rows).reset_index(drop=True)

logp_target = df_clean["logP"].astype(float)
pic50_target = df_clean["pIC50"].astype(float)


# ======================================================
# TRAIN LOGP
# ======================================================

print("=== LOGP MODEL TRAINING ===")
print("Valid rows:", len(df_clean))
print("Removed bad SMILES:", len(bad_rows))
print("LOGP Feature Shape:", X_logp_raw.shape)

logp_model, l_rmse, l_r2, ytest_l, pred_l = train_rf(X_logp_raw, logp_target)

dump(logp_model, LOGP_MODEL)
print("\nLOGP Results:")
print("RMSE:", l_rmse)
print("R2  :", l_r2)
print("Saved:", LOGP_MODEL)

show_plot(
    ytest_l,
    pred_l,
    f"LOGP Prediction Accuracy\nR2={l_r2:.3f}, RMSE={l_rmse:.3f}"
)


# ======================================================
# TRAIN PIC50
# ======================================================

print("\n=== pIC50 MODEL TRAINING ===")
print("pIC50 Feature Shape:", X_pic_raw.shape)

pic_model, p_rmse, p_r2, ytest_p, pred_p = train_rf(X_pic_raw, pic50_target)

dump(pic_model, PIC50_MODEL)
print("\npIC50 Results:")
print("RMSE:", p_rmse)
print("R2  :", p_r2)
print("Saved:", PIC50_MODEL)

show_plot(
    ytest_p,
    pred_p,
    f"pIC50 Prediction Accuracy\nR2={p_r2:.3f}, RMSE={p_rmse:.3f}"
)
# Convert continuous LogP → classes
def class_logp(v):
    if v < 2: return 0
    elif v < 4: return 1
    return 2

ytest_l_class = np.array([class_logp(v) for v in ytest_l])
pred_l_class = np.array([class_logp(v) for v in pred_l])

cm_logp = confusion_matrix(ytest_l_class, pred_l_class)

print("\n===== LOGP CLASSIFICATION METRICS =====")
print("Accuracy :", accuracy_score(ytest_l_class, pred_l_class))
print("Precision:", precision_score(ytest_l_class, pred_l_class, average='weighted'))
print("Recall   :", recall_score(ytest_l_class, pred_l_class, average='weighted'))
print("F1 Score :", f1_score(ytest_l_class, pred_l_class, average='weighted'))

# Plot Confusion Matrix (LOGP)
plt.figure(figsize=(6,5))
sns.heatmap(cm_logp, annot=True, cmap="Blues", fmt="d")
plt.title("LOGP Confusion Matrix")
plt.xlabel("Predicted")
plt.ylabel("Actual")
plt.show()

# Metric Bar Plot (LOGP)
metrics_l = ["Accuracy", "Precision", "Recall", "F1 Score"]
values_l = [
    accuracy_score(ytest_l_class, pred_l_class),
    precision_score(ytest_l_class, pred_l_class, average='weighted'),
    recall_score(ytest_l_class, pred_l_class, average='weighted'),
    f1_score(ytest_l_class, pred_l_class, average='weighted')
]
plt.figure(figsize=(7,5))
sns.barplot(x=metrics_l, y=values_l, palette="Blues")
plt.title("LOGP Classification Metrics")
plt.ylim(0, 1)
plt.show()

# -------------- PIC50 CLASSIFICATION --------------

# Convert continuous pIC50 → binary class: Active ≥ 6
ytest_p_class = np.array([1 if v >= 6 else 0 for v in ytest_p])
pred_p_class = np.array([1 if v >= 6 else 0 for v in pred_p])

cm_pic = confusion_matrix(ytest_p_class, pred_p_class)

print("\n===== pIC50 CLASSIFICATION METRICS =====")
print("Accuracy :", accuracy_score(ytest_p_class, pred_p_class))
print("Precision:", precision_score(ytest_p_class, pred_p_class))
print("Recall   :", recall_score(ytest_p_class, pred_p_class))
print("F1 Score :", f1_score(ytest_p_class, pred_p_class))

# Plot Confusion Matrix (pIC50)
plt.figure(figsize=(6,5))
sns.heatmap(cm_pic, annot=True, cmap="Greens", fmt="d")
plt.title("pIC50 Confusion Matrix")
plt.xlabel("Predicted")
plt.ylabel("Actual")
plt.show()

# Metric Bar Plot (pIC50)
metrics_p = ["Accuracy", "Precision", "Recall", "F1 Score"]
values_p = [
    accuracy_score(ytest_p_class, pred_p_class),
    precision_score(ytest_p_class, pred_p_class),
    recall_score(ytest_p_class, pred_p_class),
    f1_score(ytest_p_class, pred_p_class)
]
plt.figure(figsize=(7,5))
sns.barplot(x=metrics_p, y=values_p, palette="Greens")
plt.title("pIC50 Classification Metrics")
plt.ylim(0, 1)
plt.show()