import os
import subprocess

folder = "models"
proteins = ["thrombin", "protease", "kinase"]

OBABEL = r"C:\Program Files\OpenBabel-3.1.1\obabel.exe"

for name in proteins:
    pdb = os.path.join(folder, f"{name}.pdb")
    pdbqt = os.path.join(folder, f"{name}.pdbqt")

    print(f"Converting {pdb} → {pdbqt}")

    cmd = [
        OBABEL,
        "-ipdb", pdb,
        "-opdbqt",
        "-O", pdbqt,
        "--gen3d"
    ]

    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        print(f"❌ Failed for {name}: {result.stderr}")
    else:
        print(f"✔ Success for {name}")

print("🔥 All conversions complete")
