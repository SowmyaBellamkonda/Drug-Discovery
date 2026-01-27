import os
import uuid
import subprocess
from celery_app import celery

@celery.task
def run_docking_task(protein, ligand):
    receptor_dir = "receptor"
    results_dir = "results"
    docked_dir = "docked_outputs"

    receptor = os.path.join(receptor_dir, f"{protein}.pdbqt")
    ligand_path = os.path.join(results_dir, ligand)

    os.makedirs(docked_dir, exist_ok=True)
    pose_file = os.path.join(docked_dir, f"{uuid.uuid4()}_pose.pdbqt")

    # ⬇️ Correct Vina executable path
    VINA = r"C:\vina\vina.exe"

    vina_cmd = [
        VINA,
        "--receptor", receptor,
        "--ligand", ligand_path,
        "--center_x", "0", "--center_y", "0", "--center_z", "0",
        "--size_x", "20", "--size_y", "20", "--size_z", "20",
        "--exhaustiveness", "3",
        "--num_modes", "1",
        "--cpu", "4",
        "--out", pose_file
    ]

    result = subprocess.run(vina_cmd, capture_output=True, text=True)
    output = result.stdout

    affinity = None
    for line in output.splitlines():
        if line.strip().startswith("1 "):
            affinity = float(line.split()[1])
            break

    if affinity is None:
        return {"error": "Docking failed", "output": output}

    summary = (
        "Excellent inhibitor" if affinity < -9 else
        "Good binding" if affinity < -7 else
        "Moderate binding" if affinity < -5 else
        "Weak binding"
    )

    return {
        "affinity": affinity,
        "summary": summary,
        "pose_file": os.path.basename(pose_file)
    }
