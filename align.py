import os
import re
import subprocess
from typing import Any, Dict


def check_alignment_tool(tool_name: str) -> str:
    """Check if the specified alignment tool exists and is executable."""
    tool_path = os.path.join(os.path.dirname(__file__), "bin", 
                              "TMalign" if tool_name == "TMAlign" else "USalign")
    if not os.path.exists(tool_path) or not os.access(tool_path, os.X_OK):
        raise FileNotFoundError(f"{tool_name} executable not found at {tool_path}")
    return tool_path


def run_alignment(
    pdb_f1: str, pdb_f2: str, alignment_tool: str, use_tmscore: bool, use_rmsd: bool
) -> Dict[str, Any]:
    """Run structural alignment and extract features."""
    features = {}
    cmd = [
        check_alignment_tool(alignment_tool),
        pdb_f1,
        pdb_f2,
    ]
    try:
        output = subprocess.check_output(cmd, text=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        print(f"Alignment failed: {str(e)}")
        return features
    except Exception as e:
        print(f"Unexpected error: {str(e)}")
        return features

    if use_tmscore:
        if m := re.search(r"TM-score=\s*([0-9.]+)", output):
            features["tma_score"] = float(m.group(1))
    if use_rmsd:
        if m := re.search(r"RMSD=\s*([0-9.]+)", output):
            features["rmsd"] = float(m.group(1))

    return features
