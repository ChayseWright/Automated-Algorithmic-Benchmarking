"""Run the available independent checks; does not execute MATLAB."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

root = Path(__file__).resolve().parents[1]
names = ["check_structure.py", "check_numerics.py", "check_study_numerics.py",
         "check_revision3_numerics.py", "check_motion_numerics.py"]
results = []
for name in names:
    run = subprocess.run([sys.executable, str(root/"validation"/name)], cwd=root,
                         capture_output=True, text=True)
    record = {"check": name, "passed": run.returncode == 0,
              "kind": "lexical_structure" if name == "check_structure.py" else "independent_equations_and_source_inspection"}
    results.append(record)
    print(("PASS " if record["passed"] else "FAIL ")+name, flush=True)
    if run.returncode:
        print(run.stdout, run.stderr)
        raise SystemExit(run.returncode)
sources = list(root.rglob("*.m"))
report = {"softwareVersion": "4.0.0", "date": "2026-11-12", "checks": results,
          "matlabSourceFiles": len(sources), "MATLABNativeTestsExecuted": False,
          "MATLABNativeStatus": "runtime_unavailable_not_executed", "PatientDataAnalyzed": False,
          "motionAddendumSHA256": hashlib.sha256((root/"docs/tremor_motion_capture_benchmark_addendum.pdf").read_bytes()).hexdigest()}
(root/"validation/release_validation_summary.json").write_text(json.dumps(report, indent=2)+"\n")
print("All available checks passed. Native MATLAB tests remain unexecuted.")
