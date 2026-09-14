"""Command-line orchestration for the complete LPAB workflow."""

import argparse
import shutil
import sys
from pathlib import Path

from .alphafold import collect_results, make_jobs, run_jobs
from .external import grouped_sample, run_mmseqs, run_rosetta
from .generation import generate_candidates


STEPS = ("generate", "rosetta", "cluster", "sample", "af3-input", "af3-run", "collect")


def read_parameters(path):
    values = {}
    for number, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            raise ValueError(f"Invalid parameter line {number}: {raw_line}")
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip()
    return values


def as_bool(value):
    return value.lower() in {"1", "true", "yes", "on"}


def project_path(root, value):
    path = Path(value).expanduser()
    return path if path.is_absolute() else root / path


def require(values, name):
    if name not in values or not values[name]:
        raise ValueError(f"Missing parameter: {name}")
    return values[name]


def check_installation(root, values):
    checks = []
    try:
        import lppack  # noqa: F401
        checks.append(("lppack", True, "Python package found"))
    except ImportError:
        checks.append(("lppack", False, "install from the official GitHub repository"))

    executable_checks = [
        ("Rosetta", "rosetta_executable"),
        ("BLASTP", "blastp_executable"),
        ("MMseqs2", "mmseqs_executable"),
    ]
    if as_bool(values.get("use_clustalo", "true")):
        executable_checks.insert(0, ("Clustal Omega", "clustalo_executable"))
    for label, key in executable_checks:
        executable = values.get(key, "")
        checks.append((label, bool(shutil.which(executable)), executable or "not configured"))

    for label, key in (
        ("AlphaFold script", "alphafold_script"),
        ("AlphaFold models", "alphafold_model_dir"),
        ("AlphaFold databases", "alphafold_database_dir"),
    ):
        path = project_path(root, values.get(key, "missing"))
        checks.append((label, path.exists(), str(path)))

    for label, passed, detail in checks:
        print(f"[{'OK' if passed else 'MISSING'}] {label}: {detail}")
    return all(passed for _, passed, _ in checks)


def run_step(step, root, values):
    output = project_path(root, require(values, "output_dir"))
    output.mkdir(parents=True, exist_ok=True)
    seed_csv = project_path(root, require(values, "seed_csv"))
    antigen_fasta = project_path(root, require(values, "antigen_fasta"))
    candidates = output / "01_candidates.csv"
    rosetta_passed = output / "02_rosetta_passed.csv"
    clustered = output / "03_clustered.csv"
    selected = output / "04_selected.csv"
    jobs_dir = output / "05_af3_jobs"
    af3_output = output / "06_af3_output"
    final_csv = output / "07_final_scores.csv"

    print(f"\n=== {step} ===", flush=True)
    if step == "generate":
        clustalo = values["clustalo_executable"] if as_bool(values.get("use_clustalo", "true")) else None
        count = generate_candidates(
            seed_csv,
            candidates,
            output / "01_alignment",
            int(values["candidate_count"]),
            int(values["random_seed"]),
            int(values["max_fragment_length"]),
            clustalo,
        )
        print(f"Generated {count} candidates: {candidates}")
    elif step == "rosetta":
        count = run_rosetta(
            candidates,
            rosetta_passed,
            output / "02_rosetta",
            values["rosetta_executable"],
            values["blastp_executable"],
            int(values["rosetta_timeout_seconds"]),
            int(values["rosetta_workers"]),
        )
        print(f"Rosetta passed {count} candidates: {rosetta_passed}")
    elif step == "cluster":
        count = run_mmseqs(
            rosetta_passed,
            clustered,
            output / "03_mmseqs",
            values["mmseqs_executable"],
            float(values["mmseqs_min_identity"]),
            float(values["mmseqs_coverage"]),
            int(values["mmseqs_cov_mode"]),
        )
        print(f"MMseqs2 found {count} clusters: {clustered}")
    elif step == "sample":
        count = grouped_sample(
            clustered, selected, int(values["sample_count"]), int(values["random_seed"])
        )
        print(f"Selected {count} candidates: {selected}")
    elif step == "af3-input":
        count = make_jobs(selected, antigen_fasta, jobs_dir, int(values["random_seed"]))
        print(f"Created {count} AlphaFold 3 jobs: {jobs_dir}")
    elif step == "af3-run":
        count = run_jobs(
            jobs_dir,
            af3_output,
            values["alphafold_python"],
            str(project_path(root, values["alphafold_script"])),
            str(project_path(root, values["alphafold_model_dir"])),
            str(project_path(root, values["alphafold_database_dir"])),
            int(values["alphafold_recycles"]),
            int(values["alphafold_diffusion_samples"]),
        )
        print(f"Completed {count} AlphaFold 3 jobs: {af3_output}")
    elif step == "collect":
        count = collect_results(
            jobs_dir,
            af3_output,
            final_csv,
            root / "lpab_tools" / "ipsae.py",
            values.get("ipsae_python", sys.executable),
            float(values["ipsae_pae_cutoff"]),
            float(values["ipsae_distance_cutoff"]),
        )
        print(f"Collected {count} results: {final_csv}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parameters", default="parameters.txt")
    parser.add_argument("--step", choices=(*STEPS, "all"), default="all")
    parser.add_argument("--check", action="store_true", help="Check external software and paths.")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    parameter_path = project_path(root, args.parameters)
    values = read_parameters(parameter_path)
    if args.check:
        raise SystemExit(0 if check_installation(root, values) else 1)
    for step in STEPS if args.step == "all" else (args.step,):
        run_step(step, root, values)


if __name__ == "__main__":
    main()
