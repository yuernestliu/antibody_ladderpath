"""AlphaFold 3 input, execution, ipSAE, and summary helpers."""

import csv
import json
import random
import re
import subprocess
from pathlib import Path


OUTPUT_COLUMNS = [
    "name", "heavy", "light", "iptm", "ipsae_min", "plddt", "ptm", "pae", "ipae"
]


def read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def read_fasta(path):
    records = []
    current = None
    for raw_line in Path(path).read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current is not None:
                records.append("".join(current))
            current = []
        elif current is None:
            raise ValueError("The antigen FASTA must start with a header line.")
        else:
            current.append(line)
    if current is not None:
        records.append("".join(current))
    if len(records) != 1 or not records[0]:
        raise ValueError("The antigen FASTA must contain exactly one non-empty sequence.")
    sequence = records[0].upper()
    if not sequence.isalpha():
        raise ValueError("The antigen FASTA sequence must contain letters only.")
    return sequence


def safe_name(value):
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")
    if not name:
        raise ValueError("A candidate has an empty or invalid name.")
    return name


def make_jobs(input_csv, antigen_fasta, jobs_dir, seed):
    rows = read_csv(input_csv)
    antigen = read_fasta(antigen_fasta)
    jobs_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(seed)
    names = set()
    for row in rows:
        name = safe_name(row["name"])
        if name in names:
            raise ValueError(f"Duplicate AlphaFold job name: {name}")
        names.add(name)
        job = {
            "name": name,
            "modelSeeds": [rng.randint(1, 2_147_483_647)],
            "sequences": [
                {"protein": {"id": "H", "sequence": row["heavy"]}},
                {"protein": {"id": "L", "sequence": row["light"]}},
                {"protein": {"id": "A", "sequence": antigen}},
            ],
            "dialect": "alphafold3",
            "version": 1,
        }
        (jobs_dir / f"{name}.json").write_text(
            json.dumps(job, indent=2) + "\n", encoding="utf-8"
        )
    return len(rows)


def run_jobs(jobs_dir, output_dir, python, script, model_dir, database_dir, recycles, samples):
    output_dir.mkdir(parents=True, exist_ok=True)
    jobs = sorted(jobs_dir.glob("*.json"))
    for index, job in enumerate(jobs, 1):
        print(f"AlphaFold 3 job {index}/{len(jobs)}: {job.stem}", flush=True)
        subprocess.run(
            [
                python,
                script,
                f"--json_path={job}",
                f"--model_dir={model_dir}",
                f"--db_dir={database_dir}",
                f"--output_dir={output_dir}",
                f"--num_recycles={recycles}",
                f"--num_diffusion_samples={samples}",
            ],
            check=True,
        )
    return len(jobs)


def mean(values):
    values = list(values)
    return sum(values) / len(values) if values else None


def find_result_files(output_dir, name):
    candidates = [path for path in output_dir.glob(f"{name}*") if path.is_dir()]
    if not candidates:
        return None, None, None
    folder = sorted(candidates, key=lambda path: path.stat().st_mtime, reverse=True)[0]
    summaries = sorted(folder.glob("*summary_confidences.json"))
    full = sorted(
        path for path in folder.glob("*confidences.json")
        if "summary_confidences" not in path.name
    )
    models = sorted(folder.glob("*_model.cif"))
    return (
        summaries[0] if summaries else None,
        full[0] if full else None,
        models[0] if models else None,
    )


def parse_ipsae_max(path):
    """Read the official ipSAE text table and return H-A/L-A max scores."""
    pair_scores = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) < 6 or fields[4] != "max":
            continue
        try:
            pair_scores[frozenset((fields[0], fields[1]))] = float(fields[5])
        except ValueError:
            continue
    expected = [frozenset(("H", "A")), frozenset(("L", "A"))]
    if not all(pair in pair_scores for pair in expected):
        return None
    return min(pair_scores[pair] for pair in expected)


def run_ipsae(script, confidence_json, model_cif, python, pae_cutoff, distance_cutoff):
    subprocess.run(
        [python, str(script), str(confidence_json), str(model_cif), str(pae_cutoff), str(distance_cutoff)],
        check=True,
    )
    pae_label = f"{int(pae_cutoff):02d}"
    distance_label = f"{int(distance_cutoff):02d}"
    output = model_cif.with_suffix("").with_name(
        f"{model_cif.stem}_{pae_label}_{distance_label}.txt"
    )
    if not output.exists():
        raise FileNotFoundError(f"ipSAE output not found: {output}")
    return parse_ipsae_max(output)


def round_value(value):
    return "" if value is None else round(float(value), 6)


def collect_results(jobs_dir, output_dir, final_csv, ipsae_script, ipsae_python, pae_cutoff, distance_cutoff):
    rows = []
    incomplete = []
    job_paths = sorted(jobs_dir.glob("*.json"))
    if not job_paths:
        raise FileNotFoundError(f"No AlphaFold 3 job JSON files found in {jobs_dir}")
    for job_path in job_paths:
        job = json.loads(job_path.read_text(encoding="utf-8"))
        name = job["name"]
        summary_path, full_path, model_path = find_result_files(output_dir, name)
        if not all((summary_path, full_path, model_path)):
            incomplete.append(name)
            continue
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        full = json.loads(full_path.read_text(encoding="utf-8"))
        pae_matrix = full.get("pae", [])
        chains = full.get("token_chain_ids", [])
        pae = mean(value for row in pae_matrix for value in row)
        ipae = mean(
            pae_matrix[i][j]
            for i in range(len(chains))
            for j in range(len(chains))
            if chains[i] != chains[j]
        ) if pae_matrix and len(chains) == len(pae_matrix) else None
        ipsae_min = run_ipsae(
            ipsae_script, full_path, model_path, ipsae_python, pae_cutoff, distance_cutoff
        )
        if ipsae_min is None:
            raise RuntimeError(f"ipSAE did not produce both H-A and L-A scores for {name}")
        sequences = {
            item["protein"]["id"]: item["protein"]["sequence"]
            for item in job["sequences"] if "protein" in item
        }
        rows.append({
            "name": name,
            "heavy": sequences.get("H", ""),
            "light": sequences.get("L", ""),
            "iptm": round_value(summary.get("iptm")),
            "ipsae_min": round_value(ipsae_min),
            "plddt": round_value(mean(full.get("atom_plddts", []))),
            "ptm": round_value(summary.get("ptm")),
            "pae": round_value(pae),
            "ipae": round_value(ipae),
        })
    if incomplete:
        names = ", ".join(incomplete)
        raise RuntimeError(
            "Refusing to write a partial final CSV; incomplete AlphaFold 3 results: "
            f"{names}"
        )
    rows.sort(key=lambda row: float(row["iptm"] or -1), reverse=True)
    with open(final_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)
