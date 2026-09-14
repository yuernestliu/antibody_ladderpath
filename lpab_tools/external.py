"""Rosetta validation, MMseqs2 clustering, and grouped sampling."""

import csv
import math
import random
import re
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(rows, path, columns):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def rosetta_log_is_complete(text):
    """Require complete chain, framework, and CDR parsing by Rosetta."""
    patterns = [
        r"antibody: Heavy chain sequence:\s*[A-Z]+",
        r"antibody: Light chain sequence:\s*[A-Z]+",
        r"H1 detected:\s*[A-Z]+", r"H2 detected:\s*[A-Z]+",
        r"H3 detected:\s*[A-Z]+", r"L1 detected:\s*[A-Z]+",
        r"L2 detected:\s*[A-Z]+", r"L3 detected:\s*[A-Z]+",
    ]
    if not all(re.search(pattern, text) for pattern in patterns):
        return False
    heavy_framework = re.search(r"HeavyAntibodyFramework.*", text)
    light_framework = re.search(r"LightAntibodyFramework.*", text)
    if not heavy_framework or not light_framework:
        return False
    framework_pattern = r"fr\d\[[^\]]+,\s*[A-Z]+\]"
    return (
        len(re.findall(framework_pattern, heavy_framework.group())) >= 4
        and len(re.findall(framework_pattern, light_framework.group())) >= 4
    )


def _run_rosetta_job(row, executable, blastp, fasta_dir, log_dir, model_dir, timeout):
    name = row["name"]
    fasta = fasta_dir / f"{name}.fasta"
    fasta.write_text(
        f">heavy\n{row['heavy']}\n>light\n{row['light']}\n",
        encoding="utf-8",
    )
    prefix = model_dir / name
    prefix.mkdir(parents=True, exist_ok=True)
    command = [executable, "-fasta", str(fasta), "-blastp", blastp, "-prefix", str(prefix)]
    try:
        completed = subprocess.run(
            command,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=timeout,
            check=False,
        )
        output = completed.stdout or ""
    except subprocess.TimeoutExpired as error:
        output = error.stdout or ""
        if isinstance(output, bytes):
            output = output.decode(errors="replace")
    (log_dir / f"{name}.log").write_text(output, encoding="utf-8")
    return row if rosetta_log_is_complete(output) else None


def run_rosetta(input_csv, output_csv, work_dir, executable, blastp, timeout, workers):
    rows = read_csv(input_csv)
    fasta_dir, log_dir, model_dir = (
        work_dir / "fasta", work_dir / "logs", work_dir / "models"
    )
    for folder in (fasta_dir, log_dir, model_dir):
        folder.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [
            pool.submit(
                _run_rosetta_job, row, executable, blastp, fasta_dir, log_dir, model_dir, timeout
            )
            for row in rows
        ]
    passed = [future.result() for future in futures]
    passed = [row for row in passed if row is not None]
    write_csv(passed, output_csv, ["name", "heavy", "light"])
    if not passed:
        raise RuntimeError("No candidate passed the Rosetta antibody log checks.")
    return len(passed)


def run_mmseqs(input_csv, output_csv, work_dir, executable, min_identity, coverage, cov_mode):
    rows = read_csv(input_csv)
    work_dir.mkdir(parents=True, exist_ok=True)
    fasta = work_dir / "heavy_light.fasta"
    with open(fasta, "w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(f">{row['name']}\n{row['heavy']}{row['light']}\n")

    database = work_dir / "sequence_db"
    clusters = work_dir / "clusters"
    temporary = work_dir / "tmp"
    mapping = work_dir / "cluster_map.tsv"
    subprocess.run([executable, "createdb", str(fasta), str(database)], check=True)
    subprocess.run(
        [
            executable, "cluster", str(database), str(clusters), str(temporary),
            "--min-seq-id", str(min_identity), "-c", str(coverage),
            "--cov-mode", str(cov_mode),
        ],
        check=True,
    )
    subprocess.run(
        [executable, "createtsv", str(database), str(database), str(clusters), str(mapping)],
        check=True,
    )
    member_to_cluster = {}
    with open(mapping, encoding="utf-8") as handle:
        for line in handle:
            representative, member = line.rstrip("\n").split("\t")[:2]
            member_to_cluster[member] = representative
    for row in rows:
        row["cluster"] = member_to_cluster[row["name"]]
    write_csv(rows, output_csv, ["name", "heavy", "light", "cluster"])
    return len(set(member_to_cluster.values()))


def grouped_sample(input_csv, output_csv, count, seed):
    rows = read_csv(input_csv)
    if count > len(rows):
        raise ValueError("sample_count is larger than the clustered table.")
    groups = defaultdict(list)
    for row in rows:
        groups[row["cluster"]].append(row)
    weights = {name: math.sqrt(len(members)) for name, members in groups.items()}
    rng = random.Random(seed)
    selected = []
    while len(selected) < count:
        active = [name for name, members in groups.items() if members]
        cluster = rng.choices(active, weights=[weights[name] for name in active], k=1)[0]
        selected.append(groups[cluster].pop(rng.randrange(len(groups[cluster]))))
    write_csv(selected, output_csv, ["name", "heavy", "light", "cluster"])
    return len(selected)
