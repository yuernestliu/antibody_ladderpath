"""Ladderpath-based antibody sequence generation."""

import csv
import random
import subprocess
from collections import defaultdict
from pathlib import Path


AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWYXBZJUO")


def clean_sequence(value):
    sequence = "".join(value.upper().split()).replace("-", "")
    unexpected = set(sequence) - AMINO_ACIDS
    if not sequence or unexpected:
        raise ValueError(f"Invalid amino-acid sequence: {sorted(unexpected)}")
    return sequence


def read_seed_table(path):
    """Read common column spellings and return a normalized seed table."""
    with open(path, newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError("The seed CSV is empty.")

    aliases = {
        "name": ("name", "Name"),
        "heavy": ("heavy", "H", "VH"),
        "light": ("light", "L", "VL"),
        "score": ("score", "Score"),
    }
    columns = {}
    for standard, choices in aliases.items():
        columns[standard] = next((name for name in choices if name in rows[0]), None)
        if columns[standard] is None:
            raise ValueError(f"Missing {standard} column. Accepted names: {choices}")

    normalized = []
    for row in rows:
        normalized.append({
            "name": row[columns["name"]].strip(),
            "heavy": clean_sequence(row[columns["heavy"]]),
            "light": clean_sequence(row[columns["light"]]),
            "score": max(float(row[columns["score"]]), 0.001),
        })
    return normalized


def write_fasta(records, path):
    with open(path, "w", encoding="utf-8") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n{sequence}\n")


def read_fasta(path):
    records = []
    name = None
    parts = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(parts)))
            name, parts = line[1:].split()[0], []
        elif line.strip():
            parts.append(line.strip())
    if name is not None:
        records.append((name, "".join(parts)))
    return records


def align_sequences(names, sequences, executable, work_dir, label):
    input_path = work_dir / f"{label}_input.fasta"
    output_path = work_dir / f"{label}_aligned.fasta"
    write_fasta(zip(names, sequences), input_path)
    subprocess.run(
        [executable, "-i", str(input_path), "-o", str(output_path), "--force", "--outfmt=fasta"],
        check=True,
    )
    aligned = dict(read_fasta(output_path))
    return [aligned[name] for name in names]


def pad_sequences(sequences):
    """Validate sequences are already aligned when Clustal Omega is disabled."""
    lengths = {len(sequence) for sequence in sequences}
    if len(lengths) != 1:
        raise ValueError(
            "Clustal Omega is disabled, so all seed sequences must have equal length."
        )
    return list(sequences)


def find_positions(sequence, fragment):
    result = []
    start = 0
    while True:
        start = sequence.find(fragment, start)
        if start < 0:
            return result
        result.append(start)
        start += 1


def prepare_pool(sequences, scores, max_fragment_length):
    try:
        import lppack as lp
    except ImportError as error:
        raise RuntimeError(
            "lppack is not installed. Run: python -m pip install "
            "git+https://github.com/yuernestliu/lppack.git"
        ) from error

    lpjson = lp.get_ladderpath(
        sequences,
        fill_ladderons_STR=False,
        save_file_name=None,
        show_version=False,
        mode="cpu",
    )
    pom, _ = lp.POM_from_JSON(lpjson, display_str=False)
    multiplicities = {}
    for level in pom.values():
        for fragment, count in level.items():
            multiplicities[fragment] = max(count, multiplicities.get(fragment, 0))

    positions = defaultdict(list)
    fragments, weights = [], []
    for fragment, multiplicity in multiplicities.items():
        biological_length = len(fragment.replace("-", ""))
        if biological_length == 0 or biological_length > max_fragment_length:
            continue
        observed_scores = []
        for sequence, score in zip(sequences, scores):
            found = find_positions(sequence, fragment)
            positions[fragment].extend(found)
            observed_scores.extend([score] * len(found))
        if observed_scores:
            fragments.append(fragment)
            weights.append(
                (sum(observed_scores) / len(observed_scores))
                * max(multiplicity, 1)
                * biological_length
            )
    if not fragments:
        raise RuntimeError("Ladderpath did not return usable fragments.")
    return sequences, scores, fragments, weights, positions


def generate_chain(pool, rng):
    sequences, scores, fragments, weights, positions = pool
    length = len(sequences[0])
    result = [None] * length

    for _ in range(length * 20):
        fragment = rng.choices(fragments, weights=weights, k=1)[0]
        valid = [start for start in positions[fragment] if start + len(fragment) <= length]
        if not valid:
            continue
        start = rng.choice(valid)
        if all(result[start + i] in (None, letter) for i, letter in enumerate(fragment)):
            result[start : start + len(fragment)] = fragment
        if all(value is not None for value in result):
            break

    for position, value in enumerate(result):
        if value is None:
            residues = [sequence[position] for sequence in sequences]
            result[position] = rng.choices(residues, weights=scores, k=1)[0]
    return "".join(result).replace("-", "")


def generate_candidates(seed_csv, output_csv, work_dir, count, seed, max_fragment_length, clustalo=None):
    rows = read_seed_table(seed_csv)
    names = [row["name"] for row in rows]
    scores = [row["score"] for row in rows]
    heavy = [row["heavy"] for row in rows]
    light = [row["light"] for row in rows]
    work_dir.mkdir(parents=True, exist_ok=True)

    if clustalo:
        heavy = align_sequences(names, heavy, clustalo, work_dir, "heavy")
        light = align_sequences(names, light, clustalo, work_dir, "light")
    else:
        heavy, light = pad_sequences(heavy), pad_sequences(light)

    heavy_pool = prepare_pool(heavy, scores, max_fragment_length)
    light_pool = prepare_pool(light, scores, max_fragment_length)
    rng = random.Random(seed)
    original = {(row["heavy"], row["light"]) for row in rows}
    generated, seen = [], set(original)
    for _ in range(count * 100):
        pair = (generate_chain(heavy_pool, rng), generate_chain(light_pool, rng))
        if pair not in seen:
            seen.add(pair)
            generated.append(pair)
        if len(generated) == count:
            break
    if len(generated) < count:
        raise RuntimeError(
            f"Generated only {len(generated)} unique pairs. Add more diverse seeds or reduce candidate_count."
        )

    with open(output_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "heavy", "light"])
        writer.writeheader()
        for index, (heavy_sequence, light_sequence) in enumerate(generated, 1):
            writer.writerow({
                "name": f"candidate_{index:05d}",
                "heavy": heavy_sequence,
                "light": light_sequence,
            })
    return len(generated)
