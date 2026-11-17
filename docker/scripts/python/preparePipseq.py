#!/usr/bin/env python3

import argparse
import gzip
import shlex
import shutil
import subprocess
import sys
import tempfile
from itertools import zip_longest
from pathlib import Path

import yaml

GEX_ADAPTER = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
BASES = set("ACGT")
LIBRARY_TYPES = ("GEX", "CaTCHseq")
STAGED_CONFIG = "config.yaml"
PIPSPEAK_DIR = "pipspeak"


def fail(message):
    sys.stderr.write(f"[preparePipseq] ERROR: {message}\n")
    sys.exit(1)


def log(message):
    print(f"[preparePipseq] {message}", flush=True)


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def check_gzip(path):
    with open(path, "rb") as handle:
        magic = handle.read(2)
    if magic != b"\x1f\x8b":
        fail(f"{path} is not gzip compressed")


def read_records(path):
    with open_text(path) as handle:
        number = 0
        while True:
            header = handle.readline()
            if not header:
                return
            number += 1
            if not header.startswith("@"):
                fail(f"{path}:{number} is not a FASTQ header")
            sequence = handle.readline()
            separator = handle.readline()
            quality = handle.readline()
            if not sequence or not separator or not quality:
                fail(f"{path}:{number} is a truncated FASTQ record")
            sequence = sequence.rstrip("\n")
            quality = quality.rstrip("\n")
            if not separator.startswith("+"):
                fail(f"{path}:{number + 2} is not a FASTQ separator")
            if len(sequence) != len(quality):
                fail(
                    f"{path}:{number} has {len(sequence)} bases but "
                    f"{len(quality)} quality values"
                )
            yield header[1:].strip(), sequence
            number += 3


def record_name(header):
    token = header.split()[0] if header.split() else ""
    if token.endswith("/1") or token.endswith("/2"):
        return token[:-2]
    return token


def validate_pair(r1_path, r2_path):
    count = 0
    pairs = zip_longest(read_records(r1_path), read_records(r2_path))
    for index, (left, right) in enumerate(pairs, start=1):
        if left is None or right is None:
            fail(
                f"lane {r1_path.name}/{r2_path.name} is not paired: "
                f"record {index} has no counterpart"
            )
        if record_name(left[0]) != record_name(right[0]):
            fail(
                f"lane {r1_path.name}/{r2_path.name} is not paired: record {index} "
                f"is '{left[0]}' in R1 and '{right[0]}' in R2"
            )
        count += 1
    return count


def validate_inputs(r1_files, r2_files):
    if not r1_files:
        fail("no input files given")
    if len(r1_files) != len(r2_files):
        fail(
            f"got {len(r1_files)} R1 and {len(r2_files)} R2 files, "
            "every lane needs a paired R1/R2 file"
        )
    for r1, r2 in zip(r1_files, r2_files):
        for path in (r1, r2):
            if not path.is_file():
                fail(f"input file {path} does not exist")
            check_gzip(path)
        count = validate_pair(r1, r2)
        if count == 0:
            fail(f"lane {r1.name}/{r2.name} contains no reads")
        log(f"lane {r1.name}/{r2.name}: {count} read pairs")


def check_outputs(outputs, inputs):
    resolved_inputs = {path.resolve() for path in inputs}
    for path in outputs:
        if path.is_symlink():
            fail(f"output {path} is a symlink, refusing to write through it")
        if path.exists():
            fail(f"output {path} already exists, refusing to overwrite it")
        if path.resolve() in resolved_inputs:
            fail(f"output {path} aliases the input file {path.resolve()}")


def check_r1_length(path, expected):
    for index, (_, sequence) in enumerate(read_records(path), start=1):
        if len(sequence) != expected:
            fail(
                f"{path}:{index} has a {len(sequence)} nt R1, expected {expected} "
                "(barcode plus UMI)"
            )


def resolve_reference(yaml_path, relative, what):
    if Path(relative).is_absolute():
        fail(f"{what} '{relative}' must be relative to the config file {yaml_path}")
    candidate = (yaml_path.parent / relative).resolve()
    if not candidate.is_file():
        fail(f"{what} '{relative}' does not resolve to a file ({candidate})")
    return candidate


def read_barcode_list(path):
    sequences = []
    with path.open() as handle:
        for number, line in enumerate(handle, start=1):
            sequence = line.strip()
            if not sequence:
                fail(f"{path}:{number} is blank, barcode lists must not have empty lines")
            if sequence != sequence.upper():
                fail(f"{path}:{number} barcode '{sequence}' is not upper case")
            if not set(sequence) <= BASES:
                fail(f"{path}:{number} barcode '{sequence}' is not an ACGT sequence")
            sequences.append(sequence)
    if not sequences:
        fail(f"{path} contains no barcodes")
    return sequences


def load_config(config_path, cb_len, umi_len, barcode_files=None):
    if not config_path.is_file():
        fail(f"config YAML {config_path} does not exist")
    with config_path.open() as handle:
        raw = yaml.safe_load(handle)
    if not isinstance(raw, dict):
        fail(f"{config_path} does not contain a YAML mapping")
    barcodes = raw.get("barcodes")
    spacers = raw.get("spacers") or {}
    if not isinstance(barcodes, dict) or not barcodes:
        fail(f"{config_path} has no 'barcodes' mapping")
    if not isinstance(spacers, dict):
        fail(f"{config_path} 'spacers' must be a mapping")
    if len(spacers) != len(barcodes) - 1:
        fail(
            f"{config_path} defines {len(spacers)} spacers for {len(barcodes)} "
            f"barcode sets, expected {len(barcodes) - 1}"
        )
    overrides = None
    if barcode_files is not None:
        overrides = {}
        for name, path in barcode_files:
            if name in overrides:
                fail(f"barcode list '{name}' is given more than once")
            overrides[name] = Path(path)
        if set(overrides) != set(barcodes):
            fail(
                f"--barcode-file does not match the barcode sets of {config_path}: "
                f"missing {sorted(set(barcodes) - set(overrides))}, "
                f"unexpected {sorted(set(overrides) - set(barcodes))}"
            )
    entries = []
    lengths = []
    for name, relative in barcodes.items():
        if not isinstance(name, str) or not name.strip():
            fail(f"{config_path} has a barcode set without a name")
        if not isinstance(relative, str) or not relative.strip():
            fail(f"{config_path} barcode list '{name}' has no path")
        if overrides is not None:
            path = overrides[name]
            if not path.is_file():
                fail(f"barcode list '{name}' file {path} does not exist")
        else:
            path = resolve_reference(
                config_path, relative, f"barcode list '{name}'"
            )
        sequences = read_barcode_list(path)
        set_lengths = {len(sequence) for sequence in sequences}
        if len(set_lengths) != 1:
            fail(
                f"barcode list '{name}' ({path}) has non uniform lengths "
                f"{sorted(set_lengths)}"
            )
        lengths.append(set_lengths.pop())
        entries.append((name, path))
    total = sum(lengths)
    if total != cb_len:
        fail(
            f"barcode lengths {lengths} sum to {total}, "
            f"expected {cb_len} (--cb-len)"
        )
    spacer_entries = []
    for name, sequence in spacers.items():
        text = str(sequence).strip()
        if not text or text != text.upper() or not set(text) <= BASES:
            fail(f"{config_path} spacer '{name}' is not an upper case ACGT sequence")
        spacer_entries.append((name, text))
    parameters = raw.get("parameters") or {}
    if not isinstance(parameters, dict):
        fail(f"{config_path} 'parameters' must be a mapping")
    yaml_umi = parameters.get("umi_len")
    if yaml_umi is not None and int(yaml_umi) != umi_len:
        fail(
            f"{config_path} sets parameters.umi_len={yaml_umi}, which overrides "
            f"--umi-len {umi_len}; use the same value in both places"
        )
    log(
        f"config {config_path}: barcode lengths {lengths} (sum {total}), "
        f"spacers {[sequence for _, sequence in spacer_entries]}, umi_len {umi_len}"
    )
    return entries, spacer_entries


def stage_config(entries, spacer_entries, staged_dir, umi_len):
    staged_dir.mkdir(parents=True, exist_ok=True)
    barcodes = {}
    for index, (name, path) in enumerate(entries, start=1):
        safe = "".join(
            character if character.isalnum() or character in "-_." else "_"
            for character in name
        )
        target = staged_dir / f"{index:02d}_{safe}.txt"
        shutil.copyfile(path, target)
        barcodes[name] = str(target.resolve())
    spacers = {name: sequence for name, sequence in spacer_entries}
    payload = {"barcodes": barcodes, "spacers": spacers}
    if umi_len:
        payload["parameters"] = {"umi_len": umi_len}
    yaml_path = staged_dir / STAGED_CONFIG
    with yaml_path.open("w") as handle:
        yaml.safe_dump(payload, handle, sort_keys=False)
    log(f"staged config {yaml_path}")
    return yaml_path


def merge_lanes(inputs, out_path):
    with open(out_path, "wb") as out:
        for path in inputs:
            with open(path, "rb") as handle:
                shutil.copyfileobj(handle, out, 1 << 20)
    log(f"merged {len(inputs)} lane(s) into {out_path}")


def run_command(command, what, log_path=None):
    log(f"running {what}: {shlex.join(str(part) for part in command)}")
    if log_path is None:
        result = subprocess.run(command, check=False)
    else:
        with log_path.open("w") as handle:
            result = subprocess.run(command, check=False, stdout=handle)
    if result.returncode != 0:
        fail(f"{what} failed with exit code {result.returncode}")


def read_log_stats(log_path):
    if not log_path.is_file():
        return {}
    with log_path.open() as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        return {}
    stats = data.get("statistics") or {}
    return {
        "total_reads": stats.get("total_reads"),
        "passing_reads": stats.get("passing_reads"),
        "fraction_passing": stats.get("fraction_passing"),
        "whitelist_size": stats.get("whitelist_size"),
    }


def write_cb_whitelist(r1_path, cb_len, out_path):
    seen = set()
    for _, sequence in read_records(r1_path):
        seen.add(sequence[:cb_len])
    with out_path.open("w") as handle:
        for barcode in sorted(seen):
            handle.write(f"{barcode}\n")
    log(f"wrote {len(seen)} CB-only barcodes to {out_path}")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Run pipspeak on staged PIPseq lanes and prepare mapper input"
    )
    parser.add_argument("--r1", nargs="+", required=True, help="lane R1 FASTQ files in order")
    parser.add_argument("--r2", nargs="+", required=True, help="lane R2 FASTQ files in order")
    parser.add_argument("--prefix", required=True, help="output prefix")
    parser.add_argument("--config", required=True, help="barcode config YAML")
    parser.add_argument(
        "--barcode-file",
        action="append",
        nargs=2,
        default=None,
        metavar=("NAME", "PATH"),
        help="barcode list NAME at PATH, overriding the path in --config",
    )
    parser.add_argument("--cb-len", type=int, default=24, help="concatenated barcode length")
    parser.add_argument("--umi-len", type=int, default=12, help="UMI length")
    parser.add_argument("--offset", type=int, default=8, help="barcode offset from the R1 start")
    parser.add_argument("--umi-offset", type=int, default=10, help="UMI offset from the barcode end")
    parser.add_argument("--threads", type=int, default=1, help="threads for pipspeak and cutadapt")
    parser.add_argument("--exact", action="store_true", help="exact barcode matching")
    parser.add_argument(
        "--library-type",
        default="GEX",
        choices=LIBRARY_TYPES,
        help="GEX libraries are adapter trimmed, CaTCHseq libraries are not",
    )
    parser.add_argument("--min-length", type=int, default=20, help="cutadapt minimum length")
    parser.add_argument("--pipspeak-bin", default="pipspeak", help="pipspeak executable")
    parser.add_argument("--cutadapt-bin", default="cutadapt", help="cutadapt executable")
    parser.add_argument("--cb-whitelist-out", default=None, help="optional CB-only whitelist output")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.cb_len <= 0:
        fail("--cb-len must be positive")
    if args.umi_len < 1 or args.umi_len > 16:
        fail("--umi-len must be between 1 and 16")
    if args.threads <= 0:
        fail("--threads must be positive")
    if args.offset < 0 or args.umi_offset < 0:
        fail("--offset and --umi-offset must not be negative")
    if args.min_length < 0:
        fail("--min-length must not be negative")

    r1_files = [Path(path) for path in args.r1]
    r2_files = [Path(path) for path in args.r2]
    validate_inputs(r1_files, r2_files)

    entries, spacer_entries = load_config(
        Path(args.config), args.cb_len, args.umi_len, args.barcode_file
    )

    prefix = Path(args.prefix)
    out_r1 = Path(f"{prefix}_R1.fastq.gz")
    out_r2 = Path(f"{prefix}_R2.fastq.gz")
    outputs = [out_r1, out_r2]
    if args.cb_whitelist_out:
        outputs.append(Path(args.cb_whitelist_out))
    check_outputs(outputs, r1_files + r2_files)

    pipspeak_dir = Path(PIPSPEAK_DIR)
    if pipspeak_dir.exists():
        fail(f"output directory {pipspeak_dir} already exists, refusing to reuse it")
    pipspeak_dir.mkdir()

    pipspeak_prefix = pipspeak_dir / prefix.name
    converted_r1 = Path(f"{pipspeak_prefix}_R1.fq.gz")
    converted_r2 = Path(f"{pipspeak_prefix}_R2.fq.gz")

    with tempfile.TemporaryDirectory(dir=".", prefix="pipseq_work_") as work:
        work_dir = Path(work)
        staged_yaml = stage_config(
            entries, spacer_entries, work_dir / "config", args.umi_len
        )
        merged_r1 = work_dir / "merged_R1.fastq.gz"
        merged_r2 = work_dir / "merged_R2.fastq.gz"
        merge_lanes(r1_files, merged_r1)
        merge_lanes(r2_files, merged_r2)
        command = [
            args.pipspeak_bin,
            "--r1",
            str(merged_r1),
            "--r2",
            str(merged_r2),
            "--prefix",
            str(pipspeak_prefix),
            "--config",
            str(staged_yaml),
            "--offset",
            str(args.offset),
            "--umi-offset",
            str(args.umi_offset),
            "--umi-len",
            str(args.umi_len),
            "--threads",
            str(args.threads),
        ]
        if args.exact:
            command.append("--exact")
        run_command(command, "pipspeak")

    for path in (converted_r1, converted_r2):
        if not path.is_file():
            fail(f"pipspeak did not produce {path}")
    n1 = validate_pair(converted_r1, converted_r2)
    if n1 == 0:
        fail(f"pipspeak produced no passing read pairs from {len(r1_files)} lane(s)")
    check_r1_length(converted_r1, args.cb_len + args.umi_len)
    stats = read_log_stats(Path(f"{pipspeak_prefix}_log.yaml"))
    log(
        f"pipspeak kept {n1} read pairs "
        f"(total {stats.get('total_reads')}, passing {stats.get('passing_reads')}, "
        f"fraction {stats.get('fraction_passing')})"
    )

    if args.library_type == "GEX":
        cutadapt_log = Path(f"{pipspeak_prefix}_cutadapt.log")
        command = [
            args.cutadapt_bin,
            "--cores",
            str(args.threads),
            "-A",
            GEX_ADAPTER,
            "--minimum-length",
            str(args.min_length),
            "-o",
            str(out_r1),
            "-p",
            str(out_r2),
            str(converted_r1),
            str(converted_r2),
        ]
        run_command(command, "cutadapt", cutadapt_log)
        for path in (out_r1, out_r2):
            if not path.is_file() or path.stat().st_size == 0:
                fail(f"cutadapt did not produce a non empty {path}")
        log(f"cutadapt report written to {cutadapt_log}")
    else:
        shutil.move(str(converted_r1), str(out_r1))
        shutil.move(str(converted_r2), str(out_r2))
        log(f"kept {n1} untrimmed read pairs for {args.library_type}")

    t1 = validate_pair(out_r1, out_r2)
    check_r1_length(out_r1, args.cb_len + args.umi_len)
    log(f"wrote {t1} read pairs to {out_r1} and {out_r2}")

    if args.cb_whitelist_out:
        write_cb_whitelist(out_r1, args.cb_len, Path(args.cb_whitelist_out))

    return 0


if __name__ == "__main__":
    sys.exit(main())
