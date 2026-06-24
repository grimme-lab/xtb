#!/usr/bin/env python3
"""End-to-end numerical gate for CUDA-backed GFN1/GFN2 SCF and --opt."""

from __future__ import annotations

import argparse
import math
import os
import re
import subprocess
import tempfile
from pathlib import Path


ENERGY_RE = re.compile(r"TOTAL ENERGY\s+([-+0-9.Ee]+)\s+Eh")
NUMBER_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+)(?:[Ee][-+]?\d+)?")


def run(binary: Path, fixture: Path, method: int, args: list[str], work: Path) -> str:
    cmd = [str(binary), str(fixture), "--gfn", str(method), *args]
    proc = subprocess.run(
        cmd,
        cwd=work,
        env={**os.environ, "OMP_NUM_THREADS": "1"},
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    text = proc.stdout.replace(b"\0", b"").decode(errors="replace")
    if proc.returncode != 0 or "abnormal termination" in text:
        raise RuntimeError(f"{' '.join(cmd)} failed ({proc.returncode})\n{text[-4000:]}")
    return text


def energy(output: str) -> float:
    match = ENERGY_RE.search(output)
    if not match:
        raise RuntimeError("TOTAL ENERGY not found")
    return float(match.group(1))


def numeric_file(path: Path) -> list[float]:
    text = path.read_bytes().replace(b"\0", b"").decode(errors="replace")
    return [float(value) for value in NUMBER_RE.findall(text)]


def coordinates(path: Path) -> list[float]:
    lines = path.read_bytes().replace(b"\0", b"").decode(errors="replace").splitlines()
    values: list[float] = []
    for line in lines[2:]:
        numbers = NUMBER_RE.findall(line)
        if len(numbers) >= 3:
            values.extend(float(value) for value in numbers[-3:])
    return values


def max_delta(left: list[float], right: list[float]) -> float:
    if len(left) != len(right) or not left:
        raise RuntimeError(f"incompatible vectors: {len(left)} vs {len(right)}")
    return max(abs(a - b) for a, b in zip(left, right))


def gate_method(cpu: Path, gpu: Path, fixture: Path, method: int, root: Path) -> None:
    cpu_sp = root / f"cpu-sp-{method}"
    gpu_sp = root / f"gpu-sp-{method}"
    cpu_sp.mkdir()
    gpu_sp.mkdir()
    cpu_out = run(cpu, fixture, method, ["--sp"], cpu_sp)
    gpu_out = run(gpu, fixture, method, ["--gpu", "--sp"], gpu_sp)
    if "GPU cuSolver SCC enabled" not in gpu_out:
        raise RuntimeError(f"GFN{method}: CUDA SCC marker missing")
    de = abs(energy(cpu_out) - energy(gpu_out))
    if de > 1.0e-8:
        raise RuntimeError(f"GFN{method}: single-point |dE|={de:.3e} Eh")

    cpu_grad = root / f"cpu-grad-{method}"
    gpu_grad = root / f"gpu-grad-{method}"
    cpu_grad.mkdir()
    gpu_grad.mkdir()
    run(cpu, fixture, method, ["--grad"], cpu_grad)
    run(gpu, fixture, method, ["--gpu", "--grad"], gpu_grad)
    dg = max_delta(numeric_file(cpu_grad / "gradient"),
                   numeric_file(gpu_grad / "gradient"))
    if dg > 1.0e-6:
        raise RuntimeError(f"GFN{method}: gradient max |dG|={dg:.3e} Eh/a0")

    cpu_opt = root / f"cpu-opt-{method}"
    gpu_opt = root / f"gpu-opt-{method}"
    cpu_opt.mkdir()
    gpu_opt.mkdir()
    cpu_opt_out = run(cpu, fixture, method, ["--opt", "loose"], cpu_opt)
    gpu_opt_out = run(gpu, fixture, method, ["--gpu", "--opt", "loose"], gpu_opt)
    if "GEOMETRY OPTIMIZATION CONVERGED" not in cpu_opt_out:
        raise RuntimeError(f"GFN{method}: CPU optimization did not converge")
    if "GEOMETRY OPTIMIZATION CONVERGED" not in gpu_opt_out:
        raise RuntimeError(f"GFN{method}: GPU optimization did not converge")
    de_opt = abs(energy(cpu_opt_out) - energy(gpu_opt_out))
    dx = max_delta(coordinates(cpu_opt / "xtbopt.xyz"),
                   coordinates(gpu_opt / "xtbopt.xyz"))
    if de_opt > 1.0e-8:
        raise RuntimeError(f"GFN{method}: optimized |dE|={de_opt:.3e} Eh")
    if dx > 1.0e-7:
        raise RuntimeError(f"GFN{method}: optimized geometry max |dx|={dx:.3e} A")

    batch = root / f"gpu-batch-{method}"
    batch.mkdir()
    batch_out = run(
        gpu,
        fixture,
        method,
        ["--gpu", "--gpu-batch", str(fixture)],
        batch,
    )
    if "SCF diagonalization + density routed through CUDA" not in batch_out:
        raise RuntimeError(f"GFN{method}: CUDA batch marker missing")

    print(
        f"GFN{method}: PASS  "
        f"|dE_sp|={de:.3e} Eh  max|dG|={dg:.3e} Eh/a0  "
        f"|dE_opt|={de_opt:.3e} Eh  max|dx|={dx:.3e} A"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cpu", type=Path, required=True)
    parser.add_argument("--gpu", type=Path, required=True)
    parser.add_argument("--fixture", type=Path, required=True)
    args = parser.parse_args()

    cpu = args.cpu.resolve()
    gpu = args.gpu.resolve()
    fixture = args.fixture.resolve()
    with tempfile.TemporaryDirectory(prefix="xtb-gfn12-gpu-gate-") as tmp:
        root = Path(tmp)
        for method in (1, 2):
            gate_method(cpu, gpu, fixture, method, root)


if __name__ == "__main__":
    main()
