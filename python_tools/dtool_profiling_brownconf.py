#!/usr/bin/env python3
"""
A script to profile the code of an interacting brownian particles project.

This script copies the source code from a specified directory, compiles it with
profiling flags, runs the executable to generate a `gmon.out` file, and then
uses gprof to create a profiling report.
"""

import subprocess
from pathlib import Path
import sys
from typing import Dict, Any

# Get the directory of the script to make paths relative to it
SCRIPT_DIR = Path(__file__).parent.resolve()

# --- Configuration ---
CONFIG: Dict[str, Any] = {
    "executable_name": "main_brownconf",
    "confinement": "sept",
    "interaction": "lennardjones",
    "ext_force": "1",
    "setnumb": "5",
    "numb_mpi_tasks": "0",
    "source_dir": SCRIPT_DIR.parent / "src",
    "build_dir": SCRIPT_DIR.parent / "profiling",
}


def copy_source_files(source_dir: Path, target_dir: Path) -> None:
    """Copies all .c and .h files from the source to the target directory.

    Args:
        source_dir: The directory containing the source files.
        target_dir: The directory to copy the files to.
    """
    print(f"Copying source files from {source_dir} to {target_dir}...")
    for extension in ["*.c", "*.h"]:
        for file in source_dir.glob(extension):
            subprocess.run(["cp", str(file), str(target_dir)], check=True)


def compile_code(
    executable_name: str, confinement: str, interaction: str, build_dir: Path
) -> None:
    """Compiles the C source code with profiling flags.

    Args:
        executable_name: The name of the output executable.
        confinement: The type of confinement to use.
        interaction: The type of inter-particle interaction to use.
        build_dir: The directory containing the source files for compilation.
    """
    print("Compiling source code...")
    compile_cmd = [
        "gcc",
        "-w",
        "-Wall",
        "-Wextra",
        "-pedantic",
        "main_brownconf.c",
        "sim_config.c",
        "results_transport.c",
        "code_handling.c",
        "random_numb_gen.c",
        "simulation_core.c",
        "print_routines.c",
        "array_utils.c",
        f"conf_{confinement}.c",
        f"int_{interaction}.c",
        "-lgsl",
        "-lgslcblas",
        "-lm",
        "-o",
        executable_name,
        "-p",
        "-pg",
    ]
    subprocess.run(compile_cmd, check=True, cwd=build_dir)


def run_simulation(
    executable_name: str,
    ext_force: str,
    setnumb: str,
    numb_mpi_tasks: str,
    build_dir: Path,
) -> None:
    """Runs the compiled simulation to generate profiling data.

    Args:
        executable_name: The name of the executable to run.
        ext_force: The value of the external force.
        setnumb: The size of the interacting samples.
        numb_mpi_tasks: The number of parallel threads.
        build_dir: The directory containing the executable.
    """
    print("Running simulation...")
    execute_cmd = [
        f"./{executable_name}",
        ext_force,
        setnumb,
        numb_mpi_tasks,
    ]
    subprocess.run(execute_cmd, check=True, cwd=build_dir)


def generate_profile_report(
    executable_name: str,
    confinement: str,
    interaction: str,
    ext_force: str,
    setnumb: str,
    build_dir: Path,
) -> None:
    """Generates a gprof profiling report.

    Args:
        executable_name: The name of the profiled executable.
        confinement: The confinement type used in the simulation.
        interaction: The interaction type used in the simulation.
        ext_force: The external force value used.
        setnumb: The sample size value used.
        build_dir: The directory containing the executable and gmon.out.
    """
    print("Generating profiling report...")
    gmon_path = build_dir / "gmon.out"
    if not gmon_path.is_file():
        print(f"Error: Profiling output file '{gmon_path}' not found.", file=sys.stderr)
        sys.exit(1)

    output_filename = (
        f"gprof_result_conf_{confinement}_int_{interaction}_"
        f"force_{ext_force}_setnumb_{setnumb}.txt"
    )
    output_path = build_dir / output_filename

    with open(output_path, "w") as f:
        profile_cmd = ["gprof", f"./{executable_name}", "gmon.out"]
        subprocess.run(profile_cmd, check=True, cwd=build_dir, stdout=f)

    print(f"Profiling report saved to {output_path}")


def main() -> None:
    """Main function to run the profiling script."""
    try:
        copy_source_files(CONFIG["source_dir"], CONFIG["build_dir"])
        compile_code(
            CONFIG["executable_name"],
            CONFIG["confinement"],
            CONFIG["interaction"],
            CONFIG["build_dir"],
        )
        run_simulation(
            CONFIG["executable_name"],
            CONFIG["ext_force"],
            CONFIG["setnumb"],
            CONFIG["numb_mpi_tasks"],
            CONFIG["build_dir"],
        )
        generate_profile_report(
            CONFIG["executable_name"],
            CONFIG["confinement"],
            CONFIG["interaction"],
            CONFIG["ext_force"],
            CONFIG["setnumb"],
            CONFIG["build_dir"],
        )
    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError as e:
        print(f"File not found: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
