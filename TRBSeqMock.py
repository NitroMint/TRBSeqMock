import argparse
import os
import subprocess
import sys
from typing import Dict, Callable

MODULE_CONFIG: Dict[str, str] = {
    "TRB": "TRBSeqMock.py",
    "IGH": "IGHSeqMock.py"
    # Add new modules here: "ABC": "ABCSeqMock.py"
}


def run_module(module_type: str, number: int, output_dir: str,
                shm_params: dict = None) -> None:
    """
    Unified function to run specified module

    Args:
        module_type: Type of module to run (TRB/IGH/other configured modules)
        number: Number of sequences to generate
        output_dir: Path to output directory
        shm_params: Optional dict of SHM parameter overrides

    Raises:
        FileNotFoundError: If the module script does not exist
        subprocess.CalledProcessError: If the script execution fails
    """
    if shm_params is None:
        shm_params = {}

    # Get script name and construct full path
    script_name = MODULE_CONFIG[module_type]
    script_dir = os.path.join(os.getcwd(), f"Source_{module_type}")
    script_path = os.path.join(script_dir, script_name)

    # Check if script exists
    if not os.path.exists(script_path):
        raise FileNotFoundError(
            f"{module_type} module script not found: {script_path}\n"
            f"Please confirm {script_name} exists in Source_{module_type} directory"
        )

    project_root = os.getcwd()
    source_module_dir = os.path.join(project_root, f"Source_{module_type}")

    cmd = [
        "python3",
        "-c",
        f'import sys; sys.path.insert(0, "{project_root}"); '
        f'sys.path.insert(0, "{source_module_dir}"); '
        f'exec(open("{script_path}").read())',
        "-n", str(number),
        "-o", output_dir,
    ]

    # Append SHM arguments as CLI args (so parse_args() in write.py can read them)
    if shm_params:
        if shm_params.get("shm_rate") is not None:
            cmd.extend(["--shm-rate", str(shm_params["shm_rate"])])
        if shm_params.get("shm_hotspot_rr") is not None:
            cmd.extend(["--shm-hotspot-rr", str(shm_params["shm_hotspot_rr"])])
        if shm_params.get("shm_cdr_bias") is not None:
            cmd.extend(["--shm-cdr-bias", str(shm_params["shm_cdr_bias"])])
        if shm_params.get("shm_cold_spot_rate") is not None:
            cmd.extend(["--shm-cold-spot-rate", str(shm_params["shm_cold_spot_rate"])])
        if shm_params.get("shm_seed") is not None:
            cmd.extend(["--shm-seed", str(shm_params["shm_seed"])])

    result = subprocess.run(
        cmd,
        check=True,
        text=True,
        capture_output=True,
        cwd=project_root,
    )
    if result.stdout:
        print(result.stdout)


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments

    Returns:
        argparse.Namespace: Parsed arguments object
    """
    parser = argparse.ArgumentParser(description="Immunorepertoire sequence simulation tool")
    parser.add_argument("-n", "--number", type=int, required=True,
                        help="Number of output sequences (positive integer)")
    parser.add_argument("-t", "--type", type=str, required=True,
                        choices=MODULE_CONFIG.keys(),
                        help=f"Module type, available options: {', '.join(MODULE_CONFIG.keys())}")
    parser.add_argument("-o", "--output", type=str, default="output",
                        help="Output directory path (default: output)")

    # ---- SHM parameters (applied to IGH module only) ----
    parser.add_argument("--shm-rate", type=float, default=0.05,
                        help="Global SHM rate per base (default: 0.05 = 5%%) [IGH only]")
    parser.add_argument("--shm-hotspot-rr", type=float, default=4.0,
                        help="Hotspot relative risk multiplier (default: 4.0) [IGH only]")
    parser.add_argument("--shm-cdr-bias", type=float, default=1.5,
                        help="CDR mutation bias vs FR (default: 1.5) [IGH only]")
    parser.add_argument("--shm-cold-spot-rate", type=float, default=0.002,
                        help="Cold spot residual rate (default: 0.002) [IGH only]")
    parser.add_argument("--shm-seed", type=int, default=None,
                        help="Random seed for SHM reproducibility [IGH only]")

    return parser.parse_args()


def ensure_output_dir(output_dir: str) -> None:
    """
    Ensure output directory exists, create if not

    Args:
        output_dir: Path to output directory
    """
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory prepared: {os.path.abspath(output_dir)}")


def main():
    """Main execution function"""
    # Parse command line arguments
    args = parse_args()

    # Validate number is positive integer
    if args.number <= 0:
        print("Error: -n/--number must be a positive integer")
        exit(1)

    # Ensure output directory exists
    ensure_output_dir(args.output)

    # Build SHM params dict (only passed to IGH module)
    shm_params = {}
    if args.type == "IGH":
        shm_params["shm_rate"] = args.shm_rate
        shm_params["shm_hotspot_rr"] = args.shm_hotspot_rr
        shm_params["shm_cdr_bias"] = args.shm_cdr_bias
        shm_params["shm_cold_spot_rate"] = args.shm_cold_spot_rate
        shm_params["shm_seed"] = args.shm_seed

    try:
        # Run specified module
        print(f"Starting {args.type} module, generating {args.number} sequences...")
        if args.type == "IGH":
            print(f"  SHM: rate={args.shm_rate}, hotspot_rr={args.shm_hotspot_rr}, "
                  f"cdr_bias={args.shm_cdr_bias}")
        run_module(args.type, args.number, args.output, shm_params)
        print(f"\n{args.type} module executed successfully!")
        print(f"Output files saved to: {os.path.abspath(args.output)}")

    except FileNotFoundError as e:
        print(f"\nFile error: {str(e)}")
        exit(1)

    except subprocess.CalledProcessError as e:
        print(f"\n{args.type} module execution failed:")
        print(f"Command return code: {e.returncode}")
        print(f"Standard error output: {e.stderr}")
        exit(1)

    except Exception as e:
        print(f"\nUnknown error occurred while running {args.type} module: {str(e)}")
        print("\nFor assistance, please contact:")
        print("- Email: 2210240103@csu.edu.cn")
        print("- GitHub: https://github.com/NitroMint/TRBSeqMock")
        exit(1)


if __name__ == "__main__":
    main()