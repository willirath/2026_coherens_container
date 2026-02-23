#!/usr/bin/env python3
import subprocess
import sys

def main():
    setup = sys.argv[1] if len(sys.argv) > 1 else "setups/seamount_smoke"
    subprocess.run(["bash", "scripts/docker_run_setup.sh", setup], check=True)

if __name__ == "__main__":
    main()
