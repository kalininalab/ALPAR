#!/usr/bin/env python3

import os
import pathlib
import shutil
import time
import warnings
import subprocess
import json
import sys

warnings.filterwarnings("ignore")

# Get the path of the script
PATH_OF_SCRIPT = pathlib.Path(__file__).parent.resolve()

def ensure_conda_env(env_name, python_version="3.12"):
    # Normalize env_name
    env_name = env_name.replace("cd-hit", "cdhit")
    
    # Get list of existing environments in JSON format
    try:
        result = subprocess.run(['conda', 'env', 'list', '--json'], capture_output=True, text=True)
        if result.returncode != 0:
            print("Warning: Could not list conda environments. Conda might not be installed or in PATH.")
            return
        envs = json.loads(result.stdout).get('envs', [])
    except Exception as e:
        print(f"Warning: Failed to check conda environments: {e}")
        return
    
    # Conda returns full paths; we check if the env name matches the base name of any path
    exists = any(os.path.basename(env) == env_name for env in envs)
    
    if not exists:
        print(f"Creating environment '{env_name}'...")
        # envs folder is in the same directory as this file (sr_amr/)
        env_file_path = PATH_OF_SCRIPT / "envs" / f"{env_name}.yaml"
        
        if not env_file_path.exists():
            # Try .yml if .yaml doesn't exist
            env_file_path = PATH_OF_SCRIPT / "envs" / f"{env_name}.yml"
            
        if not env_file_path.exists():
             print(f"Error: Conda environment file for {env_name} not found at {env_file_path}.")
             sys.exit(1)
             
        subprocess.run(['conda', 'env' ,'create', '-f', str(env_file_path) ,'-y'])
    else:
        # print(f"Environment '{env_name}' exists.")
        pass

def conda_env_wrapper(env_name):
    """
    Decorator to run a function in a specific conda environment.
    If the current environment is not the target environment, it re-runs
    the entire CLI command using 'conda run -n env_name'.
    """
    def decorator(func):
        def wrapper(args):
            # Check if we are already in the target environment or in a subprocess
            current_env = os.environ.get("CONDA_DEFAULT_ENV")
            if current_env == env_name or os.environ.get("ALPAR_SUBPROCESS") == "1":
                return func(args)
            
            ensure_conda_env(env_name)
            
            # Construct the command to re-run using python -m
            # This ensures it works even if the 'alpar' command is not installed in the target env
            root_dir = PATH_OF_SCRIPT.parent
            cmd = ["conda", "run", "-n", env_name, "--no-capture-output", "python", "-m", "sr_amr.amr"] + sys.argv[1:]
            
            # Set environment variable to avoid infinite loop and ensure project root is in PYTHONPATH
            new_env = os.environ.copy()
            new_env["ALPAR_SUBPROCESS"] = "1"
            if "PYTHONPATH" in new_env:
                new_env["PYTHONPATH"] = f"{root_dir}{os.pathsep}{new_env['PYTHONPATH']}"
            else:
                new_env["PYTHONPATH"] = str(root_dir)
            
            try:
                result = subprocess.run(cmd, env=new_env)
                sys.exit(result.returncode)
            except Exception as e:
                print(f"Error re-running in environment {env_name}: {e}")
                sys.exit(1)
        return wrapper
    return decorator


def is_tool_installed(tool_name):
    # 'which' for Unix-based systems, 'where' for Windows
    command = f"which {tool_name}" if os.name != 'nt' else f"where {tool_name}"
    return os.system(command) == 0


def temp_folder_remover(path):
    for i in range(5):
        try:
            shutil.rmtree(path, ignore_errors=True)
            break
        except OSError:
            time.sleep(1)
    else:
        print(f"Failed to delete {path}")


def time_function(start_time, end_time):

    elapsed_time = end_time - start_time
    hours, rem = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(rem, 60)
    return f"Elapsed time: {int(hours)} hours, {int(minutes)} minutes, {seconds:.2f} seconds"

def copy_and_zip_file(src_file, dst_dir, zip_name):
    # Ensure the destination directory exists
    os.makedirs(dst_dir, exist_ok=True)
    
    # Copy the file to the destination directory
    dst_file = shutil.copy2(src_file, dst_dir)
    
    # Create a zip archive of the destination directory
    shutil.make_archive(zip_name, 'zip', dst_dir)
    
    # Remove the copied file from the destination directory
    os.remove(dst_file)
    
    return zip_name + ".zip"