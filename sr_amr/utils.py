#!/usr/bin/env python3

import os
import pathlib
import shutil
import time
import warnings

warnings.filterwarnings("ignore")

# Get the path of the script
PATH_OF_SCRIPT = pathlib.Path(__file__).parent.resolve()


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