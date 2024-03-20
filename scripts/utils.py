#!/usr/bin/env python3

import os
import pathlib
import shutil
import time

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
