import subprocess
import logging

# Set up logging
logging.basicConfig(filename='/scratch/SCRATCH_SAS/alper/SR-AMR/example/example_output/output.txt', level=logging.INFO)

# Define the command to run
command = "prokka --listdb"

process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

# Capture output line by line
for line in iter(process.stdout.readline, b''):
    line = line.decode('utf-8').strip()  # Decode bytes to string
    logging.info(line)  # Write the line to the log file

# Wait for the process to finish
process.communicate()


print("Output redirected to", "/scratch/SCRATCH_SAS/alper/SR-AMR/example/example_output/temp/output.txt")