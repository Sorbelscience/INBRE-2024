import os
import subprocess
import datetime

# Define the directory to back up
backup_dir = "/home/storm/INBRE-2024"

# Define the repository URL
repo_url = "https://github.com/Sorbelscience/INBRE-2024.git"

# Navigate to the backup directory
os.chdir(backup_dir)

# Initialize Git repository if not already initialized
if not os.path.exists(os.path.join(backup_dir, ".git")):
    subprocess.run(["git", "init"], check=True)
    subprocess.run(["git", "remote", "add", "origin", repo_url], check=True)

# Stage all changes
subprocess.run(["git", "add", "."], check=True)

# Commit the changes
commit_message = f"Backup on {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
subprocess.run(["git", "commit", "-m", commit_message], check=True)

# Push the changes to the remote repository
subprocess.run(["git", "push", "-u", "origin", "master"], check=True)