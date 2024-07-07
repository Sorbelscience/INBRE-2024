# check_env.py
import sys
import os

print("Python executable: ", sys.executable)
print("Python version: ", sys.version)
print("Current working directory: ", os.getcwd())
print("Environment variables: ", os.environ)
