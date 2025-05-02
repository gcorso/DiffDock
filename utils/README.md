The key changes I made to torus.py

1. Added `SCRIPT_DIR` to get the absolute path of the directory where torus.py is located
2. Created absolute paths for the cache files using `os.path.join(SCRIPT_DIR, '.p.npy')`
3. Added print statements for debugging to show when and where files are being created
4. Improved the condition to check for both cache files before loading them

This should solve the issue because:

1. Now it will look for and save the cache files in the same directory as the script itself, regardless of where the script is being run from
2. The added debugging prints will help you see if/when the cache files are being created
3. The script will explicitly check that both required cache files exist before trying to load them
