import subprocess

# List of required packages
required_packages = [
    'os',
    'sys',
    'vtk',
    'vtk.numpy_interface',
    'numpy',
    'scipy',
    'h5py',
    'multiprocessing'
]

# Check and install packages
for package in required_packages:
    try:
        __import__(package)
        print(f"{package} is already installed.")
    except ImportError:
        #print(f"{package} is not installed. Installing...")
        subprocess.run(['pip', 'install', package])

print("All required packages are installed.")
