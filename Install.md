# Installation Tips
Installing the morphometrics pipeline is usually relatively easy, but there are a few edge cases that require special attention.

## M1 Macs
**Thanks to Genevieve Buckley for this solution: (https://github.com/GrotjahnLab/surface_morphometrics/issues/8)**
This is a temporary workaround and needs additional testing!

0. Build the conda environment like normal, but don't pip install pymeshlab - just remove it from the pip_requirements.txt file
1. Install pyqt into your conda environment:
For M1 Macs, the andfoy conda channel has a working pyqt build:
```bash
conda install pyqt -c andfoy
```

2. Build pymeshlab
https://github.com/cnr-isti-vclab/PyMeshLab/blob/main/src/README.md

```bash
brew install libomp cgal xerces-c
git clone --recursive https://github.com/cnr-isti-vclab/PyMeshLab.git
cd PyMeshLab

# Build pymeshlab
mkdir src/build
cd src/build
cmake ..
make
make install
```
<!-- 
This is not what worked - we need the original build instructions
modify build script (scripts/macOS/1_build.sh) line 77, with custom qt dir and remove -DBUILD_DUMMY_BIN_MAC_DEPLOY flag.`
sh scripts/macOS/1_build.sh
sh scripts/macOS/2_deploy.sh
 -->

3. Export the `KMP_DUPLICATE_LIB_OK` environment variable
This prevents libomp errors.
You can use conda to permanently set an environment variable in the morphometrics conda environment:
```bash=
conda activate morphometrics
conda env config vars KMP_DUPLICATE_LIB_OK=TRUE

# Must reactivate conda environment
conda deactivate
conda activate morphometrics
conda env config vars list
# You should now see KMP_DUPLICATE_LIB_OK = TRUE in this list
```

4. Check the build worked:
```bash
python -m pip install pytest
python -m pytest --pyargs pymeshlab
```

4. Install pymeshlab to the conda environment:
```bash
pip install .
```

5. You're ready to go!
Try running the surface-morphometrics scripts on the example data. See the README for more details.

## Container Implementation
For older linux environments (such as CentOS7) Qt5 will not behave well. A container implementation is in progress (https://github.com/GrotjahnLab/surface_morphometrics/issues/10)