Bootstrap: docker
# Use miniconda build, wait to see if this causes issues with pymeshlab...
From:  continuumio/miniconda3 

%post
git clone https://github.com/grotjahnlab/surface_morphometrics.git
cd surface_morphometrics
conda env create -f environment.yml
source activate morphometrics
pip install -r pip_requirements.txt

