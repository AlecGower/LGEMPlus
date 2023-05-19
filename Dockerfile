FROM python:3.7-slim

LABEL maintainer.name="Alexander Gower" 
LABEL maintainer.email="gower@chalmers.se" 

# Set up root directory
WORKDIR /usr/src/logmet

# Copy requirements over to working directory
COPY requirements.txt .

RUN apt-get update && apt-get -y install cmake
RUN pip install --no-cache-dir -r requirements.txt
# # Install DFBA dependencies
# RUN pip install --no-cache-dir "pybind11[global]" glpk
# # Install DFBA
# RUN pip install --no-cache-dir dfba

# Copy other files over to working directory
COPY src/ .

CMD ["python", "./tests/load_SBML_models_test.py"]