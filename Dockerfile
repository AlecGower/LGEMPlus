FROM debian AS iproverBuilder

RUN apt-get update && apt-get install -y -q --no-install-recommends git

RUN git clone -b 2021_Nov_Bio https://gitlab.com/korovin/iprover.git

WORKDIR /iprover

RUN apt-get install -y -q --no-install-recommends opam bubblewrap libgmp-dev
# Can disable sandboxing as we are running within a container here
RUN opam init --disable-sandboxing && opam update -y \
&& opam switch create 4.14.1+flambda --package=ocaml-variants.4.14.1+options,ocaml-option-flambda \
&& opam switch 4.14.1+flambda \
&& opam install ocamlfind ocamlgraph zarith yojson z3 -y \
&& eval $(opam env) \
&& ./configure && make \
RUN ./iproveropt --help

FROM python:3.7-slim

LABEL maintainer.name="Alexander Gower" 
LABEL maintainer.email="gower@chalmers.se" 

# Set up root directory
WORKDIR /usr/src/logmet

# Copy over statically compiled iProver image

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
COPY --from=iproverBuilder ./iprover /usr/src/

CMD ["python", "./tests/load_SBML_models_test.py"]
# CMD /usr/src/iprover/iproveropt --help