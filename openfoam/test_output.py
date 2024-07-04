#!/usr/bin/env python3

# first run the container as
#
# docker run -it -p 4242:4242 <image name>  
#
# then run this script as python3 test_output.py http://localhost:4242

import argparse
import umbridge
import pytest

parser = argparse.ArgumentParser(description='Minimal HTTP model demo.')
parser.add_argument('url', metavar='url', type=str, help='the ULR on which the model is running, for example http://localhost:4242')
args = parser.parse_args()
print(f"Connecting to host URL {args.url}")

# Set up a model by connecting to URL
model = umbridge.HTTPModel(args.url, "forward")

#test get methods
output = model.get_input_sizes()
print("get_input_sizes() returns "+str(output[0]))
assert pytest.approx(output[0]) == 2, "get_input_sizes() returns wrong value"


output = model.get_output_sizes()
print("get_output_sizes() returns "+str(output[0]))
assert pytest.approx(output[0]) == 1, "get input sizes returns wrong value"

param = [[24, 0]]

#test output for another config
txt = "model output (quantity of interest) = "
# txt = ""
output = model(param,{"Fidelity": 1, "res_tol": 1e-6})
print(txt+str(output[0][0]))
output = model(param,{"Fidelity": 2, "res_tol": 1e-6})
print(txt+str(output[0][0]))
output = model(param,{"Fidelity": 3, "res_tol": 1e-6})
print(txt+str(output[0][0]))
output = model(param,{"Fidelity": 4, "res_tol": 1e-6})
print(txt+str(output[0][0]))
