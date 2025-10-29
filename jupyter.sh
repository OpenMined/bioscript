#!/bin/bash
uv venv --clear
uv pip install -U jupyter pandas cyvcf2 pytest
uv pip install -e ./python
source .venv/bin/activate
jupyter lab
