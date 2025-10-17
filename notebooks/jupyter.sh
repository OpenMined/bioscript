#!/bin/bash
uv venv
uv pip install -U jupyter
uv pip install -e ./python
source .venv/bin/activate
jupyter lab
