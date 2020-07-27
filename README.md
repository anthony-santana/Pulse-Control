# Pulse Control 

Pulse Control is a repository for developing optimized quantum controls out of Discrete Prolate Spheroidal Sequences,
through the use of the Deep Reinforcement Learning algorithm, Proximal Policy Optimization.

## Using this repository

This repository is built as an extension of Oak Ridge National Laboratory's open-source quantum framwork, XACC. 
It leverages XACC's pulse-level control functionality 

## Quick Start with Docker

Pulse Control is currently set up for use as a Docker Image. From a scratch directory in your terminal, run

```
docker run --security-opt seccomp=unconfined --init -it -p 3000:3000 xacc/xacc-drl
```
This will serve an Eclipse Theia IDE on port 3000 with XACC and all of its dependent packages already configured. Users can navigate
to the examples folder and become comfortable with the syntax and formating of files, or they can run their own custom
experiments, ensuring to prepend each file with:


```python

import sys
sys.path.insert(1, '/home/cades/dev/Pulse_Control/gym_pulsecontrol/')
import xacc_drl
import numpy as np

```
