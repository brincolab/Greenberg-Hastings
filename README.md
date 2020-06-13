# Greenberg-Hastings
Greenberg-Hasting model implementation in python

Connectome data both in matlab format and numpy .txt files are in `Data/` folder. 
Both Haimovici and Rocha dynamics are implemented in this code.

```
usage: Green_hastings.py [-h] --data-location DATA_LOCATION [--total TOTAL]
                         [--r1 R1] [--r2 R2] [--temp-min TEMP_MIN]
                         [--temp-max TEMP_MAX] [--temp-step TEMP_STEP]
                         [--dynamics DYNAMICS] [--output-dir OUTPUT_DIR]
                         [--output-prefix OUTPUT_PREFIX]
                         [--matlab-matrix-name MATLAB_MATRIX_NAME]

optional arguments:
  -h, --help            show this help message and exit
  --data-location DATA_LOCATION
                        .txt or .mat file containing the Conectome Data (NxN)
                        Matrix
  --total TOTAL         [Default: 1000] Length of the simulation
  --r1 R1               [Default: 0.001] Spontaneous activation probability
                        (it set a very small rate of background activity)
  --r2 R2               [Default: 0.2] Refractory to Quiescent probability (it
                        determine the lengh of time in which is active)
  --temp-min TEMP_MIN   [Default: 0.01] Minimum value for temperature range
                        (activation threshold)
  --temp-max TEMP_MAX   [Default: 0.1] Maximum value for temperature range
                        (activation threshold)
  --temp-step TEMP_STEP
                        [Default: 0.01] Increment for temperature range
                        (activation threshold)
  --dynamics DYNAMICS   [Default: Haimovici] Set the dynamics between
                        'Haimovici' and 'Maritan'
  --output-dir OUTPUT_DIR
                        [Default: './'] Output directory to export the
                        simulations
  --output-prefix OUTPUT_PREFIX
                        [Default: 'Sim'] prefix for the exported files
  --matlab-matrix-name MATLAB_MATRIX_NAME
                        If --data-location is a .mat file, this argument is to
                        parse the matlab cariable to python format
```

Example of running Green-hastings.py  with matlab .mat file

```
python Green_hastings.py --data-location Data/DSI_enhanced.mat --matlab-matrix-name CIJ_fbden_average
```

Usage of this code as a module in python, is shown on Example_as_module.ipynb

This code is based of Jeremi K. Ochab original matlab code, https://github.com/remolek/