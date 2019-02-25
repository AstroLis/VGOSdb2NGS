# VGOSdb2NGS
converter from VGOSdb to NGS card files

Version 0.2

Usage:

```
python read.py --archives-dir <input VGOSdb root directory>
python read.py --archive <input VGOSdb tar.gz file>
```

Example for Windows:
```
python.exe read.py --archive 18JAN25XE.tar.gz
```

# Dependences

Uses `numpy` and `scipy.io.netcdf` for reading netcdf files.
