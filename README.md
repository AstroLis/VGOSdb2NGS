# VGOSdb2NGS
converter from VGOSdb to NGS card files

Version 0.2

Usage:

```
python read.py <input VGOSdb tar.gz file> <output NGS card file name>
```

Example for Windows:
```
python.exe read.py 18JAN25XE.tar.gz 18JAN25XE.ngs
```

# Dependences

Uses `numpy` and `scipy.io.netcdf` for reading netcdf files.
