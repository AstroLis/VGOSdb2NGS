# VGOSdb2NGS
converter from VGOSdb to NGS card files

Version 0.1

Usage:

```
python read.py --archives-dir <input VGOSdb root directory>
```

Example for Windows:
```
python.exe read.py "d:/forms/2018/vgosDB/18MAR02XU" "18MAR02XU.ngs"
```

# Dependences

Uses `numpy` and `scipy.io.netcdf` for reading netcdf files.
