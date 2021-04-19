# AutoBend
Automated functions for vertebral bending in Maya Python

Python files for conducting automated bending analysis in Maya

AutoBend.py - files for setting up models and conducting a single set of bending experiments (one joint)

wrapper_vert_bend.py - Mayapy file for running AutoBend on multiple joints, models, and with various parameters
Models should be set up and arranged within a single folder, then the bending wrapper can be run directly from command line using Mayapy.exe
This allows multiple models with sensitivity repeats to be run in the background (or on a cluster).
