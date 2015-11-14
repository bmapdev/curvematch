To install curvematch for execution
===================================

Install [shapeio](https://github.com/bmapdev/shapeio) first.

Install [matplotlib](http://matplotlib.org/) by doing 

&nbsp;&nbsp;&nbsp;&nbsp;`conda install matplotlib`

Then install curvematch as: 

&nbsp;&nbsp;&nbsp;&nbsp;`pip install https://github.com/bmapdev/curvematch/archive/master.zip`


To install curvematch for development
=====================================
1. Clone the repository [curvematch](https://github.com/bmapdev/curvematch) as:
    `git clone git@github.com:bmapdev/curvematch.git`
 
     This will create a directory `curvematch` with the source code.
 
2. Change directory to DPmatch by doing `cd curvematch/curvematch/DPmatch`
 
3. Install DPmatch in the developement mode in the DPmatch directory as:

    `python setup.py develop`
    
    This will build the DPmatch C++ source and create a shared library in the `DPmatch` directory.
    

# Start developing curvematch

 