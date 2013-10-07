For testing only:

----------------------------------------------------------------------
Create a Virtual Environment using Enthought/Canopy python
----------------------------------------------------------------------
Note: You want to create a virtual environment for the specific python
interpreter for Enthought (EPD) or Canopy
For e.g. if you are in ~/sandbox, and the path to your EPD python is
/usr/local/epd/bin/python, then create the virtual environment as

virtualenv -p /usr/local/epd/bin/python ~/sandbox/myenv --system-site-packages

----------------------------------------------------------------------
Install shapeio from https://shjoshi@bitbucket.org/shjoshi/shapeio.git as follows:
----------------------------------------------------------------------
1. For e.g. if you are in ~/sandbox, type 
   git clone https://shjoshi@bitbucket.org/shjoshi/shapeio.git
   This will make a new directory ~/sandbox/shapeio and clone the repository
2. ~/sandbox/myenv/bin/pip install ~/sandbox/shapeio
   This will install the shapeio module in your virtual environment
   
----------------------------------------------------------------------
Install curvematch from https://shjoshi@bitbucket.org/shjoshi/curvematch.git as follows:
----------------------------------------------------------------------
1. For e.g. if you are in ~/sandbox, type 
   git clone https://shjoshi@bitbucket.org/shjoshi/curvematch.git
   This will make a new directory ~/sandbox/curvematch and clone the repository

----------------------------------------------------------------------
Install nose in the virtualenv for the testing environment as follows:
----------------------------------------------------------------------
~/sandbox/myenv/bin/pip install nose -I

Note the "-I" above.

----------------------------------------------------------------------
Run the nosetests for curvematch
----------------------------------------------------------------------
1. Go to the curvematch directory
cd ~/sandbox/curvematch

2. Run nosetests
~/sandbox/myenv/bin/nosetests -v -s

----------------------------------------------------------------------
Develop curvematch and/or shapeio
----------------------------------------------------------------------
Start developing curvematch and shapeio