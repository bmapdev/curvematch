This is for testing purposes only
=================================

The process below outlines the process for installing shapeio and curvematch.
This will be replaced by another documentation for those who just want to run the programs.


Quick Commands to create virtual environment, installs, and tests
----------------------------------------------------------
Assume one is in `~/sandbox`, and the path to your EPD python is
`/usr/local/epd/bin/python`


```
virtualenv -p /usr/local/epd/bin/python ~/sandbox/myenv --system-site-packages

#Clone shapeio
git clone https://shjoshi@bitbucket.org/shjoshi/shapeio.git 

#Install shapeio
~/sandbox/myenv/bin/pip install ~/sandbox/shapeio 

#Clone curvematch
git clone https://shjoshi@bitbucket.org/shjoshi/curvematch.git 

#Install Nose
~/sandbox/myenv/bin/pip install nose -I

cd ~/sandbox/curvematch

#Run nosetests
~/sandbox/myenv/bin/nosetests -v -s 
```


###For a detailed explanation of the above commands read on below.###


Create a Virtual Environment using Enthought/Canopy python
----------------------------------------------------------

>**Note:** You want to create a virtual environment for the specific python
>interpreter for Enthought (EPD) or Canopy

For e.g. if you are in `~/sandbox`, and the path to your EPD python is
`/usr/local/epd/bin/python`, then create the virtual environment as

```virtualenv -p /usr/local/epd/bin/python ~/sandbox/myenv --system-site-packages```


Clone and Install shapeio from bitbucket
----------------------------------------
For e.g. if you are in `~/sandbox`, type

```git clone https://shjoshi@bitbucket.org/shjoshi/shapeio.git```
   
This will make a new directory `~/sandbox/shapeio` and clone the repository.

Then type
   
```~/sandbox/myenv/bin/pip install ~/sandbox/shapeio```
   
This will install the shapeio module in your virtual environment


Clone curvematch from bitbucket
-------------------------------
For e.g. if you are in `~/sandbox`, type

```git clone https://shjoshi@bitbucket.org/shjoshi/curvematch.git```
  
This will make a new directory ~/sandbox/curvematch and clone the repository


Install nose in the virtualenv for the testing environment
----------------------------------------------------------

```~/sandbox/myenv/bin/pip install nose -I```

>**Note** the `-I` above.


Run the nosetests for curvematch
--------------------------------
1. Go to the curvematch directory
```cd ~/sandbox/curvematch```
2. Run nosetests
```~/sandbox/myenv/bin/nosetests -v -s```


Develop curvematch and/or shapeio
---------------------------------
Start developing curvematch and shapeio