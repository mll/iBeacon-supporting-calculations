# iBeacon-supporting-calculations

This repo contains useful tools for work with iBeacon. There are two versions of kalman filter for predicting current distance 
to the beacon from several measurements (kalman.jl). distance.jl contains functions that compute current location based on those 
filtered signals coming from multiple beacons. It uses DBSCAN-Rtree clusterizer for speed and has computational complexity of O((N logN)^2).

Just download julia http://julialang.org/downloads/ and go!

