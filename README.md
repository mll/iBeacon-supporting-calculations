# iBeacon-supporting-calculations

This repo contains some useful tools for iBeacon technology implementation. There are two versions of the Kalman Filter for predicting current distance to the beacon from several measurements (kalman.jl). distance.jl contains functions that compute current location based on those filtered signals coming from multiple beacons. It uses DBSCAN-Rtree clusteriser for speed and has computational complexity of O(N^2 log N).

Just download julia http://julialang.org/downloads/ and go!

