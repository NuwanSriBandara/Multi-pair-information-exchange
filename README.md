# Multi-pair-information-exchange-controlled-by-a-base-station
In this project, a single isolated cell for a multi-pair information exchange of 2 clusters (each with 3 users) was implemented using Mathworks MATLAB software, in regard to the controlling mechanism from a base station. The power control coordination is assumed to be generated by the base station and the distance is the factor to uniformly determine the power allocation for each user in each cluster.

Transmitted signal from base station to all clusters,

![first equation](https://latex.codecogs.com/gif.latex?y%20%3D%20%5Csum_%7Ba%3D1%7D%5E%7B2%7D%5Cleft%20%5B%20%5Csum_%7Bi%3D1%7D%5E%7B3%7D%5C%28h_%7Ba%2Ci%7D%29*%5C%5Csqrt%7B%5Cleft%20%28%20P_%7Bi%7D%5Cright%20%29%7D*%5C%28x_%7Ba%2Ci%7D%29%20%5Cright%20%5D%20&plus;n)

The relationship between transmit power (P_i) and total power (P_t),

![second equation](https://latex.codecogs.com/gif.latex?P%7B_i%7D%20%3D%20%5Cepsilon%20%7B_i%7D*%5Calpha%20%7B_i%7D*P%7B_t%7D%3B%20i%20%3D%201%2C2%2C3)

Here ![third equation](https://latex.codecogs.com/gif.latex?%5Csum_%7Bi%3D1%7D%5E%7B3%7D%20a%7B_i%7D%20%3D%201) are ![fourth equation](https://latex.codecogs.com/gif.latex?%5Cepsilon%20%7B_i%7D%20%3D) the inter pair power allocation coefficients where as ε_i= intra pair power allocation coefficient

The system is easily expandable to any cluster with any number of users for more wider implementation.
