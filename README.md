# SDC - Path planning writeup

Henry Yau
This project is part of the path planning module for Udacity's Self Driving Car nanodegree


## Model Documentation
The goal of the project is to create a path planner for a vehicle which is the car is able to drive without incident around a three lane windy road loop with other vehicles travelling at various speeds. Frenet coordinates are used to represent the path of the vehicles and the road. With $d$ representing lateral displacement and $s$ representing longitudinal displacement.

The path generated is passed as a list of 50 $(x,y)$ coordinates in map space to the controller which the simulator applies to the vehicle every 0.02 seconds. When the path planner is called, the previous list of coordinates is passed back without the elements which have already been processed, thus the objective is the fill out the remaining coordinates to extend the path back to 50 elements.

The requirements include maximum jerk and acceleration bounds. To help satisfy these requirements, a cubic spline is used to represent the path. An open source spline implementation is used here:
http://kluge.in-chemnitz.de/opensource/spline/

The starting control points of the spline are the two last points of the previous trajectory. Three additional points are created based on the desired destination. In Frenet coordinates this means setting $s$ for three regular intervals ahead. If a lane change or any other lateral movements are needed $d$ is changed.

This spline can then be interpolated, getting the required coordinates to ensure the list has 50 points. As the points are regularly spaced in time, to account for speed, we can interpolate the spline at varying distances to change the speed. 
$$x[k+1]=x[k]+v[k]\Delta t$$
where the velocity at time step k, $v[k]$ is given in m/s.

## Determining speed and lane changing
This path planner uses a simple finite state machine to decide which actions to take based on the current situation around the vehicle. The area surrounding the vehicle is split into 7 segments shown below with 0 representing the vehicle.


|       |       |   |                   
|-----|------|----------|
|a |b|c    |
|d |0 |e   |
|f | |g|

If there is a car in cell $b$, the vehicle looks at cells $a$ and $c$. If no vehicles are in one or both of those two cells and the other surrounding cells, The vehicle attempts to pass $b$, with preference to passing on the left. If there are vehicles in $a$ or $c$, the speeds of the vehicles are compared to $b$. If the speeds are faster and the other cells are free from vehicles, a lane change is made. If none of this conditions are true, the vehicle is stuck behind the vehicle in $b$ and therefore needs to match its speed.

If there is no vehicle in $b$, we attempt to drive up to the speed limit.
