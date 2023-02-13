# Particle-Filter
### Background
<p>
Particle filters are a sequential Monte Carlo methodology where we recursively compute probability distrubutions in hopes of accurately tracking a seemingly random signal or position of an object. The resultant goal is to always estimate the final state of a variable at each recursive loop through the algorithm.
</p>
The particle filter is similar to the Kalman filter in structure, however instead of tracking a model utilizing one state matrix at a time the particle filter implements a set of weighted matricies for tracking a multitude of non-linear models referred to as "particles".

### How the program works
This program reads in a file of data of three columns: actual position, actual velocity, and sensor reading.

The system outlined consists of a 1 Dimensional position moving in a line, where the data is recorded from two magnets. Utilizing a series of predict-update equations, we average the estimation of all the particles to determine the next state derivied from probability distributions.

### Results
Below represents a graphical representation of the output of the particle filter. The red line shows the estimation of our particle filter agaist the actual data plotted in gray. We can that, over time, the filter's use of particle sets provides a very accurate estimation of varied sinusoidal measurement data in a given system.
![image](https://user-images.githubusercontent.com/57545505/218597071-ba134cc4-b189-4881-9f7e-1587ea4ae79d.png)
