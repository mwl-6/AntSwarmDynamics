# Simulating Ant Swarm Aggregations Dynamics
![image](https://github.com/user-attachments/assets/a7dbf875-42b4-48c9-b2f8-934dc945e7b8)


## Repository Contents:

antsim.cpp - Contains the ant simulation

settings.txt - contains configuration details for various settings such as viscosity, adhesive strength, friction, etc. File contains instructions for how to modify values for each scene

Scenario models stored as .glb files

Other files part of Raylib library

## Compilation Instructions:

### Windows:

#### Option 1:
Download Raylib from: https://www.raylib.com/ via Windows installer

Alternatively, go to https://github.com/raysan5/raylib for repository. Links to setup instructions for various platforms available

After installing, go to installation folder and open w64devkit

EX: C:\raylib\w64devkit

Run "w64devkit.exe" and cd to repository folder

Compilation command with GCC/G++:

g++ antsim.cpp -lraylib -lopengl32 -lgdi32 -lwinmm -fopenmp -g -O2 -o antsim.exe

### Mac:

Download raylib-4.5.0_macos.tar.gz from https://github.com/raysan5/raylib/releases/tag/4.5.0

Copy libraylib.a from /lib to repository

Comment out #include <omp.h> in antsim.cpp on line 8

Compilation command:

clang++ -framework CoreVideo -framework IOKit -framework Cocoa -framework GLUT -framework OpenGL libraylib.a antsim.cpp -g -O2 -o antsim -std=c++11

## Simulation Controls:

Toggle between 3D cursor and camera using SHIFT

### 3D cursor:

WASD + QE to move ball around scene environment

Allows interactions/collisions with particles

1, 2, 3 to adjust cursor size

For plate scene, QE to move plate up and down instead of cursor

### Camera:

WASD + Mouse, scroll to zoom

IKJL + OU to change environment walls size

### Simulation:

4 to switch scenes

5 to reset simulation

6/7 to start/stop recording frame data to text file

8/9 to start/stop recording FPS/performance data to csv

Frame data file stored in following format:

index,x,y,z

Each line contains all of the particles for one frame

Mesh data file stored in the following format:

Model Name

Render Scale

List of x,y,z values

## How to Adjust Parameters:
See settings.txt

## Current Scenes:
Funnel, Bracket, Plate, Dispersal, Coin Drop, Ant Tower, Crushed Ants
