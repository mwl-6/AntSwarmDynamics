#include "raylib.h"
#include <iostream>
#include <vector>
#include <math.h>
#include "raymath.h"
#include "rlgl.h"
#include <algorithm>
#include <omp.h>
#include <map>
#include <string>
#include <fstream>
#include <iomanip>
#include <chrono>
#define RLIGHTS_IMPLEMENTATION
#include "rlights.h"

//Spatial Hash
#define BUCKET_WIDTH 15.0

//Particle radius
#define RADIUS 7

//Used to render scene at a smaller scale
#define GLOBAL_RENDER_SCALE 0.01f

//Limits the number of simulation frames that can be recorded on a text file
//Sizes will vary, but 1000 frames of simulation for 5500 particles is about 116 MB
#define FRAME_DATA_LIMIT 28000

//Timestep
float timestep = 0.01f;
//Gravity
float gravity = 280.f;

//Adjustable Controls (also in settings.txt)
double globalViscSetting = 0.0;
double globalBondDistance = 30.0;
double globalSplitSetting = 0.0;
int globalBondLimit = 2; 
double globalFrictionSetting = 0.0;
double globalParticleFriction = 0.0;
double globalWallBondStrength = 1.0;
//Measure of stiffness
float distanceConstraintIntensity = timestep;
//Maximum absolute velocity of a particle
const double cRange = 10/timestep;
//Floor friction values
float floorFriction[] = {0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8};
//How much influence does P1 have compared to P2
float distanceConstraintBalance = 0.5f;
//Default number of particles
int numParts = 5500;
//Number of particles by scene
int sceneMaxParticles[] = {8500, 8500, 8500, 8500, 8500, 8500, 8500};
//Solver iterations - usually 4
const int solverIter = 4;


//Boundaries of environment
float X_MIN = 0.0f;
float X_MAX = 2000.0f;
float Y_MIN = 0.0f;
float Y_MAX = 2380.0f;
float Z_MIN = -1000;
float Z_MAX = 1000;

//Bracket Scenario Variables
//Indicates locations of where there are holes that need to be filled
const int decisionRegionCount = 3;
const int fillThreshold[] = {300, 300, 300};
Vector2 decisionRegions[] = { {320, 580}, {770, 1030}, {1230, 1490}};
float holeXThreshold = 720;
int holeCounts[] = {0, 0, 0};



using namespace std;

const int screenWidth = 1280;
const int screenHeight = 720;



//Default scene settings stored in settings.txt
void loadSettingsFromFile(string fileName, int sceneNum){
    ifstream userFile;
	userFile.open(fileName.c_str());

	if (userFile.is_open()) {
		string curr;
        int line = 0;
        int fileScene = -1;
        //Ignore directions
        for(int i = 0; i < 12; i++){
            getline(userFile, curr);
        }
        
		while (getline(userFile, curr)) {
            //Get current scene from line number
            if(line % 10 == 0){
                line++;
                fileScene++;
                continue;
            }
            
            if(fileScene != sceneNum){
                line++;
                continue;
            }
            
            //Validate that line is either an int or float
            bool isFloat = false;
            bool pointFound = false;

            for(int i = 0; i < curr.length(); i++){
                if(isdigit(curr[i]) || (curr[i] == '-' && i == 0)){
                    continue;
                }
                else if(curr[i] == '.' && !pointFound){
                    pointFound = true;
                    isFloat = true;
                }
                else if(curr[i] == '.' && pointFound){
                    cout << "Invalid Number on Line: " << line+11 << endl;
                    exit(0);
                }
                else {
                    cout << "Invalid Number on Line: " << line+11 << endl;
                    exit(0);
                }
            }

            //Make sure max bonds is an integer
            if(isFloat && line % 10 == 4){
                cout << "Max Bonds must be an integer (see line " << line+11 << ")" << endl; 
                exit(0);
            }

            //Edit settings
            int setting = line % 10;
            switch(setting){
                case 1:
                    globalViscSetting = stof(curr);
                break;
                case 2:
                    distanceConstraintIntensity = stof(curr);
                break;
                case 3:
                    globalSplitSetting = stof(curr);
                break;
                case 4:
                    globalBondLimit = stoi(curr);
                break;
                case 5:
                    globalFrictionSetting = stof(curr);
                break;
                case 6:
                    globalWallBondStrength = stof(curr);
                break;
                case 7:
                    globalParticleFriction = stof(curr);
                break;
                case 8:
                    timestep = stof(curr);
                break;
                case 9:
                    gravity = stof(curr);
                break;
            }
            line++;
		}

	}
	else {
		cout << "Unable to open file: " + fileName << endl;
        exit(0);
	}
    userFile.close();
};

typedef struct Vector {
    double x;
    double y;
    double z;
} Vector;

//Used for spatial hash
typedef struct SpatialKey {
    int index;
    int bucket;
} SpatialKey;


//Hash Function by Matthias Müller
//See: https://github.com/matthias-research/pages/blob/master/tenMinutePhysics/11-hashing.pdf
int hashCoords(int x, int y, int z) { 
  int hashVal = (x * 92837111) ^ (y * 689287499) ^ (z * 283923481);
  return abs(hashVal);
}
//Used when sorting spatial hash array
int compareKeys(const void *a, const void *b){
  int l = ((struct SpatialKey*)a)->bucket;
  int r = ((struct SpatialKey*)b)->bucket;
  return l-r;
}

//Utility functions
float randRange(int min, int max){
	return (float)(rand() % (max-min) + min);
}
double distance(double x1, double y1, double z1, double x2, double y2, double z2){
	return sqrt((y2-y1)*(y2-y1) + (x2-x1)*(x2-x1) + (z2-z1)*(z2-z1));
}
double constrain(double val, double min, double max){
    if(val < min){
        return min;
    }
    else if(val > max){
        return max;
    }
    return val;
}


class Particle {
    private:
        
        float mass = 1.0;
        float invmass = 1/mass;
        
        Color color = {255, 0, 0};
        vector<int> neighbors;
        vector<Vector> gradients;
        
        //Current density
        double density = 0;
        
        //Target density
        double idealDensity = 0.00001;
        //Density constraint
        double Ci = density / idealDensity;
        //Lambda 
        double lambda = 0;

        double neighborSearchDistance = 60;

        int minNeighbors = 15;

        //Viscosity Force
        Vector viscosityForce = {0.0, 0.0, 0.0};
        
        //Is the ant touching the floor or on top of the pile
        bool onSurface = true;

        
    public:
        Vector3 pos;
        Vector3 opos;
        Vector3 vel = {0.0, 0.0, 0.0};
        double r;
        bool isWater;
        double viscosityMu = globalViscSetting;
        int activeDistConstraints = 0;

        //Flags used that affect final velocity calculations
        bool activateFriction = false;
        bool activateFloorFriction = false;
        bool particleFriction = false;

        //Is particle adhering to environment walls/floor
        bool wallBond = false;
        //Set true if particle is adhering to an object that moves
        bool objBond = false;
        //Location of adhesive bond
        float wallBondX;
        float wallBondY;
        float wallBondZ;

        //Used for testing
        bool setForce = false;
        Vector3 force = {0.0, 0.0, 0.0};

        //Used for random movement
        float forwardAngle = 0.0f;

        //Used for erasing distance constraints after a certain number of frames (dispersal scenario)
        float timeToSeverBonds;
        int severCounter = 0;
        bool severedBonds = false;
        bool activateRandomMovement = false;

        //For scene-specific behaviors
        int sceneNum = 0;

        //Used for freely moving ants traversing the surface of a mesh (bracket scenario)
        //Keeps track of the "gravity" vectors that push them towards the closest surface
        Vector3 closestSurface = {-100000, 0, 0};
        Vector3 oldSurface = closestSurface;

        //A vector parallel to the surface of a mesh triangle for random movement
        Vector3 randomForward;
        
        Particle(float x, float y, float z, double radius, bool water, double massV, int scene){
            pos = {x, y, z};
            opos = {x, y, z};
            r = radius;
            isWater = water;
            mass = massV;
            invmass = 1/mass;
            forwardAngle = GetRandomValue(0, 360);
            sceneNum = scene;
            
            timeToSeverBonds = GetRandomValue(720, 1460);
        }
        void draw(){
            Vector3 drawPos = {(float)pos.x, (float)pos.y, (float)pos.z};
            DrawSphere(drawPos, r, color);
        }
        void setColor(unsigned char r, unsigned char g, unsigned char b){
            color = {r, g, b, 255};
        }
        Color getColor(){
            return color;
        }
        double getDensity(){
            return density;
        }

        void applyForces(){
            //Break bonds after certain amount of time (only used in dispersal/bracket scenarios)
            if(activateRandomMovement){
                severCounter++;

                if(severCounter > timeToSeverBonds){
                    if(sceneNum != 1){
                        severedBonds = true;
                    }
                    else if(sceneNum == 1 && pos.x > holeXThreshold+5 && pos.y < decisionRegions[2].y){
                        severedBonds = true;
                        isWater = false;
                    }
                }
            }
            //Dispersal + bracket scene ground behavior random movement
            if(activateRandomMovement && (activateFloorFriction || (onSurface && sceneNum != 1)) && sceneNum != 7){
                
                if(rand() % 5 == 0){
                    if(GetRandomValue(0, 1) == 0){
                        forwardAngle -= 10;
                    }
                    else {
                        forwardAngle += 10;
                    }
                }
                float speed = 120.0f;

                Vector3 forward = {cos(forwardAngle * PI / 180.0f) * speed, 0, sin(forwardAngle * PI / 180.0f)*speed};
                
                
                vel = {forward.x, vel.y - invmass*gravity*timestep, forward.z};
            }//Bracket main mass behavior (filling bracket gaps)
            else if(setForce && !activateRandomMovement){
                //Bracket scenario
                //apply a rightward force to fill bracket spaces
                

                //Should the particle move right
                bool applyRightForce = false;
                //Is the particle inside a filled hole
                bool applyCorrectiveForce = false;

                //Is the particle within the hole space
                bool withinBounds = false;
                int region = -1;
                for(int i = 0; i < decisionRegionCount; i++){
                    if(pos.y > decisionRegions[i].x && pos.y < decisionRegions[i].y){
                        if(pos.x > holeXThreshold-10){
                            withinBounds = true;
                            region = i;
                        }

                        if(holeCounts[i] < fillThreshold[i]){
                            applyRightForce = true;
                            region = i;
                            break;
                        }
                        else if(withinBounds){
                            applyCorrectiveForce = true;
                            break;
                        }
                        
                        
                    }
                }

                //The hole count check is to account for extraneous particles that fall past the general mass of ants
                //Those particles shouldn't have forces applied or it will look unnatural
                if(applyRightForce && (region == 2 || (region == 1 && holeCounts[2] > fillThreshold[2]) || (region == 0 && holeCounts[1] > fillThreshold[1])) ){
                    //Fill gap by moving along ceiling
                    if(withinBounds){
                        vel = {vel.x + 4.f, vel.y, vel.z};
                    }
                    else if(pos.x > 500){
                        //If not in gap yet move right
                        vel = {vel.x + 4.f, vel.y - invmass*gravity*timestep, vel.z};
                    }
                }
                else if(applyCorrectiveForce){
                    //Optional additional forces to keep particles inside gap
                    vel = {vel.x, vel.y, vel.z};
                    
                }
                else {
                    //Continue as usual
                    vel = {vel.x, vel.y - invmass*gravity*timestep, vel.z};
                }
            }//Free ants scenario (bracket)
            else if(!isWater && closestSurface.x > -9999.0f){

                //"Gravity" force towards closest surface intended to keep ants attached to mesh
                float gFactor = 280.f;
                Vector3 gravityForce = Vector3Scale(closestSurface, gFactor);
                
                //Affects speed for walking along the walls of a mesh
                float speed = 10.0f;

                //First frame, perform set up
                if(oldSurface.x < -9999){
                    //0 degrees is the vector perpendicular to the gravity force vector
                    //It is then rotated around the gravity axis according to the ant's rotation
                    //That forward force is then stored in randomForward
                    Vector3 normal = Vector3Perpendicular(closestSurface);
                    normal = Vector3RotateByAxisAngle(normal, closestSurface, forwardAngle * PI / 180.0f);
                    randomForward = Vector3Scale(normal, speed);
                }
                else if(closestSurface.x != oldSurface.x || closestSurface.y != oldSurface.y || closestSurface.z != oldSurface.z){
                    /*
                    In the event the ant moves to a new triangle that has a different rotation than the former, the forward
                    vector must be adjusted accordingly in a way that feels consistent to "forward" movement

                    The idea here is to find the rotation between the old triangle normal and the new triangle normal,
                    then apply that rotation to the randomForward vector

                    This should only happen once when the ant changes triangles
                    */
                    
                    //Perform a rotation
                    Quaternion q = QuaternionFromVector3ToVector3(oldSurface, closestSurface);
                    randomForward = Vector3RotateByQuaternion(randomForward, q);
                    
                }

                //Change directions by rotating randomForward using the surface normal as an axis.
                if(rand() % 5 == 0){
                    if(GetRandomValue(0, 1) == 0){
                        randomForward = Vector3RotateByAxisAngle(randomForward, closestSurface, 0.03f);
                        
                    }
                    else {
                        randomForward = Vector3RotateByAxisAngle(randomForward, closestSurface, -0.03f);
                    }
                }
                
                

                //Apply this force if the ant is on the surface of a mesh
                if(activateFriction){
                    vel = Vector3Add(vel, randomForward);
                }
                
                vel = Vector3Add(vel, Vector3Scale(gravityForce, timestep));
            }//Regular ant forces
            else{
                vel = {vel.x, vel.y - invmass*gravity*timestep, vel.z};
            }
            
            
            activateFloorFriction = false;
            activateFriction = false;
            onSurface = true;
        }
        void predictPosition(){
            pos = {pos.x + vel.x*timestep, pos.y + vel.y*timestep, pos.z + vel.z*timestep};
        }
        void setVel(float x, float y, float z){
            vel = {x, y, z};
        }
        void unboundedPerformCollisions(int index, SpatialKey *keyTable, int *startingIndexTable, int tableSize, vector<Particle*> &particles){
            int bX = floor(particles[index]->pos.x / BUCKET_WIDTH);
            int bY = floor(particles[index]->pos.y / BUCKET_WIDTH);
            int bZ = floor(particles[index]->pos.z / BUCKET_WIDTH);
            if(!isWater){
                return;
            }

            //Check for collisions in current square and adjacent
            for(int i = -1; i <= 1; i++){
                for(int j = -1; j <=1; j++){
                    for(int k = -1; k <= 1; k++){
                        //Get bucket
                        int bucket = hashCoords(bX+i, bY+j, bZ+k) % tableSize;

                        //Find starting index of index table given bucket and iterate through its contents
                        int startingIndex = startingIndexTable[bucket];

                        while(startingIndex != -1 && startingIndex < particles.size() && keyTable[startingIndex].bucket == bucket){
                            int pIndex = keyTable[startingIndex].index;
                            
                            if(pIndex == index){
                                startingIndex++;
                                continue;
                                
                            }

                            
                            //Check for intersection of particles, if so move apart
                            const double pX = particles[pIndex]->pos.x;
                            const double pY = particles[pIndex]->pos.y;
                            const double pZ = particles[pIndex]->pos.z;
                            double dX = (pos.x - pX);
                            double dY = (pos.y - pY);
                            double dZ = (pos.z - pZ);
                            double d = sqrt(dX*dX + dY*dY + dZ*dZ);
                            
                            double normaldX = dX/d;
                            double normaldY = dY/d;
                            double normaldZ = dZ/d;
                            double cd = d-(2*r);
                            if(cd < 0 && d > 0.0001){
                                //Move particles away from each other
                                if(!particles[pIndex]->isWater || !isWater){
                                    startingIndex++;
                                    continue;
                                }

                                pos.x -= normaldX * cd * 0.5;
                                pos.y -= normaldY * cd * 0.5;
                                pos.z -= normaldZ * cd * 0.5;
                                particles[pIndex]->pos.x += normaldX * cd * 0.5;
                                particles[pIndex]->pos.y += normaldY * cd * 0.5;
                                particles[pIndex]->pos.z += normaldZ * cd * 0.5;


                                particleFriction = true;
                                particles[pIndex]->particleFriction = true;
                                
                                

                                //Friction Implementation
                                Vector3 deltaPI = {(float)(pos.x - opos.x), (float)(pos.y - opos.y), (float)(pos.z - opos.z)};
                                Vector3 deltaPJ = {(float)(particles[pIndex]->pos.x - particles[pIndex]->opos.x), (float)(particles[pIndex]->pos.y - particles[pIndex]->opos.y), (float)(particles[pIndex]->pos.z - particles[pIndex]->opos.z)};
                                Vector3 netdelta = Vector3Subtract(deltaPI, deltaPJ);
                                float netdeltaMag = distance(0, 0, 0, netdelta.x, netdelta.y, netdelta.z);
                                
                                Vector3 normal = {(float)(pos.x - particles[pIndex]->pos.x), (float)(pos.y - particles[pIndex]->pos.y), (float)(pos.z - particles[pIndex]->pos.z)};
                                normal = Vector3Normalize(normal);
                                Vector3 perpendicularNormal = Vector3Perpendicular(normal);

                                //Find components perpendicular to the normal vector
                                //Projection of netdelta into perpendicularNormal
                                Vector3 tangentialDisplacement = Vector3Scale(perpendicularNormal, Vector3DotProduct(netdelta, perpendicularNormal));
                                float frictionMu = globalParticleFriction;
                                float tangentialScale = distance(0, 0, 0, tangentialDisplacement.x, tangentialDisplacement.y, tangentialDisplacement.z);
                                
                                if(tangentialScale == 0 || isnan(tangentialDisplacement.x) || isnan(tangentialDisplacement.y) || isnan(tangentialDisplacement.z)){
                                    startingIndex++;
                                    continue;
                                }
                                
                                pos.x += -tangentialDisplacement.x * frictionMu;
                                pos.y += -tangentialDisplacement.y * frictionMu;
                                pos.z += -tangentialDisplacement.z * frictionMu;
                                particles[pIndex]->pos.x -= -tangentialDisplacement.x * frictionMu;
                                particles[pIndex]->pos.y -= -tangentialDisplacement.y * frictionMu;
                                particles[pIndex]->pos.z -= -tangentialDisplacement.z * frictionMu;


                                

                            }
                            
                            startingIndex++;
                        }
                    }
                }
            }

            if(isnan(pos.x) || isnan(pos.y) || isnan(pos.z)){
                cout << "Error at Perform Collisions" << endl;
                exit(1);
            }

        }
        //Cursor collision
        void ballCollisions(float &x, float &y, float &z, float nr){
            const float pX = x;
            const float pY = y;
            const float pZ = z;
            float dX = (pos.x - pX);
            float dY = (pos.y - pY);
            float dZ = (pos.z - pZ);
            float d = sqrt(dX*dX + dY*dY + dZ*dZ);
            
            float normaldX = dX/d;
            float normaldY = dY/d;
            float normaldZ = dZ/d;
            float cd = d-(nr+r);
            
            
            if(cd < 0 && d > 0.0001){
                pos.x -= normaldX * cd;
                pos.y -= normaldY * cd;
                pos.z -= normaldZ * cd;
            }
            
            
        }
        //Make sure particle stays within walls
        //Adhesive bonds may also be activated if wall/floor collisions are detected
        void keepInWalls(){
            if(pos.y + r / 2 > Y_MAX){
                pos.y += (Y_MAX - r/2 - pos.y);
                activateFriction = true;

            }
            if(pos.y - r / 2 < Y_MIN){
                pos.y += (Y_MIN- (pos.y-r/2));
                activateFriction = true;
                activateFloorFriction = true;
                onSurface = true;

                if(!setForce){
                    wallBondX = pos.x;
                    wallBondY = pos.y;
                    wallBondZ = pos.z;

                    
                    wallBond = true;
                    
                    if(sceneNum != 2 && sceneNum != 6)
                    objBond = true;
                }
            }
            
            if(pos.x+r  > X_MAX){
                pos.x += (X_MAX - (pos.x+r));
                activateFriction = true;

                if(!setForce){
                    wallBondX = pos.x;
                    wallBondY = pos.y;
                    wallBondZ = pos.z;
                    wallBond = true;
                }
                
            }
            if(pos.x-r  < X_MIN){
                pos.x += (X_MIN - (pos.x-r));
                activateFriction = true;

                if(!setForce){
                    wallBondX = pos.x;
                    wallBondY = pos.y;
                    wallBondZ = pos.z;
                    wallBond = true;
                }
            }
            if(pos.z+r  > Z_MAX){
                pos.z += (Z_MAX - (pos.z+r));
                activateFriction = true;

                if(!setForce){
                    wallBondX = pos.x;
                    wallBondY = pos.y;
                    wallBondZ = pos.z;
                    wallBond = true;
                }
            }
            if(pos.z-r  < Z_MIN){
                pos.z += (Z_MIN - (pos.z-r));
                activateFriction = true;

                if(!setForce){
                    wallBondX = pos.x;
                    wallBondY = pos.y;
                    wallBondZ = pos.z;
                    wallBond = true;
                }
            }



            if(setForce && isWater){
                if(pos.y > 1900 && pos.x < 700){
                    pos.x = 700;
                }
                if(pos.y > 1700 && pos.x > 900){
                    pos.x = 900;
                }

                //Wall boundaries defined here for bracket scenario
                //so that it no longer applies once particle touches ground
                if(pos.x < 480.0f && !activateRandomMovement){
                    pos.x = 480.0f;
                }
                if(pos.x > 1250.0f && !activateRandomMovement){
                    pos.x = 1250.0f;
                }
                if(pos.y < 0.0f && !activateRandomMovement){
                    pos.y = 0.0f;
                }
                if(pos.y > 2380.0f && !activateRandomMovement){
                    pos.y = 2380.0f;
                }
                if(pos.z < -350.0f && !activateRandomMovement){
                    pos.z = -350.0f;
                }
                if(pos.z > -100.0f && !activateRandomMovement){
                    pos.z = -100.0f;
                }

                if(pos.y < 20 && !activateRandomMovement){
                    activateRandomMovement = true;
                    timeToSeverBonds = 1;
                }

            }

            //Funnel scene, have ants disperse on ground
            if(sceneNum == 0 && pos.y < 20 && !activateRandomMovement){
                activateRandomMovement = true;
                timeToSeverBonds = 1;
            }

            if(isnan(pos.x) || isnan(pos.y) || isnan(pos.z)){
                cout << "Error at Keep In Walls" << endl;
                exit(1);
            }

        }
        
        //Find neighboring particles and calculate density lambda
        void findNeighbors(int j, SpatialKey *keyTable, int *startingIndexTable, int tableSize, vector<Particle*> &particles, map<pair<int, int>, float> &restDistance){
            neighbors.clear();
            
            
            const double poly6h = neighborSearchDistance;
            
            double pi = 0;
            //Each gradient calculation is slightly different depending on whether or not its a neighbor
            Vector localGradient = {0.0, 0.0, 0.0};
            double localGradientMag = 0.0;


            gradients.clear();
            viscosityForce = {0, 0, 0};

            double gradMag = 0;
            double polyConst = 315/(64*PI*pow(poly6h, 9));
            double spikyConst = -45/(PI*pow(poly6h, 6));
            double viscosityConst = 45/(PI*pow(poly6h, 6));
            
            
            //Find how many adjacent spatial hash squares to check
            int viewerDepth = floor(neighborSearchDistance / BUCKET_WIDTH);
            int bX = floor(particles[j]->pos.x / BUCKET_WIDTH);
            int bY = floor(particles[j]->pos.y / BUCKET_WIDTH);
            int bZ = floor(particles[j]->pos.z / BUCKET_WIDTH);
            
            
            
            for(int i = -viewerDepth; i <= viewerDepth; i++){
                for(int k = -viewerDepth; k <= viewerDepth; k++){
                    for(int l = -viewerDepth; l <= viewerDepth; l++){

                        //Search through range of neighbors after finding particle bucket
                        int bucket = hashCoords(bX+i, bY+k, bZ+l) % tableSize;
                        int startingIndex = startingIndexTable[bucket];

                        while(startingIndex != -1 && startingIndex < particles.size() && keyTable[startingIndex].bucket == bucket){
                            int pIndex = keyTable[startingIndex].index;
                            
                            if(pIndex == j){
                                
                                startingIndex++;
                                continue;
                            }

                            const double pX = particles[pIndex]->pos.x;
                            const double pY = particles[pIndex]->pos.y;
                            const double pZ = particles[pIndex]->pos.z;

                            if(pY-pos.y > 40){
                                onSurface = false;
                            }

                            //Make sure neighbor particle is in range to be considered a neighbor
                            if(distance(pos.x, pos.y, pos.z, pX, pY, pZ) < neighborSearchDistance){
                                
                                neighbors.push_back(pIndex);
                                double relX = pos.x - pX;
                                double relY = pos.y - pY;
                                double relZ = pos.z - pZ;
                                double relMag = distance(0, 0, 0, relX, relY, relZ);

                                if(relX == 0 && relY == 0 && relZ == 0){
                                    startingIndex++;
                                    continue;
                                }

                                if(relMag < globalBondDistance && activeDistConstraints <= globalBondLimit && !severedBonds){
                                    
                                    #pragma omp critical
                                    {
                                        //Create distance constraint bonds if bond space avaliable
                                        if(particles[pIndex]->activeDistConstraints <= globalBondLimit && !particles[pIndex]->severedBonds){
                                            if(j < pIndex && restDistance.count(make_pair(j, pIndex)) == 0){
                                                restDistance[make_pair(j, pIndex)] = globalBondDistance;
                                                activeDistConstraints++;
                                                particles[pIndex]->activeDistConstraints++;
                                            
                                            }
                                            else if (restDistance.count(make_pair(pIndex, j)) == 0){
                                                restDistance[make_pair(pIndex, j)] = globalBondDistance;
                                                activeDistConstraints++;
                                                particles[pIndex]->activeDistConstraints++;
                                            
                                            }
                                        }
                                        
                                    }
                                    
                                }
                                
                                
                                if(relMag < poly6h){
                                    
                                    //Poly6
                                    double p6val = polyConst * pow((pow(poly6h, 2) - pow(relMag, 2)), 3);
                                    pi += p6val;
                                    //Local implementation of spikygrad
                                    double factor = spikyConst * pow(poly6h-relMag, 2) / relMag;
                                    
                                    
                                    //Compute gradient of neighbors
                                    relX/=relMag; 
                                    relY/=relMag;
                                    relZ/=relMag;
                                    
                                    Vector spikygradMag = {(relX*factor), (relY*factor), (relZ*factor)};

                                    //k = i case for lambda calculation
                                    localGradientMag += factor*factor;
                                    
                                    //Gradient of neighbor 
                                    gradients.push_back({spikygradMag.x, spikygradMag.y, spikygradMag.z});

                                    //Instead do the gradient sum here rather than later
                                    //k = j case for lambda calculation
                                    gradMag += factor*factor;


                                    //Everything below relates to viscosity force correction
                                    double viscosityFactor = viscosityConst * (poly6h - relMag);
                                    Vector relVel = {(particles[pIndex]->vel.x - vel.x) / density, (particles[pIndex]->vel.y - vel.y) /density, (particles[pIndex]->vel.z - vel.z) /density};
                                    viscosityForce.x += relVel.x * viscosityFactor;
                                    viscosityForce.y += relVel.y * viscosityFactor;
                                    viscosityForce.z += relVel.z * viscosityFactor;
                                    
                                }
                            }

                            startingIndex++;

                        }
                    }
                }

            }
            
            
            density = pi;
            

            Ci = density/idealDensity - 1;
            
            //Square k=i case
            gradMag += localGradientMag;
            gradMag/=idealDensity;
            
            if(gradMag == 0){
                lambda = 0;
                return;
            }

            //Now that a list of gradients for all particles are gathered we can solve for lambda
            const double e = 1;
            //e is a relaxation parameter that softens the constraint
            double nlambda = -Ci / (gradMag + e);

            
            lambda = nlambda;

            if(isnan(pos.x) || isnan(pos.y) || isnan(pos.z)){
                cout << "Error at Find Neighbors" << endl;
                exit(1);
            }
        
        }
        
        //Position corrections from fluid constraints
        void applyLambdaPositionUpdate(vector<Particle*> &particles){

            double updateX = 0;
            double updateY = 0;
            double updateZ = 0;

            
            if(neighbors.size() < minNeighbors){
                return;
            }
            
            for(int i = 0; i < neighbors.size(); i++){
                int index = neighbors[i];
                Vector spikygradMag = gradients[i];
                double neighbor_lambda = particles[index]->lambda;

                updateX += spikygradMag.x * (lambda + neighbor_lambda);
                updateY += spikygradMag.y * (lambda + neighbor_lambda);
                updateZ += spikygradMag.z * (lambda + neighbor_lambda);

                
            
            }
            
            updateX /= idealDensity;
            updateY /= idealDensity;
            updateZ /= idealDensity;
            

            //Apply position corrections according to density constraint and viscosity constraint
            pos.x += constrain(updateX, -20, 20);
            pos.y += constrain(updateY, -20, 20);
            pos.z += constrain(updateZ, -20, 20);
            
            if(!isnan(viscosityForce.x) && !isnan(viscosityForce.y) && !isnan(viscosityForce.z)){
                pos.x += constrain(viscosityForce.x * viscosityMu, -10, 10);
                pos.y += constrain(viscosityForce.y * viscosityMu, -10, 10);
                pos.z += constrain(viscosityForce.z * viscosityMu, -10, 10);
            }

        

            if(isnan(pos.x) || isnan(pos.y) || isnan(pos.z)){
                cout << "Error at Apply Lambda Position Update" << endl;
                exit(1);
            }
        }

        
};

//Used for storing mesh triangle data
typedef struct Triangle {
    Vector3 a;
    Vector3 b;
    Vector3 c;
} Triangle;

//Used in collision detection
typedef struct pLine {
    Vector3 point;
    float dist;
} pLine;

//Gets closest point on a line segment to any other point
//see https://paulbourke.net/geometry/pointlineplane/
pLine pointToLine(Vector3 p1, Vector3 p2, Vector3 p3){
    //p = p1 + u(p2-p1)

    pLine res;
    res.dist = 999999;
    res.point = {0, 0, 0};
    float d2 = Vector3DistanceSqr(p1, p2);
    
    float u = ( ( ( p3.x - p1.x ) * ( p2.x - p1.x ) ) +
        ( ( p3.y - p1.y ) * ( p2.y - p1.y ) ) +
        ( ( p3.z - p1.z ) * ( p2.z - p1.z ) ) ) / d2;
    if(u < 0){
        u = 0;
    }
    else if(u > 1){
        u = 1;
    }
    res.point.x = p1.x + u*(p2.x-p1.x);
    res.point.y = p1.y + u*(p2.y-p1.y);
    res.point.z = p1.z + u*(p2.z-p1.z);
    res.dist = Vector3Distance(p3, res.point);
    return res;
}

//Counts the number of times a ray intersects with a triangle on a mesh
int countIntersections(Mesh *mesh, Ray testRay, Vector3 meshPosition){
    int intersections = 0;
    for(int i = 0; i < mesh->triangleCount * 3; i+=3){
        //Extract triangle information from indices and vertices arrays
        int index = mesh->indices[i];
        Triangle t;
        t.a = {(mesh->vertices[index * 3] + meshPosition.x) / GLOBAL_RENDER_SCALE, (mesh->vertices[index * 3 + 1] + meshPosition.y) / GLOBAL_RENDER_SCALE, (mesh->vertices[index * 3 + 2] + meshPosition.z) / GLOBAL_RENDER_SCALE};
        index = mesh->indices[i+1];
        t.b = {(mesh->vertices[index * 3] + meshPosition.x) / GLOBAL_RENDER_SCALE, (mesh->vertices[index * 3 + 1] + meshPosition.y) / GLOBAL_RENDER_SCALE, (mesh->vertices[index * 3 + 2] + meshPosition.z) / GLOBAL_RENDER_SCALE};
        index = mesh->indices[i+2];
        t.c = {(mesh->vertices[index * 3] + meshPosition.x) / GLOBAL_RENDER_SCALE, (mesh->vertices[index * 3 + 1] + meshPosition.y) / GLOBAL_RENDER_SCALE, (mesh->vertices[index * 3 + 2] + meshPosition.z) / GLOBAL_RENDER_SCALE};

        //Check if ray intersects with triangle on mesh
        RayCollision triangleTest = GetRayCollisionTriangle(testRay, t.a, t.b, t.c);

        if(triangleTest.hit){
            intersections++;
        }
    }
    return intersections;
}
Vector3 findClosestTrianglePoint(Mesh *mesh, Vector3 pos, Vector3 opos, Vector3 meshPosition){
    Vector3 closestPoint = {0, 0, 0};
    float closestDistance = 1000000.f;
    
    
    for(int i = 0; i < mesh->triangleCount * 3; i+=3){
        //Extract triangle information from indices and vertices arrays
        int index = mesh->indices[i];
        Triangle t;
        t.a = {(mesh->vertices[index * 3] + meshPosition.x) / GLOBAL_RENDER_SCALE, (mesh->vertices[index * 3 + 1] + meshPosition.y) / GLOBAL_RENDER_SCALE, (mesh->vertices[index * 3 + 2] + meshPosition.z) / GLOBAL_RENDER_SCALE};
        index = mesh->indices[i+1];
        t.b = {(mesh->vertices[index * 3] + meshPosition.x) / GLOBAL_RENDER_SCALE, (mesh->vertices[index * 3 + 1] + meshPosition.y) / GLOBAL_RENDER_SCALE, (mesh->vertices[index * 3 + 2] + meshPosition.z) / GLOBAL_RENDER_SCALE};
        index = mesh->indices[i+2];
        t.c = {(mesh->vertices[index * 3] + meshPosition.x) / GLOBAL_RENDER_SCALE, (mesh->vertices[index * 3 + 1] + meshPosition.y) / GLOBAL_RENDER_SCALE, (mesh->vertices[index * 3 + 2] + meshPosition.z) / GLOBAL_RENDER_SCALE};

        
        //Use triangle normal to compute perpendicular distance to triangle
        Vector3 N = Vector3Normalize(Vector3CrossProduct(Vector3Subtract(t.b, t.a), Vector3Subtract(t.c, t.a)) );
        
        
        Ray testRay;
        testRay.position = pos;
        testRay.direction = N;
       
        
        /*
        Check to see if there is a valid path from the particle's current position in the direction of the mesh triangle's
        normal that leads to the surface of that triangle
        */
        RayCollision triangleTest = GetRayCollisionTriangle(testRay, t.a, t.b, t.c);
        
        /*
        If successful and the distance is the shortest current distance to the mesh's surface,
        record the point
        */
        if(triangleTest.hit){
            float d = Vector3Distance(triangleTest.point, pos);
            if(triangleTest.distance < closestDistance){
                closestDistance = triangleTest.distance;
                closestPoint = triangleTest.point;
                
            }
        }
        
        
        //See if closest edge or vertex is applicable
        pLine al;
        pLine bl;
        pLine cl;

        //Check if moving particle to a triangle vertex is the best option
        al.point = t.a;
        al.dist = Vector3Distance(t.a, pos);
        bl.point = t.b;
        bl.dist = Vector3Distance(t.b, pos);
        cl.point = t.c;
        cl.dist = Vector3Distance(t.c, pos);

        pLine ideal = al;
        if(bl.dist < ideal.dist){
            ideal = bl;
        }
        if(cl.dist < ideal.dist){
            ideal = cl;
        }
        if(ideal.dist < closestDistance){
            closestDistance = ideal.dist;
            closestPoint = ideal.point;
        }
        
        //Check if moving particle to a triangle line segment is the best option
        pLine ab = pointToLine(t.a, t.b, pos);
        pLine bc = pointToLine(t.b, t.c, pos);
        pLine ac = pointToLine(t.a, t.c, pos);

        pLine ideal2 = ab;
        if((bc.dist < ab.dist && bc.dist >= 0) || ab.dist < 0){
            ideal2 = bc;
        }
        if((ac.dist < ideal2.dist && ac.dist >= 0) || ideal2.dist < 0){
            ideal2 = ac;
        }

        if(ideal2.dist >= 0 && ideal2.dist < closestDistance && !triangleTest.hit){
            closestDistance = ideal2.dist;
            closestPoint = ideal2.point;
        }
                

        
    }
    
    return closestPoint;
}

//This function returns a vector towards the closest triangle point
Vector3 findClosestTrianglePointFromOutside(Mesh *mesh, Vector3 pos, Vector3 opos, Vector3 meshPosition){
    Vector3 closestPoint = {0, 0, 0};
    float closestDistance = 1000000.f;
    Vector3 closestForce = {0, 0, 0};
    
    
    for(int i = 0; i < mesh->triangleCount * 3; i+=3){
        //Extract triangle information from indices and vertices arrays
        int index = mesh->indices[i];
        Triangle t;
        t.a = {(mesh->vertices[index * 3] + meshPosition.x) / GLOBAL_RENDER_SCALE, (mesh->vertices[index * 3 + 1] + meshPosition.y) / GLOBAL_RENDER_SCALE, (mesh->vertices[index * 3 + 2] + meshPosition.z) / GLOBAL_RENDER_SCALE};
        index = mesh->indices[i+1];
        t.b = {(mesh->vertices[index * 3] + meshPosition.x) / GLOBAL_RENDER_SCALE, (mesh->vertices[index * 3 + 1] + meshPosition.y) / GLOBAL_RENDER_SCALE, (mesh->vertices[index * 3 + 2] + meshPosition.z) / GLOBAL_RENDER_SCALE};
        index = mesh->indices[i+2];
        t.c = {(mesh->vertices[index * 3] + meshPosition.x) / GLOBAL_RENDER_SCALE, (mesh->vertices[index * 3 + 1] + meshPosition.y) / GLOBAL_RENDER_SCALE, (mesh->vertices[index * 3 + 2] + meshPosition.z) / GLOBAL_RENDER_SCALE};

        
        //Use triangle normal to compute perpendicular distance to triangle
        Vector3 N = Vector3Normalize(Vector3CrossProduct(Vector3Subtract(t.b, t.a), Vector3Subtract(t.c, t.a)) );
        Vector3 IN = Vector3Normalize(Vector3CrossProduct(Vector3Subtract(t.a, t.b), Vector3Subtract(t.c, t.a)) );
   
        Vector3 P = Vector3Perpendicular(N);  
        
        Ray testRay;
        testRay.position = Vector3Add(pos, Vector3Scale(IN, 1));
        testRay.direction = N;
       
        Ray testRay2;
        testRay2.position = Vector3Add(pos, Vector3Scale(N, 1));
        testRay2.direction = IN;
        
        /*
        Check to see if there is a valid path from the particle's current position in the direction of the mesh triangle's
        normal that leads to the surface of that triangle
        */
        RayCollision triangleTest;
        RayCollision triangleTestOUT;
        RayCollision triangleTestIN;

        
        //triangleTestOUT = GetRayCollisionTriangle(testRay, t.a, t.b, t.c);
        triangleTestIN = GetRayCollisionTriangle(testRay2, t.a, t.b, t.c);
        
        /*
        If successful and the distance is the shortest current distance to the mesh's surface,
        record the point
        */
        
        if(triangleTestIN.hit){
            
            if(triangleTestIN.distance < closestDistance){
                closestDistance = triangleTestIN.distance;
                closestPoint = triangleTestIN.point;
                closestForce = Vector3Normalize(Vector3Subtract(closestPoint, pos));
            }
        }
        
        //Check if moving particle to a triangle line segment is the best option
        pLine ab = pointToLine(t.a, t.b, pos);
        pLine bc = pointToLine(t.b, t.c, pos);
        pLine ac = pointToLine(t.a, t.c, pos);



        pLine ideal2 = ab;
        if((bc.dist < ab.dist && bc.dist >= 0) || ab.dist < 0){
            ideal2 = bc;
        }
        if((ac.dist < ideal2.dist && ac.dist >= 0) || ideal2.dist < 0){
            ideal2 = ac;
        }

        if(ideal2.dist >= 0 && ideal2.dist < closestDistance && !triangleTestIN.hit){
            closestDistance = ideal2.dist;
            closestPoint = ideal2.point;
            closestForce = IN;
        }
                

        
    }
    
    return closestForce;
}

bool isObjectOnMesh(Mesh *mesh, Vector3 pos, Vector3 meshPosition){
    bool isOn = false;

    for(int i = 0; i < mesh->triangleCount * 3; i+=3){
        //Extract triangle information from indices and vertices arrays
        int index = mesh->indices[i];
        Triangle t;
        t.a = {(mesh->vertices[index * 3] + meshPosition.x) / GLOBAL_RENDER_SCALE, (mesh->vertices[index * 3 + 1] + meshPosition.y) / GLOBAL_RENDER_SCALE, (mesh->vertices[index * 3 + 2] + meshPosition.z) / GLOBAL_RENDER_SCALE};
        index = mesh->indices[i+1];
        t.b = {(mesh->vertices[index * 3] + meshPosition.x) / GLOBAL_RENDER_SCALE, (mesh->vertices[index * 3 + 1] + meshPosition.y) / GLOBAL_RENDER_SCALE, (mesh->vertices[index * 3 + 2] + meshPosition.z) / GLOBAL_RENDER_SCALE};
        index = mesh->indices[i+2];
        t.c = {(mesh->vertices[index * 3] + meshPosition.x) / GLOBAL_RENDER_SCALE, (mesh->vertices[index * 3 + 1] + meshPosition.y) / GLOBAL_RENDER_SCALE, (mesh->vertices[index * 3 + 2] + meshPosition.z) / GLOBAL_RENDER_SCALE};

        
        //Use triangle normal to compute perpendicular distance to triangle
        Vector3 N = Vector3Normalize(Vector3CrossProduct(Vector3Subtract(t.b, t.a), Vector3Subtract(t.c, t.a)) );
        Vector3 IN = Vector3Normalize(Vector3CrossProduct(Vector3Subtract(t.a, t.b), Vector3Subtract(t.c, t.a)) );

        Ray testRay;
        testRay.position = Vector3Add(pos, Vector3Scale(IN, 1));
        testRay.direction = N;
       
        Ray testRay2;
        testRay2.position = Vector3Add(pos, Vector3Scale(N, 1));
        testRay2.direction = IN;
        
        RayCollision triangleTestOUT;
        RayCollision triangleTestIN;

        
        //Send two rays to check: one going out and one going in
        
        triangleTestOUT = GetRayCollisionTriangle(testRay, t.a, t.b, t.c);
        triangleTestIN = GetRayCollisionTriangle(testRay2, t.a, t.b, t.c);
        
        //If the raycast hits fall within particle radius, the particle is considered to be on the surface
        //Right now this applies to both sides of the mesh triangles
        if(triangleTestOUT.hit && triangleTestOUT.distance <= 7.f){
            isOn = true;
            break;
        }
        if(triangleTestIN.hit && triangleTestIN.distance <= 7.f){
            isOn = true;
            break;
        }
        
    }
    
    return isOn;
}

void ApplyMeshCollision(Mesh *mesh, Particle *particle, Vector3 meshPosition, bool saferCollisions, BoundingBox meshBounds){
    
    //Mesh bounding box
    Vector3 minB = meshBounds.min;
    Vector3 maxB = meshBounds.max;


    //Simple bounding box check
    if(!(particle->pos.x > minB.x && particle->pos.x < maxB.x 
    && particle->pos.y > minB.y && particle->pos.y < maxB.y
    && particle->pos.z > minB.z && particle->pos.z < maxB.z)  ){
        return;
    }

    /*
    Take a ray from the particle's position to some direction
    If the number of times the ray passes through the mesh is odd, it is inside the mesh and a collision is taking place
    */
    bool hasCollided = false;
    int intersections = 0;

    Ray testRay;
    testRay.position = particle->pos;
    testRay.direction = {1, 0, 0};

    intersections = countIntersections(mesh, testRay, meshPosition);
    
    /*
    On the rare chance the test ray runs along a vertex or edge, intersection counts might increase,
    leading to a false assessment of whether point is within mesh.

    If safercollisions enabled, check with another ray to confirm there are no discrepancies
    If so, proceed as usual. If not, use third check to classify.
    */
    if(saferCollisions){
        testRay.direction = {0, 1, 0};
        int secondIntersections = countIntersections(mesh, testRay, meshPosition);
        if(intersections % 2 == secondIntersections % 2){
            if(intersections % 2 == 1){
                hasCollided = true;
            }
        }
        else {
            testRay.direction = {0, 0, 1};
            int thirdIntersections = countIntersections(mesh, testRay, meshPosition);
            if(thirdIntersections % 2 == 1){
                hasCollided = true;
            }
        }
    }
    else {
        if(intersections % 2 == 1){
            hasCollided = true;
        }
        else {
            return;
        }
        
        
    }
    

    if(hasCollided){
        //Move the particle to the closest available point on the mesh surface
        Vector3 result = findClosestTrianglePoint(mesh, particle->pos, particle->opos, meshPosition);
        if(!(result.x == 0 && result.y == 0 && result.z == 0)){
            
            
            particle->pos = result;

            

            //Establish an adhesive bond
            if(!particle->wallBond){
                particle->wallBondX = result.x;
                particle->wallBondY = result.y;
                particle->wallBondZ = result.z;
                particle->wallBond = true;
                particle->objBond = true;
            }
            
            particle->activateFriction = true;
        }
        
        
    }
    
}

//Cylinder Collision
void ApplyCoinCollision(Particle *particle, Vector3 &coinPosition, float coinRadius, float coinDepth){
    
    float dX = (particle->pos.x - coinPosition.x);
    float dY = (particle->pos.y - coinPosition.y);
    float dZ = (particle->pos.z - coinPosition.z);
    float d = sqrt(dX*dX + dY*dY + dZ*dZ);
    
    float normaldX = dX/d;
    float normaldY = dY/d;
    float normaldZ = dZ/d;
    float cd = d-(particle->r+coinRadius);

    float weight = 0.99978f;
    if(cd < 0 && d > 0.0001 && particle->pos.z >= coinPosition.z - coinDepth && particle->pos.z <= coinPosition.z + coinDepth){
        
        particle->pos.x -= normaldX * cd * weight;
        particle->pos.y -= normaldY * cd * weight;
        particle->pos.z -= normaldZ * cd * weight;

        coinPosition.x += normaldX * cd * (1-weight);
        coinPosition.y += normaldY * cd * (1-weight);
        coinPosition.z += normaldZ * cd * (1-weight);
        
    }
    
}

//Scene configuration
void initializeParticles(vector<Particle*> &particles, int numParts, int sceneNum, map<pair<int, int>, float> &restDistance, Matrix *transforms, bool withPreset){
    particles.clear();
    restDistance.clear();

    //Maintain default boundaries where needed
    if(sceneNum != 4 && sceneNum != 1){
        X_MIN = 0.0f;
        X_MAX = 2000.0f;
        Y_MIN = 0.0f;
        Y_MAX = 2380.0f;
        Z_MIN = -1000;
        Z_MAX = 1000;
    }
    else if(sceneNum == 1){
        X_MIN = 0.0f;
        X_MAX = 2000.0f;
        Y_MIN = 0.0f;
        Y_MAX = 2380.0f;
        Z_MIN = -1000;
        Z_MAX = 1000;
    }
    else {
        X_MIN = 0.0f;
        X_MAX = 750.0f;
        Y_MIN = 0.0f;
        Y_MAX = 1450.0f;
        Z_MIN = -25;
        Z_MAX = 25;
    }

    if(withPreset){
        loadSettingsFromFile("settings.txt", sceneNum);
    }

    //Funnel
    if(sceneNum == 0){
        //Spawn a cube block of particles
        
        int cnt = 0;

        int cubeWidth = (int)floor(sqrt(numParts));

        int xv = 800;
        int yv = Y_MAX-10;
        int zv = 0;
        int spacing = 20;
        while(cnt < numParts){
            particles.push_back(new Particle(xv, yv, zv, RADIUS, true, 1, sceneNum));
            
            
            if(xv < 1200){
                xv+=spacing;
            }
            else if(zv < 600){
                zv+=spacing;
                xv = 800;
            }
            else {
                zv = 0;
                yv-=spacing;
            }
            cnt++;
        }
    }
    //Bracket
    //Plate
    if(sceneNum == 2){
        //Spawn a spherical configuration of particles
        
        int cnt = 0;
        while(cnt < numParts){
            
            Vector3 newloc = {randRange(800, 1300), randRange(100, 600), randRange(-250, 250)};
            
            if(Vector3Distance(newloc, {1050, 350, 0}) < 350){
                particles.push_back(new Particle(newloc.x, newloc.y, newloc.z, RADIUS, true, 1, sceneNum));
                
                cnt++;
            }
        }

    }
    //None
    if(sceneNum == 3){
        //Spawn a spherical configuration of particles
        
        int cnt = 0;

        while(cnt < numParts){
            
            Vector3 newloc = {randRange(800, 1300), randRange(0, 500), randRange(-250, 250)};
            
            if(Vector3Distance(newloc, {1050, 250, 0}) < 250){
                particles.push_back(new Particle(newloc.x, newloc.y, newloc.z, RADIUS, true, 1, sceneNum));
                
                particles[cnt]->activateRandomMovement = true;
                cnt++;
            }
        }
        
        
    }//Coin drop
    if(sceneNum == 4){
        
        //Fill the coin drop chamber randomly
        int cnt = 0;

        while(cnt < numParts){
            
            Vector3 newloc = {randRange(0, 750), randRange(0, 1390), randRange(-25, 25)};
            particles.push_back(new Particle(newloc.x, newloc.y, newloc.z, RADIUS, true, 1, sceneNum));
            
            cnt++;
            
        }
        
    }//Ant tower
    if(sceneNum == 5){
        
        int cnt = 0;

        while(cnt < numParts){
            
            float h = 1400;
            float nY = randRange(0, h);
            float nX = randRange(950, 1050);
            float nZ = randRange(0, 100) - 50;

            if(cnt > 3000){
                nY = randRange(0, 600);
                nX = randRange(600, 1400);
                nZ = randRange(0, 800) - 400;
                if(distance(nX, nY, nZ, 1000, 200, 0) < 200){
                        particles.push_back(new Particle(nX, nY, nZ, RADIUS, true, 1, sceneNum));
                        
                        cnt++;
                }
            }
            else {
            
                particles.push_back(new Particle(nX, nY, nZ, RADIUS, true, 1, sceneNum));
                
                
                cnt++;
            }
        }
        
    }
    if(sceneNum == 6){
        //Spawn a spherical configuration of particles
        
        int cnt = 0;
        while(cnt < numParts){
            
            Vector3 newloc = {randRange(900, 1200), randRange(0, 700), randRange(-150, 150)};
            
            if(Vector3Distance(newloc, {1050, 0, 0}) < 650){
                particles.push_back(new Particle(newloc.x, newloc.y, newloc.z, RADIUS, true, 1, sceneNum));
                
                cnt++;
            }
        }
    }

    
    //Record the transform matrix of the particles for rendering
    for(int i = 0; i < numParts; i++){
        transforms[i] = MatrixTranslate((particles[i]->pos.x)*GLOBAL_RENDER_SCALE - (X_MAX-X_MIN)*GLOBAL_RENDER_SCALE/2.0f, (-particles[i]->pos.y)*GLOBAL_RENDER_SCALE + Y_MAX*GLOBAL_RENDER_SCALE, (particles[i]->pos.z)*GLOBAL_RENDER_SCALE - (Z_MAX-Z_MIN)*GLOBAL_RENDER_SCALE/2.0f);
            
    }
       
}

void resolveDistanceConstraints(vector<Particle*> &particles, int numParts, map<pair<int, int>, float> &restDistance, float balance, float intensity){
    
    //Wall bond resolution
    for(int i = 0; i < numParts; i++){
        if(particles[i]->wallBond && particles[i]->isWater){
            
            float restDist = 3;
            if(restDist <= 0 || isnan(particles[i]->wallBondX) || isnan(particles[i]->wallBondY) || isnan(particles[i]->wallBondZ)){
                continue;
            }

            float xa = particles[i]->pos.x;
            float xb = particles[i]->wallBondX;
            float ya = particles[i]->pos.y; 
            float yb = particles[i]->wallBondY;
            float za = particles[i]->pos.z; 
            float zb = particles[i]->wallBondZ;
            float normalDirX = (xa-xb);
            float normalDirY = (ya-yb);
            float normalDirZ = (za-zb);  
            float distSq = normalDirX*normalDirX + normalDirY*normalDirY + normalDirZ*normalDirZ;
            float dist = sqrt(distSq);

            //Allow bonds to break
            if(dist > globalSplitSetting*restDist*5 && globalSplitSetting > 0.01 && !particles[i]->setForce){
                particles[i]->wallBond = false;
                particles[i]->objBond = false;
               
                continue;
            }
            if(dist == 0){
                continue;
            }
            
            float constraintDistance = dist - restDist;
            normalDirX /= dist;
            normalDirY /= dist;
            normalDirZ /= dist;
            
            float w_idxA = 0;
            float w_idxB = 1;
            

            
            float localBondStrength = globalWallBondStrength;
            

            //Move particles towards adhesive bonds
            particles[i]->pos.x += -w_idxB / (w_idxA + w_idxB) * constraintDistance / 2 * normalDirX * localBondStrength;
            particles[i]->pos.y += -w_idxB / (w_idxA + w_idxB) * constraintDistance / 2 * normalDirY * localBondStrength;
            particles[i]->pos.z += -w_idxB / (w_idxA + w_idxB) * constraintDistance / 2 * normalDirZ * localBondStrength;


            if( isnan(particles[i]->pos.x) ||  isnan(particles[i]->pos.y) || isnan(particles[i]->pos.z)){
                cout << particles[i]->wallBondX << ", " << particles[i]->wallBondY << ", " << particles[i]->wallBondZ << endl;
                cout << "Error at ResolveDistConstraints (wallBonds)" << endl;
                exit(1);
            }

        }
    }
    
    if(intensity <= 0){
        return;
    }
    
    vector<int> pairsToErase;
    
    //Particle distance constraint bonds
    for(auto it = restDistance.begin(); it != restDistance.end(); it++){
        float restDist = it->second;
        int i = it->first.first;
        int j = it->first.second;
        //Mark bonds for deletion if necessary
        if(restDist <= 0){
            pairsToErase.push_back(i);
            pairsToErase.push_back(j);
            continue;
        }

        float xa = particles[i]->pos.x;
        float xb = particles[j]->pos.x;
        float ya = particles[i]->pos.y; 
        float yb = particles[j]->pos.y;
        float za = particles[i]->pos.z; 
        float zb = particles[j]->pos.z;
        float normalDirX = (xa-xb);
        float normalDirY = (ya-yb);  
        float normalDirZ = (za-zb);  
        float distSq = normalDirX*normalDirX + normalDirY*normalDirY + normalDirZ*normalDirZ;
        float dist = sqrt(distSq);

        //Allow bonds to break
        if((dist > globalSplitSetting*globalBondDistance && globalSplitSetting > 1.5) || (particles[i]->severedBonds || particles[j]->severedBonds)){
            it->second = -1;
            particles[i]->activeDistConstraints--;
            particles[j]->activeDistConstraints--;
            continue;
        }
        
        float constraintDistance = dist - restDist;
        normalDirX /= dist;
        normalDirY /= dist;
        normalDirZ /= dist;
        
        float w_idxA = balance;
        float w_idxB = 1-balance;
        
        //Apply distance constraint adjustments
        particles[j]->pos.x += w_idxA / (w_idxA + w_idxB) * constraintDistance / 2 * normalDirX * intensity;
        particles[j]->pos.y += w_idxA / (w_idxA + w_idxB) * constraintDistance / 2 * normalDirY * intensity;
        particles[j]->pos.z += w_idxA / (w_idxA + w_idxB) * constraintDistance / 2 * normalDirZ * intensity;
        particles[i]->pos.x += -w_idxB / (w_idxA + w_idxB) * constraintDistance / 2 * normalDirX * intensity;
        particles[i]->pos.y += -w_idxB / (w_idxA + w_idxB) * constraintDistance / 2 * normalDirY * intensity;
        particles[i]->pos.z += -w_idxB / (w_idxA + w_idxB) * constraintDistance / 2 * normalDirZ * intensity;

        if( isnan(particles[i]->pos.x) || isnan(particles[i]->pos.y) || isnan(particles[i]->pos.z)){
            cout << "Error at ResolveDistConstraints (particle bonds)" << endl;
            exit(1);
        }
        
    }

    //Delete marked distance constraints
    for(int i = 0; i < pairsToErase.size(); i+=2){
        int a = pairsToErase[i];
        int b = pairsToErase[i+1];
        restDistance.erase(make_pair(a, b));
    }
}

//Get bounding box info for mesh (used for collisions)
BoundingBox AdjustedBoundingBox(Mesh mesh, Vector3 position){
    BoundingBox meshBounds = GetMeshBoundingBox(mesh);
    meshBounds.min = Vector3Add(meshBounds.min, position);
    meshBounds.max = Vector3Add(meshBounds.max, position);
    meshBounds.min = Vector3Scale(meshBounds.min, 1/GLOBAL_RENDER_SCALE);
    meshBounds.max = Vector3Scale(meshBounds.max, 1/GLOBAL_RENDER_SCALE);
    return meshBounds;
}

//Used to record scenes in a text file for future rendering
void writeFrameToFile(string path, vector<Particle*> &particles, int numParts){
    
    ofstream UserFile(path.c_str(), ios::app);

    if(UserFile.is_open()){
        
        //Particle coordinates reduced to 2 decimal places to decrease file size
        UserFile << fixed << setprecision(2);
        for(int i = 0; i < numParts-1; i++){
            if(particles[i]->isWater){
                UserFile << i << "," << particles[i]->pos.x << "," << particles[i]->pos.y << "," << particles[i]->pos.z << ",";
            }
            else {
                UserFile << -i << "," << particles[i]->pos.x << "," << particles[i]->pos.y << "," << particles[i]->pos.z << ",";
            }
        }
        if(particles[numParts-1]->isWater){
            UserFile << numParts-1 << "," << particles[numParts-1]->pos.x << "," << particles[numParts-1]->pos.y << "," << particles[numParts-1]->pos.z << endl;
        }
        else {
            UserFile << -(numParts-1) << "," << particles[numParts-1]->pos.x << "," << particles[numParts-1]->pos.y << "," << particles[numParts-1]->pos.z << endl;
        }

        UserFile.close();
    }
    else {
        cout << "Failed to open recording file\n";
        exit(0);
    }
}
void writeMeshHeaderToFile(string path, string modelName){
    ofstream UserFile(path.c_str(), ios::app);

    if(UserFile.is_open()){
        
        //Write Location to File
        UserFile << modelName << "\n" << 1.0/GLOBAL_RENDER_SCALE << "\n";
        
        
        UserFile.close();
    }
    else {
        cout << "Failed to open recording file\n";
        exit(0);
    }
}
void writeMeshLocToFile(string path, Vector3 pos, bool scale){
    ofstream UserFile(path.c_str(), ios::app);

    if(UserFile.is_open()){
        
        //Write Location to File
        
        if(scale){
            //Convert to particle coords
            UserFile << pos.x / GLOBAL_RENDER_SCALE <<"," << pos.y / GLOBAL_RENDER_SCALE << "," << pos.z / GLOBAL_RENDER_SCALE << "\n";
        }
        else {
            UserFile << pos.x <<"," << pos.y << "," << pos.z << "\n";
        }
        
        UserFile.close();
    }
    else {
        cout << "Failed to open recording file\n";
        exit(0);
    }
}
void clearFrameFile(string path){
    //Delete file if exists
    remove(path.c_str());
    cout << "Deleted " << path << " if already exists" << endl;
    cout << "Writing to " << path << "\n";
}

int main(void){
	
	
    SetConfigFlags(FLAG_MSAA_4X_HINT);
	//Protection from screen tearing (may limit FPS)
	//SetConfigFlags(FLAG_VSYNC_HINT);
    SetTargetFPS(150);
	InitWindow(screenWidth, screenHeight, "Water Sim");

    //Recording variables:
    vector<float> benchmark_apply_forces_predict_pos;
    vector<float> benchmark_find_neighbors;
    vector<float> benchmark_apply_lambda;
    vector<float> benchmark_resolve_dist;
    vector<float> benchmark_wall_self;
    vector<float> benchmark_mesh;
    vector<float> benchmark_velocities;
    vector<float> benchmark_total;
    vector<float> frametime;
    vector<float> framerate;
    //Used for testing
    vector<float> raylibFramerate;
    bool startBenchmark = false;
    
    
    //Scene management
    int currentScene = 0;
    vector<string> sceneNames = {"Funnel", "Bracket", "Plate", "Dispersal", "Coin Drop", "Ant Tower", "Crushed Ants"};

    vector<Particle*> particles;
    
    //Distance constraints
    map<pair<int, int>, float> restDistance;
    
    //For rendering
    Matrix *transforms = (Matrix *)calloc(numParts*2, sizeof(Matrix));
    
    //Debugging
    Matrix *adhesionTransforms = (Matrix *)calloc(numParts*2, sizeof(Matrix));

    

    //Set up particle positions using scene information
    numParts = sceneMaxParticles[currentScene];
    initializeParticles(particles, numParts, currentScene, restDistance, transforms, true);

    
    
    //Spatial hash
    int tableSize = numParts * 10;
    SpatialKey *keyTable = new SpatialKey[numParts];
    int *startingIndexTable = new int[tableSize];

    //3D cursor size
    int handRadius = 60;
    
    //Camera setup
    Camera camera = { 0 };
    camera.position = { 0.0f, 5.5f, 20.0f };
    camera.target = { 0.0f, 5.5f, 0.0f };
    camera.up = { 0.0f, 1.0f, 0.0f };
    camera.fovy = 45.0f;
    camera.projection = CAMERA_PERSPECTIVE;

    //How many particles are in the bracket scene (increases slowly over time until max count is reached)
    int testFmarker = 0;

    Shader shader = LoadShader(TextFormat("lighting.vs"),
                               TextFormat("lighting.fs"));
    shader.locs[SHADER_LOC_MATRIX_MVP] = GetShaderLocation(shader, "mvp");
    shader.locs[SHADER_LOC_VECTOR_VIEW] = GetShaderLocation(shader, "viewPos");
	int ambientLoc = GetShaderLocation(shader, "ambient");

    Shader shader2 = LoadShader(TextFormat("lighting_instancing.vs"),
                               TextFormat("lighting.fs"));
    shader2.locs[SHADER_LOC_MATRIX_MVP] = GetShaderLocation(shader2, "mvp");
    shader2.locs[SHADER_LOC_VECTOR_VIEW] = GetShaderLocation(shader2, "viewPos");
    shader2.locs[SHADER_LOC_MATRIX_MODEL] = GetShaderLocationAttrib(shader2, "instanceTransform");

	int ambientLoc2 = GetShaderLocation(shader2, "ambient");
    float v[4] = { 0.4f, 0.4f, 0.4f, 1.0f };
    SetShaderValue(shader2, ambientLoc2, v, SHADER_UNIFORM_VEC4);
    SetShaderValue(shader, ambientLoc, v, SHADER_UNIFORM_VEC4);

    //Create light
	CreateLight(LIGHT_DIRECTIONAL, { 250.0f, 250.0f, 20.0f }, Vector3Zero(), WHITE, shader2);
	CreateLight(LIGHT_DIRECTIONAL, { 250.0f, 250.0f, 20.0f }, Vector3Zero(), WHITE, shader);
	

    //Particle material
	Material defaultMat = LoadMaterialDefault();
    defaultMat.maps[MATERIAL_MAP_DIFFUSE].color = BLUE;
    defaultMat.shader = shader2;

    //Solid object material
    Material defaultMat2 = LoadMaterialDefault();
    defaultMat2.maps[MATERIAL_MAP_DIFFUSE].color = GRAY;
    defaultMat2.shader = shader;

    Material redMat = LoadMaterialDefault();
    redMat.maps[MATERIAL_MAP_DIFFUSE].color = RED;
    redMat.shader = shader2;

    //Cursor mesh
    Mesh sphere = GenMeshSphere(RADIUS*GLOBAL_RENDER_SCALE, 8, 8);

    //If true, camera can be controlled
    bool viewerMode = false;

    //3D cursor
    float cursorX = -250;
    float cursorY = 0;
    float cursorZ = (Z_MAX+Z_MIN)/2.0f;

    //3D Coin
    Vector3 coinPos = {375, 1090, 0};
    Vector3 oCoinPos = coinPos;
    Vector3 coinVel = {0, 0, 0};

    //Load in funnel
    string funnelPath = "models/funnel.glb";
    Model funnel = LoadModel(funnelPath.c_str());
    Mesh funnelMesh = funnel.meshes[0];
    Vector3 funnelPos = {10.0f, 6.0f, 0.0f};
    BoundingBox funnelMeshBounds = AdjustedBoundingBox(funnelMesh, funnelPos);
    

    //Load in bracket
    string bracketPath = "models/bracket.glb";
    Model bracket = LoadModel(bracketPath.c_str());
    Mesh bracketMesh = bracket.meshes[0];
    Vector3 bracketPos = {10.0f, 1.f, 0.0f};
    BoundingBox bracketMeshBounds = AdjustedBoundingBox(bracketMesh, bracketPos);

    
    //Load in plate
    string platePath = "models/plate.glb";
    Model plate = LoadModel(platePath.c_str());
    Mesh plateMesh = plate.meshes[0];
    Vector3 platePos = {10.0f, 12.f, 0.0f};
    BoundingBox plateMeshBounds = AdjustedBoundingBox(plateMesh, platePos);

    //Load in tower
    string towerPath = "models/tower.glb";
    Model tower = LoadModel(towerPath.c_str());
    Mesh towerMesh = tower.meshes[0];
    Vector3 towerPos = {10.0f, 10.f, 0.0f};
    BoundingBox towerMeshBounds = AdjustedBoundingBox(towerMesh, towerPos);
    
    //Load in coin
    string coinPath = "models/coin.glb";
    Model coin = LoadModel(coinPath.c_str());
    Mesh coinMesh = coin.meshes[0];

    //Testing Box
    Mesh cubeTest = GenMeshCube(5, 5, 5);
    Model cubeTestModel;
    cubeTestModel = LoadModelFromMesh(cubeTest);
    Vector3 cubeTestPos = {10, 10, 0};
    BoundingBox cubeTestBounds = AdjustedBoundingBox(cubeTest, cubeTestPos);

    //Scene Recording Tools
    bool isRecording = false;
    int currentFrame = 0;

	
	while(!WindowShouldClose()){
        //Record start time of frame
        auto fullFrameStart = std::chrono::high_resolution_clock::now();


        currentFrame++;
        //If viewmode active the user can control the camera
        if(viewerMode && !IsKeyDown(KEY_Q) && !IsKeyDown(KEY_E)){
            UpdateCamera(&camera, CAMERA_THIRD_PERSON);
        }
        //Switch between viewing modes
        if(IsKeyPressed(KEY_LEFT_SHIFT)){
            if(!viewerMode){
                viewerMode = true;
            }
            else {
                viewerMode = false;
            }
        }
        
        //Reset the buckets
        for(int i = 0; i < tableSize; i++){
            startingIndexTable[i] = -1;
        }
       
        //Adjust cursor size
        if(IsKeyPressed(49)){
            handRadius = 5;
        }
        if(IsKeyPressed(50)){
            handRadius = 100;
        }
        if(IsKeyPressed(51)){
            handRadius = 200;
        }
        //Change scene
        if(IsKeyPressed(KEY_FOUR)){
            currentScene+=1;
            if(currentScene >= sceneNames.size()){
                currentScene = 0;
                testFmarker = 0;
                coinPos.x = 375;
                coinPos.y = 1090;
                coinPos.z = 0;
                coinVel.y = 0;
            }
            numParts = sceneMaxParticles[currentScene];
            initializeParticles(particles, numParts, currentScene, restDistance, transforms, true);
            
            
    
        }
        if(IsKeyPressed(KEY_EIGHT)){
            startBenchmark = true;
            cout << "Starting Benchmark Recording" << endl;
        }
        
        //Change environment size
        if(IsKeyDown(76)){
            X_MIN -= 570 * GetFrameTime();
            X_MAX += 570*GetFrameTime();
        }
        if(IsKeyDown(74)){
            X_MIN += 570 * GetFrameTime();
            X_MAX -= 570*GetFrameTime();
        }
        if(IsKeyDown(75)){
            Y_MIN -= 570 * GetFrameTime();
            Y_MAX += 570*GetFrameTime();
        }
        if(IsKeyDown(73)){
            Y_MIN += 570 * GetFrameTime();
            Y_MAX -= 570*GetFrameTime();
        }
        if(IsKeyDown(KEY_U)){
            Z_MIN -= 570 * GetFrameTime();
            Z_MAX += 570*GetFrameTime();
        }
        if(IsKeyDown(KEY_O)){
            Z_MIN += 570 * GetFrameTime();
            Z_MAX -= 570*GetFrameTime();
        }
        //Allow scene to be reset
        if(IsKeyPressed(KEY_FIVE)){
            testFmarker = 0;
            initializeParticles(particles, numParts, currentScene, restDistance, transforms, false);
            if(currentScene == 4){
                coinPos.x = 375;
                coinPos.y = 1330;
                coinPos.z = 0;
                coinVel.y = 0;
            }
        }
        //Print data to console
        if(IsKeyPressed(KEY_M)){
            for(int i = 0; i < numParts; i++){
                cout << particles[i]->pos.x << ", " << particles[i]->pos.y << ", "<< particles[i]->pos.z << endl;
            }
            cout << endl;
        }
        if(IsKeyPressed(KEY_MINUS)){
            for(int i = 0; i < decisionRegionCount; i++){
                cout << "Gap " << i << ": " << holeCounts[i] << "\n";
            }
        }

        //Scene 2 - particles pour in over time instead of spawning instantly
        
        if(testFmarker < numParts && currentScene == 1 && currentFrame % 55 == 0){
            
            for(int i = 0; i < 5; i++){
                for(int j = 0; j < 5; j++){
                    for(int k = 0; k < 5; k++){
                        particles.push_back(new Particle(800+15*i, 2000+20*j, -240+15*k, RADIUS, true, 1, 1));
                        
                        particles[testFmarker]->setForce = true;
                        if(rand() % 6 == 0){
                            particles[testFmarker]->activateRandomMovement = true;
                        }

                        testFmarker++;
                    }

                }
            }
            
        }



        

        //Start recording
        if(IsKeyPressed(KEY_SIX)){
            currentFrame = 1;
            isRecording = true;
            clearFrameFile("frames.txt");
            clearFrameFile("mesh.txt");

            switch(currentScene){
                case 0:
                    writeMeshHeaderToFile("mesh.txt", funnelPath);
                    writeMeshLocToFile("mesh.txt", funnelPos, true);
                break;
                case 1:
                    writeMeshHeaderToFile("mesh.txt", bracketPath);
                    writeMeshLocToFile("mesh.txt", bracketPos, true);
                break;
                case 2:
                    writeMeshHeaderToFile("mesh.txt", platePath);
                break;
                case 4:
                    writeMeshHeaderToFile("mesh.txt", coinPath);
                break;
                case 5:
                    writeMeshHeaderToFile("mesh.txt", towerPath);
                    writeMeshLocToFile("mesh.txt", towerPos, true);
                break;
            }
            
            

        }
        //Cancel Recording Operation Entirely
        if(IsKeyPressed(KEY_SEVEN) && isRecording){
            isRecording = false;
            cout << "Done writing to frames.txt, " << currentFrame << " frames recorded" << endl;
            currentFrame = 1;
            
            
        }

        //Move 3D cursor with WASD/QE
        if(currentScene != 2 && currentScene != 6){
            if(IsKeyDown(KEY_A) && !viewerMode){
                cursorX -= 400.0f * GetFrameTime();
            }
            if(IsKeyDown(KEY_D) && !viewerMode){
                cursorX += 400.0f * GetFrameTime();
            }
            if(IsKeyDown(KEY_W) && !viewerMode){
                cursorZ -= 400.0f * GetFrameTime();
            }
            if(IsKeyDown(KEY_S) && !viewerMode){
                cursorZ += 400.0f * GetFrameTime();
            }
            if(IsKeyDown(KEY_Q) && !viewerMode){
                cursorY -= 400.0f * GetFrameTime();
            }
            if(IsKeyDown(KEY_E) && !viewerMode){
                cursorY += 400.0f * GetFrameTime();
            }
        }
        else {
            //Control plate that crushes particles
            if(IsKeyDown(KEY_Q) && !viewerMode){
                platePos.y -= 0.02f;
                for(int i = 0; i < particles.size(); i++){
                    if(particles[i]->objBond){
                        
                        particles[i]->wallBondY -= 0.02f / GLOBAL_RENDER_SCALE;
                        
                    }
                }
                
            }
            if(IsKeyDown(KEY_E) && !viewerMode){
                platePos.y += 0.02f;
                for(int i = 0; i < particles.size(); i++){
                    if(particles[i]->objBond){
                        
                        particles[i]->wallBondY += 0.02 / GLOBAL_RENDER_SCALE;
                        
                    }
                }
            } 
        }

        

		BeginDrawing();
        ClearBackground({230, 230, 230});
        
        BeginMode3D(camera);
            DrawGrid(100, 1.0f);
            //Draw particles as spheres
            DrawMeshInstanced(sphere, defaultMat, transforms, particles.size()-1);

            
            //Visualize adhesion bonds
            //DrawMeshInstanced(sphere, redMat, adhesionTransforms, particles.size()-1);

        
            //Visualize distance constraints
            /*
            for(auto it = restDistance.begin(); it != restDistance.end(); it++){
                float restDist = it->second;
                int i = it->first.first;
                int j = it->first.second;
                DrawLine3D(Vector3Scale({particles[i]->pos.x, particles[i]->pos.y, particles[i]->pos.z}, GLOBAL_RENDER_SCALE), Vector3Scale({particles[j]->pos.x, particles[j]->pos.y, particles[j]->pos.z}, GLOBAL_RENDER_SCALE), BLACK);
            }
            */
            

            //Visualize Decision Regions
            if(currentScene == 1){
                for(int i = 0; i < decisionRegionCount; i++){
                    BoundingBox b;
                    b.min = {holeXThreshold * GLOBAL_RENDER_SCALE, decisionRegions[i].x * GLOBAL_RENDER_SCALE, -800 * GLOBAL_RENDER_SCALE};
                    b.max = {(holeXThreshold+1400)*GLOBAL_RENDER_SCALE, decisionRegions[i].y*GLOBAL_RENDER_SCALE, 10*GLOBAL_RENDER_SCALE};
                    DrawBoundingBox(b, RED);
                }
            }

            //Environment boundaries
            DrawCubeWires({(X_MAX+X_MIN)*GLOBAL_RENDER_SCALE / 2.0f, (Y_MAX-Y_MIN)*GLOBAL_RENDER_SCALE/2.0f, (Z_MAX+Z_MIN)/2.0f*GLOBAL_RENDER_SCALE}, (X_MAX-X_MIN)*GLOBAL_RENDER_SCALE, (Y_MAX-Y_MIN)*GLOBAL_RENDER_SCALE, (Z_MAX-Z_MIN)*GLOBAL_RENDER_SCALE, BLACK);
            
            //Draw either cursor or coin
            if(currentScene != 2 && currentScene != 4)
                DrawSphere({cursorX * GLOBAL_RENDER_SCALE, cursorY * GLOBAL_RENDER_SCALE, cursorZ*GLOBAL_RENDER_SCALE}, handRadius*GLOBAL_RENDER_SCALE, {255, 100, 0, 200});
            if(currentScene == 4){
                DrawMesh(coinMesh, defaultMat2, MatrixTranslate(coinPos.x * GLOBAL_RENDER_SCALE, coinPos.y * GLOBAL_RENDER_SCALE, coinPos.z * GLOBAL_RENDER_SCALE));
            }

            //Draw solid object models (funnel, bracket, etc.)
            if(currentScene == 0){
                DrawMesh(funnelMesh, defaultMat2, MatrixTranslate(funnelPos.x, funnelPos.y, funnelPos.z));
                DrawModelWires(funnel, funnelPos, 1, BLACK);
            }
            else if(currentScene == 1){
                DrawMesh(bracketMesh, defaultMat2, MatrixTranslate(bracketPos.x, bracketPos.y, bracketPos.z));
                DrawModelWires(bracket, bracketPos, 1, BLACK);
            }
            else if(currentScene == 2 || currentScene == 6){
                DrawMesh(plateMesh, defaultMat2, MatrixTranslate(platePos.x, platePos.y, platePos.z));
                DrawModelWires(plate, platePos, 1, BLACK);
            }
            else if(currentScene == 5){
                DrawMesh(towerMesh, defaultMat2, MatrixTranslate(towerPos.x, towerPos.y, towerPos.z));
            }
        
        //Move coin in coin drop scene
        if(currentScene == 4){
            oCoinPos = coinPos;
            coinVel.y -= gravity*timestep;
            coinPos.y += coinVel.y*timestep;
            coinPos.x += coinVel.x*timestep;
            coinPos.z += coinVel.z*timestep;
            
        }
        
        
        //Simulation steps
        auto start = std::chrono::high_resolution_clock::now();

        for(int i = 0; i < particles.size(); i++){
            const float oldX = particles[i]->pos.x;
            const float oldY = particles[i]->pos.y;
            const float oldZ = particles[i]->pos.z;
            particles[i]->opos = {oldX, oldY, oldZ};
            
            

            particles[i]->applyForces();
            particles[i]->predictPosition();
            particles[i]->viscosityMu = globalViscSetting;
            
            int bX = floor(particles[i]->pos.x / BUCKET_WIDTH);
            int bY = floor(particles[i]->pos.y / BUCKET_WIDTH);
            int bZ = floor(particles[i]->pos.z / BUCKET_WIDTH);
            
            //Get bucket
            int bucket = hashCoords(bX, bY, bZ) % tableSize;
            
            //Track index of particle and associated bucket
            keyTable[i].bucket = bucket;
            keyTable[i].index = i;

            
            
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = end - start;
        
        if(startBenchmark){
            benchmark_apply_forces_predict_pos.push_back(elapsed.count());
            benchmark_total.push_back(elapsed.count());
        }

        /*
        Sort keytable by bucket values so that particles of similar buckets end up next to one another
        See: https://youtu.be/rSKMYc1CQHE?t=1421
        
        Other info on spatial hashing:
        https://github.com/matthias-research/pages/blob/master/tenMinutePhysics/11-hashing.pdf
        */
        qsort((void*) keyTable, numParts, sizeof(SpatialKey), compareKeys);

        //Mark the starting index of each bucket
        for(int i = 0; i < numParts; i++){
            int currentBucket = keyTable[i].bucket;
            if(startingIndexTable[currentBucket] > i || startingIndexTable[currentBucket] == -1){
                startingIndexTable[currentBucket] = i;
            }
        }
        start = std::chrono::high_resolution_clock::now();
        //Find adjacent particles
        #pragma omp parallel for
        for(int i = 0; i < particles.size(); i++){
            if(particles[i]->isWater){
                particles[i]->findNeighbors(i, keyTable, startingIndexTable, tableSize, particles, restDistance);
                
            }
        }
        
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        
        if(startBenchmark){
            benchmark_find_neighbors.push_back(elapsed.count());
            benchmark_total[benchmark_total.size() - 1] += elapsed.count();
        }
        
        int iter = 0;
        while(iter < solverIter){
            
            start = std::chrono::high_resolution_clock::now();

            #pragma omp parallel for
            for(int i = 0; i < particles.size(); i++){
            //Calc lambda
                if(particles[i]->isWater){
                    particles[i]->applyLambdaPositionUpdate(particles);
                }
            }

            end = std::chrono::high_resolution_clock::now();
            elapsed = end - start;
            
            if(startBenchmark){
                if(iter == 0){
                    benchmark_apply_lambda.push_back(elapsed.count());
                }
                else {
                    benchmark_apply_lambda[benchmark_apply_lambda.size()-1] += elapsed.count();
                }
                benchmark_total[benchmark_total.size() - 1] += elapsed.count();
            }

            start = std::chrono::high_resolution_clock::now();
            //Dist constraints / Adhesive bonds
            resolveDistanceConstraints(particles, numParts, restDistance, distanceConstraintBalance, distanceConstraintIntensity);
        
            end = std::chrono::high_resolution_clock::now();
            elapsed = end - start;
            
            if(startBenchmark){
                if(iter == 0){
                    benchmark_resolve_dist.push_back(elapsed.count());
                }
                else {
                    benchmark_resolve_dist[benchmark_resolve_dist.size()-1] += elapsed.count();
                }
                benchmark_total[benchmark_total.size() - 1] += elapsed.count();
            }
            start = std::chrono::high_resolution_clock::now();
            //#pragma omp parallel for
            for(int i = 0; i < particles.size(); i++){
                
                particles[i]->unboundedPerformCollisions(i, keyTable, startingIndexTable, tableSize, particles);
                if(currentScene != 4){
                    particles[i]->ballCollisions(cursorX, cursorY, cursorZ, handRadius);
                }
                
                
                particles[i]->keepInWalls();
                
            }

            end = std::chrono::high_resolution_clock::now();
            elapsed = end - start;
            
            if(startBenchmark){
                if(iter == 0){
                    benchmark_wall_self.push_back(elapsed.count());
                }
                else {
                    benchmark_wall_self[benchmark_wall_self.size()-1] += elapsed.count();
                }
                benchmark_total[benchmark_total.size() - 1] += elapsed.count();
            }
            
            iter++;
        }

        start = std::chrono::high_resolution_clock::now();

        //Solid object collisions
        if(currentScene == 0){
            #pragma omp parallel for
            for(int i = 0; i < particles.size(); i++){
                ApplyMeshCollision(&funnelMesh, particles[i], funnelPos, false, funnelMeshBounds);
            }
        }
        else if(currentScene == 1){
            #pragma omp parallel for
            for(int i = 0; i < particles.size(); i++){
                ApplyMeshCollision(&bracketMesh, particles[i], bracketPos, false, bracketMeshBounds);
            }

            for(int i = 0; i < particles.size(); i++){

                if(!particles[i]->isWater){
                    //Store original surface particle is attracted to
                    particles[i]->oldSurface.x = particles[i]->closestSurface.x;
                    particles[i]->oldSurface.y = particles[i]->closestSurface.y;
                    particles[i]->oldSurface.z = particles[i]->closestSurface.z;

                    //If the object is not connected to the mesh find the closest available triangle to move towards
                    if(!isObjectOnMesh(&bracketMesh, particles[i]->pos, bracketPos)){
                        particles[i]->closestSurface = findClosestTrianglePointFromOutside(&bracketMesh, particles[i]->pos, particles[i]->opos, bracketPos);
                    }
                
                }
            }
        }
        else if(currentScene == 2 || currentScene == 6){
            plateMeshBounds = AdjustedBoundingBox(plateMesh, platePos);
            
            #pragma omp parallel for
            for(int i = 0; i < particles.size(); i++){
                ApplyMeshCollision(&plateMesh, particles[i], platePos, false, plateMeshBounds);
            }
        }
        else if(currentScene == 4){
            
            #pragma omp parallel for
            for(int i = 0; i < particles.size(); i++){
                ApplyCoinCollision(particles[i], coinPos, 115, 50);
            }
        }
        else if(currentScene == 5){
            #pragma omp parallel for
            for(int i = 0; i < particles.size(); i++){
                ApplyMeshCollision(&towerMesh, particles[i], towerPos, false, towerMeshBounds);
            }
        }
        


        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        
        if(startBenchmark){
            benchmark_mesh.push_back(elapsed.count());
            benchmark_total[benchmark_total.size() - 1] += elapsed.count();
        }

        for(int i = 0; i < decisionRegionCount; i++){
            holeCounts[i] = 0;
        }

        start = std::chrono::high_resolution_clock::now();


        #pragma omp parallel for
        for(int i = 0; i < particles.size(); i++){
            //Calc new velocities
            const double rX = particles[i]->pos.x - particles[i]->opos.x;
            const double rY = particles[i]->pos.y - particles[i]->opos.y;
            const double rZ = particles[i]->pos.z - particles[i]->opos.z;
            

            float drag = 0;

            //If friction/drag flags are active, velocity will be dampened more severely
            if(particles[i]->activateFriction){
                drag = globalFrictionSetting;
            }
            if(particles[i]->activateFloorFriction && globalFrictionSetting > 0){
                //Floor friction of different scenes
                drag = floorFriction[currentScene];


                //Special plate rule
                if(currentScene == 2){
                    particles[i]->pos = particles[i]->opos;
                }
            }
            //Special plate rule
            if(particles[i]->wallBond && (currentScene == 2 || currentScene == 6)){
                drag = 0.998;
            }
            
            if(particles[i]->particleFriction){
                particles[i]->particleFriction = false;
            }

            //Set particle velocities
            particles[i]->setVel(constrain(rX/timestep*(0.999 - drag), -cRange, cRange), constrain(rY/timestep*(0.999 - drag), -cRange, cRange), constrain(rZ/timestep*(0.999 - drag), -cRange, cRange));
            
            //Record matrix data for particles (for instanced rendering)
            transforms[i] = MatrixTranslate((particles[i]->pos.x)*GLOBAL_RENDER_SCALE, (particles[i]->pos.y)*GLOBAL_RENDER_SCALE, (particles[i]->pos.z)*GLOBAL_RENDER_SCALE);
            
            if(particles[i]->wallBond){
                adhesionTransforms[i] = MatrixTranslate(particles[i]->wallBondX * GLOBAL_RENDER_SCALE, particles[i]->wallBondY * GLOBAL_RENDER_SCALE, particles[i]->wallBondZ * GLOBAL_RENDER_SCALE);

            }
            else {
                adhesionTransforms[i] = MatrixTranslate(-10000 * GLOBAL_RENDER_SCALE, -10000 * GLOBAL_RENDER_SCALE, -10000 * GLOBAL_RENDER_SCALE);

            }


            //Bracket scenario
            //Keep track of how many particles filled each gap
            if(currentScene == 1){
                for(int j = 0; j < 3; j++){
                    if(particles[i]->pos.x > holeXThreshold && particles[i]->pos.y > decisionRegions[j].x && particles[i]->pos.y < decisionRegions[j].y){
                        holeCounts[j]++;
                    }
                }
            }
            
        }
        if(currentScene == 4){
            //Keep coin in wall boundaries
            if(coinPos.y < 0){
                coinPos.y = 0;
                coinVel.y = 0;
            }
            if(coinPos.y + handRadius / 2.0 > Y_MAX){
                coinPos.y += (Y_MAX - handRadius/2.0 - coinPos.y);
                
            }
            if(coinPos.y - handRadius / 2.0 < Y_MIN){
                coinPos.y += (Y_MIN- (coinPos.y-handRadius/2.0));
            }
            
            if(coinPos.x+handRadius/2.0  > X_MAX){
                coinPos.x += (X_MAX - (coinPos.x+handRadius));
                
            }
            if(coinPos.x-handRadius  < X_MIN){
                coinPos.x += (X_MIN - (coinPos.x-handRadius));
            }
            if(coinPos.z+handRadius  > Z_MAX){
                coinPos.z = 0;
            }
            if(coinPos.z-handRadius < Z_MIN){
                coinPos.z = 0;
            }

            //Compute coin velocity like particles
            Vector3 delta = Vector3Subtract(coinPos, oCoinPos);
            
            coinVel.x = constrain(delta.x/timestep*(0.9999), -cRange, cRange);
            coinVel.y = constrain(delta.y/timestep*(0.9999), -cRange, cRange);
            coinVel.z = constrain(delta.z/timestep*(0.9999), -cRange, cRange);
            
        }

        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        
        if(startBenchmark){
            benchmark_velocities.push_back(elapsed.count());
            benchmark_total[benchmark_total.size() - 1] += elapsed.count();
        }

        //Recording frame data
        if(isRecording){
            if(currentFrame < FRAME_DATA_LIMIT){
                if(currentScene != 1){
                writeFrameToFile("frames.txt", particles, numParts);
                }
                else {
                    writeFrameToFile("frames.txt", particles, testFmarker);
                }
                if(currentScene == 2 || currentScene == 6){
                    writeMeshLocToFile("mesh.txt", platePos, true);
                }
                else if(currentScene == 4){
                    writeMeshLocToFile("mesh.txt", coinPos, false);
                }
                
                
            }
            else {
                isRecording = false;
                cout << "Done writing to frames.txt, " << currentFrame << " frames recorded" << endl;
                currentFrame = 1;

            }
        }

        EndMode3D();
        
		//Directions/Information Text
		DrawFPS(10, 20);
        DrawText("IJKLOU to adjust environment size\n4 to Cycle Through Scenes\n5 to reset simulation\n6 to START recording\n7 to STOP recording\n8/9 to START/STOP recording FPS data\n\nRecordings saved to \"frames.txt\" + \"mesh.txt\"", 10, 380, 20, DARKGREEN);
        DrawText(TextFormat("Total Particles: %d", numParts), 10, 60, 20, DARKGREEN);
        DrawText(TextFormat("Scene: %s", sceneNames[currentScene].c_str()), screenWidth/2.0f - 100, 30, 20, DARKGREEN);
        if(!viewerMode){
            DrawText("Controlling: 3D CURSOR\nWASD + QE to move, 1-3 to change size\n\nSHIFT to switch to CAMERA", 10, 100, 20, DARKGREEN);
        }
        else {
            DrawText("Controlling: CAMERA\nWASD + Mouse to Move, Scroll to Zoom\n\nSHIFT to switch to 3D CURSOR", 10, 100, 20, DARKGREEN);
        }

        if(isRecording){
            DrawText(TextFormat("RECORDING - Frame: %d", currentFrame), 10, 360, 20, DARKGREEN);
        }
        
        
        EndDrawing();

        
        
        

        auto fullFrameEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> fullTimeElapsed = fullFrameEnd - fullFrameStart;

        if(startBenchmark){
            frametime.push_back(fullTimeElapsed.count());
            framerate.push_back(1000.0f / fullTimeElapsed.count());
            raylibFramerate.push_back(GetFPS());
        }

        //Write benchmark time data to csv
        if(IsKeyPressed(KEY_NINE)){
            startBenchmark = false;
            cout << "Ending Benchmark Recording" << endl;
            string csvPath = sceneNames[currentScene];
            
            csvPath.append(".csv");
            for(int i = 0; i < csvPath.length(); i++){
                if(csvPath[i] == ' '){
                    csvPath[i] = '_';
                }
                csvPath[i] = tolower(csvPath[i]);
            }
            remove(csvPath.c_str());
            ofstream UserFile(csvPath.c_str(), ios::app);
            
            if(UserFile.is_open()){
                
                //Particle coordinates reduced to 2 decimal places to decrease file size
                UserFile << fixed << setprecision(2);
                UserFile << "scene,force-predict,calculate-density-gradients,apply-density-constraint,resolve-distance-constraints,ant-ant+ant-wall-collisions,mesh-collisions,calculate-velocities,simulation-time,total-frame-time,frames-per-second" << endl;
                for(int i = 0; i < benchmark_apply_forces_predict_pos.size(); i++){
                    UserFile << sceneNames[currentScene] << ","
                    << benchmark_apply_forces_predict_pos[i] << "," 
                    << benchmark_find_neighbors[i] << ","
                    << benchmark_apply_lambda[i] << ","
                    << benchmark_resolve_dist[i] << ","
                    << benchmark_wall_self[i] << ","
                    << benchmark_mesh[i] << ","
                    << benchmark_velocities[i] << ","
                    << benchmark_total[i] << ","
                    << frametime[i] << ","
                    << framerate[i]
                    << endl;
                }


                UserFile.close();
            }
            else {
                cout << "Failed to open recording file\n";
                exit(0);
            }
        }
        
        
        
        
        
	}
	
    //Clear Objects
    free(transforms);
    free(adhesionTransforms);
    for(int i = 0; i < particles.size(); i++){
        delete particles[i];
    }
    delete[] keyTable;
    delete[] startingIndexTable;

	return 0;
}
