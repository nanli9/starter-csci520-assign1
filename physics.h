/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#ifndef _PHYSICS_H_
#define _PHYSICS_H_

void computeAcceleration(struct world * jello, struct point a[8][8][8]);

// perform one step of Euler and Runge-Kutta-4th-order integrators
// updates the jello structure accordingly
void Euler(struct world * jello);
void RK4(struct world * jello);
point findNetForce(struct world* jello, int i, int j, int k);
double distanceBetweenPoints(struct world* jello,int x1, int y1, int z1, int x2, int y2, int z2);
point getL(struct world* jello, int x1, int y1, int z1, int x2, int y2, int z2);
point calculateDamping(point Va, point Vb, double dElastic, point L);
void nomralize(point &p);
double dotProduct(point a, point b);
point findExternalForce(struct world* jello, int i, int j, int k);
point getGrid(struct world* jello, int i,int j,int k);
point trilinearInterpolation(point grid, struct world* jello);
point collisionHandler(struct world* jello,double a,double b,double c,double d,int i,int j,int k);
point userInputForce(struct world* jello);
#endif

