/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"

/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(struct world * jello, struct point a[8][8][8])
{
  /* for you to implement ... */
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            for (int k = 0; k < 8; k++) {
                point netForce = findNetForce(jello, i, j, k);
                a[i][j][k].x = netForce.x/jello->mass;
                a[i][j][k].y = netForce.y/jello->mass;
                a[i][j][k].z = netForce.z/jello->mass;
            }
        }
    }

}
point findNetForce(struct world* jello,int i,int j,int k) {
    //initialize the forces
    point structuralForce = { 0.0, 0.0, 0.0 };
    point shearForce = { 0.0, 0.0, 0.0 };
    point bendForce = { 0.0, 0.0, 0.0 };
    point netForce = { 0.0, 0.0, 0.0 };
    point dampForce = { 0.0, 0.0, 0.0 };
    point externalForce = { 0.0, 0.0, 0.0 };
    point collisionForce = { 0.0, 0.0, 0.0 };
    //go over the structural springs
    for (int x = 0; x < jello->neighborsInfo[i][j][k].structuralNeighborList.size(); x++) {
        int i2, j2, k2;
        i2 = jello->neighborsInfo[i][j][k].structuralNeighborList[x].p.x;
        j2 = jello->neighborsInfo[i][j][k].structuralNeighborList[x].p.y;
        k2 = jello->neighborsInfo[i][j][k].structuralNeighborList[x].p.z;
        double restLength = jello->neighborsInfo[i][j][k].structuralNeighborList[x].length;
        double distance = distanceBetweenPoints(jello, i, j, k, i2, j2, k2);
        point L = getL(jello, i, j, k, i2, j2, k2);
        nomralize(L);
        structuralForce.x += -jello->kElastic * (distance - restLength)* L.x;
        structuralForce.y += -jello->kElastic * (distance - restLength)* L.y;
        structuralForce.z += -jello->kElastic * (distance - restLength)* L.z;
        dampForce = calculateDamping(jello->v[i][j][k], jello->v[i2][j2][k2],jello->dElastic, L);
        pSUM(structuralForce, dampForce, structuralForce);
    };
    //go over the bend springs
    for (int x = 0; x < jello->neighborsInfo[i][j][k].bendNeighborList.size(); x++) {
        int i2, j2, k2;
        i2 = jello->neighborsInfo[i][j][k].bendNeighborList[x].p.x;
        j2 = jello->neighborsInfo[i][j][k].bendNeighborList[x].p.y;
        k2 = jello->neighborsInfo[i][j][k].bendNeighborList[x].p.z;
        double restLength = jello->neighborsInfo[i][j][k].bendNeighborList[x].length;
        double distance = distanceBetweenPoints(jello, i, j, k, i2, j2, k2);
        point L = getL(jello, i, j, k, i2, j2, k2);
        nomralize(L);
        bendForce.x += -jello->kElastic * (distance - restLength) * L.x;
        bendForce.y += -jello->kElastic * (distance - restLength) * L.y;
        bendForce.z += -jello->kElastic * (distance - restLength) * L.z;
        dampForce = calculateDamping(jello->v[i][j][k], jello->v[i2][j2][k2], jello->dElastic, L);
        pSUM(bendForce, dampForce, bendForce);
    };
    //go over the shear springs
    for (int x = 0; x < jello->neighborsInfo[i][j][k].shearNeighborList.size(); x++) {
        int i2, j2, k2;
        i2 = jello->neighborsInfo[i][j][k].shearNeighborList[x].p.x;
        j2 = jello->neighborsInfo[i][j][k].shearNeighborList[x].p.y;
        k2 = jello->neighborsInfo[i][j][k].shearNeighborList[x].p.z;
        double restLength = jello->neighborsInfo[i][j][k].shearNeighborList[x].length;
        double distance = distanceBetweenPoints(jello, i, j, k, i2, j2, k2);
        point L = getL(jello, i, j, k, i2, j2, k2);
        nomralize(L);
        shearForce.x += -jello->kElastic * (distance - restLength) * L.x;
        shearForce.y += -jello->kElastic * (distance - restLength) * L.y;
        shearForce.z += -jello->kElastic * (distance - restLength) * L.z;
        dampForce = calculateDamping(jello->v[i][j][k], jello->v[i2][j2][k2], jello->dElastic, L);
        pSUM(shearForce, dampForce, shearForce);
    };
    externalForce = findExternalForce(jello,i,j,k);
    //collision detection for six plane top bottom left right far near
    collisionForce = collisionForce+collisionHandler(jello,0.0,0.0,-1.0,2.0,i,j,k);
    collisionForce = collisionForce + collisionHandler(jello,0.0,0.0,1.0,2.0,i,j,k);
    collisionForce = collisionForce + collisionHandler(jello,0.0,-1.0,0.0,2.0,i,j,k);
    collisionForce = collisionForce + collisionHandler(jello,0.0,1.0,0.0,2.0,i,j,k);
    collisionForce = collisionForce + collisionHandler(jello,-1.0,0.0,0.0,2.0,i,j,k);
    collisionForce = collisionForce + collisionHandler(jello,1.0,0.0,0.0,2.0,i,j,k);
    
    //collision check for incline plane
    if(jello->incPlanePresent==1)
        collisionForce = collisionForce + collisionHandler(jello, jello->a, jello->b, jello->c, jello->d, i, j, k);
    //add up the forces without user input
    netForce.x = structuralForce.x + bendForce.x + shearForce.x + externalForce.x + collisionForce.x;
    netForce.y = structuralForce.y + bendForce.y+ shearForce.y+ externalForce.y + collisionForce.y;
    netForce.z = structuralForce.z + bendForce.z+ shearForce.z+ externalForce.z + collisionForce.z;
    //add user input forces for single points or all points
    if (mark)
    {
        if (i == pickedPoint[0] && j == pickedPoint[1]&& k == pickedPoint[2])
        {
            netForce.x = netForce.x + jello->userInputForce.x;
            netForce.y = netForce.y + jello->userInputForce.y;
        }
    }
    else
    {
        if (jello->userInputForce.x != 0 || jello->userInputForce.y != 0) {
            netForce.x = netForce.x + jello->userInputForce.x;
            netForce.y = netForce.y + jello->userInputForce.y;
        }
    }
    return netForce;
}
double distanceBetweenPoints(struct world* jello,int x1, int y1, int z1,int x2,int y2,int z2) {
    double x, y, z;
    x = pow((jello->p[x1][y1][z1].x - jello->p[x2][y2][z2].x),2);
    y = pow((jello->p[x1][y1][z1].y - jello->p[x2][y2][z2].y),2);
    z = pow((jello->p[x1][y1][z1].z - jello->p[x2][y2][z2].z),2);
    return sqrt(x+y+z);
}
point getL(struct world* jello, int x1, int y1, int z1, int x2, int y2, int z2) {
    double x, y, z;
    x = jello->p[x1][y1][z1].x - jello->p[x2][y2][z2].x;
    y = jello->p[x1][y1][z1].y - jello->p[x2][y2][z2].y;
    z = jello->p[x1][y1][z1].z - jello->p[x2][y2][z2].z;
    return { x,y,z };
}
void nomralize(point &p) {
    double length = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
    p.x /= length;
    p.y /= length;
    p.z /= length;
}
point calculateDamping(point Va,point Vb,double dElastic,point L) {
    point result = { 0.0, 0.0, 0.0 };
    result.x = -dElastic * (Va.x - Vb.x) * L.x*L.x;
    result.y = -dElastic * (Va.y - Vb.y) * L.y * L.y;
    result.z = -dElastic * (Va.z - Vb.z) * L.z * L.z;
    return result;
}
double dotProduct(point a, point b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
point findExternalForce(struct world * jello, int i, int j, int k) {
    point result = { 0.0, 0.0, 0.0 };
    if(jello->resolution==0)
        return result;
    point grid = getGrid(jello,i,j,k);
    result = trilinearInterpolation(grid,jello);
    return result;

}
point getGrid(struct world* jello, int i, int j, int k) {
    point pos = jello->p[i][j][k];
    //double h = 4 / jello->resolution;
    double inverseh = jello->resolution/4;
    double x, y, z;
    x = (2 + pos.x) * inverseh;
    y =(2 + pos.y) * inverseh;
    z = (2 + pos.z) * inverseh;
    //return the x,y,z integer part for grid number, decimal part for interpolation
    return {x,y,z};
}
point trilinearInterpolation(point grid, struct world* jello) {
    point result = { 0.0, 0.0, 0.0 };
    //may avoid type casting to make code more efficient
    int x, y, z;
    double alpha, beta, gamma;
    x = (int)(grid.x);
    y = (int)(grid.y);
    z = (int)(grid.z);
    alpha = grid.x - x;
    beta = grid.y - y;
    gamma = grid.z - z;
    if (x<jello->resolution-1&&y < jello->resolution-1&&z < jello->resolution-1&&x>=0&&y>=0&&z>=0)
    {
        point A000 = jello->forceField[x * jello->resolution * jello->resolution + y * jello->resolution + z];
        point A010 = jello->forceField[x * jello->resolution * jello->resolution + (y + 1) * jello->resolution + z];
        point A100 = jello->forceField[(x + 1) * jello->resolution * jello->resolution + y * jello->resolution + z];
        point A110 = jello->forceField[(x + 1) * jello->resolution * jello->resolution + (y + 1) * jello->resolution + z];

        point A001 = jello->forceField[x * jello->resolution * jello->resolution + y * jello->resolution + z + 1];
        point A011 = jello->forceField[x * jello->resolution * jello->resolution + (y + 1) * jello->resolution + z + 1];
        point A101 = jello->forceField[(x + 1) * jello->resolution * jello->resolution + y * jello->resolution + z + 1];
        point A111 = jello->forceField[(x + 1) * jello->resolution * jello->resolution + (y + 1) * jello->resolution + z + 1];
        //do trilinear Interpolation here
            result = ((1 - alpha) * (1 - beta) * (1 - gamma)) * A000 + ((1 - alpha) * beta * (1 - gamma)) * A010 +
            (alpha * (1 - beta) * (1 - gamma)) * A100 + (alpha * beta * (1 - gamma)) * A110 +
            ((1 - alpha) * (1 - beta) * gamma) * A001 + ((1 - alpha) * beta * gamma) * A011 +
            (alpha * (1 - beta) * gamma) * A101 + (alpha * beta * gamma) * A111;
    }
    return result;
}
point collisionHandler(struct world* jello, double a, double b, double c, double d, int i, int j, int k) {
    point collisionForce = {0.0,0.0,0.0};
    double x, y, z;
    x = jello->p[i][j][k].x;
    y = jello->p[i][j][k].y;
    z = jello->p[i][j][k].z;
    point n = { a,b,c };
    //check normal direction
    if (jello->incPlanePresent == 1)
    {
        double divisor = sqrt(a * a + b * b + c * c);
        point vector = { -d * a / divisor ,-d * b / divisor ,-d * c / divisor };
        //if the normal is point outward the origin, flip the sign
        if (vector.x * n.x + vector.y * n.y + vector.z * n.z > 0)
        {
            a = -a;
            b = -b;
            c = -c;
            d = -d;
            point newNormal = { a,b,c };
            n = newNormal;
        }
    }
    double fx = a * x + b * y + c * z + d;
    //define f(x,y,z) > 0 no collision, the normal should be all point inside the cube
    if (fx > 0)
        return collisionForce;
    //do collision handler here
    nomralize(n);
    double distance = (fabs(fx)/sqrt(a*a+b*b+c*c));
    //rest length is 0, the direction should be same as nomral
    point elasticityForce = jello->kCollision* distance *n;
    //damp force always depends on v
    point dampingForce = -jello->dCollision*jello->v[i][j][k]*n;
    collisionForce = elasticityForce + dampingForce;
    return collisionForce;
}
/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
  int i,j,k;
  point a[8][8][8];

  computeAcceleration(jello, a);
  
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
        jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
        jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
        jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
        jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
        jello->v[i][j][k].z += jello->dt * a[i][j][k].z;

      }
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
  point F1p[8][8][8], F1v[8][8][8], 
        F2p[8][8][8], F2v[8][8][8],
        F3p[8][8][8], F3v[8][8][8],
        F4p[8][8][8], F4v[8][8][8];

  point a[8][8][8];


  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(jello, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);
         pMULTIPLY(a[i][j][k],jello->dt,F1v[i][j][k]);
         pMULTIPLY(F1p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F1v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);
         // F2v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F2v[i][j][k]);
         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],1.0,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],1.0,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }
         
  computeAcceleration(&buffer, a);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
      }

  return;  
}
