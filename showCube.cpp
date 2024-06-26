/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "showCube.h"
#include "physics.h"
#include "stb_image.h"

int pointMap(int side, int i, int j)
{
  int r;

  switch (side)
  {
  case 1: //[i][j][0] bottom face
    r = 64 * i + 8 * j;
    break;
  case 6: //[i][j][7] top face
    r = 64 * i + 8 * j + 7;
    break;
  case 2: //[i][0][j] front face
    r = 64 * i + j;
    break;
  case 5: //[i][7][j] back face
    r = 64 * i + 56 + j;
    break;
  case 3: //[0][i][j] left face
    r = 8 * i + j;
    break;
  case 4: //[7][i][j] right face
    r = 448 + 8 * i + j;
    break;
  }

  return r;
}

void showCube(struct world * jello, GLenum mode)
{
  int i,j,k,ip,jp,kp;
  point r1,r2,r3; // aux variables
  
  /* normals buffer and counter for Gourad shading*/
  struct point normal[8][8];
  int counter[8][8];

  int face;
  double faceFactor, length;

  //RK4(jello);
  if (fabs(jello->p[0][0][0].x) > 10)
  {
    printf ("Your cube somehow escaped way out of the box.\n");
    exit(0);
  }

  
  #define NODE(face,i,j) (*((struct point * )(jello->p) + pointMap((face),(i),(j))))

  
  #define PROCESS_NEIGHBOUR(di,dj,dk) \
    ip=i+(di);\
    jp=j+(dj);\
    kp=k+(dk);\
    if\
    (!( (ip>7) || (ip<0) ||\
      (jp>7) || (jp<0) ||\
    (kp>7) || (kp<0) ) && ((i==0) || (i==7) || (j==0) || (j==7) || (k==0) || (k==7))\
       && ((ip==0) || (ip==7) || (jp==0) || (jp==7) || (kp==0) || (kp==7))) \
    {\
      glVertex3f(jello->p[i][j][k].x,jello->p[i][j][k].y,jello->p[i][j][k].z);\
      glVertex3f(jello->p[ip][jp][kp].x,jello->p[ip][jp][kp].y,jello->p[ip][jp][kp].z);\
    }\

 
  if (viewingMode==0) // render wireframe
  {
    glLineWidth(1);
    glPointSize(5);
    glDisable(GL_LIGHTING);
    for (i=0; i<=7; i++)
      for (j=0; j<=7; j++)
        for (k=0; k<=7; k++)
        {
          if (i*j*k*(7-i)*(7-j)*(7-k) != 0) // not surface point
            continue;
          //add name stack here
          if (mode == GL_SELECT)
              glLoadName(i*8*8+8*j+k);
          glBegin(GL_POINTS); // draw point
            glColor4f(0,0,0,0.5);
            glVertex3f(jello->p[i][j][k].x,jello->p[i][j][k].y,jello->p[i][j][k].z);        
          glEnd();

          //
          //if ((i!=7) || (j!=7) || (k!=7))
          //  continue;

          glBegin(GL_LINES);      
          // structural
          if (structural == 1)
          {
            glColor4f(0,0,1,1);
            PROCESS_NEIGHBOUR(1,0,0);
            PROCESS_NEIGHBOUR(0,1,0);
            PROCESS_NEIGHBOUR(0,0,1);
            PROCESS_NEIGHBOUR(-1,0,0);
            PROCESS_NEIGHBOUR(0,-1,0);
            PROCESS_NEIGHBOUR(0,0,-1);
          }
          
          // shear
          if (shear == 1)
          {
            glColor4f(0,1,0,1);
            PROCESS_NEIGHBOUR(1,1,0);
            PROCESS_NEIGHBOUR(-1,1,0);
            PROCESS_NEIGHBOUR(-1,-1,0);
            PROCESS_NEIGHBOUR(1,-1,0);
            PROCESS_NEIGHBOUR(0,1,1);
            PROCESS_NEIGHBOUR(0,-1,1);
            PROCESS_NEIGHBOUR(0,-1,-1);
            PROCESS_NEIGHBOUR(0,1,-1);
            PROCESS_NEIGHBOUR(1,0,1);
            PROCESS_NEIGHBOUR(-1,0,1);
            PROCESS_NEIGHBOUR(-1,0,-1);
            PROCESS_NEIGHBOUR(1,0,-1);

            PROCESS_NEIGHBOUR(1,1,1)
            PROCESS_NEIGHBOUR(-1,1,1)
            PROCESS_NEIGHBOUR(-1,-1,1)
            PROCESS_NEIGHBOUR(1,-1,1)
            PROCESS_NEIGHBOUR(1,1,-1)
            PROCESS_NEIGHBOUR(-1,1,-1)
            PROCESS_NEIGHBOUR(-1,-1,-1)
            PROCESS_NEIGHBOUR(1,-1,-1)
          }
          
          // bend
          if (bend == 1)
          {
            glColor4f(1,0,0,1);
            PROCESS_NEIGHBOUR(2,0,0);
            PROCESS_NEIGHBOUR(0,2,0);
            PROCESS_NEIGHBOUR(0,0,2);
            PROCESS_NEIGHBOUR(-2,0,0);
            PROCESS_NEIGHBOUR(0,-2,0);
            PROCESS_NEIGHBOUR(0,0,-2);
          }           
          glEnd();
        }
    glEnable(GL_LIGHTING);
  }
  
  else
  {
    glPolygonMode(GL_FRONT, GL_FILL); 
    
    for (face=1; face <= 6; face++) 
      // face == face of a cube
      // 1 = bottom, 2 = front, 3 = left, 4 = right, 5 = far, 6 = top
    {
      
      if ((face==1) || (face==3) || (face==5))
        faceFactor=-1; // flip orientation
      else
        faceFactor=1;
      

      for (i=0; i <= 7; i++) // reset buffers
        for (j=0; j <= 7; j++)
        {
          normal[i][j].x=0;normal[i][j].y=0;normal[i][j].z=0;
          counter[i][j]=0;
        }

      /* process triangles, accumulate normals for Gourad shading */
  
      for (i=0; i <= 6; i++)
        for (j=0; j <= 6; j++) // process block (i,j)
        {
          pDIFFERENCE(NODE(face,i+1,j),NODE(face,i,j),r1); // first triangle
          pDIFFERENCE(NODE(face,i,j+1),NODE(face,i,j),r2);
          CROSSPRODUCTp(r1,r2,r3); pMULTIPLY(r3,faceFactor,r3);
          pNORMALIZE(r3);
          pSUM(normal[i+1][j],r3,normal[i+1][j]);
          counter[i+1][j]++;
          pSUM(normal[i][j+1],r3,normal[i][j+1]);
          counter[i][j+1]++;
          pSUM(normal[i][j],r3,normal[i][j]);
          counter[i][j]++;

          pDIFFERENCE(NODE(face,i,j+1),NODE(face,i+1,j+1),r1); // second triangle
          pDIFFERENCE(NODE(face,i+1,j),NODE(face,i+1,j+1),r2);
          CROSSPRODUCTp(r1,r2,r3); pMULTIPLY(r3,faceFactor,r3);
          pNORMALIZE(r3);
          pSUM(normal[i+1][j],r3,normal[i+1][j]);
          counter[i+1][j]++;
          pSUM(normal[i][j+1],r3,normal[i][j+1]);
          counter[i][j+1]++;
          pSUM(normal[i+1][j+1],r3,normal[i+1][j+1]);
          counter[i+1][j+1]++;
        }

        /* the actual rendering */
        for (j=1; j<=7; j++) 
        {
          if (faceFactor  > 0)
            glFrontFace(GL_CCW); // the usual definition of front face
          else
            glFrontFace(GL_CW); // flip definition of orientation
         
          glBegin(GL_TRIANGLE_STRIP);
          for (i=0; i<=7; i++)
          {
            glNormal3f(normal[i][j].x / counter[i][j],normal[i][j].y / counter[i][j],
              normal[i][j].z / counter[i][j]);
            glVertex3f(NODE(face,i,j).x, NODE(face,i,j).y, NODE(face,i,j).z);
            glNormal3f(normal[i][j-1].x / counter[i][j-1],normal[i][j-1].y/ counter[i][j-1],
              normal[i][j-1].z / counter[i][j-1]);
            glVertex3f(NODE(face,i,j-1).x, NODE(face,i,j-1).y, NODE(face,i,j-1).z);
          }
          glEnd();
        }
       
    }  
  } // end for loop over faces
  glFrontFace(GL_CCW);
}
void loadTexture() {

    unsigned int texture;
    glGenTextures(1, &texture);
    // set the texture wrapping parameters
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
    glBindTexture(GL_TEXTURE_2D, texture); // all upcoming GL_TEXTURE_2D operations now have effect on this texture object

    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    // load image, create texture and generate mipmaps
    int width, height, nrChannels;
    // The FileSystem::getPath(...) is part of the GitHub repository so we can find files on any IDE/platform; replace it with your own image path.
    unsigned char* data = stbi_load("glass.png", &width, &height, &nrChannels, 0);
    if (data)
    {
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
    }
    else
    {
        printf("Failed to load texture");
    }
    stbi_image_free(data);
    glEnable(GL_TEXTURE_2D);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    //use blend to make the texture transparent
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
    glDisable(GL_CULL_FACE);
    glColor4f(0.6, 0.6, 0.6, 0.2);
    glBegin(GL_QUADS);
    //front face texture
    glTexCoord2f(0.0, 0.0);
    glVertex3f(2, -2, 2);
    glTexCoord2f(0.0, 1.0);
    glVertex3f(2, -2, -2);
    glTexCoord2f(1.0, 1.0);
    glVertex3f(2, 2, -2);
    glTexCoord2f(1.0, 0.0);
    glVertex3f(2, 2, 2);
    //back face texture
    glTexCoord2f(0.0, 0.0);
    glVertex3f(-2, -2, 2);
    glTexCoord2f(0.0, 1.0);
    glVertex3f(-2, -2, -2);
    glTexCoord2f(1.0, 1.0);
    glVertex3f(-2, 2, -2);
    glTexCoord2f(1.0, 0.0);
    glVertex3f(-2, 2, 2);
    //right face texture
    glTexCoord2f(0.0, 0.0);
    glVertex3f(-2, 2, 2);
    glTexCoord2f(0.0, 1.0);
    glVertex3f(-2, 2, -2);
    glTexCoord2f(1.0, 1.0);
    glVertex3f(2, 2, -2);
    glTexCoord2f(1.0, 0.0);
    glVertex3f(2, 2, 2);
    //left face texture
    glTexCoord2f(0.0, 0.0);
    glVertex3f(-2, -2, 2);
    glTexCoord2f(0.0, 1.0);
    glVertex3f(-2, -2, -2);
    glTexCoord2f(1.0, 1.0);
    glVertex3f(2, -2, -2);
    glTexCoord2f(1.0, 0.0);
    glVertex3f(2, -2, 2);
    //bottom face texture
    glTexCoord2f(0.0, 0.0);
    glVertex3f(-2, 2, -2);
    glTexCoord2f(0.0, 1.0);
    glVertex3f(-2, -2, -2);
    glTexCoord2f(1.0, 1.0);
    glVertex3f(2, -2, -2);
    glTexCoord2f(1.0, 0.0);
    glVertex3f(2, 2, -2);
    //top face texture
    glTexCoord2f(0.0, 0.0);
    glVertex3f(-2, 2, 2);
    glTexCoord2f(0.0, 1.0);
    glVertex3f(-2, -2, 2);
    glTexCoord2f(1.0, 1.0);
    glVertex3f(2, -2, 2);
    glTexCoord2f(1.0, 0.0);
    glVertex3f(2, 2, 2);
    glEnd();
    glDisable(GL_TEXTURE_2D);
    glEnable(GL_CULL_FACE);
    glDisable(GL_BLEND);

}
void showBoundingBox()
{
  int i,j;
  if (addTexture == 1)
      loadTexture();
  else 
  {
      glColor4f(0.6, 0.6, 0.6, 1.0);

      glBegin(GL_LINES);
      // front face
      for (i = -2; i <= 2; i++)
      {
          glVertex3f(i, -2, -2);
          glVertex3f(i, -2, 2);
      }
      for (j = -2; j <= 2; j++)
      {
          glVertex3f(-2, -2, j);
          glVertex3f(2, -2, j);
      }

      // back face
      for (i = -2; i <= 2; i++)
      {
          glVertex3f(i, 2, -2);
          glVertex3f(i, 2, 2);
      }
      for (j = -2; j <= 2; j++)
      {
          glVertex3f(-2, 2, j);
          glVertex3f(2, 2, j);
      }

      // left face
      for (i = -2; i <= 2; i++)
      {
          glVertex3f(-2, i, -2);
          glVertex3f(-2, i, 2);
      }
      for (j = -2; j <= 2; j++)
      {
          glVertex3f(-2, -2, j);
          glVertex3f(-2, 2, j);
      }

      // right face
      for (i = -2; i <= 2; i++)
      {
          glVertex3f(2, i, -2);
          glVertex3f(2, i, 2);
      }
      for (j = -2; j <= 2; j++)
      {
          glVertex3f(2, -2, j);
          glVertex3f(2, 2, j);
      }

      glEnd();
  }
  return;
}
void showInclinePlane(struct world* jello) {
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
    glColor4f(0.1, 0, 0, 1.0);
    double a = jello->a;
    double b = jello->b;
    double c = jello->c;
    double d = jello->d;
    double z1, z2, z3, z4;
    //check if c is 0
    if (c != 0)
    {
         z1 = (-d - a * (2) - b * (2)) / c;
         z2 = (-d - a * (2) - b * (-2)) / c;
         z3 = (-d - a * (-2) - b * (2)) / c;
         z4 = (-d - a * (-2) - b * (-2)) / c;
         glBegin(GL_QUADS);
         glVertex3f(2, 2, z1);
         glVertex3f(2, -2, z2);
         glVertex3f(-2, -2, z4);
         glVertex3f(-2, 2, z3);
         glEnd();
         glDisable(GL_BLEND);
    }
    else
    {
        glBegin(GL_QUADS);
        glVertex3f(2, (-d-(a*2))/b, 2);
        glVertex3f(2, ( - d - (a * -2)) / b, 2);
        glVertex3f(-2, (-d - (a * -2)) / b, -2);
        glVertex3f(-2, (-d - (a * 2)) / b, -2);
        glEnd();
        glDisable(GL_BLEND);
    }
    
}

