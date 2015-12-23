//copyright: (c) 2006 Justin Gardner, Jonas Larsson (GPL see mgl/COPYING)
//$Id: mglStrokeText.c 71 2006-10-20 21:19:00Z jonas $

/////////////////////////
//   include section   //
/////////////////////////
#include "mgl.h"
#include <string.h>
#include <stdlib.h>

////////////////////////
//   define section   //
////////////////////////

///////////////////////////////
//   function declarations   //
///////////////////////////////

void drawStrokeCharacter( float x, float y, float * strokesX, float * strokesY, float scaleX, float scaleY, int numstrokes, int hasdots, float rot, int rectwidth ) {
  // draw the strokes as a series of lines
  int i;
  float *tx,*ty;
  float globalscale=0.2; // scales to width==1
  float pointSize;
  scaleX=scaleX*globalscale;
  scaleY=scaleY*globalscale;
  tx=(float*)malloc(sizeof(float)*2*numstrokes);
  ty=(float*)malloc(sizeof(float)*2*numstrokes);

  if (rot!=0.0) {
    for (i=0;i<2*numstrokes;i++) {
      tx[i]=cos(rot)*(scaleX*strokesX[i]+x)+sin(rot)*(scaleY*strokesY[i]+y);
      ty[i]=-sin(rot)*(scaleX*strokesX[i]+x)+cos(rot)*(scaleY*strokesY[i]+y);
    }
  } else {
     for (i=0;i<2*numstrokes;i++) {
       tx[i]=scaleX*strokesX[i]+x;
       ty[i]=scaleY*strokesY[i]+y;
     }
  }

  glBegin(GL_LINES);
  for (i=0;i<2*numstrokes-1;i+=2){
    if (!hasdots || !(strokesX[i]==strokesX[i+1] && strokesY[i]==strokesY[i+1])) {      
      glVertex2f(tx[i],ty[i]);
      glVertex2f(tx[i+1],ty[i+1]);
    }
  }
  glEnd();

  if (hasdots) {
    if (rectwidth) {
      pointSize=(float)rectwidth/2;
      if (!(mglGetGlobalDouble("screenCoordinates")>0)) {
	pointSize=pointSize*mglGetGlobalDouble("xPixelsToDevice");
      }
      for(i=0;i<2*numstrokes-1;i+=2) {
	if (strokesX[i]==strokesX[i+1] && strokesY[i]==strokesY[i+1]) {
	  glRectd(tx[i]-pointSize,ty[i]-pointSize, tx[i]+pointSize,ty[i]+pointSize);
	}
      }

    } else {
      glBegin(GL_POINTS);
      for (i=0;i<2*numstrokes-1;i+=2) {
	if (strokesX[i]==strokesX[i+1] && strokesY[i]==strokesY[i+1]) {
	  glVertex2f(tx[i],ty[i]);
	}
      }
      glEnd();    
    }
  }
  free(tx);
  free(ty);
}

// Define static variables: 
// Characters: string of uppercase characters in ASCII order
#define MGLNUMCHARS 84
#define MGLMAXSTROKES 30 
// 15 max strokes per char; two points per stroke
//static char MGLChars[]="ABCDEFGHIJKLMNOPQRSTUVWXYZ-+:;!?*#/\\.,abcdefghijklmnopqrstuvwxyz";
static char MGLChars[]="ABCDEFGHIJKLMNOPQRSTUVWXYZ-+:;!?*#/\\.,()[]&=<>abcdefghijklmnopqrstuvwxyz\"'1234567890";
static int MGLCharsWithDots[]={0,0,0,0,0,0,0,0,0,0,\
			       0,0,0,0,0,0,0,0,0,0,\
			       0,0,0,0,0,0,0,0,1,1,\
			       1,1,0,0,0,0,1,0,0,0,\
			       0,0,0,1,1,0,0,0,0,0,\
			       0,0,0,0,1,1,0,0,0,0,\
			       0,0,0,0,0,0,0,0,0,0,\
			       0,0,0,0,0,0,0,0,0,0,\
			       0,0,0,0};
// CharacterStrokesX,CharacterStrokesY: matrices of character vertices in alphabetical (ASCII) order
// CharacterNumStrokes: vector of number of strokes per character
static int MGLCharNumStrokes[]={\
				3,10,7,6,4,3,9,3,1,6,\
				3,2,4,3,8,6,9,7,9,2,\
				5,2,4,2,3,3,1,2,2,2,\
				2,9,3,4,1,1,1,2,5,5,\
				3,3,9,2,2,2,7,7,7,7,\
				9,5,10,5,2,3,3,1,7,5,\
				8,7,7,4,7,5,5,2,4,2,\
				3,3,4,2,3,7,8,3,9,11,\
				3,15,11,9};

static float MGLCharStrokesX[MGLNUMCHARS][MGLMAXSTROKES]={\
							  {-2,0,0,2,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,-2,-2,1,1,2,2,2,2,1,1,2,2,2,2,1,1,-2,-2,1,0,0,0,0,0,0,0,0,0,0},\
							  {2,1,1,-1,-1,-2,-2,-2,-2,-1,-1,1,1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,-2,-2,1,1,2,2,2,2,1,1,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {2,-2,-2,-2,-2,2,-2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,2,-2,-2,-2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {2,1,1,-1,-1,-2,-2,-2,-2,-1,-1,1,1,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,-2,2,2,-2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,2,2,2,2,1,1,-1,-1,-2,-2,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,-2,-2,2,-1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,-2,-2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,-2,-2,0,0,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,-2,-2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-1,1,1,2,2,2,2,1,1,-1,-1,-2,-2,-2,-2,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,-2,-2,1,1,2,2,2,2,1,1,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-1,1,1,2,2,2,2,1,1,-1,-1,-2,-2,-2,-2,-1,1,2,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,-2,-2,1,1,2,2,2,2,1,1,-2,-1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {2,1,1,-1,-1,-2,-2,-2,-2,2,2,2,2,1,1,-1,-1,-2,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,-2,-2,-1,-1,1,1,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,-1,-1,0,0,1,1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,2,-2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,2,2,-2,-2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,-2,-2,-1,-1,1,1,2,2,2,2,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,2,-1,1,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,2,-2,2,-1,-1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {2,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {1,0,0,-1,-1,-1,-1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-1,0,0,1,1,1,1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {1,-1,-1,-1,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-1,1,1,1,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {2,1,1,-1,-1,-2,-2,-2,-2,1,1,0,0,-1,-1,-2,-2,2,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,2,2,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {2,-2,-2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,2,2,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {2,1,1,-1,-1,-2,-2,-2,-2,-1,-1,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,-2,-2,-1,-1,1,1,2,2,2,2,1,1,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {1,-1,-1,-2,-2,-2,-2,-1,-1,1,1,2,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {2,2,2,-1,-1,-2,-2,-2,-2,-1,-1,1,1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,2,2,1,1,-1,-1,-2,-2,-2,-2,-1,-1,1,1,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-1,-1,-1,0,0,1,1,2,-2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {2,1,1,-1,-1,-2,-2,-2,-2,-1,-1,2,2,2,2,1,1,-1,-1,-2,0,0,0,0,0,0,0,0,0,0},\
							  {-2,-2,-2,-1,-1,1,1,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,-2,-2,1,-1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,-2,-2,-1,0,-1,0,1,1,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,-2,-2,-1,-1,1,1,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,-1,-1,1,1,2,2,2,2,1,1,-1,-1,-2,-2,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,-2,-2,-1,-1,1,1,2,2,2,2,1,1,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {2,2,2,1,1,-1,-1,-2,-2,-2,-2,-1,-1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,-2,-2,-1,-1,1,1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {2,1,1,-1,-1,-2,-2,-1,-1,1,1,2,2,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-1,-1,-1,0,0,1,1,2,-2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,-2,-2,-1,-1,1,1,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,-1,-1,0,0,1,1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,2,-2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-1,2,0,-2,-1,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,2,2,-2,-2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {0,0,0,-1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-1,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,2,2,2,-2,2,2,1,1,-1,-1,-2,-2,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,2,-2,-2,-2,-1,-1,1,1,2,2,2,2,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {1,-2,-2,2,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,2,-2,-1,-1,1,1,2,2,2,2,1,1,-1,-1,-2,-2,-2,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,-1,-1,1,1,2,2,2,2,1,1,-1,-1,-2,-2,-2,-2,-1,-1,1,1,2,0,0,0,0,0,0,0,0},\
							  {-2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-1,1,1,2,2,2,2,1,1,-1,-1,-2,-2,-2,-2,-1,-1,-2,-2,-2,-2,-1,-1,1,1,2,2,2,2,1},\
							  {-2,-1,-1,1,1,2,-2,-2,-2,-1,-1,1,1,2,2,2,2,1,1,-1,-1,-2,0,0,0,0,0,0,0,0},\
							  {-2,-1,-1,1,1,2,2,2,2,1,1,-1,-1,-2,-2,-2,-2,2,0,0,0,0,0,0,0,0,0,0,0,0}};

static float MGLCharStrokesY[MGLNUMCHARS][MGLMAXSTROKES]={\
							  {-4,4,4,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-4,4,4,4,4,3,3,1,1,0,0,-1,-1,-3,-3,-4,-4,-4,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {3,4,4,4,4,3,3,-3,-3,-4,-4,-4,-4,-3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,-4,-4,-4,-4,-3,-3,3,3,4,4,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,4,4,-4,-4,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,4,4,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {3,4,4,4,4,3,3,-3,-3,-4,-4,-4,-4,-3,-1,-3,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,-4,-4,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,4,4,-3,-3,-4,-4,-4,-4,-3,-2,-3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-4,4,0,4,1,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,-4,-4,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-4,4,4,-1,-1,4,4,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-4,4,4,-4,-4,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,4,4,3,3,-3,-3,-4,-4,-4,-4,-3,-3,3,3,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,-4,4,4,4,3,3,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,4,4,3,3,-3,-3,-4,-4,-4,-4,-3,-3,3,3,4,-3,-4,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,-4,4,4,4,3,3,1,1,0,0,0,0,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {3,4,4,4,4,3,3,1,1,-1,-1,-3,-3,-4,-4,-4,-4,-3,0,0,0,0,0,0,0,0,0,0,0,0}, \
							  {4,4,4,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,-3,-3,-4,-4,-4,-4,-3,-3,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,-4,-4,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,-4,-4,1,1,-4,-4,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,-4,-4,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,0,0,4,0,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,4,4,-4,-4,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {0,0,2,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {1.5,1.5,-3,-3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {2,2,-3,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,-2,-4,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {2,3,3,4,4,4,4,3,1,3,1,0,0,-1,-1,-3,-4,-4,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {0,0,2,-2,2,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {1,1,-1,-1,2,-2,2,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {3,-3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {3,-3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-4,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-3,-4,-4,-5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0},\
							  {4,3,3,1,1,-1,-1,-3,-3,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,3,3,1,1,-1,-1,-3,-3,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,4,4,-4,-4,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,4,4,-4,-4,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-1,-4,-4,-4,-4,-3,-3,-1,-1,3,3,4,4,4,4,3,3,-4,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {1,1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {2,0,0,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {2,0,0,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-3,-4,-4,-4,-4,-3,-3,-1,-1,0,0,0,0,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,-4,-1,0,0,0,0,-1,-1,-3,-3,-4,-4,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {0,0,0,-1,-1,-3,-3,-4,-4,-4,-4,-3,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,-4,0,0,0,-1,-1,-3,-3,-4,-4,-4,-4,-3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-2,-2,-3,-4,-4,-4,-4,-3,-3,-1,-1,0,0,0,0,-1,-1,-2,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-4,2,2,3,3,3,3,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-3,-4,-4,-4,-4,-3,-3,-1,-1,0,0,0,0,-6,-6,-7,-7,-7,-7,-6,0,0,0,0,0,0,0,0,0,0},\
							  {4,-4,-1,0,0,0,0,-1,-1,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {0,-4,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {0,-5,-5,-6,1,1,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,0,0,0,0,0,0,0,0,0,0},\
							  {4,-4,-2,1,-1,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-4,0,-1,0,-1,0,-1,0,0,-1,-1,-4,-4,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {0,-4,-1,0,0,0,0,-1,-1,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-1,0,0,0,0,-1,-1,-3,-3,-4,-4,-4,-4,-3,-3,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {0,-7,-1,0,0,0,0,-1,-1,-3,-3,-4,-4,-4,-4,-4,-4,-4,-4,-4,0,0,0,0,0,0,0,0,0,0},\
							  {0,-7,-1,0,0,0,0,-1,-1,-3,-3,-4,-4,-4,-4,-4,-4,-4,-4,-4,0,0,0,0,0,0,0,0,0,0},\
							  {-4,0,-1,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-1,0,0,0,0,-1,-3,-4,-4,-4,-4,-3,-3,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {2,-3,-3,-4,-4,-4,-4,-3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {0,-3,-3,-4,-4,-4,-4,-3,0,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {0,-4,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {0,-4,-4,-2,-2,-4,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {0,-4,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-6,0,-4,0,-6,-7,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,0,0,0,0,0,0,0,0,0,0},\
							  {0,0,0,-4,-4,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,3,3,2,4,3,3,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,3,3,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {3,4,4,-4,-4,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-4,-4,3,1,-4,1,3,4,4,4,4,3,3,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,4,-2,-3,-3,-4,-4,-4,-4,-3,-3,0,0,1,1,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,-1,-1,-1,4,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,4,-3,-4,-4,-4,-4,-3,-3,0,0,1,1,1,1,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {-3,-4,-4,-4,-4,-3,-3,0,0,1,1,1,1,0,-3,3,3,4,4,4,4,3,0,0,0,0,0,0,0,0},\
							  {4,4,4,0,0,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},\
							  {4,4,4,3,3,1,1,0,0,0,0,1,1,3,3,4,0,-1,-1,-3,-3,-4,-4,-4,-4,-3,-3,-1,-1,0},\
							  {3,4,4,4,4,3,3,0,0,-1,-1,-1,-1,0,3,-3,-3,-4,-4,-4,-4,-3,0,0,0,0,0,0,0,0},\
							  {3,4,4,4,4,3,3,-3,-3,-4,-4,-4,-4,-3,-3,3,-3,3,0,0,0,0,0,0,0,0,0,0,0,0}\
};


/////////////
//   main   //
//////////////
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  // mglStrokeText( string, x, y, scalex, scaley, linewidth, color, rotation )
  if ( nrhs<3 )
    return;


  // Set character string variable
  // Set scale
  // Set starting position
  int buflen=mxGetN( prhs[0] )*mxGetM( prhs[0] )+1;
  char * inputString=(char*) malloc(buflen);
  mxGetString( prhs[0], inputString, buflen);
  float x,y,scaleX,scaleY;
  float linewidth,rot;
  float color[3];
  
  x=(float) *mxGetPr( prhs[1] );
  y=(float) *mxGetPr( prhs[2] );
  if (nrhs>3) 
    scaleX=(float) *mxGetPr( prhs[3] );
  else
    scaleX=1.0;
  if (nrhs>4) 
    scaleY=(float) *mxGetPr( prhs[4] );
  else
    scaleY=scaleX;
  if (nrhs>5) 
    linewidth=(float) *mxGetPr( prhs[5] );
  else
    linewidth=1;
  if (nrhs>6) {
    color[0]=(float)*mxGetPr( prhs[6] );
    color[1]=(float)*(mxGetPr( prhs[6] )+1);
    color[2]=(float)*(mxGetPr( prhs[6] )+2);
  } else {
    color[0]=1.0;
    color[1]=1.0;
    color[2]=1.0;
  }
  if (nrhs>7) 
    rot=(float) *mxGetPr( prhs[7] );
  else
    rot=0.0;

  glLineWidth(linewidth);

  int rectWidth=0;
  GLint range[2];
  glGetIntegerv(GL_POINT_SIZE_RANGE,range);
  if (range[1]<linewidth) {
    rectWidth=linewidth+1;
  } else {
    glPointSize(linewidth+1);
  }

  char space[]=" ";
  int iS,iC,numstrokes;
  float * strokesX;
  float * strokesY;
  glColor3f(color[0],color[1],color[2]);
  glEnable(GL_LINE_SMOOTH);
  bool charFound;
  
  // Select each character in string
  for (iS=0; iS<buflen-1; iS++) {
    // note that buflen includes null \0 character which we ignore
    iC=0;
    charFound=false;
    if (inputString[iS]==(*space)) {
      // advance x
      x+=(scaleX*1.1);
      continue;
    } else {
      // Find index of corresponding character      
      for (iC=0; iC<MGLNUMCHARS; iC++) {
	if (inputString[iS]==MGLChars[iC]) {
	  //	  mexPrintf("%c %c\n",inputString[iS],MGLChars[iC]);
	  numstrokes=MGLCharNumStrokes[iC];
	  strokesX=MGLCharStrokesX[iC];
	  strokesY=MGLCharStrokesY[iC];
	  // Draw character
	  //	  mexPrintf("%i %i\n",iC, iS);
	  drawStrokeCharacter( x, y, strokesX, strokesY, scaleX, scaleY, numstrokes, MGLCharsWithDots[iC], rot, rectWidth );	
	  // Update x,y
	  x+=(scaleX*1.1); // letter spacing
	  charFound=true;
	  break;
	}	
      }
      if (!charFound) {
	// draw hash
	iC=33;
	numstrokes=MGLCharNumStrokes[iC];
	strokesX=MGLCharStrokesX[iC];
	strokesY=MGLCharStrokesY[iC];
	// Draw character
	//	mexPrintf("%i %i\n",iC, iS);
	drawStrokeCharacter( x, y, strokesX, strokesY, scaleX, scaleY, numstrokes, MGLCharsWithDots[iC], rot, rectWidth );	
	// Update x,y
	x+=(scaleX*1.1); // letter spacing
      }
    }
  }

  if (nlhs==2) {
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(plhs[0]) = x;
    *mxGetPr(plhs[1]) = y;
  }
}

