/*
 * auxfuncs.hpp
 *
 *  Created on: Jan 29, 2013
 *      Author: luzdora
 */

#ifndef AUXFUNCS_HPP_
#define AUXFUNCS_HPP_

#include "../../../KLT/include/boost_incl.hpp"


#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>

namespace cluster{

using namespace boost_incl;
using namespace std;

#define PI 3.14159265
#define NUM_REG 20  // Number of test regions
//#define CARD_REG NUM_REG*NUM_REG
#define CARD_REG NUM_REG*NUM_REG*NUM_REG*NUM_REG
/** Define node (the head) in the binary tree */

typedef struct node {
	float *v; // caracterization vector (x, y, vx, vy)
	float dist; // Euclidean distance between s1 and s2 elements
	double nfa; // Number of False Alarms
	double nfag; // Number of False Alarms between a group pair
	double p; // probability of group
	int n; // number of points in the node
	int id; // number of node in the binary tree
	node *s1; //left element of node
	node *s2;  //right element of node
}node;

typedef node *nodeptr;

// Auxiliar list: Save elements that: NFA <1 and both childs are groups

typedef struct grupl {
	nodeptr g;
	grupl *prev;
	grupl *next;
}grupl;

typedef grupl *glptr;

/* Define a double linked list to save clusters resulted of singlelinkage and final meaningful clusters */

typedef struct list {
	float *v;
	list *prev;
	list *next;
}list;

typedef list *listptr;

// Some functions prototypes
listptr findheadlist(listptr _lc);
nodeptr findnodeaddress(nodeptr &, int _id);
glptr findheadgrouplist(glptr _gl);
int getnumpointspernode(nodeptr &);
void getnodepoints(nodeptr & , listptr &);
int minimalValue(float* eucdis, int ndis);
double factorial (double num);

double binomial(int m, int k, double p);
double trinomi(int Mp, int k1, int k2, double p1, double p2);
float intersectRegion(float sR1, float iR1, float sR2, float iR2);
void pointsInIntersectionRegion(glptr & eG, listptr *lc, int & p1, int & p2, int nP);
int indivisible(nodeptr &_p, int nP);
int checkMaximalEGroupsChilds(nodeptr p, nodeptr ch, int nP);

void printpoints(listptr &_lc);
void freeMatrix(float **m);
void freeArrayList( listptr *_alc, int nP);
void deletree(nodeptr &);
void deleteGroupList (glptr & _eG);
float **matrixReserve(int nr, int nc);
double **matrixDoubleReserve(int nr, int nc);
float *memorySpace(double numfeat);

} //end namespace
#endif /* AUXFUNCS_HPP_ */
