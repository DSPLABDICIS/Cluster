/*
 * clustersFuncs.hpp
 *
 *  Created on: Jan 29, 2013
 *      Author: luzdora
 */

#ifndef CLUSTERSFUNCS_HPP_
#define CLUSTERSFUNCS_HPP_

#include "auxfuncs.hpp"
#include "../../../KLT/include/boost_incl.hpp"
#include <opencv/cv.h>
#include <opencv/highgui.h>


 namespace cluster {

 using namespace boost_incl;
/** This class let to contruct the binary tree */
    class stack
    {
    private:
      struct snode
      {
	nodeptr ele;
	snode *next;
      };
      snode *top;
    public:
      stack()
      {
	top=NULL;
      }
      void push(nodeptr p)
      {
	snode *temp;
	temp = new snode;
	temp->ele = p;
	temp->next = top;
	top=temp;
      }

      void pop(int _id)
      {
	if (top != NULL)
	  {
	    snode *temp1;
	    snode *temp2;
	    temp2 = top;

	    do {
	      if (temp2->ele->id == _id)
		{
		  if (temp2 == top)
		    top = temp2->next;
		  else
		    temp1->next = temp2->next;
		  delete temp2;
		}
	      temp1 = temp2;
	      temp2 = temp1->next;
	    }while (temp2 != NULL);
	  }
      }

      nodeptr getele(int _id)
      {
	snode* temp = top;
	while(temp != NULL)
	  {
	    if(temp->ele->id == _id)
	      return temp->ele;
	    temp=temp->next;

	  }
	return NULL;
      }
      nodeptr topele()
      {
	if (top !=NULL)
	  return top->ele;
	else
	  return NULL;
      }

      void printeles()
      {
	snode *atemp = top;
	while(atemp != NULL)
	  {
	    std::cout << atemp->ele->id << " " ;
	    atemp= atemp->next;
	  }
	std::cout << "\n" ;
      }

      int isempty()
      {
	return ((top == NULL) ? 1 : 0);
      }

    };

/** class to use the binary tree */
    class treeb {
    //	using namespace boost_incl;
    public:
      // void inserte(mat & _v, float , int , int , int , int , int , nodeptr &, stack &);
      void inserte( mat & _v , float _dist, int _id1, int _id2, int _id3,  int _n, int _nP, nodeptr &_p, stack &_lt);
      void saveFinalRegionNode(double *ind, nodeptr &, glptr &_eG, int nP);
      void saveTemporalEMeaninHeadNode(nodeptr &_aux, glptr &_eG);
      int detectMeaningfulGroup(nodeptr &_p, int nP);
      int valideNodeGroups(glptr & _eG, int &i);
      int isChildNode(nodeptr nodP, nodeptr nodCh);
      void printree(nodeptr &);

      void printeG(glptr &_eG);
      treeb(){}
      ~treeb(){}

    };

    /** Main functions of clustering method */
    class clusFuncs {
      double **R;   // R = Region sizes;
      double **auxR; /* pZ = proba, Maximal and minimal dimension region, points per group, centre and NFA; auxR = auxiliar that saves regions limits and proba */
      float maxVel;
      float minVel;
    public:

      /**
	  @brief default constructor and destructor
      */

      clusFuncs()
      {
	R = NULL;
	auxR = NULL;
      }

      ~clusFuncs()
      {
	if (R != NULL)
	  freeMatrixclass(R);
	if (auxR != NULL)
	  freeMatrixclass(auxR);
	std::cout << "bye bye clusFuncs" << std::endl;
      }

      void matrixInitialization( int dim, int nP);
      void minValueVel (mat & _Xpol);
      void maxValueVel (mat & _Xpol);
      void calcRegionSizes(int ndim, int imaC, int imaR);
      void calcEuclideanDistance( mat & Xp_, int dim, float *pdistEuc, int tPoin);
      void singleLinkage(mat & Xp_, float* eucd, double sEucdis, int numPT, int dim, nodeptr &_nod);

      void bestRegion(listptr &_lc,  int c, int dim, int maxd, int mind);
      double bestRegionVel(listptr &_lc,  int c, mat & _Xp);
      void saveRegionTemp(double *ind);
      void findNFAg(glptr & , listptr *_lc, int nP, int iC, int iR);
      void freeMatrixclass(double **m);
      void printNZ(int k, nodeptr & );
      void printRegions(int );

    }; //end class

    void rectangularToPolar(mat& _Xp, mat& _Xpol);
    void getnodepoints_inx(nodeptr &_p, mat &_xc, int l);
    void  WriteClustersToImg(mat & _Xp, cv::Mat & imgR, int nc);
    void drawPointToImg(cv::Mat &jRGBimg, int x, int y, int R, int G, int B);

  } //end namespace cluster

#endif /* CLUSTERSFUNCS_HPP_ */
