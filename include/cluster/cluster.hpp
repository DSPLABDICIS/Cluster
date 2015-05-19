#ifndef _CLUSTER_HPP_
#define _CLUSTER_HPP_


#include <cmath>
#include <cstdlib>
#include "clustersFuncs.hpp"
//#include <iostream>
//#include <cv.h>
//#include <highgui.h>


namespace cluster {

/** @ingroup cluster
	@brief Clustering method for interest points in an image plane. First, euclidean distance among all the points is obtained, then that is used to carry out the single linkage method. Next the Number of Falses Alarms is obtained for each group resulted of single linkage using the a contrario theory. Finally the resulted clusters are delivered selon the acontrario desition. (see. Thomas Veit, et. al. Space time A contrario clustering for Detecting Coherent Motions) 
	@pre This function receives a matrix with the feature points. Each feature point is characterized by 4 parameters: X and Y position and velocity in X and Y (velocity vector). Output of this function is a jblas vector with the number of cluster corresponding to each point of input feature vector 
	format of the data: x y dx dy val with a tabular space between them
	x = x position
	y = y position
	dx= (x2 - x1) where= x1: x position at image t-1; x2: x position at image t (current image) x2 is the same value as x position.
	dy= (y2 - y1) where= y1: y position at image t-1; y2: y position at image t (current image) y2 is the same value as y position. 
	p = Position in the global table.
 */
//  glptr makeglptrdata();


// PARA CAMBIAR POR JAFAR IMAGE  cv::Mat * imgR,


mat auxiliopelicanmurio();
void probandoMatrices(mat & mA, int b, int c);
int clusteringTrackingPoints(mat & Xp,
		int nC,
		int nR);
int clusteringMethod(mat & featl, int imaC, int imaR, nodeptr &bt, glptr & emg);

void saveclusters(mat & _xp,
		cv::Mat * imgR,
		std::string path2FeatureFile,
		int c,
		int fima);

}//end namespace cluster

#endif
