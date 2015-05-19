/*
 * clustersFuncs.cpp
 *
 *  Created on: Jan 29, 2013
 *      Author: luzdora
 */
#include "../include/cluster/clustersFuncs.hpp"


namespace cluster {

   // using namespace cluster;
	  using namespace boost::numeric::ublas;

    void drawPointToImg(cv::Mat &jRGBimg, int x, int y, int R, int G, int B)
    {
     // int ncols = jRGBimg.width();
     // int nrows = jRGBimg.height();
      int ncols = jRGBimg.cols;
      int nrows = jRGBimg.rows;
      int width = jRGBimg.step/3;
      int offset;
      unsigned char* RGBimg = (unsigned char *)jRGBimg.data;

      for (int yy = y - 1 ; yy <= y + 1 ; yy++)
	for (int xx = x - 1 ; xx <= x + 1 ; xx++)
	  if (xx >= 0 && yy >= 0 && xx < ncols && yy < nrows)
	    {
	      offset = (yy * width + xx)*3;
	      RGBimg[offset+0] = B;
	      RGBimg[offset+1] = G;
	      RGBimg[offset+2] = R;
	    }
    }




    /****************************************************************************************
	     #########################################################################################
		WriteClustersToImg
	     #########################################################################################
       ****************************************************************************************/
void  WriteClustersToImg(mat & _Xp, cv::Mat & jRGBimg, int nc)
{
  int i, j;

 // JFR_PRECOND(jRGBimg.colorSpace()==JfrImage_CS_BGR || jRGBimg.colorSpace()==JfrImage_CS_RGB, "jRGBimg is not in the RGB color space !");
 // JFR_PRECOND(jRGBimg.depth()==IPL_DEPTH_8U, "jRGBimg has not u8 depth !");
 // JFR_PRECOND(jRGBimg.channels()==3, "jRGBimg has not 3 channels !");
  if (nc == 0)
    { //que imprima los puntos contenidos en el arreglo de datos
      for (i=0;i<(int)_Xp.size1() ; i++)
        drawPointToImg(jRGBimg, (int)(_Xp(i,0) + 0.5), (int)(_Xp(i,1) + 0.5), 0,0,255);
    }
  else{
    if (nc==1)
      {
	for (i=0;i<(int)_Xp.size1() ; i++)
	  if (_Xp(i,4)==1)
	    drawPointToImg(jRGBimg, (int)(_Xp(i,0) + 0.5), (int)(_Xp(i,1) + 0.5), 255,150,0);

      }
    else{
      int r, g, b;
      long int aux;
      // Overlay features in red
      for ( i = 0 ; i < nc ; i++)
	{
	  r = (int)((255.0*i)/(nc-1));
	  aux = (65025.0*i)/(nc-1);
	  g = (int)(aux)%255;
	  aux = (16581375.0*i)/(nc-1);
	  b = (int)(aux)%65025;

	  for (j=0;j<(int)_Xp.size1() ; j++)
	    if (_Xp(j,4)==i+1)
	      drawPointToImg(jRGBimg, (int)(_Xp(j,0) + 0.5), (int)(_Xp(j,1) + 0.5), r,g,b);
	  //	  std::cout << "Cluster " << i+1 << " En color r= " << r << " g= " << g << " b= "<< b << std::endl;
	}
    }
  }
}

void rectangularToPolar (mat& _Xp, mat& _Xpol )
{
  int i;
  double aux;

  for (i=0; i <(int)_Xp.size1() ; i++)
    {
      _Xpol(i,2) = sqrt(_Xp(i,2)*_Xp(i,2) + _Xp(i,3)*_Xp(i,3)); //magnitude
      aux = _Xp(i,3)/_Xp(i,2);

      if ( (_Xp(i,2)>0) && (_Xp(i,3) >0) )  // first quad
	_Xpol(i,3) = atan(aux)*180/PI;
      else if (((_Xp(i,2) < 0) && (_Xp(i,3) < 0)) || ((_Xp(i,2) < 0)&&(_Xp(i,3) > 0)) )  // second, third quad
	_Xpol(i,3) = 180 + atan(aux)*180/PI;
      else
	_Xpol(i,3) = 360 + atan(aux)*180/PI;  // fourth quad

     std::cout << _Xp(i,0)  << "\t" << _Xp(i,1) << "\t" << _Xpol(i,2) << "\t" << _Xpol(i,3);
     std::cout << "\n";
    }
}

void getnodepoints_inx(nodeptr &_p, mat &_xc, int l)
{
  if(_p != NULL)
    {
      if ( _p->s1==NULL && _p->s2==NULL)
	{
	  _xc(_p->id,4)=l+1;
	}
      else
	{
	  getnodepoints_inx(_p->s1, _xc, l);
	  getnodepoints_inx(_p->s2, _xc, l);
	}
    }
}

void clusFuncs::matrixInitialization( int dim, int nP)
{
  // memory reservation of the matrix elements in the class
  R = matrixDoubleReserve( dim , NUM_REG );
  auxR = matrixDoubleReserve(dim, 3);   // save the intern results of group test regions
}

void clusFuncs::freeMatrixclass(double **m)
{
  if (m != NULL){
    free(m[0]);
    free(m);
  }
}
/*Obtain 20 different sizes of regions based on the geometric rule a*r^i
 where i = 0 .. 19, a= 1. The length maximal is set to image dimensions, maximal and
minimal velocity and maximal orientation (360) */

void clusFuncs::calcRegionSizes( int nDim, int imaC, int imaR)
{
  double a = 1.0, r[4], aux;
  int i, j;

  aux = log(imaC)/(NUM_REG-1);
  r[0]  = exp(aux);
  aux = log(imaR)/(NUM_REG-1);
  r[1]  = exp(aux);

  for (j = 0; j<nDim; j++){
     if (j == 2){
      a = minVel;
      aux = log((maxVel/minVel))/(NUM_REG-1);
      r[j] = exp(aux);
    }
    if (j == 3){
      a = 0.5;
      aux = log((360/a))/(NUM_REG-1);
      r[j] = exp(aux);
    }
    for (i=0; i<NUM_REG; i++)
      R[j][i] = a*pow(r[j],i);
  }
}

void clusFuncs::minValueVel(mat& _Xpol)
{
  int i;
  minVel = _Xpol(0,2);
  for (i=1; i<(int)_Xpol.size1(); i++)
    if (_Xpol(i,2) < minVel)
      minVel = _Xpol(i,2);
}

void clusFuncs::maxValueVel(mat& _Xpol)
{
  int i;
  maxVel = _Xpol(0,2);
  for (i=1; i<(int)_Xpol.size1(); i++)
    if (_Xpol(i,2) > maxVel)
      maxVel = _Xpol(i,2);
}

// Find the minimal one dimesion region which contains all the points of the group and its probability
// Region could be not simetrical in the center point

void clusFuncs::bestRegion(listptr &_lc, int c, int dim, int maxd, int mind)
{
  int i, flag;
  listptr cen=_lc;
  listptr aux=_lc;
  for (i=0; i<c; i++)
    cen=cen->next;

  for (i=0; i<NUM_REG; i++)
    {
      aux=_lc;
      flag=1;
      while ( aux != NULL){
	if ( ((cen->v[dim]+R[dim][i])< (aux->v[dim])) || ((cen->v[dim]-R[dim][i]) > (aux->v[dim])) ){
	  flag = 0;
	  break;
	}
	aux=aux->next;
      }
      if (flag == 1)
	break; // no more regions will be tested
    } //end for
  auxR[dim][0] = cen->v[dim]+R[dim][i];  // upper limit of the region
  auxR[dim][1] = cen->v[dim]-R[dim][i]; // lower limit of the region

  if ( auxR[dim][0] > maxd )
    auxR[dim][0]=maxd;
  if ( auxR[dim][1] < mind )
    auxR[dim][1]=mind;

  auxR[dim][2] = (auxR[dim][0] - auxR[dim][1])/maxd; //probability of the region
  // printf("\n probab de la region %d = %f", dim, auxR[dim][2]);
}

double clusFuncs::bestRegionVel(listptr &_lc, int c, mat & _Xp )
{
  int i, flag, dim=2, pIn=0;
  listptr cen=_lc;
  listptr aux=_lc;
  double piR;
  for (i=0; i<c; i++)
    cen=cen->next;

  for (i=0; i<NUM_REG; i++)
    {
      aux=_lc;
      flag=1;
      while ( aux != NULL){
	if ( ((cen->v[dim]+R[dim][i])< (aux->v[dim])) || ((cen->v[dim]-R[dim][i]) > (aux->v[dim])) ){
	  flag = 0;
	  break;
	}
	aux=aux->next;
      }
      if (flag == 1)
	break; // no more regions will be tested
    }
  auxR[dim][0] = cen->v[dim]+R[dim][i];  // upper limit of the region
  auxR[dim][1] = cen->v[dim]-R[dim][i]; // lower limit of the region

  if ( auxR[dim][0] > maxVel )
    auxR[dim][0]=maxVel;
  if ( auxR[dim][1] < minVel )
    auxR[dim][1]=minVel;

  auxR[dim][2] = (auxR[dim][0] - auxR[dim][1])/maxVel; //probability of the region
  //  printf("\n probab de la region %d = %f", dim, auxR[dim][2]);

  //Find probability density function of the region
  for (i=0; i < (int)_Xp.size1(); i++) // for each point of the data set
    {
      if ( (auxR[0][0] >= _Xp(i,0)) && (auxR[0][1] <=_Xp(i,0))   &&   (auxR[1][0] >= _Xp(i,1)) && (auxR[1][1] <=_Xp(i,1))   &&  (auxR[2][0] >= _Xp(i,2)) && (auxR[2][1] <=_Xp(i,2))    &&   (auxR[3][0] >= _Xp(i,3)) && (auxR[3][1] <=_Xp(i,3)) )
	pIn++;
    }
  auxR[dim][2]*= (((double)pIn)/((double)_Xp.size1()));
  //  printf("\n puntos en la region = %d, segunda probab = %f", pIn, auxR[dim][2]);

  piR = auxR[0][2]*auxR[1][2]*auxR[2][2]*auxR[3][2];
  return piR;
}

void clusFuncs::saveRegionTemp(double *ind)
{
  ind[0] = auxR[0][0];
  ind[1] = auxR[0][1];
  ind[2] = auxR[1][0];
  ind[3] = auxR[1][1];
  ind[4] = auxR[2][0];
  ind[5] = auxR[2][1];
  ind[6] = auxR[3][0];
  ind[7] = auxR[3][1];
}

void clusFuncs::findNFAg(glptr & _eG, listptr *_lc, int nP, int iC, int iR)
{
  float iX, iY, iV, iTet;
  double aux, piR1, piR2;
  int n1, n2;

  n1 = _eG->g->s1->n;
  n2 = _eG->g->s2->n;
  piR1 = _eG->g->s1->p;
  piR2 = _eG->g->s2->p;

  iX = intersectRegion(_eG->g->s1->v[0],_eG->g->s1->v[1],_eG->g->s2->v[0],_eG->g->s2->v[1]);
  iY = intersectRegion(_eG->g->s1->v[2],_eG->g->s1->v[3],_eG->g->s2->v[2],_eG->g->s2->v[3]);
  iV = intersectRegion(_eG->g->s1->v[4],_eG->g->s1->v[5],_eG->g->s2->v[4],_eG->g->s2->v[5]);
  iTet = intersectRegion(_eG->g->s1->v[6],_eG->g->s1->v[7],_eG->g->s2->v[6],_eG->g->s2->v[7]);

  if  ((iX*iY*iV*iTet) > 0 )
    { // There is an intersection region
      int ni1, ni2, intP;
      double piR1R2;
      pointsInIntersectionRegion(_eG, _lc, ni1, ni2, nP);
      intP = ni1 + ni2;
      if ( intP != 0)
	{
	  piR1R2 = (iX/iC)*(iY/iR)*(iV/maxVel)*(iTet/360)*(intP/nP);
	  piR1 = piR1 - piR1R2;
	  piR2 = piR2 - piR1R2;
	  n1 = n1 - ni1;
	  n2 = n2 - ni2;
	}
      /*  else
	{
	  piR1R2 = (iX/iC)*(iY/iR)*(iV/maxVel)*(iTet/360);
	  if (piR1R2 < piR1)
	    piR1 = piR1 - piR1R2;
	  if (piR1R2 < piR2)
	    piR2 = piR2 - piR1R2;
	  std::cout << "No hay puntos en la region de intersection pero si se intersectan " << std::endl;
	  }*/
    }
  //  else
  //  std::cout << "Not intersection region In node " << _eG->g->id << std::endl;

  if ((n1>0) && (n2>0) && (piR1>0) && (piR2>0)){
    aux = trinomi(nP-2,n1-1,n2-1,piR1,piR2);
    //  _eG->g->nfag = pow(nP,4)*CARD_REG*aux;
    _eG->g->nfag = pow(nP,4)*CARD_REG*CARD_REG*aux;
  }
  else
    _eG->g->nfag = 1.0; // Groups non interesting for merging

  // std::cout << " NFAg de " << _eG->g->id  <<" = " <<_eG->g->nfag << "NFA del mismo " <<_eG->g->nfa << " Con hijos " << _eG->g->s1->id << " y " << _eG->g->s2->id << std::endl;
}

/****************** AUXILIAR PRINT FUNCTIONS ******************/

void clusFuncs::printRegions(int _nDim)
{
  for (int j = 0; j<_nDim; j++){
    for (int i=0; i<NUM_REG; i++)
      std::cout << " R" << j << "= " << R[j][i] ;
    std::cout << std::endl;
  }
}

void clusFuncs::printNZ(int k, nodeptr &_aux)
{
  std::cout << "node " << k << std::endl;
  for (int j=0; j<8; j++)
    std::cout << "   " << _aux->v[j] << "\t";
}

// Obtain euclidean distance among the features, and save the result in pdistEuc
void clusFuncs::calcEuclideanDistance(mat & Xp_, int dim, float *pdistEuc, int tPoin)
{
  int i,j, k=0;
  vec Diff(dim);

  mat aux (Xp_.size1(), Xp_.size2());
  aux= Xp_;
  aux.resize(Xp_.size1(), 4, 1);

  for ( i=0; i < tPoin-1; i++ ) //number of features
    {
      for ( j=i+1; j<tPoin; j++) //residual points
	{
	  //  Diff= ublas::row(Xp_,i) - ublas::row(Xp_,j);
	  Diff= matrix_row< matrix<double> > (aux,i) - matrix_row< matrix<double> > (aux,j);
	  Diff= element_prod(Diff, Diff);
     	  *(pdistEuc+k) = inner_prod(Diff, trans(Diff));
	  *(pdistEuc+k) = sqrt(*(pdistEuc+k));

	  //	  std::cout << " Dist  " << *(pdistEuc+k) << " \t ";
	  k++;
	} //end for residual
    } //end for
}

void clusFuncs::singleLinkage(mat& _Xp, float* eucd, double sEucdis, int numPT, int dim, nodeptr &_nod )
{
  float **Z, *Y=NULL, *eucdis=NULL;
  int *npClus, *R, *I1=NULL, *I2=NULL, *I3=NULL, *I, *J;
  int indV, aux, i, j, indx, indY, a;
  int s, m, sI1, sI2, sI3, sI;   //sY;
  double sY;
//  nodeptr nod = NULL ;
  stack sav;
  treeb c;
  sY = sEucdis; //save the number of original distances
  m = numPT; //save the number of originals points

  //create memory for auxiliars vectors and matrix
  Z = matrixReserve(m-1, 3);
  Y = memorySpace(sEucdis);

  npClus = (int *) calloc (2*m-1, sizeof(int)); //points in each cluster
  R = (int *) calloc (m, sizeof(int)); //R is an index vector mapping the original index to the current

  //Copy of euclidean distances
  for ( i=0 ; i<sEucdis; i++ )
      *(Y+i) = *(eucd+i);

  //At the beginig there are number of clusters as numbers of points,
  //And R is at the begining C1 to Cm
  for (i=0; i<m; i++){
    *(npClus+i)=1;
    *(R+i)=i+1;
  }

  //init the single linkage algorithmic
  for (s=0; s<numPT-1; s++){

    // New data distances
    if (eucdis != NULL)
      free(eucdis);

    eucdis = memorySpace(sEucdis);
    if (eucdis == NULL || Y == NULL )
      {
	std::cout << " Tant pis no memory " << std::endl;
	exit(0);
      }

    for ( i=0 ; i<sEucdis; i++ )
      *(eucdis+i) = *(Y+i);

    //Find the minimal distance in euclidean distance vector
    indV = minimalValue(eucdis, sEucdis);

    //find the i and j index clusters of minimal value

    aux=indV+1;

    for (i=1; i<=m; i++){   //all the possible clusters
      if (aux>(m-i))
	aux= aux - m + i;
      else
	break;
    }

    j = i + aux;

    //save the index in vector result
    if ( *(R+i-1) < *(R+j-1) ){
      Z[s][0] = *(R+i-1)-1;
      Z[s][1] = *(R+j-1)-1;
    }
    else{
      Z[s][0] = *(R+j-1)-1;
      Z[s][1] = *(R+i-1)-1;
    }
    Z[s][2] = *(eucdis+indV);

/* save binary tree */

    //   std::cout << " a salvar Z[s][2] " << Z[s][2] << " Z[s][0] " << Z[s][0] <<  "Z[s][1] " <<  Z[s][1] << "\t";
    c.inserte(_Xp,  Z[s][2], (int) Z[s][0], (int) Z[s][1], numPT+s, dim , numPT, _nod, sav);

    //Searching the following values to be erased into eucdis Vector
    // Auxiliar vectors generated to each iteration
    sI1 = i-1;
    sI2 = j-i-1;
    sI3 = m-j;
    sI = sI1+sI2+sI3;

    if (sI1 > 0){
      I1 = (int *) calloc (sI1, sizeof(int));
      for (indx=0; indx<sI1; indx++)
	*(I1+indx)=indx+1;
    }

    if (sI2 > 0){
      I2 = (int *) calloc (sI2, sizeof(int));
      for (indx=0; indx<sI2; indx++)
	*(I2+indx)=i+1+indx;
    }

    if (sI3 > 0){
      I3 = (int *) calloc (sI3, sizeof(int));
      for (indx=0; indx<sI3; indx++)
	*(I3+indx)=j+1+indx;
    }
    if (sI>0)
      {
      	//All memory positions to be compared in eucdis because i cluster
	I = (int *) calloc (sI, sizeof(int));

	for(indx=0; indx<sI1; indx++)
	  *(I+indx)= (*(I1+indx)*m - ( ( *(I1+indx) + 1 )*( *(I1+indx))/2 ) -m + i) - 1;

	for ( a=0; indx < sI1+sI2; indx++, a++)
	  *(I+indx)= ( ( i*m - (i*(i+1))/2) - m + *(I2+a) ) - 1 ;

	for ( a=0; indx < sI; indx++, a++)
	  *(I+indx)= (( i*m -(i*(i+1))/2) - m + *(I3+a) ) - 1 ;

	//All memory positions to be compared in eucdis because j cluster

	J = (int *) calloc (sI+1, sizeof(int));

	for (indx=0; indx<sI1; indx++)
	  *(J+indx)= (*(I1+indx)*m - ( ( *(I1+indx) + 1 )*( *(I1+indx))/2 ) -m + j) - 1;

	for ( a=0 ; indx < sI1+sI2; indx++, a++)
	   *(J+indx)= (*(I2+a)*m - ( ( *(I2+a) + 1 )*( *(I2+a))/2 ) -m + j) - 1;

	for ( a=0; indx < sI; indx++, a++)
	  *(J+indx)= (( j*m - (j*(j+1))/2 ) - m + *(I3+a)) - 1 ;

	*(J+sI) = indV;

	//Find the minimal value between index I and J  and save them into eucdis vector
	for (indx=0; indx<sI; indx++)
	  {
	    if ( *(eucdis + (*(I+indx)) ) > *(eucdis + (*(J+indx))) )
	      *(eucdis + (*(I+indx))) = *(eucdis + (*(J+indx)));
	    *(eucdis+(*(J+indx))) = -1;
	  }
	*( eucdis + (*(J+sI)) ) = -1;

	// free auxiliar memory vectors

	if (sI1 != 0)
	  free(I1);

	if (sI2 != 0)
	  free(I2);

	if (sI3 != 0)
	  free(I3);

	free(I);
	free(J);
      }

    //update variables and distance vector
    m=m-1;
    *(npClus + numPT + s) =  *(npClus + *(R+i-1)) + *(npClus + *(R+j-1));
    *(R+i-1) = numPT+s+1;
    for (indx=j-1; indx<numPT-1; indx++)
      *(R+indx) = *(R+indx+1);

    sY=sEucdis-(sI+1); //minus J vector size

    if (Y != NULL){
      free(Y);
      Y = memorySpace(sY);
    }

    for (indx=0, indY=0; indx<sEucdis; indx++){
      if (*(eucdis+indx) == -1)
	continue;
      *(Y+indY)= *(eucdis+indx);
      indY++;
    }
    sEucdis=sY; //new size of euclidean vector
  }//end for linkage

  if (Y != NULL)
    free(Y);
  if (Z != NULL)
    freeMatrix(Z);
  if (R != NULL)
    free(R);
  if (npClus != NULL)
    free(npClus);
}


}//end namespace cluster




