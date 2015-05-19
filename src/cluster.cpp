
#include "cluster.hpp"


namespace cluster {

using namespace std;
using namespace boost_incl;
//using namespace cluster;


int clusteringTrackingPoints(mat & Xp,
				 int          nC,
				 int          nR)
   {

     int l, c;
     //  double NFA;
     // To save cluster results
     nodeptr btree = NULL;  // binary tree
     glptr eGroup = NULL; // head of groups in the binary tree that are clusters
     // Call clustering method which returns the number of cluster founded

     c = clusteringMethod(Xp, nC, nR, btree, eGroup);
     std::cout << "Numero de clusters " << c << std::endl;

     //c = 0; //todo: this is temporary remember to remove

     if(c>0)   // save and display cluster points
	{
	  l=0;
	  do
	    { //node aux is meaningful, then save it
	      getnodepoints_inx(eGroup->g, Xp,l);

	      //    std::cout << " Number of group to save " << eGroup->g->id << " num points = " << eGroup->g->n << "NFA : " <<  -log(eGroup->g->nfa) <<std::endl;
	      l++;
	      eGroup = eGroup->next;
	    }while(eGroup != NULL);

	  deletree(btree);
	  //	  int fc = secondFusion(Xp, c);
	  return c;
	}
     else
	{
	  std::cout<< "There is not clusters" << std::endl;
	  return 0; //there is not clusters in dataset
	}
   }



mat auxiliopelicanmurio()

{
  /*
  jblas::mat _Xp(11,4) ;
  _Xp(0,0)=40.4476; _Xp(0,1)=206.506; _Xp(0,2)=-1.83746; _Xp(0,3)=0.532051;
  _Xp(1,0)=50.4344; _Xp(1,1)=181.449; _Xp(1,2)=-1.84017; _Xp(1,3)=0.227005;
  _Xp(2,0)=45.4587; _Xp(2,1)=148.387; _Xp(2,2)= -1.85561; _Xp(2,3)=0.0834351;
  _Xp(3,0)= 68.6357 ; _Xp(3,1)= 182.511; _Xp(3,2)=-1.81659 ; _Xp(3,3)=0.233932 ;
  _Xp(4,0)= 52.3116 ; _Xp(4,1)=264.735 ; _Xp(4,2)=  -1.87646 ; _Xp(4,3)= 0.398346;
  _Xp(5,0)= 44.4737 ; _Xp(5,1)=171.398 ; _Xp(5,2)=-1.90638 ; _Xp(5,3)= 0.202469 ;
  _Xp(6,0)= 49.3264  ; _Xp(6,1)= 191.475; _Xp(6,2)= -1.90369; _Xp(6,3)=0.212616 ;
  _Xp(7,0)= 55.4727 ; _Xp(7,1)=149.378 ; _Xp(7,2)=-1.90636 ; _Xp(7,3)= 0.134079 ;
  _Xp(8,0)= 97.6935 ; _Xp(8,1)= 169.627 ; _Xp(8,2)=-1.86367  ; _Xp(8,3)=0.279007 ;
  _Xp(9,0)= 56.6327 ; _Xp(9,1)=137.394 ; _Xp(9,2)= -1.84868; _Xp(9,3)= 0.10701 ;
  _Xp(10,0)=  76.7403 ; _Xp(10,1)=116.307  ; _Xp(10,2)= -1.81574; _Xp(10,3)= 0.00523376;

  */

  mat X(35,6);
 
  for (int i=0; i<35; i++)
    for(int j=0; j<6; j++)
      X(i,j)=0;

  X(0,0)= 52.4861; X(0,1) = 257.923; X(0,2) = 124.861; X(0,3) = -0.770874;       
  X(1,0)= 114.339; X(1,1)= 317.485 ; X(1,2)= 93.3878 ; X(1,3)= 34.8529;     
  X(2,0)= 70.6485; X(2,1)= 290.536 ; X(2,2)= 116.485 ; X(2,3)= -4.63501;   
  X(3,0)= 81.6085; X(3,1)= 282.534 ; X(3,2)= 116.085 ; X(3,3)= -4.65851;    
  X(4,0)= 61.1525; X(4,1)= 274.6 ; X(4,2)= 121.525 ;   X(4,3)= -4.00177;    
  X(5,0)= 57.4335; X(5,1)= 289.523 ; X(5,2)= 124.335 ; X(5,3)= -4.76715;   
  X(6,0)= 119.898; X(6,1)= 308.48 ; X(6,2)= 118.978 ;  X(6,3)= 44.8038 ;   
  X(7,0)= 82.7006; X(7,1)= 307.981 ; X(7,2)= 127.006 ; X(7,3)=-60.1898;    
  X(8,0)= 64.2541; X(8,1)= 257.31 ; X(8,2)= 117.68 ;       X(8,3)= -6.12518 ;
  X(9,0)= 122.58;  X(9,1)= 320.963; X(9,2)= 82.4085 ;      X(9,3)= 34.7791; 
  X(10,0)= 82.0365 ; X(10,1)= 289.667 ; X(10,2)= 113.88 ;  X(10,3)= -8.69965;
  X(11,0)= 92.5519 ; X(11,1)= 281.737 ; X(11,2)= 109.434 ; X(11,3)= -7.97546;
  X(12,0)= 54.0814 ; X(12,1)= 257.512 ; X(12,2)= 120.814 ; X(12,3)= -4.88403;
  X(13,0)= 72.737 ; X(13,1) =273.668 ; X(13,2)= 115.844 ;  X(13,3)= -9.32068;
  X(14,0)= 69.0632 ; X(14,1)= 288.597 ; X(14,2)= 116.297 ; X(14,3)= -9.26697;
  X(15,0)= 129.944 ; X(15,1)= 310.777 ; X(15,2)= 100.46 ;  X(15,3)= 22.9645 ;
  X(16,0)= 96.3861 ; X(16,1)= 303.031 ; X(16,2)= 136.855 ; X(16,3)= -49.4986;
  X(17,0)= 75.4449 ; X(17,1)= 257.369 ; X(17,2)= 111.908 ; X(17,3)= 0.58197; 
  X(18,0)= 128.094 ; X(18,1)= 324.424 ; X(18,2)= 55.1404 ; X(18,3)= 34.6057;
  X(19,0)= 92.8351 ; X(19,1)= 289.524 ; X(19,2)= 107.986 ; X(19,3)= -1.42792 ;
  X(20,0)= 103.168 ; X(20,1)= 281.584 ; X(20,2)= 106.158 ; X(20,3)= -1.52344; 
  X(21,0)= 65.9226 ; X(21,1)= 257.682 ; X(21,2)= 118.412 ; X(21,3)= 1.70624; 
  X(22,0)= 83.8475 ; X(22,1)= 273.513 ; X(22,2)= 111.105 ; X(22,3)= -1.54694 ;
  X(23,0)= 80.2008 ; X(23,1)= 288.427 ; X(23,2)= 111.376 ; X(23,3)= -1.70044; 
  X(24,0)= 138.9 ; X(24,1)= 314.007 ; X(24,2)= 89.5651 ;   X(24,3)= 32.3065 ;
  X(25,0)= 110.939 ; X(25,1)= 298.569 ; X(25,2)= 145.526 ; X(25,3)= -44.6188;
  X(26,0)= 86.6825 ; X(26,1)= 258.731 ; X(26,2)= 112.376; X(26,3)= 13.6264 ;
  X(27,0)= 132.65 ; X(27,1)= 327.602 ; X(27,2)= 45.5661;  X(27,3)= 31.7844; 
  X(28,0)= 103.582 ; X(28,1)= 290.388 ; X(28,2)= 107.467; X(28,3)= 8.64258 ;
  X(29,0)= 113.687 ; X(29,1)= 282.499 ; X(29,2)= 105.19;  X(29,3)= 9.15161; 
  X(30,0)= 77.4565 ; X(30,1)= 259.064 ; X(30,2)= 115.339; X(30,3)= 13.8141 ;
  X(31,0)= 94.8209 ; X(31,1)= 274.496 ; X(31,2)= 109.735; X(31,3)= 9.83002; 
  X(32,0)= 91.2858 ; X(32,1)= 289.322 ; X(32,2)= 110.85 ; X(32,3)= 8.95172; 
  X(33,0)= 147.228 ; X(33,1)= 317.462 ; X(33,2)= 83.2735; X(33,3)= 34.5468 ;
  X(34,0)= 125.604 ; X(34,1)= 296.132 ; X(34,2)= 146.649; X(34,3)= -24.3738;

  return (X);
}

int clusteringMethod(mat &_Xp, int imaC, int imaR, nodeptr & nod, glptr & meaninG)
{
  using namespace cluster;

  //Declaration
  clusFuncs node;
  int  dataDimension, i, l, sons, trakP; // sVector;
  float *pdistEuclidean;
  double ind[12], piR, nz, sVector;
  nodeptr aux = NULL;
  listptr lc = NULL ;
  listptr *alc;
  treeb c;

      // Some initialization variables
      trakP = _Xp.size1(); //number of tracking and moving points
      dataDimension = _Xp.size2()-2;
      //  std::cout << "dataDimension " << dataDimension << std::endl;

      node.matrixInitialization(dataDimension, trakP-1);
      nz  = trakP*trakP;
      alc = (listptr *) calloc(trakP-1, sizeof(listptr));
      mat Xpol (trakP, dataDimension);
      Xpol = _Xp;

    /**********************************************************************************/
      /**** FIRST STEP: Preparing data set and form the initial groups by singlelinkage */
      /**********************************************************************************/

      //Transform rectangular coodinates to polar
      rectangularToPolar(_Xp, Xpol);
      //   _Xp = Xpol;
      node.minValueVel(Xpol);
      node.maxValueVel(Xpol);
      /* Determine the size of euclidean distance vector
	 Get size of the vector distance by the equation= ((n-1)n)/2  where n=number of tracked features*/
      sVector = (trakP*(trakP-1))/2;

      //  std::cout << " sVector " << sVector << " nz " << nz << std::endl;

      /* Reserve the distance vector just for the tracked features */
      pdistEuclidean = (float *) calloc(sVector, sizeof(float));


      /*Obtain region sizes, 20 regions per coordinate */
      node.calcRegionSizes(dataDimension, imaC, imaR);
      node.calcEuclideanDistance(_Xp, dataDimension, pdistEuclidean, trakP); // rectan coord
      node.singleLinkage(Xpol, pdistEuclidean, sVector, trakP, dataDimension, nod);
      //     std::cout << "End First Step " << std::endl;
      /* Para debug */
      //   node.printRegions(dataDimension);
      // c.printree(nod);
      /* Fin zona debug */

      //   _Xp = Xpol;


 /**********************************************************************************/
 /*** SECOND STEP: Validity Step. Calcule the NFA for each node in the binary tree */
 /**********************************************************************************/

      // Start analysis for each group in the binary tree
  for (i=0; i < trakP-1; i++ ) {
	aux = findnodeaddress(nod, i+trakP);
	sons = getnumpointspernode(aux);
	getnodepoints(aux, lc);
	if (lc!=NULL){
	  alc[i] = findheadlist(lc);
	  //	 printpoints(alc[i]);
	  lc=NULL;
	}
	  //initialisation
	ind[8] = 2.0;
	// Analyse of each point in the group
	for (l=0; l<sons; l++){
	  node.bestRegion(alc[i],l,0,imaC,0);  //Find minimal region in X coordinate
	  node.bestRegion(alc[i],l,1,imaR,0);  //Find minimal region in Y coordinate
	  node.bestRegion(alc[i],l,3,360,0);   //Find minimal region in orientation coordinate
	  piR = node.bestRegionVel(alc[i],l,Xpol);   //Find minimal region in Velocity coordinate and density function
	  if (piR < ind[8]){
	    ind[8]=piR;       // proba of region
	    //   std::cout << " piR temporal " << piR << "\n";
	    ind[9]=l;        //centre
	    ind[10]=sons;   // num of points
	    //	     std::cout << " hijos " << sons << "\n";
	    node.saveRegionTemp(ind);
	  }
	} //end group points
	// Calcule NFA per node

       	ind[11] = CARD_REG*nz*binomial(trakP-1, ind[10]-1, ind[8]);

	//	std::cout << " NFA para el nodo " << i+trakP << " " << ind[11] << "\n";
	c.saveFinalRegionNode(ind, aux, meaninG, trakP);

	/*std::cout << " z= " << z[i][0] << " " << z[i][1] << " dist=  " << z[i][2] ;
	  std::cout << "\n"; */

	/*	if (aux != NULL)
	std::cout << " hijos de " << aux->id <<" = " << aux->s1->id << " y " << aux->s2->id << " NFA =  " << aux->nfa;
	std::cout << "\n";

	std::cout << std::endl;
	node.printNZ(i+trakP, aux);
*/
      }  //end binary tree

      if (meaninG != NULL )
	{
	  meaninG = findheadgrouplist(meaninG);
	 // c.printeG(meaninG);
	}
      // std::cout << "End Second Step " << std::endl;

      /**********************************************************************************/
      /***  THIRD STEP : MERGING STEP. Calcule the NFAg of the sibling groups */
      /**********************************************************************************/
    while(meaninG != NULL)
	{
	  node.findNFAg(meaninG, alc, trakP, imaC, imaR);
	  if (meaninG->next == NULL)
	    break;
	  meaninG = meaninG->next;
	}
      // std::cout << "End Third Step " << std::endl;

      /**********************************************************************************/
      /***  FOURTH STEP : Find the e-Meaningful groups. *********************************/
      /**********************************************************************************/
      // meaninG will save the temporal groups e-meaningful
	if (meaninG != NULL)
	   deleteGroupList (meaninG);
      // alc will save the groups resulted of the clustering process
       freeArrayList( alc, trakP);
      // Top to down analysis for each node in the binary tree

      for (i=2*trakP-2, l=0; i>=trakP; i--) {
	//	std::cout << " Grupo Analizado " << i << " valor de l " << l << std::endl;
	aux = findnodeaddress(nod, i);
	if(c.detectMeaningfulGroup(aux, trakP))
	  { //node aux is meaningful, then save the node address
	    c.saveTemporalEMeaninHeadNode(aux, meaninG);
	    //   std::cout << " Numero de group salvado como meaninFul " << i << std::endl;
	    l++;
	  }
      }

      // Find the maximal e-meaningful groups in the list meaninG: some groups could be included in bigger node groups
      if (l>0){
	meaninG = findheadgrouplist(meaninG);
	i=0;
	c.valideNodeGroups(meaninG, i);
	meaninG = findheadgrouplist(meaninG);
	//std::cout << " en i (interno, los q c van)= " << i << "l (e-meanGr)es igual a= " << l << std::endl;

	//     bt=nod;
	//    emg = meaninG;   // estaba comentado
      }
      else
	  {
	  nod = NULL;
	  meaninG = NULL;
	  l=i;
	  }
    //  std::cout << "End of clustering " << std::endl;

	/*
       alc = (listptr *) calloc(l-i, sizeof(listptr));
       l=0;
       do
	 { //node aux is meaningful, then save it
	   getnodepoints(emg->g, lc);
	   if (lc!=NULL){
	     alc[l] = findheadlist(lc);
	     printpoints(alc[l]);
	     lc=NULL;
	     std::cout << " Numero de grupo guardado " << emg->g->id << std::endl;
	     l++;
	   }
	   emg = emg->next;
	 }while(emg != NULL);
	*/
       //  deletree(nod);
      free(pdistEuclidean);
      return l-i; //numero de clusters
    }// end of clustering Method


void probandoMatrices(mat & mA, int b, int c)
{

	for (uint i=0; i<mA.size1();i++)
		{
		for(uint j=0; i<mA.size2(); j++)
			std::cout << " " << mA(i,j);
		}
	std::cout << " " << std::endl;
	std::cout << "Valor mA(3,1)*b " << mA(34,1)*b << std::endl;
	std::cout << "Valor mA(3,1)*c " << mA(34,1)*c << std::endl;
}
} //end namespace cluster

/*
int main()
{
	using namespace cluster;
	using namespace boost_incl;

	mat XP;
	int clusters;
	XP = cluster::auxiliopelicanmurio();
	// cluster::probandoMatrices(XP, 1, 2);
	clusters=clusteringTrackingPoints(XP,640,480);
	std::cout << "Numero de clusters " << clusters << std::endl;
	printf("fin");
	return 0;
}
*/
