/*
 * auxfuncs.cpp
 *
 *  Created on: Jan 29, 2013
 *      Author: luzdora
 */
#include "auxfuncs.hpp"


namespace cluster{

listptr findheadlist(listptr _lc)
{
	if (_lc->prev != NULL)
		return findheadlist(_lc->prev);
	else
		return (_lc);
}

glptr findheadgrouplist(glptr _gl)
{
	if (_gl->prev != NULL)
		return findheadgrouplist(_gl->prev);
	else
		return (_gl);

}

void printpoints(listptr &_lc)
{
	std::cout << "puntos= " << _lc->v[0] <<"  " << _lc->v[1] << "  " << _lc->v[2] << " " << _lc->v[3] << "\n";
	if (_lc->next != NULL)
		printpoints(_lc->next);
}

nodeptr findnodeaddress(nodeptr &_p, int _id)
{
	nodeptr aux;
	if (_p != NULL)
	{
		if (_p->id == _id)
			return(_p);
		else
		{
			aux = findnodeaddress(_p->s1, _id);
			if (aux != NULL)
				return(aux);
			aux = findnodeaddress(_p->s2, _id);
			if (aux != NULL)
				return(aux);
			else
				return NULL;
		}
	}
	else
		return NULL;
}

int getnumpointspernode(nodeptr &_p)
{
	if (_p != NULL)
	{
		if (_p->s1 == NULL && _p->s2 == NULL)
			return 1;
		else
			return ( getnumpointspernode(_p->s1) + getnumpointspernode(_p->s2));

	}
	else
		return 0;
}

void getnodepoints(nodeptr &_p, listptr &_lc)
{
	if(_p != NULL)
	{
		if ( _p->s1==NULL && _p->s2==NULL)
		{
			if (_lc == NULL)
			{
				_lc = (listptr) calloc(1,sizeof(listptr));
				_lc->v = _p->v;
				_lc->prev = NULL;
				_lc->next = NULL;
			}
			else
			{
				_lc->next =(listptr) calloc(1,sizeof(listptr));
				_lc->next->v = _p->v;
				_lc->next->prev = _lc;
				_lc->next->next = NULL;
				_lc = _lc->next;
			}
		}
		else
		{
			getnodepoints(_p->s1, _lc);
			getnodepoints(_p->s2, _lc);
		}
	}
}


int minimalValue(float* eucdis, int ndis)
{
	int min_pos=0, i;

	for(i = 1; i < ndis; i++)
	{
		if(*(eucdis+i) < *(eucdis+min_pos))
			min_pos = i;
	}

	return min_pos;
}

double factorial (double num)
{
	if (num==1.0 || num==0.0)
		return 1.0;
	return factorial(num-1.0)*num;
}

/*
double binomial(int m, int k, double p)
{
  double bino=0.0, aux=0.0, factom;
  int j;

  factom = factorial((double)m);
  for (j=k; j<=m; j++)
    {
      aux = (factom/(factorial((double)(m-j))*factorial((double)j)))*pow(p,j)*pow((1.0-p),(m-j));
      bino = bino + aux;
    }
  //     std::cout << "binomial =  " << bino << " proba = "<< p <<std::endl ;
  return bino;
}


double trinomi(int Mp, int k1, int k2, double p1, double p2)
{
  double tripo = 0.0, aux, pro, expo, factoMp, p1i, facti;
  int i, j;
  factoMp = factorial((double)Mp);
  pro = (double) (1.0 - p1 - p2);

  for (i = k1; i<=Mp; i++){
    p1i = pow(p1,i);
    facti = factorial((double)i);

    for(j = k2; j<=Mp-k1; j++){
      expo =(double) (Mp-i-j);
      if ( expo >= 0){
	aux = (factoMp/( factorial(expo)*factorial((double)j)*facti))*(p1i*pow(p2,j)*pow(pro,expo));
	tripo = tripo + aux;
	// debuging
      }
    }
    //  std::cout <<"factoMp= " << factoMp << " expo= " << expo << " factoi= " << facti << " i= "<< i <<" powP2j "<< pow(p2,j) << "aux = " << aux << " trinomial =  " << tripo <<std::endl ;
  }
    std::cout <<"Mp= " << Mp << " k1= " << k1 << " k2= " << k2 << " proba 1= "<< p1 <<" proba 2= "<< p2 << "facMp = " << factoMp << " trinomial =  " << tripo <<std::endl ;
  return tripo;
}

 */


double trinomi(int Mp, int k1, int k2, double p1, double p2)
{
	double tripo = 0.0, aux, pro, expo, factoMp, p1i;
	//double facti;
	int i, j, k;
	// factoMp = factorial((double)Mp);
	pro = (double) (1.0 - p1 - p2);

	for (i = k1; i<=Mp; i++){
		p1i = pow(p1,i);
		//  facti = factorial((double)i);

		for(j = k2; j<=Mp-k1; j++){
			expo =(double) (Mp-i-j);
			if ( expo >= 0){
				factoMp = 1.0;

				for (k=1; k<=i ; k++)
					factoMp *= ((1.0*(expo+k))/k) ;
				for (k=1; k<=j ; k++)
					factoMp *= ((1.0*(Mp-j+k))/k) ;

				aux = (factoMp)*(p1i*pow(p2,j)*pow(pro,expo));
				tripo = tripo + aux;
				// debuging
			}
		}
		//  std::cout <<"factoMp= " << factoMp << " expo= " << expo << " factoi= " << facti << " i= "<< i <<" powP2j "<< pow(p2,j) << "aux = " << aux << " trinomial =  " << tripo <<std::endl ;
	}
	//   std::cout <<"Mp= " << Mp << " k1= " << k1 << " k2= " << k2 << " proba 1= "<< p1 <<" proba 2= "<< p2 << "facMp = " << factoMp << " trinomial =  " << tripo <<std::endl ;
	return tripo;
}



double binomial(int m, int k, double p)
{
	double bino=0.0, aux=0.0, factom;
	int j, i;

	for (j=k; j<=m; j++)
	{
		factom = 1.0;
		for (i = 1 ; i <=j; i++)
			factom *= ((1.0*(m-j+i))/i) ;

		aux = (factom)*pow(p,j)*pow((1.0-p),(m-j));
		bino = bino + aux;
	}
	//     std::cout << "binomial =  " << bino << " proba = "<< p <<std::endl ;
	return bino;
}



// looking for intersection in one dimension
float intersectRegion(float sR1, float iR1, float sR2, float iR2)
{
	float inteR;

	// std::cout << " sup1, sup2 " << sR1, sR2 << " inf1, inf2 " << iR1, iR2 << std::endl;

	if ((sR1 <= sR2) && (sR1>=iR2))  //case 1
	{
		if (iR2>iR1)
			inteR = sR1-iR2;
		else  // R2 has all the region R1
			inteR=sR1-iR1;
	}
	else if ((sR2 <= sR1) && (sR2 >= iR1)){  // case 2
		if (iR1>iR2)
			inteR= sR2-iR1;
		else // R1 has all the region R2
			inteR= sR2-iR2;
	}
	else
	{ //there is not an intersection between regions
		inteR=0;
		//    std::cout "There is not an intersection between regions desde la funcion " << std::endl;
	}
	return inteR;
}



//Count the points in the intersection
void pointsInIntersectionRegion(glptr & eG, listptr *lc, int & p1, int & p2, int nP)
{
	int i;
	listptr aux;

	// group 1
	i = (eG->g->s1->id)-nP; // index in the linked list of first child
	p1=0;
	aux =  lc[i];

	while (aux != NULL ){ //each point in the group 1
		if  ( (aux->v[0] >= eG->g->s2->v[1]) && (aux->v[0] <= eG->g->s2->v[0]) && (aux->v[1] >= eG->g->s2->v[3]) && (aux->v[1] <= eG->g->s2->v[2] ) && (aux->v[2] >= eG->g->s2->v[5]) && (aux->v[2] <= eG->g->s2->v[4]) && (aux->v[3] >= eG->g->s2->v[7]) && (aux->v[3] <= eG->g->s2->v[6]) ) // there is a point in the intersection region (X-Y-Vel-Teta dimension)
			p1++;
		aux = aux->next;
	}

	// group 2
	i = (eG->g->s2->id)-nP; // index in the linked list of second child
	aux =  lc[i];

	p2=0;

	while (aux != NULL ){ //each point in the group 2
		if  ( (aux->v[0] >= eG->g->s1->v[1]) && (aux->v[0] <= eG->g->s1->v[0]) && (aux->v[1] >= eG->g->s1->v[3]) && (aux->v[1] <= eG->g->s1->v[2] ) && (aux->v[2] >= eG->g->s1->v[5]) && (aux->v[2] <= eG->g->s1->v[4]) && (aux->v[3] >= eG->g->s1->v[7]) && (aux->v[3] <= eG->g->s1->v[6]) ) // there is a point in the intersection region (X-Y-Vel-Teta dimension)
			p2++;
		aux = aux->next;
	}
	//  std::cout << "puntos en la intersection: " << p1 << " y " << p2 << std::endl;

}

int indivisible(nodeptr &_p, int nP)
{
	// Check if childs are groups
	if (_p->s1->id >= nP && _p->s2->id >= nP)
	{
		if ((_p->nfa < _p->nfag) && (_p->nfag < 0.1)) // Second test:  the group is indivisible
		{
			return 1;
		}
		else
			return 0; // Group divisible
	}
	else
		return 0; // Childs are one point or both points
}

int checkMaximalEGroupsChilds(nodeptr p, nodeptr ch, int nP)
{
	if (indivisible(ch, nP)){
		if (p->nfa < ch->nfa)  //3rd test: more meaningful than all its indivisible childs
			return checkMaximalEGroupsChilds(p,ch->s1, nP)*checkMaximalEGroupsChilds(p,ch->s2, nP);
		else
			return 0;
	}
	else
	{
		if (ch->s2->id < nP)
			return 1;
		else if (ch->s1->id < nP)
			return checkMaximalEGroupsChilds(p,ch->s2, nP);
		else
			return checkMaximalEGroupsChilds(p,ch->s1, nP)*checkMaximalEGroupsChilds(p,ch->s2, nP);
	}
}

float *memorySpace(double numfeat){

	float *vectorf;

	vectorf = (float *) calloc (numfeat, sizeof(float));

	return vectorf;
}

float **matrixReserve(int nr, int nc)
{
	float **m;
	int a;
	m = (float **) calloc(nr, sizeof(float*));
	m[0] = (float *) calloc((nr*nc), sizeof(float));
	for(a = 1; a < nr; a++) m[a] = m[a-1]+nc;
	return m;
}

double **matrixDoubleReserve(int nr, int nc)
{
	double **m;
	int a;
	m = (double **) calloc(nr, sizeof(double*));
	m[0] = (double *) calloc((nr*nc), sizeof(double));
	for(a = 1; a < nr; a++) m[a] = m[a-1]+nc;
	return m;
}

void deleteGroupList (glptr & _eG)
{
	_eG = findheadgrouplist(_eG);

	while(_eG != NULL)
	{
		if (_eG->next != NULL){
			_eG = _eG->next;
			free(_eG->prev);
		}
		else
		{
			free(_eG);
			_eG = NULL;
		}
	}
}

void deletree(nodeptr &_p)
{
	if (_p != NULL)
	{
		if (_p->s1 != NULL)
			deletree(_p->s1);
		if (_p->s2 != NULL)
			deletree(_p->s2);
		free(_p->v);
		delete(_p);
	}
}


void freeMatrix(float **m)
{
	free(m[0]);
	free(m);
}

void freeArrayList( listptr *_alc, int nP)
{
	int i;
	for (i=0; i < nP-1; i++){
		while(_alc[i] != NULL)
		{
			if (_alc[i]->next != NULL){
				_alc[i] = _alc[i]->next;
				free(_alc[i]->prev);
			}
			else
			{
				free(_alc[i]);
				_alc[i] = NULL;
			}
		}
	}
	free(_alc);

}

} //end namespace
