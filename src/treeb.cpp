/*
 * treeb.cpp
 *
 *  Created on: Jan 29, 2013
 *      Author: luzdora
 */

#include "../include/cluster/clustersFuncs.hpp"

namespace cluster {

    //using namespace cluster;

int treeb::isChildNode(nodeptr nodP, nodeptr nodCh)
{
  if (nodP->s1 != NULL && nodP->s2 != NULL)
    {
      if ((nodP->s1 == nodCh) || (nodP->s2 == nodCh))
	{
	  return 1;
	}
      else
	return isChildNode(nodP->s1, nodCh)+isChildNode(nodP->s2, nodCh);
    }
  else
    return 0;

}

int treeb::valideNodeGroups(glptr & _eG, int &i)
{
  glptr aux1, aux2;

  if (_eG->next != NULL)
    {
      aux1 = _eG->next;
      while (aux1 != NULL)
	{
	  if (isChildNode( _eG->g, aux1->g))
	    {
	      i++; //one group eliminated
	      aux1->prev->next = aux1->next;
	      if (aux1->next != NULL)
		  aux1->next->prev = aux1->prev;
	      aux2 = aux1->next;
	      free(aux1);
	      aux1 = aux2;
	    }
	  else
	    aux1 = aux1->next;
	} //end of while
      if (_eG->next != NULL)
	return valideNodeGroups(_eG->next, i);
      else
	return 0;
    }
  else
    return 0;
}

void treeb::saveTemporalEMeaninHeadNode(nodeptr &_aux, glptr &_eG)
{
  if (_eG==NULL) // First group listed
    {
      _eG = (glptr) calloc(1,sizeof(glptr));
      _eG->g = _aux;
      _eG->prev = NULL;
      _eG->next = NULL;
    }
  else
    {
      _eG->next =(glptr) calloc(1,sizeof(glptr));
      _eG->next->g = _aux;
      _eG->next->prev = _eG;
      _eG->next->next = NULL;
      _eG = _eG->next;
    }
}
void treeb::saveFinalRegionNode(double *ind, nodeptr &_aux, glptr &_eG, int nP)
{
  _aux->v = (float *)calloc(8, sizeof(float));
  for (int i=0; i<8; i++) {
    _aux->v[i]=ind[i];
    //    std::cout << " Coord: " << _aux->v[i] ;
  }
  //  std::cout << " Group " << _aux->id << " number of points " << ind[10] <<std::endl;
  _aux->p = ind[8]; // probability of the region
  _aux->n = ind[10]; // number of points in the group
  _aux->nfa = ind[11];
  if ((_aux->nfa < 0.1) && (_aux->s1->id >= nP) && (_aux->s2->id >= nP))
    {
      // E-meaningful group with sibling groups as childs
      saveTemporalEMeaninHeadNode(_aux,_eG);
      //     std::cout << "Grupo salvado por hijos y significancia " << _aux->id <<std::endl;
    }
}



int treeb::detectMeaningfulGroup(nodeptr &_p, int nP)
{
  if (_p->nfa < 0.1) // First test:  the group is e-Meaningful
    {
      // Second test:  the group is indivisible
      if ( indivisible(_p, nP) ){
	return checkMaximalEGroupsChilds(_p, _p->s1, nP)*checkMaximalEGroupsChilds(_p, _p->s2, nP);
      }
      else
	return 0; // group divisible
    }
  else
    return 0; // Not e-Meaningful
}


void treeb::inserte( mat & _v , float _dist, int _id1, int _id2, int _id3,  int _n, int _nP, nodeptr &_p, stack &_lt)
{
  nodeptr aux;

  if ((_id1<_nP) && (_id2<_nP))   // two points
    {
      if (_p != NULL)
	_lt.push(_p);
      // create a father node
      _p = new node;
      _p->dist =_dist;
      _p->nfa = 1.0;
      _p->nfag= 1.0;
      _p->n = 2;  // father node with two points
      _p->id = _id3;
      _p->v = NULL;
      _p->s1 = new node;
      _p->s2 = new node;

      if (_id1 < _id2)
	{
	  _p->s1->id = _id1;
	  _p->s2->id = _id2;
	}
      else
	{
	  _p->s1->id = _id2;
	  _p->s2->id = _id1;
	}
      _p->s1->s1 = NULL;
      _p->s1->s2 = NULL;
      _p->s1->dist = 0;
      _p->s1->nfa = 1.0;
      _p->s1->nfag= 1.0;
      _p->s1->n = 1; // first point
      _p->s2->s1 = NULL;
      _p->s2->s2 = NULL;
      _p->s2->dist = 0;
      _p->s2->nfa = 1.0;
      _p->s2->nfag= 1.0;
      _p->s2->n = 1; // second point
      _p->s1->v =(float *) calloc (_n, sizeof(float));
      _p->s2->v =(float *) calloc (_n, sizeof(float));
      for (int i=0; i<_n ; i++){
	_p->s1->v[i] = _v(_p->s1->id , i);
	_p->s2->v[i] = _v(_p->s2->id , i);

      }
    }
  else if ((_id1<_nP) || (_id2<_nP))  // one group and one point
    {
      _lt.push(_p);
      _p = new node;
      _p->dist = _dist;
      _p->nfa = 1.0;
      _p->nfag= 1.0;
      // _p->n = _n; unknown yet
      _p->id = _id3;
      _p->v = NULL;
      _p->s1 = new node;

      if (_id1 < _id2)
	{
	  _p->s1->id = _id1;
	  aux = _lt.getele(_id2);
	  _lt.pop(_id2);
	}
      else
	{
	  _p->s1->id = _id2;
	  aux = _lt.getele(_id1);
	  _lt.pop(_id1);
	}

      _p->s1->s1 = NULL;
      _p->s1->s2 = NULL;
      _p->s1->dist = 0;
      _p->s1->nfa = 1.0;
      _p->s1->nfag= 1.0;
      _p->s1->n = 1; // first point
      _p->s1->v = (float *)calloc (_n, sizeof(float));

      for (int i=0; i<_n ; i++)
	_p->s1->v[i] = _v(_p->s1->id , i);
      if (aux != NULL)
	_p->s2 = aux;
    }
  else   // two groups
    {
      _lt.push(_p);
      _p = new node;
      _p->dist = _dist;
      _p->nfa = 1.0;
      _p->nfag= 1.0;
      //  _p->n = _n; unknown yet
      _p->id = _id3;
      _p->v = NULL;
      if (_id1 < _id2 )
	{
	  aux = _lt.getele(_id1);
	  _lt.pop(_id1);
	  if (aux != NULL)
	    _p->s1 = aux;
	  aux = _lt.getele(_id2);
	  _lt.pop(_id2);
	  if (aux != NULL)
	    _p->s2 = aux;
	}
      else
	{
	  aux =_lt.getele(_id2);
	  _lt.pop(_id2);
	  if (aux != NULL)
	    _p->s1 = aux;
	  aux = _lt.getele(_id1);
	  _lt.pop(_id1);
	  if (aux != NULL)
	    _p->s2 = aux;
	}
    }
} //end function


/**********************************************/
/************  Print Functions  **************/
/********************************************/

void treeb::printree(nodeptr &_p)
{
  std::cout << _p->id << " " ;
  if (_p->s1 != NULL)
    printree(_p->s1);
  if (_p->s2 != NULL)
    printree(_p->s2);
}


void treeb::printeG(glptr &_eG)
{
  std::cout << "grupo = " << _eG->g->id << " hijos " << _eG->g->n << " : " <<_eG->g->s1->id << " y " << " " << _eG->g->s2->id << "\n";
  if (_eG->next != NULL)
    printeG(_eG->next);
}

} //end cluster
