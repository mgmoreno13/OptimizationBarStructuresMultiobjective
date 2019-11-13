/*==========================================================================
// //  Author: Carlos Segura, Joel Chacón 
//     Description: 
//
// ===========================================================================*/


#ifndef _PROBLEM_H
#define _PROBLEM_H
#include "OPTBAR/Optimizador.h"



void Bar_Opt(std::vector< double >& F_original, std::vector< int >& X)
{


  for(unsigned int i=1; i<=X.size(); i++) 
    mata[i]=X[i-1];


  std::vector<double> v;


  v=Funcion_fitness(isho, npno,nele, npre,  ncar, ntipo, ngd, nmat,
                 ndim, iwri, isal, nres, ninc, nfam, ncas,
                  nnod, mnod, nprp, ngau, nten, neva, inds, ntot,
                 lnod, inpr, narea, naref, narec, narcon,
                 lres, linc, ntip, nodp, matn, mata, matm,
                 iffi, indf, node, maxa,mhig, leq,
                 coor, cargaa,carp, pres, prop,
                 fuem, fuep, rig, xmas, vmat,
                 gir, girt, tem, fuef, fueb,
                 srm, smm, cArma, cAceRF,
                 cAceRC,cCon, res, xinc,
                 fuec, aslo, desp, fix, reac,
                 veca, vec1, des, xlon, valr,
                 valm, angl, anglg, wcar, vecr,
                 asti, vecp, vecu, listaCata, ubiC,
                 vectores_f, npes, efi_mx, stif, neq,
                 nwkt);//@mg'

  F_original[0]=v[0];//with normalization
  F_original[1]=v[1];//with normalization


}


#endif
