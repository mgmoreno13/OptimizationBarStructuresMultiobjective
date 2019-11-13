#ifndef __INDIVIDUAL_H_
#define __INDIVIDUAL_H_

#include "global.h"
#include "problem.h"


class CIndividual{
public:
	CIndividual();
	virtual ~CIndividual();

	std::vector <int> x_var;//@opt_bar
	std::vector <double> y_obj;
	std::vector <double> y_obj_original;
	std::vector<double> obj_max; //objs without normalization
	std::vector<double> obj_min; //objs without normalization
	std::vector <CIndividual *> ptr_dominate;
	int    rank;
	std::vector <double> valores;//@opt_bar
	double nearest_variable_distance;//will be double because distance is normalized
	double neares_objective_distance;
	double times_dominated;
	void   rnd_init();
	void   obj_eval();
	void   show_objective();
	void   show_variable();

    bool   operator<(const CIndividual &ind2);
    bool   operator<<(const CIndividual &ind2);
    bool   operator==(const CIndividual &ind2);
    void   operator=(const CIndividual &ind2);
};

CIndividual::CIndividual()
{
    x_var = std::vector<int>(nvar, 0);
    y_obj = std::vector<double>(nobj, 0); //save objs normalization + objs
    y_obj_original = std::vector<double>(nobj, 0); //save objs normalization + objs
	rank = 0;
}

CIndividual::~CIndividual()
{

}
void CIndividual::rnd_init()
{
    for(int n=0;n<nvar;n++){
        x_var[n] = vlowBound[n] + (int)(rnd_uni(&rnd_uni_init)*(vuppBound[n] - vlowBound[n]));    
    }

}

void CIndividual::obj_eval()
{


    if(!strcmp("bar", strTestInstance)) Bar_Opt(y_obj_original, x_var);


}

void CIndividual::show_objective()
{
    for(int n=0; n<nobj; n++)
		printf("%f ",y_obj[n]);
	printf("\n");
}

void CIndividual::show_variable()
{
    for(int n=0; n<nvar; n++)
		printf("%f ",x_var[n]);
	printf("\n");
}

void CIndividual::operator=(const CIndividual &ind2)
{
	x_var = ind2.x_var;
	y_obj = ind2.y_obj;
	y_obj_original = ind2.y_obj_original;
	rank  = ind2.rank;
}

/*bool CIndividual::operator<(const CIndividual &ind2)
{
	bool dominated = true;
    for(int n=0; n<nobj; n++)
	{
		if(ind2.y_obj[n]<y_obj[n]) return false;
	}
	if(ind2.y_obj==y_obj) return false;
	return dominated;
}*/

bool CIndividual::operator<(const CIndividual &ind2)
{
	bool  similar = true;
    for(int n=0; n<nobj; n++)
	{
		if(ind2.y_obj[n]<y_obj[n]) return false;
	        if(fabs(ind2.y_obj[n]-y_obj[n])>1e-5) similar=false;
	}
	if(similar) return false;
	return true;
}


bool CIndividual::operator<<(const CIndividual &ind2)
{
	bool dominated = true;
    for(int n=0; n<nobj; n++)
	{
		if(ind2.y_obj[n]<y_obj[n]  - 0.0001) return false;
	}
	if(ind2.y_obj==y_obj) return false;
	return dominated;
}

bool CIndividual::operator==(const CIndividual &ind2)
{
	if(ind2.y_obj==y_obj) return true;
	else return false;
}
#endif

