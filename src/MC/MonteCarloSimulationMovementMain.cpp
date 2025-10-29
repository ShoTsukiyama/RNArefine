
#include "MonteCarloSimulation.hpp"

double MonteCarloSimulation::getrandomnewparam(double l, double r){
	return l+(r-l)*Random();
}

int MonteCarloSimulation::getrandval(vector<double> tmov,double trandnum)
{
	int i;
	for(i=0;i<tmov.size()-1;i++)
	{
		if(trandnum>=tmov[i] && trandnum<tmov[i+1])
		{
			return i;
		}
	}
	return 0;
}

bool MonteCarloSimulation::metropolisUpdate(double oldenergy, double newenergy)
{
	double v, r;

	if (temperature <= 0) {
        return (newenergy < oldenergy);
    }

	if (newenergy < oldenergy) {
        return true;
    }else {
        double v = exp(-(newenergy - oldenergy) / temperature);
        double r = Random();
        return (r < v); 
    }
}

void MonteCarloSimulation::setsimulationparams(double _temperature, int _lsweepnum)
{
	temperature=_temperature;
	lsweepnum=_lsweepnum;
}

double MonteCarloSimulation::gettemperature()
{
	return temperature;
}


bool MonteCarloSimulation::calcmovementmain(double oldenergy, double &newenergy, int &moveidx)
{
	int flag;
	double movetype=Random();
	vector<double> tmov;

	//tmov={0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0};
	tmov={0, 0.2, 0.4, 0.6, 0.8, 1.0};
	moveidx=getrandval(tmov,movetype);

	//cout<<"=========="<<endl;
	//cout<<moveidx<<endl;
	//cout<<"=========="<<endl;

	//flag=fmovementloopfrag(oldenergy, newenergy);
	//flag=fmovementloopfrag(oldenergy, newenergy);
	//flag=fmovementfrag(oldenergy, newenergy);
    //flag=fmovementbpfrag(oldenergy, newenergy);
        //flag=fmovementstackfrag(oldenergy, newenergy);
	//return flag;

	switch(moveidx)
	{
		case 0:
			flag=fmovementang(oldenergy, newenergy);
			break;
		case 1:
			flag=fmovementleng(oldenergy, newenergy);
			break;
		case 2:
			flag=fmovementtor(oldenergy, newenergy);
			break;
		//case 3:
			//flag=fmovementcombtor(oldenergy, newenergy);
			//break;
		case 3: 
			flag=fmovementLMP(oldenergy, newenergy);
			break;
		//case 5:
			//flag=fmovementrot(oldenergy, newenergy);
			//break;

		//case 4:
			//flag=fmovementstackfrag(oldenergy, newenergy);
			//break;
		//case 5:
			//flag=fmovementbpfrag(oldenergy, newenergy);
			//break;

		//case 5:
			//flag=fmovementloopfrag(oldenergy, newenergy);
			//break;
		case 4:
			flag=fmovementfrag(oldenergy, newenergy);
			break;
	}
	return flag;
}



























































