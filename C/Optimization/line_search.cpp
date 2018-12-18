#include "line_search.h"

using namespace std;

double line_search(Eigen::VectorXd& x, const Eigen::VectorXd& d, double step_size, Objective& f, double cur_energy) {
	double old_energy;
  if (cur_energy >= 0)
  {
    old_energy = cur_energy;
  }
  else
  {
    old_energy = f.obj(x); // no energy was given -> need to compute the current energy
  }
  double new_energy = old_energy;
  int cur_iter = 0; int MAX_STEP_SIZE_ITER = 22;

  while (new_energy >= old_energy && cur_iter < MAX_STEP_SIZE_ITER)
  {
    Eigen::VectorXd new_x = x + step_size * d;

    double cur_e = f.obj(new_x);
    
    //cout << "cur_e = " << cur_e << endl;
    if (cur_e >= old_energy)
    {
      step_size /= 2;
      //cout << "step_size = " << step_size << endl;
    }
    else
    {
      x = new_x;
      new_energy = cur_e;
    }
    cur_iter++;
  }
  if (cur_iter < MAX_STEP_SIZE_ITER) {
  	//cout << "ls success!" << endl;
  } else {
  	cout << "ls failure, is it a local minimum?" << endl;
  	return old_energy;
  }
  //cout << "step = " << step_size << endl;
  return new_energy;
}



double exact_l2_merit_lineserach(Eigen::VectorXd& x, const Eigen::VectorXd& d, double step_size, Objective& f, Constraints& constraints,
				const double& merit_penalty,
				double cur_energy = -1) {

	ExactL2MeritObjective meritObj(f,constraints,merit_penalty);
	return line_search(x,d,step_size, meritObj, cur_energy);
}