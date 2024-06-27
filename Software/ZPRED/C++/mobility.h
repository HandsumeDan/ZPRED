# include "solution.h"
# include "zpredUtils.h"

# ifndef mobility_H
# define mobility_H

	typedef vector<double> state_type;

	struct push_back_state_and_time
	{
	    vector<state_type>& m_states;
	    vector<double>& m_times;

	    push_back_state_and_time(vector<state_type> &states,vector<double> &times ) : m_states(states),m_times(times){}

	    void operator()(const state_type &x ,double t)
	    	{m_states.push_back(x);
		   m_times.push_back(t);}
	};

	struct Error
	{double absolute;double relative;double a_x;double a_dxdt;};

	typedef runge_kutta_dopri5<state_type> error_stepper_type;

	typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;

	// MOBILITY Functions
	double two_ion_problem(double ZP,double a,Solution S,Error E,bool VERBOSE);
	double three_ion_problem(double ZP,double a,Solution S,Error E,bool VERBOSE);
	double four_ion_problem(double ZP,double a,Solution S,Error E,bool VERBOSE);
	double five_ion_problem(double ZP,double a,Solution S,Error E,bool VERBOSE);
	double six_ion_problem(double ZP,double a,Solution S,Error E,bool VERBOSE);
	double seven_ion_problem(double ZP,double a,Solution S,Error E,bool VERBOSE);
	double eight_ion_problem(double ZP,double a,Solution S,Error E,bool VERBOSE);
	double nine_ion_problem(double ZP,double a,Solution S,Error E,bool VERBOSE);
	double calc_PBE_coefficient(double Zeta,double &r0,double rf,double &dr,Solution S,Error E,bool VERBOSE);
	//
	void homogeneous_form_2ion(const state_type &y,state_type &dy,const double r);
	void homogeneous_form_3ion(const state_type &y,state_type &dy,const double r);
	void homogeneous_form_4ion(const state_type &y,state_type &dy,const double r);
	void homogeneous_form_5ion(const state_type &y,state_type &dy,const double r);
	void homogeneous_form_6ion(const state_type &y,state_type &dy,const double r);
	void homogeneous_form_7ion(const state_type &y,state_type &dy,const double r);
	void homogeneous_form_8ion(const state_type &y,state_type &dy,const double r);
	void homogeneous_form_9ion(const state_type &y,state_type &dy,const double r);
	//
	void inhomogeneous_Prob1_2ion(const state_type &y,state_type &dy,const double r);
	void inhomogeneous_Prob1_3ion(const state_type &y,state_type &dy,const double r);
	void inhomogeneous_Prob1_4ion(const state_type &y,state_type &dy,const double r);
	void inhomogeneous_Prob1_5ion(const state_type &y,state_type &dy,const double r);
	void inhomogeneous_Prob1_6ion(const state_type &y,state_type &dy,const double r);
	void inhomogeneous_Prob1_7ion(const state_type &y,state_type &dy,const double r);
	void inhomogeneous_Prob1_8ion(const state_type &y,state_type &dy,const double r);
	void inhomogeneous_Prob1_9ion(const state_type &y,state_type &dy,const double r);
	//
	void inhomogeneous_Prob2_2ion(const state_type &y,state_type &dy,const double r);
	void inhomogeneous_Prob2_3ion(const state_type &y,state_type &dy,const double r);
	void inhomogeneous_Prob2_4ion(const state_type &y,state_type &dy,const double r);
	void inhomogeneous_Prob2_5ion(const state_type &y,state_type &dy,const double r);
	void inhomogeneous_Prob2_6ion(const state_type &y,state_type &dy,const double r);
	void inhomogeneous_Prob2_7ion(const state_type &y,state_type &dy,const double r);
	void inhomogeneous_Prob2_8ion(const state_type &y,state_type &dy,const double r);
	void inhomogeneous_Prob2_9ion(const state_type &y,state_type &dy,const double r);
	//
	void poisson_boltzmann(const state_type &y,state_type &dy,const double r);
	void initialize_PBE(state_type &y,double r0,double C,Solution S);

# endif
