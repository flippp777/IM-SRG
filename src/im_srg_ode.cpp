#include "Hamiltonian.hpp"
#include "im_srg_ode.hpp"
#include <gsl/gsl_errno.h>
#include <vector>
#include <iostream>

int im_srg_ode(double s, const double y[], double dyds[], void* H) 
{
    if (!y || !dyds || !H) {
        return GSL_EFAULT;
    }

    Hamiltonian* h = static_cast<Hamiltonian*>(H);
    
    try {
        // Set E, f and gamma from input vector
        h->set_elements(y);

        // Calculate generators and flow equations
        h->set_generators();
        h->set_flow_eq();

        // Get flow equations into output vector
        std::vector<double> flow_eqs = h->get_full_eqs();
        
        // Safety check
        if(flow_eqs.empty()) {
            return GSL_EFAULT;
        }

        std::copy(flow_eqs.begin(), flow_eqs.end(), dyds);
        return GSL_SUCCESS;
    }
    catch(const std::exception& e) {
        std::cerr << "Error in ODE function: " << e.what() << std::endl;
        return GSL_EFAULT;
    }
}
