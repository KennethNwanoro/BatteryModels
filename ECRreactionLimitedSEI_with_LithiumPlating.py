

# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 21:50:56 2021

@author: nwanoro
"""
# Class for reaction limited SEI growth with Li plating

#

import pybamm

from .base_sei import BaseModel





class ECReactionLimitedSEI_With_Lithiumplating(BaseModel):

    """

    Class for reaction limited SEI growth. This model assumes the "inner"

    SEI layer is of zero thickness and only models the "outer" SEI layer.



    Parameters

    ----------

    param : parameter class

        The parameters to use for this submodel

    domain : str

        The domain of the model either 'Negative' or 'Positive'



    **Extends:** :class:`pybamm.sei.BaseModel`

    """



    def __init__(self, param, domain):

        super().__init__(param, domain)



    def get_fundamental_variables(self):

###  L_sei = L_inner + L_outer, defined in the sei base class

        L_inner = pybamm.FullBroadcast(

            0, self.domain.lower() + " electrode", "current collector"

        )

        L_outer = pybamm.standard_variables.L_outer



        j_inner = pybamm.FullBroadcast(

            0, self.domain.lower() + " electrode", "current collector"

        )

        j_outer = pybamm.Variable(

            "Outer " + self.domain + " electrode sei interfacial current density",

            domain=self.domain.lower() + " electrode",

            auxiliary_domains={"secondary": "current collector"},

        )



        j_plating = pybamm.Variable(

            "Outer " + self.domain + " electrode lithium plating interfacial current density",

            domain=self.domain.lower() + " electrode",

            auxiliary_domains={"secondary": "current collector"},

        )





        variables = self._get_standard_thickness_variables(L_inner, L_outer)

        variables.update(self._get_standard_reaction_variables(j_inner, j_outer, j_plating)) # j_plating needs to be updated as well



        return variables



    def get_coupled_variables(self, variables):


        # Get variables related to the concentration

        variables.update(self._get_standard_concentration_variables(variables))



        # Update whole cell variables, which also updates the "sum of" variables

        if (

            "Negative electrode sei interfacial current density" in variables

            and "Positive electrode sei interfacial current density" in variables

            and "Sei interfacial current density" not in variables

        ):

            variables.update(

                self._get_standard_whole_cell_interfacial_current_variables(variables) 

            ) #This function is from interface base class or from base_kinetics if jtot is taken as a state variable



        return variables



    def set_rhs(self, variables):

        domain = self.domain.lower() + " electrode"

        L_sei = variables["Outer " + domain + " sei thickness"]

        j_sei = variables["Outer " + domain + " sei interfacial current density"]
        
        j_plating = variables["Outer " + domain + " lithium plating interfacial current density"]


        if self.domain == "Negative":
            
            Gamma_SEI = self.param.Gamma_SEI_n
            
            Gamma_Plating = self.param.Gamma_plating



        self.rhs = {L_sei: -Gamma_SEI * j_sei / 2 - Gamma_Plating * j_plating /2 }



    def set_algebraic(self, variables):

        phi_s_n = variables[self.domain + " electrode potential"]

        phi_e_n = variables[self.domain + " electrolyte potential"]

        j_sei = variables[

            "Outer "

            + self.domain.lower()

            + " electrode sei interfacial current density"

        ]
        
        j_plating = variables[

            "Outer "

            + self.domain.lower()

            + " electrode lithium plating interfacial current density"

        ]

        L_sei = variables["Outer " + self.domain.lower() + " electrode sei thickness"]

        c_ec = variables[self.domain + " electrode EC surface concentration"]



        # Look for current that contributes to the -IR drop

        # If we can't find the interfacial current density from the main reaction, j,

        # it's ok to fall back on the total interfacial current density, j_tot

        # This should only happen when the interface submodel is "InverseButlerVolmer"

        # in which case j = j_tot (uniform) anyway

        if (

            "Total "

            + self.domain.lower()

            + " electrode interfacial current density variable"

            in variables

        ):

            j = variables[

                "Total "

                + self.domain.lower()

                + " electrode interfacial current density variable"

            ]

        else:

            j = variables[

                "X-averaged "

                + self.domain.lower()

                + " electrode total interfacial current density"

            ]



        if self.domain == "Negative":

            C_sei_ec = self.param.C_sei_ec_n   #  This is the scale factor for solvent concentration c_ec

            R_sei = self.param.R_sei_n
            
            C_plating = self.param.C_lithium_plating_reaction


        # need to revise for thermal case
        if isinstance(self, pybamm.lithium_plating.Tafel):
           # Model lithium plating using Cathodic Tafel equation
  
     
           self.algebraic = {

            j_sei: j_sei

            + C_sei_ec

            * c_ec

            * pybamm.exp(-0.5 * (phi_s_n - phi_e_n - j * L_sei * R_sei)),
            
            
            j_plating: j_plating

            + C_plating

            * pybamm.exp(-0.5 * (phi_s_n - phi_e_n - j * L_sei * R_sei)), # L_sei is required here becuase R_sei is ohms.m instead of ohms.m^2


        }


        elif isinstance(self, pybamm.lithium_plating.Butler_Volver):
             
             eta = phi_s_n - phi_e_n - j * L_sei * R_sei  # Expression for lithium plating overpotential
             
             alpha_a= self.param. alpha_a_BV_plating      # anodic transfer coefficient for plating
             
             alpha_c= self.param. alpha_c_BV_plating     # cathodic transfer coefficient for plating
             
             self.algebraic = {

            j_sei: j_sei

            + C_sei_ec

            * c_ec

            * pybamm.exp(-0.5 * (phi_s_n - phi_e_n - j * L_sei * R_sei)),
            
            
             j_plating : j_plating + (1 / C_plating) * (
                 
                 pybamm.exp( alpha_a *eta)
                 
                 - pybamm.exp(-alpha_c*eta)
                 
        )
            
      

        }


        def set_initial_conditions(self, variables):

           L_sei = variables["Outer " + self.domain.lower() + " electrode sei thickness"]

           j_sei = variables[

            "Outer "

            + self.domain.lower()

            + " electrode sei interfacial current density"]

        
        
           j_plating = variables[

            "Outer "

            + self.domain.lower()

            + " electrode lithium plating interfacial current density"]



           L_sei_0 = pybamm.Scalar(1)

           j_sei_0 = pybamm.Scalar(0)
        
           j_plating_0 = pybamm.Scalar(0)
 


           self.initial_conditions = {L_sei: L_sei_0, j_sei: j_sei_0, j_plating:j_plating_0}