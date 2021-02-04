# Class for reaction limited SEI growth

#

import pybamm

from .base_sei import BaseModel





class ReactionLimitedSEI_With_Lithium_Plating(BaseModel):

    """

    Class for reaction limited SEI growth with lithium plating reaction.



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

        L_inner = pybamm.standard_variables.L_inner

        L_outer = pybamm.standard_variables.L_outer



        variables = self._get_standard_thickness_variables(L_inner, L_outer)

        variables.update(self._get_standard_concentration_variables(variables))



        return variables



    def get_coupled_variables(self, variables):

        param = self.param

        phi_s_n = variables[self.domain + " electrode potential"]

        phi_e_n = variables[self.domain + " electrolyte potential"]



        # Look for current that contributes to the -IR drop

        # If we can't find the interfacial current density from the main reaction, j,

        # it's ok to fall back on the total interfacial current density, j_tot

        # This should only happen when the interface submodel is "InverseButlerVolmer"

        # in which case j = j_tot (uniform) anyway

        if self.domain + " electrode interfacial current density" in variables:

            j = variables[self.domain + " electrode interfacial current density"]

        else:

            j = variables[

                "X-averaged "

                + self.domain.lower()

                + " electrode total interfacial current density"

            ]

        L_sei = variables["Total " + self.domain.lower() + " electrode sei thickness"]



        if self.domain == "Negative":

            R_sei = self.param.R_sei_n

        alpha = 0.5

        # alpha = param.alpha

        if self.domain == "Negative":

            C_sei = param.C_sei_reaction_n
            
            C_plating = param.C_lithium_plating_reaction


        # need to revise for thermal case
        if isinstance(self, pybamm.lithium_plating.Tafel):     # Using Tafel cathodic expression for Li plating
           j_sei = -(1 / C_sei) * pybamm.exp(

            -0.5 * (phi_s_n - phi_e_n - j * L_sei * R_sei)

        )
                   
           #C_lithium_plating_reaction = 
           #(self.J_scale_n/self.m_plating_dimensional)*pybamm.exp(-self.F*(self.U_n_ref)/(2*self.R * self.Tref)) add this scale factor in the parameters
           
           j_plating = -(1 / C_plating ) * pybamm.exp(

            -0.5 * (phi_s_n - phi_e_n - j * L_sei * R_sei)

        )
           
        elif isinstance(self, pybamm.lithium_plating.ButlerVolmer):
            
             alpha_a = param.alpha_a_BV_plating     # anodic transfer coefficient for Li plating
             
             
             alpha_c = param.alpha_c_BV_plating      # cathodic transfer coefficient for Li plating
            
             j_sei = -(1 / C_sei) * pybamm.exp(

            -0.5 * (phi_s_n - phi_e_n - j * L_sei * R_sei)

        )    
             
             
             eta = phi_s_n - phi_e_n - j * L_sei * R_sei
             
             
             j_plating = -(1 / C_plating) * (
                 
                 pybamm.exp( alpha_a *eta)
                 
                 - pybamm.exp(-alpha_c*eta)
                 
                 # Need to enforce the condition j_plating = min(0,j_plating ) if pybamm.lithium_plating.ButlerVolmer is True
                 

        )
           
            

             j_inner = alpha * j_sei
             

             j_outer = (1 - alpha) * j_sei

         # I have assumed that plating occur only on the surface of the graphite and that
         # that plated lithium contributes only to the outer layer sei growth

        variables.update(self._get_standard_reaction_variables(j_inner, j_outer,j_plating))



        # Update whole cell variables, which also updates the "sum of" variables

        if (

            "Negative electrode sei interfacial current density" in variables

            and "Positive electrode sei interfacial current density" in variables

            and "Sei interfacial current density" not in variables

        ):

            variables.update(

                self._get_standard_whole_cell_interfacial_current_variables(variables)

            )



        return variables



    def set_rhs(self, variables):

        domain = self.domain.lower() + " electrode"

        L_inner = variables["Inner " + domain + " sei thickness"]

        L_outer = variables["Outer " + domain + " sei thickness"]

        j_inner = variables["Inner " + domain + " sei interfacial current density"]

        j_outer = variables["Outer " + domain + " sei interfacial current density"]
        
        j_plating = variables["Outer " + domain + " lithium plating interfacial current density"]


        v_bar = self.param.v_bar

        if self.domain == "Negative":

            Gamma_SEI = self.param.Gamma_SEI_n
            
            #Add Gamma_plating parameters in the parameter class
            #self.Gamma_Plating = (self.M_li * self.j_scale_n * self.tau_discharge)/((self.F * self.L_sei_0_dim*self.Den_li)
            # M_li = "Lithium metal molar mass [kg.mol-3]"
            # Den_li = "Lithium metal density kg. m-3]"
            
            Gamma_Plating = self.param.Gamma_plating


        self.rhs = {

            L_inner: -Gamma_SEI * j_inner,   # i have assumed that plating contributes to only outer layer sei growth

            L_outer: (-v_bar * Gamma_SEI * j_outer)-(Gamma_Plating*j_plating/2),

        }



    def set_initial_conditions(self, variables):

        domain = self.domain.lower() + " electrode"

        L_inner = variables["Inner " + domain + " sei thickness"]

        L_outer = variables["Outer " + domain + " sei thickness"]
        



        L_inner_0 = self.param.L_inner_0

        L_outer_0 = self.param.L_outer_0



        self.initial_conditions = {L_inner: L_inner_0, L_outer: L_outer_0}