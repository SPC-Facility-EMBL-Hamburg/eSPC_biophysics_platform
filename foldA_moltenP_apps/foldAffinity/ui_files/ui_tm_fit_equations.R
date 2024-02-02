box(title = "Equations", width = 12, solidHeader = T, status = "primary", 
    fluidRow(
      withMathJax(),
      column(6, p(HTML("One Site"),
                  withMathJax("$$ T_{mObs} = \\frac{T_{m}}{1 - \\frac{RT_{m} * ln(1+\\frac{L}{Kd})}{\\Delta H}}  $$"))),
      column(6, p(HTML(""),
                  h5("where L (free ligand concentration) is supposed to be similar as the total ligand concentration,
                     R is the gas constant in kcal/mol, DH is the enthalpy in kcal/mol, Tm is the melting temperature 
                     (without ligand),and TmObs is the observed melting temperature with ligand.")))),
    fluidRow(
      column(6, p(HTML("Two Sites"),
                  withMathJax("$$ T_{mObs} = \\frac{T_{m}}{1 - \\frac{RT_{m} * ln(1+L*\\frac{L+Kd1+Kd2}{Kd1*Kd2})}{\\Delta H}}  $$"))),
      column(6, p(HTML(""), 
                  h5("where L (free ligand concentration) is supposed to be similar as the total ligand concentration,
                     R is the gas constant in kcal/mol, DH is the enthalpy in kcal/mol, Tm is the melting temperature 
                     (without ligand), TmObs is the observed melting temperature with ligand, 
                     Kd1 is the dissociation constant L*F/LF, and Kd2 is the dissociation constant L*F/FL")))
      ))

#(Tm - Tm0) / Tm=RTm0ln (1+LKd) / H 
