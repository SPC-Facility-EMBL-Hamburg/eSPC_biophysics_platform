box(title = "Equations", width = 12, solidHeader = T, status = "primary", 
    fluidRow(
      withMathJax(),
      column(4, p(HTML("<b</b>"),
                  withMathJax("$$ Kd(T) = \\frac{1}{Kb(T)} $$"))),
      column(4, p(HTML("<b</b>"),
                  withMathJax("$$ Kb(T) =exp(\\frac{\\Delta G_B(T)}{RT}) $$"))),
      column(4, p(HTML("<b</b>"),
                  withMathJax("$$ Ku(T) =exp(\\frac{\\Delta G_U(T)}{RT}) $$")))),
    
    fluidRow(
      column(4, p(HTML("<b</b>"),
                  withMathJax("$$ \\Delta G_B (T) = \\Delta H_{B,T_B}+\\Delta Cp_B*(T-(T_B+273.15)) - 
                                              T(\\Delta S_{B,T_B} + \\Delta Cp_B ln(\\frac{T}{(T_B+273.15)})) $$")))),
    fluidRow(
      column(4, p(HTML("<b</b>"),
                  withMathJax("$$ \\Delta G_U (T) = \\Delta H_{U,T_m} * (1 - \\frac{T}{(T_m+273.15)}) -
                                              \\Delta Cp_U * ((T_m+273.15) - T + T * ln(\\frac{T}{(T_m+273.15)}) )$$")))),
    fluidRow(
      column(4, p(HTML("<b</b>"),
                  withMathJax("$$ Fluorescence(T,L0,P0) = R1(T)*F(Kd(T),Ku(T),T,L0,P0)+
                                              R2(T)*FL(Kd(T),Ku(T),T,L0,P0)+
                                              R3(T)*U(Kd(T),Ku(T),T,L0,P0)$$")))),
    fluidRow(
      column(4, p(HTML("<b</b>"),
                  withMathJax("$$ R1(T)=b1+m1*T $$"))),
      column(4, p(HTML("<b</b>"),
                  withMathJax("$$ R2(T)=b2+m2*T $$"))),
      column(4, p(HTML("<b</b>"),
                  withMathJax("$$ R3(T)=b3+m3*T $$")))))