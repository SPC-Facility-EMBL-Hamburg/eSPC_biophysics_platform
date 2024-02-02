box(title = "CD data (selected 'working units')", width = 4, solidHeader = T, status = "primary", 
    fluidRow(
      column(12, p(style = "font-size: 120%",HTML(""),
                  downloadLink('download_cd_data_col_wise'          ,'All CD data       (column-wise format)'))),
      
      column(12, p(style = "font-size: 120%",HTML(""),
                  downloadLink('download_cd_data_row_wise'          ,'All CD data       (row-wise format)'))),
      
      column(12, p(style = "font-size: 120%",HTML(""),
                  downloadLink('download_generated_cd_data_col_wise','Generated CD data (column-wise format)'))),

      column(12, p(style = "font-size: 120%",HTML(""),
                  downloadLink('download_generated_cd_data_row_wise','Generated CD data (row-wise format)'))),      
            
      column(12, p(style = "font-size: 120%",HTML(""),
                  selectInput('selected_cd_exp',NULL,choices = 'None'))),
      
      column(12, p(style = "font-size: 120%",HTML(""),
                  downloadLink('download_selected_cd_exp'  ,'Specific CD curve(s)')))
))