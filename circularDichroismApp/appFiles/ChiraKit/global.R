packages <- c("shinydashboard","shinycssloaders","rhandsontable","plotly",
              "shinyalert","reticulate","DT","reshape2","tidyverse",
              "colourpicker",'signal','ggdendro','factoextra','FactoMineR')

invisible(lapply(packages, library, character.only = TRUE))

appName     <- "ChiraKit"
user        <- Sys.info()['user']

reticulate::use_python(paste0("/home/",user,"/myenv/bin/python"), required = TRUE)

# developer path
base_dir <- paste0("/home/",user,"/spc_shiny_servers/circularDichroismApp/appFiles/",appName,"/")

# path for the docker user
if (user == 'shiny') {
  base_dir <- paste0("/home/shiny/",appName,'/')
}

colorPalette9 <- c("#E41A1C","#377EB8","#4DAF4A",
                   "#984EA3","#FF7F00","#FFFF33",
                   "#A65628","#F781BF","#999999")

colorPalette12 <- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3",
                    "#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD",
                    "#CCEBC5","#FFED6F")

global_cd_units_choice_short <- c(
  'Millidegrees'       = 'millidegrees',
  'Degrees'           = 'degrees',
  'Differential absorbance'        = 'absorbance',
  'Differential milliabsorbance'    = 'milliabsorbance'
)

global_cd_units_choice <- c(
  'Millidegrees (m°)'       = 'millidegrees',
  'Degrees (°)'           = 'degrees',
  'Molar ellipticity' = 'molarEllipticity',
  'Mean unit molar ellipticity ([θ])'  = 'meanUnitMolarEllipticity',
  'Differential absorbance (\u0394A)'        = 'absorbance',
  'Differential milliabsorbance (\u0394mA)'    = 'milliabsorbance',
  'Molar extinction'  = 'molarExtinction',
  'Mean unit molar extinction (\u0394\u03B5)'  = 'meanUnitMolarExtinction'
)

myrenderer <- "function(instance, td, row, col, prop, value, cellProperties) {
                Handsontable.renderers.TextRenderer.apply(this, arguments);
                if (instance.params) {
                    hcols = instance.params.col_highlight
                    hcols = hcols instanceof Array ? hcols : [hcols]
                    hrows = instance.params.row_highlight
                    hrows = hrows instanceof Array ? hrows : [hrows]
                    
                    for (i = 0; i < hcols.length; i++) { 
                        if (hcols[i] == col && hrows[i] == row) {
                            td.style.background = instance.getDataAtCell(row, col);
                        }
                    }
                }
  }"     

myrendererBoolean <- "function(instance, td, row, col, prop, value, cellProperties) {
            Handsontable.renderers.CheckboxRenderer.apply(this, arguments);}"


