packages <- c("shinydashboard","shinycssloaders","plotly","shinyalert","reticulate","reshape2",
              "rhandsontable","colourpicker")

invisible(lapply(packages, library, character.only = TRUE))

appName     <- "PhotoMol"

user      <- Sys.info()['user']
users_dir <- paste0("/home/",user,"/data_users/")

notebook_app  <- (Sys.info()["nodename"] == "osvaldo")

if (notebook_app) {
  use_python("/home/osvaldo/miniconda3/bin/python")
} else  {
  reticulate::use_condaenv("r-reticulate",required = TRUE)
}

# developer path
base_dir <- paste0("/home/",user,"/spc_shiny_servers/refeynApp/",appName,"/")

# set the corrrect path for the docker user
if (user == 'shiny') {
  base_dir <- paste0("/home/shiny/",appName,'/')
}

factorForContrast <- 1e3 # Contrasts will be multiplied by this value 

colorPalette9 <- c("#E41A1C","#377EB8","#4DAF4A",
                  "#984EA3","#FF7F00","#FFFF33",
                  "#A65628","#F781BF","#999999")

colorPalette12 <- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3",
                    "#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD",
                    "#CCEBC5","#FFED6F")

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

