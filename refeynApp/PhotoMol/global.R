packages <- c("shinydashboard","shinycssloaders","plotly","shinyalert","reticulate","reshape2",
              "rhandsontable","colourpicker")

invisible(lapply(packages, library, character.only = TRUE))

appName     <- "PhotoMol"
user        <- Sys.info()['user']

reticulate::use_python(paste0("/home/",user,"/myenv/bin/python"), required = TRUE)

# developer path
base_dir <- paste0("/home/",user,"/spc_shiny_servers/refeynApp/",appName,"/")

# path for the docker user
if (user == 'shiny') {
  base_dir <- paste0("/home/shiny/",appName,'/')
}

cstFactorForContrast <- 1e3 # Contrasts will be multiplied by this value

colorPalette9 <- c("#E41A1C","#377EB8","#4DAF4A",
                  "#984EA3","#FF7F00","#FFFF33",
                  "#A65628","#F781BF","#999999")

colorPalette12 <- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3",
                    "#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD",
                    "#CCEBC5","#FFED6F")

colorPalette40 <- c(
  "#8DD3C7", "#ADDFC1", "#CDEBBB", "#EDF8B6", "#F6F6B8", "#E4E2C3", "#D2CFCE",
  "#BFBBD9", "#CDABBF", "#DE9AA2", "#F08A84", "#EE857B", "#CB9297", "#A9A0B2",
  "#86AECE", "#9CB1B8", "#C0B299", "#E3B379", "#F7B762", "#E2C364", "#CDCE66",
  "#B8DA68", "#C1DA82", "#D6D5A5", "#EBD0C8", "#FACDE4", "#F0D1E1", "#E6D4DD",
  "#DCD7DA", "#D3C9D3", "#CBAFCC", "#C396C4", "#BC82BD", "#C0A0BF", "#C5BFC1",
  "#C9DDC3", "#D3EBB7", "#E2EB9F", "#F0EC87", "#FFED6F"
)

histogram_palette <- c('#AEC6CF','#FFB347','#77DD77','#CFCFC4',
                       '#e02b35','#b42009','#009275','#b66e7d')


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

