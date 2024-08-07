#' Function to choose colors for cell types
#'
#' Uses Giotto::getDistinctColors to begin with. Orders colors so the most
#' common cell types get the lightest colors. Removes colors that are too light
#' (sum of rgb values > 600)
#' @param names Vector of cell type names
#' @param freqs Optional, named vector of cell type abundance (e.g. c(T = 1000,
#'   tumor = 15000...))
#' @param init_colors Optional, a named vector of cell colors. This will be used
#'   for all cell types in the "names" vector that match names(init_colors).
#'   Intended for use with the iocolors vector (found in the Ptolemy package
#'   data).
#' @param max_sum_rgb Don't return any colors with total rgb values above this
#'   level. (Removes excessively light colors.)
#' @param palette One of "tableau20", "brewers" or "earthplus".
#' @return A named color vector
#' @importFrom grDevices col2rgb colors
#' @export
#' @examples
#' data("mini_nsclc")
#' unsup <- insitutype(
#'  x = mini_nsclc$counts,
#'  neg = Matrix::rowMeans(mini_nsclc$neg),
#'  n_clusts = 8,
#'  n_phase1 = 200,
#'  n_phase2 = 500,
#'  n_phase3 = 2000,
#'  n_starts = 1,
#'  max_iters = 5,
#'  assay_type="RNA"
#' ) # choosing inadvisably low numbers to speed the vignette; using the defaults in recommended.
#' colorCellTypes(freqs = table(unsup$clust), palette = "brewers")

colorCellTypes <- function(names = NULL, freqs = NULL, init_colors = NULL, max_sum_rgb = 600, 
                           palette = "earthplus") {
  
  if (is.null(freqs) && is.null(names)) {
    stop("must specify either names or freqs")
  } 
  
  if (is.null(freqs) && palette == "earthplus") {
    warning("this palette is best used when cell frequencies are known.")
  }
  
  if (is.null(freqs)) {
    # format names into freqs, then work with freqs henceforth
    freqs <- rep(1, length(names))
    names(freqs) <- names
  }
  
  ### "brewers" version: increasingly bright Rcolorbrewer paletted:
  if (palette == "brewers") {
    # start with R colorbrewer pallettes, then add a ton of filler colors:
    cols <- c('#8DD3C7','#FFFFB3','#BEBADA','#FB8072','#80B1D3','#FDB462','#B3DE69','#FCCDE5',
              '#D9D9D9','#BC80BD','#CCEBC5','#FFED6F','#66C2A5','#FC8D62','#8DA0CB','#E78AC3',
              '#A6D854','#FFD92F','#E5C494','#B3B3B3','#E41A1C','#377EB8','#4DAF4A','#984EA3',
              '#FF7F00','#FFFF33','#A65628','#F781BF','#999999','firebrick','darkorange2','tan3',
              'magenta','wheat4','palevioletred2','dodgerblue4','tomato3','mediumspringgreen',
              'grey26','antiquewhite4','red1','blue2','olivedrab4','lightyellow1','rosybrown3',
              'lightsteelblue4','rosybrown','rosybrown2','snow1','pink4','ghostwhite','ivory4',
              'lightgoldenrod','royalblue1','deeppink1','white','violetred2','hotpink2',
              'lightblue3','chartreuse4','azure2','plum','springgreen2','lemonchiffon1',
              'goldenrod2','grey6','darkorchid','palevioletred4','green4','lightsalmon1',
              'saddlebrown','rosybrown1','antiquewhite1','whitesmoke','plum4','cyan2',
              'forestgreen','burlywood3','lightyellow4','firebrick1','khaki3','salmon3',
              'sienna2','coral1','tan1','mediumvioletred','springgreen1','lemonchiffon',
              'lightgoldenrod4','darkred','navajowhite1','lightcoral','mediumturquoise',
              'lavenderblush','mistyrose1','indianred2','darkgoldenrod4','lightgoldenrod1',
              'lightsalmon3','lavender','magenta4','tomato2','seashell3','purple','tan2',
              'palevioletred3','coral3','lightblue1','darkorange4','orange1','darkolivegreen',
              'maroon1','skyblue3','cadetblue2','mediumorchid3','gold3','violetred1',
              'ivory2','snow4','aquamarine','darkgrey','darkolivegreen3','turquoise4',
              'sienna4','springgreen4','peachpuff4','seashell','violet','turquoise',
              'bisque2','lightsteelblue2','honeydew','lightsteelblue3','lawngreen',
              'tomato4','lightsalmon4','chocolate2','black','lightpink4','deepskyblue4',
              'aquamarine3','dodgerblue1','salmon1','yellow3','wheat','skyblue4','navajowhite4',
              'purple2','lavenderblush1','darkorange1','khaki2','aquamarine1','honeydew2',
              'cornsilk','lightskyblue4','mediumpurple2','paleturquoise1','seashell1',
              'darkcyan','orchid','royalblue','darkseagreen2','seagreen4','darkmagenta',
              'lightblue','mediumblue','chocolate3','yellow','darkgoldenrod2','mediumorchid4',
              'palegreen2','olivedrab','darkslateblue','chocolate1','maroon2','grey36',
              'orangered','goldenrod1','bisque3','deeppink3','peachpuff3','darkgreen',
              'royalblue4','darkgoldenrod1','blanchedalmond','mistyrose4','turquoise2',
              'ivory3','orchid1','limegreen','mediumpurple1','darkorange3','lemonchiffon4',
              'palevioletred1','magenta2','blue4','cyan1','thistle4','peru','grey56','cornsilk4',
              'mediumorchid2','green2','lightblue4','salmon4','burlywood4','burlywood1','orange',
              'burlywood','purple4','plum1','violetred3','khaki4','lightgoldenrodyellow',
              'lavenderblush3','lightpink3','azure4','orangered4','yellow2','mistyrose2',
              'deepskyblue2','mediumaquamarine','slateblue1','orange2','coral2','darkorchid4',
              'lightsalmon','gold2','darkseagreen')
    cols <- cols[!duplicated(cols)]
    
    # remove colors that are too light:
    sum_rgb <- colSums(grDevices::col2rgb(cols))
    cols <- cols[sum_rgb < max_sum_rgb]
    # add more colors if needed:
    n_removed <- sum(sum_rgb >= max_sum_rgb)
    if (n_removed > 0) {
      newcols <-  sample(colors()[!grepl("grey", colors())], length(freqs) * 2)[length(freqs) + seq_len(length(freqs))]
      newcols <- newcols[colSums(grDevices::col2rgb(newcols)) < max_sum_rgb]
      cols <- c(cols, newcols[seq_len(n_removed)])
    }
    
    # order so the most common cells have lighter colors:
    cols <- cols[seq_along(freqs)]
    names(cols) <- names(freqs)[order(freqs, decreasing = TRUE)]
  }
  
  ### "tableau20" palette: start with the tablueau20 colors:
  if (palette == "tableau20") {
    tab20 <- c('#aec7e8','#ffbb78','#98df8a','#ff9896','#c5b0d5','#c7c7c7',
               '#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b',
               '#e377c2','#17becf','#7f7f7f',
               '#8DD3C7','#FFFFB3','#BEBADA','#FB8072','#80B1D3','#FDB462',
               '#B3DE69','#FCCDE5','#D9D9D9','#BC80BD','#CCEBC5','#FFED6F',
               sample(colors()[!grepl("grey", colors())], 200, replace = FALSE))
    cols <- tab20[seq_along(freqs)]
    names(cols) <- names(freqs)[order(freqs)]
  }
  
  ### "earthplus" palette: earthtones for common cells, radiant colors for rare cells:
  if (palette == "earthplus") {
    # step 1: top least common cells, as long as <1% freq, get "radiant" colors:
    radiantcolors <-
      c(
        "#FF0000",
        "#00CCFF",
        "#00FF00",
        "#FFFF00",
        "#FF00CC",
        "#00FFFF",
        "#FF3300",
        "#CC00FF",
        "#CCFF00",
        "#66FF33"
      )
    richcolors <- c("#660099", "#006600", "#000000", "#000066")
    nlow <- min(sum(freqs < 0.01), 14)
    lowcols <- c(radiantcolors, richcolors)[seq_len(nlow)]
    
    # step 2: most common cells, as long as >10% freq, get "earth" colors:
    earthtones <- c('#D9AF6B','#AF6458','#526A83','#68855C','#9C9C5E','#855C75')
    nhigh <- min(sum(freqs > 0.1), length(earthtones)) 
    highcols <- earthtones[seq_len(nhigh)]
    
    # step 3: remainder get mid-range colors
    moderatecolors <- c('#1D6996','#73AF48','#E17C05','#94346E','#EDAD08','#38A6A5', 
                        '#CC503E','#0F8554','#5F4690',   
                        '#8DD3C7','#FFFFB3','#BEBADA','#FB8072','#80B1D3','#FDB462',
                        '#B3DE69','#FCCDE5','#D9D9D9','#BC80BD','#CCEBC5','#FFED6F',
                        sample(colors()[!grepl("grey", colors())], 200, replace = FALSE))
    nmid <- length(freqs) - nlow - nhigh
    if (nmid < length(moderatecolors)) {
      midcols <- moderatecolors[seq_len(nmid)]
    } else {
      stop("too many cell types")
    }
    cols <- c(lowcols, midcols, highcols)
    names(cols) <- names(freqs)[order(freqs)]
  }
  
  # if init_colors are provided, use them when possible:
  if (!is.null(init_colors)) {
    overlap <- intersect(names(cols), names(init_colors))
    cols[overlap] <- init_colors[overlap]
  }
  
  return(cols)
}
