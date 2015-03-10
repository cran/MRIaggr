#### 1- accept.coords #### 
Line 18 :   \item{accept.coords}{should coordinates be accepted as parameters ? \emph{logical}.} # initParameter.Rd

#### 2- accept.index #### 
Line 20 :   \item{accept.index}{should index be accepted as a parameter ? \emph{logical}.} # initParameter.Rd

#### 3- accept.mask #### 
Line 19 :   \item{accept.mask}{should mask be accepted as a parameter ? \emph{logical}.} # initParameter.Rd

#### 4- add #### 
Line 14 :   \item{add}{should \code{value} be added to the existing clinical slot ? \emph{logical}.} # affectClinic.Rd

#### 5- angle_test #### 
Line 28 :   \item{angle_test}{the angle or the number of angle to test. \emph{numeric vector} or \emph{postive integer}.} # calcHemisphere.Rd

#### 6- arg_name #### 
Line 21 :   \item{arg_name}{a short name for the error message ? \emph{character}.} # initParameter.Rd

#### 7- array #### 
Line 12 :   \item{array}{the array that should be converted into a data.frame. \emph{array} or \emph{matrix}. REQUIRED.} # array2df.Rd
Line 12 :   \item{array}{alternative specification of the spatial coordinates using an array where the non-NA values indicates the points of interest. \emph{array} or \code{NULL} leading to consider the \code{coords} argument.} # calcGroupsCoords.Rd
Line 12 :   \item{array}{the values of the contrast parameter. \code{array}. REQUIRED.} # constCarto3D.Rd
Line 12 :   \item{array}{alternative specification of the spatial coordinates using an array where the non-NA values indicates the points of interest. \emph{array} or \code{NULL} leading to use the \code{coords} argument.} # plotOutline.Rd

#### 8- as.logical #### 
Line 27 :   \item{as.logical}{should \code{mask} be convert to logical ? \emph{logical}.} # boxplotMask.Rd
Line 18 :   \item{as.logical}{should \code{mask} be convert to logical ? \emph{logical}.} # calcDistMask.Rd
Line 17 :   \item{as.logical}{should \code{mask} be convert to logical ? \emph{logical}.} # calcGroupsMask.Rd
Line 24 :   \item{as.logical}{should \code{mask} be convert to logical ? \emph{logical}.} # calcHemisphere.Rd
Line 19 :   \item{as.logical}{should \code{mask} be convert to logical ? \emph{logical}.} # calcROCthreshold.Rd
Line 18 :   \item{as.logical}{should \code{mask} be convert to logical ? \emph{logical}.} # calcSmoothMask.Rd
Line 21 :   \item{as.logical}{should \code{mask} be converted to logical ? \code{logical}.} # calcTableHypoReperf.Rd
Line 16 :   \item{as.logical}{should \code{mask} and \code{maskN} values be converted to logical ? \emph{logical}.} # calcTableLesion.Rd
Line 32 :   \item{as.logical}{should \code{seed} be converted to logical ? \emph{logical}.} # calcThresholdMRIaggr.Rd
Line 14 :   \item{as.logical}{should \code{mask} be converted to logical ? \emph{logical}.} # constReduceMRIaggr.Rd
Line 16 :   \item{as.logical}{if a parameter is specified for the index arguments, should it be converted to logical ? \emph{logical}.} # initIndex.Rd
Line 51 :   \item{as.logical}{if a parameter is specified for the index arguments, should it be converted to logical ? \emph{logical}.} # multiplot.Rd
Line 22 :   \item{as.logical}{if a parameter is specified for the index argument, should it be converted to logical ? \emph{logical}.} # outlineMRIaggr.Rd
Line 19 :   \item{as.logical}{should \code{mask} be convert to logical ? \emph{logical}.} # plotLesion3D.Rd

#### 9- as.numeric #### 
Line 17 :   \item{as.numeric}{should the default values be converted from character to numeric. \emph{logical}.} # selectDefault_value.Rd

#### 10- asp #### 
Line 73 :   \item{asp}{the aspect ratio y/x. \emph{numeric}.} # multiplot.Rd
Line 17 :   \item{asp}{the aspect ratio y/x. \emph{numeric}.} # plotMRI.Rd

#### 11- author #### 
Line 37 :   \item{author}{the author of the latex document. \emph{character}.} # constSweave.Rd

#### 12- axes #### 
Line 65 :   \item{axes}{should the axes be plotted ? \emph{logical}.} # multiplot.Rd
Line 23 :   \item{axes}{should the axes be plotted ? \emph{logical}.} # plotMRI.Rd

#### 13- beta #### 
Line 18 :   \item{beta}{the parameter 'inverse temperature' of the Potts model. \emph{numeric}.} # calcTissueType.Rd

#### 14- bg #### 
Line 20 :   \item{bg}{the color used for the background. \emph{character}.} # initDisplayWindow.Rd
Line 74 :   \item{bg}{the color used for the background. \emph{character}.} # multiplot.Rd

#### 15- breaks #### 
Line 16 :   \item{breaks}{the break points to use to categorize the contrast distribution. \code{numeric vector}.} # calcCriteriaGR.Rd
Line 20 :   \item{breaks}{the break points or the number of break points to use to categorize the contrast distribution. \code{numeric vector} or \code{postive integer}.} # calcGR.Rd
Line 22 :   \item{breaks}{the break points or the number of break points to use to categorize the contrast distribution. \code{numeric vector} or \code{postive integer}.} # calcSigmaGR.Rd
Line 18 :   \item{breaks}{the break points to use to categorize the contrast distribution. \code{numeric vector}.} # GRalgo.Rd
Line 27 :   \item{breaks}{the break points to use to generate the color intervals. \emph{numeric vector} or \code{NULL} leading to automatic breakpoints generation.} # heatmapMRIaggr.Rd
Line 18 :   \item{breaks}{the break points or the number of breakpoints to use to generate the color intervals. \emph{numeric vector} or \emph{postive integer}.} # initCol.Rd
Line 17 :   \item{breaks}{the break points or the number of break points to use to categorize the contrast distribution. \code{numeric vector} or \code{postive integer}.} # initGR.Rd
Line 11 :   \item{breaks}{the break points to use to generate the color intervals. \emph{numeric vector}.} # legendMRI.Rd
Line 53 :   \item{breaks}{the break points or the number of breakpoints to use to generate the color intervals . \emph{numeric vector} or \emph{postive integer}.} # multiplot.Rd
Line 29 :   \item{breaks}{the break points or the number of breakpoints to use to generate the color intervals . \emph{numeric vector} or \emph{postive integer}.} # outlineMRIaggr.Rd
Line 14 :   \item{breaks}{the break points to use to generate the color intervals. \emph{numeric vector}.} # plotMRI.Rd

#### 16- bw.adjust #### 
Line 25 :   \item{bw.adjust}{the smoothing bandwidth to use. \emph{numeric}. See \code{\link{density}} for more details.} # plotDistClass.Rd

#### 17- cex #### 
Line 16 :   \item{cex}{the expansion factor for the legend labels. \emph{positive numeric}.} # legendMRI.Rd
Line 58 :   \item{cex}{the expansion factor used to plot the observations. \emph{positive numeric}.} # multiplot.Rd
Line 17 :   \item{cex}{the expansion factor used to plot the edge points, the interpolated edge points and the interior points. \emph{positive numeric vector of size 3}.} # outline.Rd
Line 39 :   \item{cex}{the expansion factor used to plot the edge points, the interpolated edge points and the interior points. \emph{positive numeric vector of size 3}.} # outlineMRIaggr.Rd
Line 22 :   \item{cex}{the expansion factor for the observation labels. \emph{positive numeric}.} # plotMRI.Rd

#### 18- cex.axis #### 
Line 35 :   \item{cex.axis}{the magnification to be used for axis annotation relative to the current setting of cex. \emph{positive numeric}.} # heatmapMRIaggr.Rd
Line 31 :   \item{cex.axis}{the magnification to be used for axis annotation relative to the current setting of cex. \emph{positive numeric}.}		 # plotTableLesion.Rd

#### 19- cex.default #### 
Line 18 :   \item{cex.default}{the default expansion factor used to plot the observations. \emph{numeric}.} # initIndex.Rd

#### 20- cex.lab #### 
Line 32 :   \item{cex.lab}{the magnification to be used for x and y labels relative to the current setting of cex. \emph{positive numeric}.}		 # plotTableLesion.Rd

#### 21- cex.legend #### 
Line 36 :   \item{cex.legend}{the expansion factor of the legend. \emph{positive numeric}.}		 # boxplotMask.Rd
Line 84 :   \item{cex.legend}{the expansion factor of the legend. \emph{positive numeric}.} # multiplot.Rd
Line 39 :   \item{cex.legend}{the expansion factor of the legend. \emph{positive numeric}.}		 # plotDistClass.Rd
Line 30 :   \item{cex.legend}{the expansion factor of the legend. \emph{positive numeric}.}		 # plotTableLesion.Rd

#### 22- cex.main #### 
Line 15 :   \item{cex.main}{the expansion factor for the legend title. \emph{positive numeric}.} # legendMRI.Rd
Line 80 :   \item{cex.main}{the expansion factor for the main title. \emph{numeric}.} # multiplot.Rd
Line 29 :   \item{cex.main}{the expansion factor for the main title. \emph{numeric}.} # plotMRI.Rd
Line 29 :   \item{cex.main}{the expansion factor for the main title. \emph{numeric}.} # plotTableLesion.Rd

#### 23- class #### 
Line 14 :   \item{class}{the parameters indicating the probabilistic membership of each observations to each cerebral structure. \emph{character vector}.}  # calcDistTissues.Rd
Line 19 :   \item{class}{the parameters indicating the probabilistic membership to the tissue classes. \emph{character vector}.} # plotDistClass.Rd

#### 24- clinic #### 
Line 12 :   \item{clinic}{the clinical data of the patient. \emph{data.frame}.}   # MRIaggr-class.Rd
Line 13 :   \item{clinic}{should detailed information be printed for the clinical attribute ? \emph{logical}.} # summary.MRIaggr.Rd

#### 25- col #### 
Line 31 :   \item{col}{the colors of the boxplots for observations inside and outside the mask(s). \emph{character vector of size 2}.} # boxplotMask.Rd
Line 30 :   \item{col}{the colors with which the correlations will be displayed. \emph{character vector}.} # heatmapMRIaggr.Rd
Line 16 :   \item{col}{the color with which the observations will be displayed. \emph{character vector}.} # initCol.Rd
Line 56 :   \item{col}{the color with which the observations will be displayed. \emph{character vector} or \code{NULL} leading to determine the colors using the \code{palette} and \code{breaks} arguments.} # multiplot.Rd
Line 15 :   \item{col}{the colors in which the user-defined edge points, the interpolated edge points and the interior points should be ploted. \emph{character vector of size 3}.} # outline.Rd
Line 28 :   \item{col}{the color to use to plot the observations. \emph{character vector} or \code{NULL} leading to automatic generation of the colors using the \code{breaks} and \code{palette} arguments.} # outlineMRIaggr.Rd
Line 31 :   \item{col}{the colors with which the distributions will be displayed. \emph{character vector}.} # plotDistClass.Rd
Line 28 :   \item{col}{the color of the core of the lesion. \emph{character}.} # plotLesion3D.Rd
Line 16 :   \item{col}{the colors with which the observations will be displayed. \emph{character vector}.} # plotMRI.Rd
Line 23 :   \item{col}{the colors with which the volumes will be displayed. \emph{character vector} or \emph{numeric vector}.} # plotTableLesion.Rd
Line 13 :   \item{col}{the color with which the mid-saggital plan will be plotted. \emph{character}.} # pointsHemisphere.Rd

#### 26- col.default #### 
Line 20 :   \item{col.default}{the default color used to plot the observations. \emph{character vector}.} # initIndex.Rd

#### 27- col.edge #### 
Line 29 :   \item{col.edge}{the color of the edge of the lesion. \emph{character}.} # plotLesion3D.Rd

#### 28- col.midplane #### 
Line 62 :   \item{col.midplane}{the color in which the mid-saggital plan should appear. \emph{character}.} # multiplot.Rd

#### 29- col.NA #### 
Line 60 :   \item{col.NA}{the color to use to plot the NAs. \emph{character}.} # multiplot.Rd
Line 24 :   \item{col.NA}{the color to use to plot the NAs. \emph{character}.} # plotMRI.Rd

#### 30- col.outline #### 
Line 37 :   \item{col.outline}{the colors in which the user-defined edge points, the interpolated edge points and the interior points should be ploted. \emph{character vector}[3].} # outlineMRIaggr.Rd

#### 31- contrast #### 
Line 12 :   \item{contrast}{the contrast value of each voxel on a given slice. \emph{matrix}.} # calcContro_cpp.Rd
Line 12 :   \item{contrast}{the contrast value of each observations. \emph{numeric vector}. REQUIRED.} # calcCriteriaGR.Rd
Line 13 :   \item{contrast}{the contrast value of each observation. \emph{numeric vector}. REQUIRED.} # calcGR.Rd
Line 15 :   \item{contrast}{the contrast value of each observation. \emph{numeric vector}. REQUIRED.} # calcSigmaGR.Rd
Line 22 :   \item{contrast}{the dataset containing the contrast parameter to be thresholded. \code{matrix}. REQUIRED.} # calcThresholdMRIaggr.Rd
Line 12 :   \item{contrast}{the contrast and the spatial coordinates of each voxel. \emph{data.frame} with four columns named \code{"i"} \code{"j"} \code{"k"} and the name of the contrast parameter.}   # Carto3D-class.Rd
Line 11 :   \item{contrast}{the dataset containing the observations in rows and the contrast parameters in columns. \emph{vector} or \emph{data.frame}. REQUIRED.} # df2array.Rd
Line 12 :   \item{contrast}{the contrast value of each observation. \emph{numeric vector}.} # GRalgo.Rd
Line 12 :   \item{contrast}{the contrast value of each observation. \emph{matrix}.} # initCol.Rd
Line 11 :   \item{contrast}{the contrast value of each observation. \emph{numeric vector}.} # initGR.Rd
Line 11 :   \item{contrast}{the value of the contrast parameters for each observation. \emph{data.frame}.}   # MRIaggr-class.Rd
Line 42 :   \item{contrast}{the intensities to display. \emph{numerical vector} or \code{NULL} leading to use the same color for all observations.} # multiplot.Rd
Line 12 :   \item{contrast}{the intensities to display. \emph{numeric vector}.} # plotMRI.Rd

#### 32- coords #### 
Line 13 :   \item{coords}{the spatial coordinates of the observations contained in \code{array}. \emph{matrix} or \code{NULL}.} # array2df.Rd
Line 11 :   \item{coords}{the spatial coordinates of the observations. \emph{data.frame}. REQUIRED.} # calcGroupsCoords.Rd
Line 11 :   \item{coords}{the spatial coordinates of the observations. \emph{matrix} with a number of rows equal to the length of \code{sample}.} # calcRadius_cpp.Rd
Line 12 :   \item{coords}{the spatial coordinates of the observations. \emph{matrix} with a number of rows equal to the number of rows of \code{data}. REQUIRED.} # df2array.Rd
Line 13 :   \item{coords}{the spatial coordinates of the observations. \emph{data.frame}.} # initCol.Rd
Line 13 :   \item{coords}{the spatial coordinates of the observations. \emph{data.frame}.} # plotMRI.Rd
Line 11 :   \item{coords}{the spatial coordinates of the observations. \emph{data.frame}.} # plotOutline.Rd
Line 23 :   \item{coords}{the coordinates that sould be extracted. \emph{logical} or any of \code{"i"} \code{"j"} \code{"k"}.} # selectContrast.Rd
Line 17 :   \item{coords}{the coordinates that sould be extracted. Any of \code{"i"} \code{"j"} \code{"k"} or \code{"index"}.} # selectCoords.Rd

#### 33- coords_max #### 
Line 14 :   \item{coords_max}{the maximum coordinate in each dimension. \emph{numerical vector}} # calcGroupsCoords_cpp.Rd

#### 34- coords_NNA #### 
Line 11 :   \item{coords_NNA}{the spatial coordinates of the observations in C version (beginning at 0). \emph{matrix}.} # calcGroupsCoords_cpp.Rd

#### 35- coords_px #### 
Line 13 :   \item{coords_px}{the coordinates of the observations. \emph{matrix}.} # calcContro_cpp.Rd

#### 36- criterion_d1 #### 
Line 18 :   \item{criterion_d1}{should the d1 criterion be computed ? \emph{logical}.  Require \code{Wweight} to be computed.}   # calcCriteriaGR.Rd
Line 27 :   \item{criterion_d1}{should the d1 criterion be computed ? \emph{logical}.  Require \code{Wweight} to be computed.}   # calcSigmaGR.Rd

#### 37- criterion_entropy #### 
Line 19 :   \item{criterion_entropy}{should the entropy criterion be computed ? \code{logical}.}   # calcCriteriaGR.Rd
Line 28 :   \item{criterion_entropy}{should the entropy criterion be computed ? \code{logical}.}   # calcSigmaGR.Rd

#### 38- criterion_Kalinsky #### 
Line 20 :   \item{criterion_Kalinsky}{should the Kalinsky criterion be computed ? \code{logical}.}   # calcCriteriaGR.Rd
Line 29 :   \item{criterion_Kalinsky}{should the Kalinsky criterion be computed ? \code{logical}.}   # calcSigmaGR.Rd

#### 39- criterion_Laboure #### 
Line 21 :   \item{criterion_Laboure}{should the Laboure criterion be computed ? \code{logical}.} # calcCriteriaGR.Rd
Line 30 :   \item{criterion_Laboure}{should the Laboure criterion be computed ? \code{logical}.} # calcSigmaGR.Rd

#### 40- d_lim #### 
Line 16 :   \item{d_lim}{the distance within which the controlateral values are considered. \emph{numeric}.} # calcContro_cpp.Rd

#### 41- date #### 
Line 36 :   \item{date}{the date on the latex document. \emph{character}.} # constSweave.Rd

#### 42- decreasing #### 
Line 27 :   \item{decreasing}{should the increasing thresholding (\code{FALSE}) or decreasing thresholding (\code{TRUE}) be used. \code{logical}.} # calcThresholdMRIaggr.Rd

#### 43- default_value #### 
Line 15 :   \item{default_value}{the reference values of the contrast parameters (e.g background values). A one row \emph{data.frame} where the column length must match the length of the \code{param} argument.} # affectContrast.Rd
Line 14 :   \item{default_value}{the reference values of the contrast parameters (e.g. the background value). \emph{character}.}    # Carto3D-class.Rd
Line 15 :   \item{default_value}{the reference values of the contrast parameters (e.g. the background value). \emph{character} or \code{NULL} leading to search the reference value in \code{array[pos_default_value]}.} # constCarto3D.Rd
Line 15 :   \item{default_value}{the reference values of the contrast parameters (e.g. the background value). \emph{character} or \code{NULL} leading to search the reference value in \code{array[pos_default_value]}.} # constMRIaggr.Rd
Line 14 :   \item{default_value}{the element used to fill the missing observations. \emph{numeric}.} # df2array.Rd
Line 14 :   \item{default_value}{the reference values of the contrast parameters. \emph{data.frame}.}   # MRIaggr-class.Rd

#### 44- descStats #### 
Line 14 :   \item{descStats}{should detailed information be printed for the \code{ls_descStats} attribute ? \emph{logical}.} # summary.MRIaggr.Rd

#### 45- diagonal #### 
Line 24 :   \item{diagonal}{should the diagonal be added to the neighborhood matrix ? \emph{logical}.} # calcRegionalIntensity.Rd

#### 46- digit #### 
Line 20 :   \item{digit}{the number of decimal places to use when generating the thresholds. \emph{positive integer}.} # calcROCthreshold.Rd
Line 20 :   \item{digit}{the number of decimal places to use for the initialization. \emph{positive integer}.} # calcTissueType.Rd
Line 26 :   \item{digit}{the number of decimal places to use for the labels. \emph{positive integer}.} # heatmapMRIaggr.Rd
Line 19 :   \item{digit}{the number of decimal places to use for the legend label. \emph{positive integer}.} # legendMRI.Rd

#### 47- digit.legend #### 
Line 83 :   \item{digit.legend}{the number of decimal places to use for the legend labels. \emph{integer}.} # multiplot.Rd

#### 48- digit.plot #### 
Line 22 :   \item{digit.plot}{the number of decimal places to use for the legend. \emph{positive integer}.} # calcROCthreshold.Rd

#### 49- dim #### 
Line 16 :   \item{dim}{the number of bytes per element in the byte stream. The default, NA_integer_, uses the natural size. Size changing is not supported for raw and complex vectors. Only active if \code{format} equals \code{rawb.gz}.} # readMRI.Rd

#### 50- dir #### 
Line 15 :   \item{dir}{the path to the root(s) directory(ies) containing the image files. \emph{character} or \emph{character vector}. REQUIRED.} # constSweave.Rd

#### 51- distband #### 
Line 18 :   \item{distband}{the distance within which the controlateral values are considered. \emph{postive numeric}.} # calcControlateral.Rd
Line 15 :   \item{distband}{the neighborhood range. \emph{postive numeric}. REQUIRED.} # calcGroupsMask.Rd
Line 21 :   \item{distband}{the neighborhood range. \emph{postive numeric}. Required only if \code{W} have to be computed.} # calcRegionalIntensity.Rd
Line 19 :   \item{distband}{only distances smaller than delta are recorded. \code{postive numeric}. REQUIRED.} # calcW.Rd
Line 12 :   \item{distband}{the distband of the kernel. \emph{numeric}.} # EDK.Rd

#### 52- distband_EDK #### 
Line 22 :   \item{distband_EDK}{the distband of the kernel. \emph{postive numeric}. REQUIRED.} # calcRegionalIntensity.Rd

#### 53- edge #### 
Line 17 :   \item{edge}{should the edges of the lesion be ploted instead of the core ? \emph{logical}.} # plotLesion3D.Rd

#### 54- erosion_th #### 
Line 27 :   \item{erosion_th}{the threshold below which the observations will be removed by the erosion. \emph{numeric between 0 and 1}.} # calcSmoothMask.Rd

#### 55- extra_text #### 
Line 19 :   \item{extra_text}{additionnal text to display. \emph{list of character vector} or \code{NULL} if there is no extra text to display.} # constSweave.Rd

#### 56- factor #### 
Line 13 :   \item{factor}{the compression factor. \emph{postive integer}. REQUIRED.} # constCompressMRIaggr.Rd

#### 57- file #### 
Line 12 :   \item{file}{the file name of the imaging file. \emph{character}. REQUIRED.} # readMRI.Rd

#### 58- filename #### 
Line 38 :   \item{filename}{the name of the file used to export the plot. \emph{character}.} # boxplotMask.Rd
Line 33 :   \item{filename}{the name of the file used to export the plot. \emph{character}.} # calcBrainMask.Rd
Line 34 :   \item{filename}{the name of the file used to export the plot. \emph{character}.} # calcHemisphere.Rd
Line 25 :   \item{filename}{the name of the file used to export the plot. \emph{character}.} # calcROCthreshold.Rd
Line 38 :   \item{filename}{the name of the file used to export the plot. \emph{character}.} # calcSigmaGR.Rd
Line 37 :   \item{filename}{the name of the file used to export the plot. \emph{character}.} # heatmapMRIaggr.Rd
Line 12 :   \item{filename}{the name of the file used to export the plot. \emph{character}.} # initDisplayWindow.Rd
Line 13 :   \item{filename}{the name of the file used to export the plot. \emph{character}.} # initWindow.Rd
Line 88 :   \item{filename}{the name of the file used to export the plot. \emph{character}.} # multiplot.Rd
Line 41 :   \item{filename}{the name of the file used to export the plot. \emph{character}.} # plotDistClass.Rd
Line 34 :   \item{filename}{the name of the file used to export the plot. \emph{character}.} # plotTableLesion.Rd

#### 59- fill #### 
Line 31 :   \item{fill}{should the spatial region be filled ? Otherwise only the edge is used. \emph{logical}.} # outlineMRIaggr.Rd

#### 60- filter #### 
Line 18 :   \item{filter}{the filter to use. Can be a \emph{matrix} or an \emph{array}, or a name indicating which filter should be used. REQUIRED.} # calcFilter.Rd
Line 11 :   \item{filter}{the filter to be initialized. \emph{character}.} # initFilter.Rd
Line 13 :   \item{filter}{the type of filter, see \code{\link{calcFilter}} for more details. \emph{character}.} # plotOutline.Rd

#### 61- filter_default #### 
Line 21 :   \item{filter_default}{the default filter used to define the neighborhood. \emph{character}.}  # initIndex.Rd

#### 62- flipud #### 
Line 19 :   \item{flipud}{is a logical variable for vertical flipping of the image (default is TRUE). Only active if \code{format} equals \code{dicom}.} # readMRI.Rd

#### 63- format #### 
Line 26 :   \item{format}{the format of the output. Could be \code{"spam"} or \code{"dgCMatrix"}.} # calcW.Rd
Line 13 :   \item{format}{the format of the output. Can be \code{"any"},\code{"matrix"},\code{"data.frame"} or \code{"list"}.} # df2array.Rd
Line 13 :   \item{format}{the format of the image file. Can be \code{"raw.gz"}, \code{"analyze"},  \code{"nifti"} or \code{"dicom"}. REQUIRED.}  # readMRI.Rd
Line 21 :   \item{format}{the format of the output. Can be \code{"matrix"}, \code{"data.frame"} or \code{"any"}.} # selectContrast.Rd
Line 23 :   \item{format}{the format of the output. Can be \code{"matrix"}, \code{"data.frame"} or \code{"any"}.} # selectCoords.Rd

#### 64- from,to #### 
Line 27 :   \item{from,to}{the left and right-most points of the grid at which the density is to be estimated. \emph{numeric} or \code{NULL} leading to automatic adjustment. See \code{\link{density}} for more details.} # plotDistClass.Rd

#### 65- GRalgo #### 
Line 29 :   \item{GRalgo}{should a Growing Region algorithm be used to clean the thresholded parameter ? \code{logical}.} # calcThresholdMRIaggr.Rd

#### 66- groups #### 
Line 13 :   \item{groups}{the indicator of group membership. \code{logical vector}. REQUIRED.} # calcCriteriaGR.Rd

#### 67- height #### 
Line 40 :   \item{height}{the height of the device used to export the plot. \emph{postive numeric}.} # boxplotMask.Rd
Line 35 :   \item{height}{the height of the device used to export the plot. \emph{postive numeric}.} # calcBrainMask.Rd
Line 36 :   \item{height}{the height of the device used to export the plot. \emph{postive numeric}.} # calcHemisphere.Rd
Line 27 :   \item{height}{the height of the device used to export the plot. \emph{postive numeric}.} # calcROCthreshold.Rd
Line 40 :   \item{height}{the height of the device used to export the plot. \emph{postive numeric}.} # calcSigmaGR.Rd
Line 39 :   \item{height}{the height of the device used to export the plot. \emph{postive numeric}.} # heatmapMRIaggr.Rd
Line 15 :   \item{height}{the height of the device used to export the plot. \emph{postive numeric}.} # initDisplayWindow.Rd
Line 16 :   \item{height}{the height of the device used to export the plot. \emph{postive numeric}.} # initWindow.Rd
Line 90 :   \item{height}{the height of the device used to export the plot. \emph{postive numeric}.} # multiplot.Rd
Line 43 :   \item{height}{the height of the device used to export the plot. \emph{postive numeric}.} # plotDistClass.Rd
Line 36 :   \item{height}{the height of the device used to export the plot. \emph{postive numeric}.} # plotTableLesion.Rd

#### 68- hemisphere #### 
Line 22 :   \item{hemisphere}{the hemisphere to consider. \emph{character}.} # boxplotMask.Rd
Line 16 :   \item{hemisphere}{the hemisphere to use. \emph{character}.} # calcDistTissues.Rd
Line 18 :   \item{hemisphere}{the hemisphere to use. \emph{character}.} # calcRegionalIntensity.Rd
Line 24 :   \item{hemisphere}{the hemisphere to consider. \emph{character} or \code{NULL}.} # calcThresholdMRIaggr.Rd
Line 23 :   \item{hemisphere}{the hemisphere to use. \emph{character}.} # calcW.Rd
Line 20 :   \item{hemisphere}{the hemisphere to use. \emph{character}. See the details section of \code{\link{selectContrast}}.} # heatmapMRIaggr.Rd
Line 15 :   \item{hemisphere}{the hemisphere to display. \emph{character}.} # initIndex.Rd
Line 48 :   \item{hemisphere}{the hemisphere to display. \emph{character}.} # multiplot.Rd
Line 21 :   \item{hemisphere}{the hemisphere to consider. \emph{character}.} # outlineMRIaggr.Rd
Line 21 :   \item{hemisphere}{the hemisphere to use. \emph{character}.} # plotDistClass.Rd
Line 24 :   \item{hemisphere}{the hemisphere to extract. \emph{character}. See the details section.} # selectContrast.Rd
Line 20 :   \item{hemisphere}{the hemisphere to extract. \emph{character}.} # selectCoords.Rd
Line 16 :   \item{hemisphere}{the hemisphere to use. \emph{character}.} # selectDescStats.Rd
Line 13 :   \item{hemisphere}{the hemisphere of interest. \emph{character}. See the details section.} # selectHemispheres.Rd
Line 14 :   \item{hemisphere}{the hemisphere to consider. \emph{character}.} # selectN.Rd
Line 18 :   \item{hemisphere}{the hemisphere to extract. \emph{character}.} # selectNormalization.Rd

#### 69- hemispheres #### 
Line 18 :   \item{hemispheres}{the presence or absence of lesion in each cerebral hemisphere. \emph{data.frame}.}   # MRIaggr-class.Rd

#### 70- history #### 
Line 16 :   \item{history}{the list of the \code{calc} or \code{const} methods that have been already applied on the \code{MRIaggr} object. \emph{data.frame}.}   # MRIaggr-class.Rd
Line 15 :   \item{history}{should the\code{calc} and \code{const} methods that have been applied to the object be listed ? \emph{logical}.} # summary.MRIaggr.Rd

#### 71- history_front #### 
Line 28 :   \item{history_front}{should the propagation front of the GR set be recorded ? \emph{logical}.} # calcGR.Rd
Line 25 :   \item{history_front}{should the propagation front of the GR set be recorded ? \emph{logical}.} # GRalgo.Rd

#### 72- history_sigma #### 
Line 26 :   \item{history_sigma}{should the values of sigma be recorded ? \emph{logical}.} # calcGR.Rd
Line 23 :   \item{history_sigma}{should the values of sigma be recorded ? \emph{logical}.} # GRalgo.Rd

#### 73- history_step #### 
Line 27 :   \item{history_step}{should the number of observations included in the growing region set be recorded ? \emph{logical}.} # calcGR.Rd
Line 24 :   \item{history_step}{should the number of observations included in the GR set be recorded ? \emph{logical}.} # GRalgo.Rd

#### 74- i_test #### 
Line 26 :   \item{i_test}{the abscissa or the number of abscissa to test. \emph{numeric vector} or \emph{positive integer}.} # calcHemisphere.Rd

#### 75- identifier #### 
Line 10 :   \item{identifier}{the patient identifier. \emph{character}.}   # Carto3D-class.Rd
Line 13 :   \item{identifier}{the identifier of the patient from which the data originated. \emph{character}. REQUIRED.} # constCarto3D.Rd
Line 13 :   \item{identifier}{the identifier of the patient to which belong the contrast parameters. \emph{character}. REQUIRED.} # constMRIaggr.Rd
Line 16 :   \item{identifier}{the identifiers of the patients for which the graphics should be displayed. \emph{character vector} or \code{NULL} leading to use all patients.} # constSweave.Rd
Line 10 :   \item{identifier}{the patient identifier. \emph{character}.}   # MRIaggr-class.Rd

#### 76- index #### 
Line 13 :   \item{index}{the coordinates of additionnal points to display. \emph{data.frame} or \emph{list}.} # initIndex.Rd

#### 77- index_data #### 
Line 13 :   \item{index_data}{index of the non NA data.} # filtrage2D_cpp.Rd
Line 13 :   \item{index_data}{index of the non NA data.} # filtrage2Dmed_cpp.Rd
Line 16 :   \item{index_data}{index of the non NA data.} # filtrage3D_cpp.Rd
Line 15 :   \item{index_data}{index of the non NA data.} # filtrage3Dmed_cpp.Rd

#### 78- index_k #### 
Line 14 :   \item{index_k}{the index of the observations on the hemisphere of interest. \emph{integer vector}.} # calcContro_cpp.Rd

#### 79- index_k_contro #### 
Line 15 :   \item{index_k_contro}{the index of the observations on the controlateral hemisphere. \emph{integer vector}.} # calcContro_cpp.Rd

#### 80- index_NNA #### 
Line 12 :   \item{index_NNA}{the index of the coordinates in a array in C version (beginning at 0). \emph{numerical vector}.} # calcGroupsCoords_cpp.Rd

#### 81- index_subsection #### 
Line 22 :   \item{index_subsection}{the position of the images in the subsections. \emph{numeric vector} or \code{NULL} leading to use a subsection for each image.} # constSweave.Rd

#### 82- index_subsubsection #### 
Line 25 :   \item{index_subsubsection}{A list of positions of the images in the subsubsections. \emph{list of numeric vector} or \code{NULL}.} # constSweave.Rd

#### 83- index1 #### 
Line 19 :   \item{index1}{the coordinates of additionnal points to display. \emph{data.frame} or \emph{list} or \code{NULL}.} # outlineMRIaggr.Rd

#### 84- index1,index2,index3 #### 
Line 44 :   \item{index1,index2,index3}{the coordinates of additionnal points to display. \emph{data.frame} or \emph{list} or \code{NULL}.} # multiplot.Rd

#### 85- indexNum #### 
Line 17 :   \item{indexNum}{the number associated to the index (for display). \emph{numeric}.} # initIndex.Rd

#### 86- init #### 
Line 18 :   \item{init}{should the slice numbers be initialized if \code{num} equals \code{NULL} ? \emph{logical}.  } # initNum.Rd
Line 17 :   \item{init}{should the parameters be initialized if \code{param} equals \code{NULL} ? \emph{logical}.} # initParameter.Rd

#### 87- iter_max #### 
Line 22 :   \item{iter_max}{the maximum number of iterations for the expansion of the growing region. \code{postive integer}.}   # calcGR.Rd
Line 24 :   \item{iter_max}{the maximum number of iterations for the expansion of the growing region. \code{postive integer}.}   # calcSigmaGR.Rd
Line 20 :   \item{iter_max}{the maximum number of iterations for the expansion of the growing region. \code{postive integer}.}   # GRalgo.Rd

#### 88- j_test #### 
Line 27 :   \item{j_test}{the ordinates or the number of ordinates to test. \emph{numeric vector} or \emph{positive integer}.} # calcHemisphere.Rd

#### 89- keep.index #### 
Line 15 :   \item{keep.index}{should the previous \code{index} parameter be saved in the \code{ls_descStats} slot ? \emph{logical}.} # constReduceMRIaggr.Rd

#### 90- kernel #### 
Line 26 :   \item{kernel}{the smoothing kernel to use. \emph{character}. See \code{\link{density}} for more details.} # plotDistClass.Rd

#### 91- kmeans.n_groups #### 
Line 27 :   \item{kmeans.n_groups}{the number of groups to use in the kmeans algorithm. \emph{postive integer vector}.} # calcBrainMask.Rd

#### 92- kmeans.Neighborhood #### 
Line 28 :   \item{kmeans.Neighborhood}{the range of the neighborhood. \emph{postive integer}.} # calcBrainMask.Rd

#### 93- lambda #### 
Line 17 :   \item{lambda}{the importance of the penalization by the distance. \emph{numeric}.} # calcContro_cpp.Rd
Line 19 :   \item{lambda}{the importance of the penalization. \emph{numeric}.} # calcControlateral.Rd

#### 94- las #### 
Line 34 :   \item{las}{the style of the axis labels. Any of \code{0}, \code{1}, \code{2} or \code{3}.} # heatmapMRIaggr.Rd

#### 95- legend #### 
Line 27 :   \item{legend}{the legend of each image. \emph{character vector} or \code{NULL} leading to use \code{param} for the legend.} # constSweave.Rd
Line 19 :   \item{legend}{how should the legend be displayed ? \emph{logical} or \code{NULL}.} # initCol.Rd
Line 68 :   \item{legend}{how should the legend be displayed ? \emph{logical} or \code{NULL}.} # multiplot.Rd
Line 26 :   \item{legend}{how should the legend be displayed ? \emph{logical} or \code{NULL}.} # outlineMRIaggr.Rd

#### 96- long_name #### 
Line 22 :   \item{long_name}{the complete name for the error message ? \emph{character}..} # initParameter.Rd

#### 97- ls.array #### 
Line 12 :   \item{ls.array}{the value of the contrast parameter(s) for each observation. \code{list of array}. REQUIRED.} # constMRIaggr.Rd

#### 98- ls.Carto3D #### 
Line 12 :   \item{ls.Carto3D}{a list of \code{\linkS4class{Carto3D}} objects. REQUIRED.} # Carto3D2MRIaggr.Rd

#### 99- ls_statDesc #### 
Line 25 :   \item{ls_statDesc}{a slot to store additional data. \emph{list}.}   # MRIaggr-class.Rd

#### 100- lty #### 
Line 24 :   \item{lty}{the line type used to represent the volume. \emph{numeric vector}.} # plotTableLesion.Rd
Line 15 :   \item{lty}{the type of line used to represent the mid-saggital plan. \emph{numeric vector}.} # pointsHemisphere.Rd

#### 101- lwd #### 
Line 36 :   \item{lwd}{the line width. \emph{postive numeric}.} # plotDistClass.Rd
Line 25 :   \item{lwd}{the line width. \emph{postive numeric}.} # plotTableLesion.Rd
Line 14 :   \item{lwd}{the line width. \emph{postive numeric}.} # pointsHemisphere.Rd

#### 102- M_data #### 
Line 11 :   \item{M_data}{matrix to which the filter will be applied.} # filtrage2D_cpp.Rd
Line 11 :   \item{M_data}{matrix to which the filter will be applied.} # filtrage2Dmed_cpp.Rd

#### 103- M_operateur #### 
Line 12 :   \item{M_operateur}{the filter to be applied.} # filtrage2D_cpp.Rd
Line 12 :   \item{M_operateur}{the filter to be applied.} # filtrage2Dmed_cpp.Rd

#### 104- main #### 
Line 32 :   \item{main}{an overall title for the plot. \emph{character}.} # boxplotMask.Rd
Line 35 :   \item{main}{an overall title for the plot. \emph{character}.} # calcSigmaGR.Rd
Line 31 :   \item{main}{an overall title for the plot. \emph{character}.} # heatmapMRIaggr.Rd
Line 17 :   \item{main}{an overall title for the legend. \emph{character}.} # legendMRI.Rd
Line 78 :   \item{main}{an overall title for the plot. \emph{character}.} # multiplot.Rd
Line 32 :   \item{main}{an overall title for the plot. \emph{character}.} # plotDistClass.Rd
Line 28 :   \item{main}{an overall title for the plot. \emph{character}.} # plotMRI.Rd
Line 26 :   \item{main}{an overall title for the plot. \emph{character}.} # plotTableLesion.Rd

#### 105- main.legend #### 
Line 86 :   \item{main.legend}{a main title for the legend. \emph{character}.} # multiplot.Rd

#### 106- mar #### 
Line 33 :   \item{mar}{the number of margin lines to be specified on the four sides of the plot. \emph{positive numeric vector of size 4}.} # calcSigmaGR.Rd
Line 33 :   \item{mar}{the number of margin lines to be specified on the four sides of the plot. \emph{positive numeric vector of size 4}.} # heatmapMRIaggr.Rd
Line 22 :   \item{mar}{the number of margin lines to be specified on the four sides of the plot. \emph{positive numeric vector of size 4}.} # initDisplayWindow.Rd
Line 14 :   \item{mar}{the number of margin lines to be specified on the four sides of the legend. \emph{positive numeric vector of size 4}.} # legendMRI.Rd
Line 70 :   \item{mar}{the number of margin lines to be specified on the four sides of the plot. \emph{positive numeric vector of size 4}.} # multiplot.Rd
Line 28 :   \item{mar}{the number of margin lines to be specified on the four sides of the plot. \emph{positive numeric vector of size 4}.} # plotTableLesion.Rd

#### 107- mar.legend #### 
Line 85 :   \item{mar.legend}{the number of margin lines to be specified on the four sides of the legend. \emph{numeric vector  of size 4}.} # multiplot.Rd

#### 108- mask #### 
Line 20 :   \item{mask}{the binary contrast parameter(s) defining the spatial group(s). \emph{character vector}. REQUIRED.} # boxplotMask.Rd
Line 15 :   \item{mask}{the binary contrast parameter(s) defining the spatial groups from which the distance will be computed. \emph{character vector}. REQUIRED.} # calcDistMask.Rd
Line 14 :   \item{mask}{the binary contrast parameter that should be used to identifying the spatial groups. \emph{character}. REQUIRED.} # calcGroupsMask.Rd
Line 23 :   \item{mask}{the binary contrast parameter(s) indicating the lesion. \emph{character vector} or \code{NULL} if no mask is available.} # calcHemisphere.Rd
Line 18 :   \item{mask}{the binary contrast parameter that will be used as the outcome in the ROC analysis. \emph{character vector}. REQUIRED.} # calcROCthreshold.Rd
Line 17 :   \item{mask}{the binary contrast parameter that should be smoothed. \emph{character}.} # calcSmoothMask.Rd
Line 20 :   \item{mask}{the binary contrast parameter indentifing the lesion at time[1]. \code{character} or \code{NULL} if no mask is available.} # calcTableHypoReperf.Rd
Line 15 :   \item{mask}{the binary contrast parameter indicating the brain. \emph{character} or \code{NULL} if no mask is available.} # calcTableLesion.Rd
Line 15 :   \item{mask}{the binary contrast parameter(s). \emph{character vector}.} # constCompressMRIaggr.Rd
Line 13 :   \item{mask}{the binary contrast parameter or a vector indicating the observations to be keeped. \emph{character} or \emph{logical vector} with length equal to the number of observations in \code{object}. REQUIRED.} # constReduceMRIaggr.Rd
Line 16 :   \item{mask}{the binary contrast parameter indicating the lesion. \emph{character}. REQUIRED.} # plotLesion3D.Rd
Line 18 :   \item{mask}{the binary contrast parameter indicating the lesion. \emph{character}. REQUIRED.} # plotTableLesion.Rd
Line 17 :   \item{mask}{should the mask be considered as an available contrast parameter ? \emph{logical}.} # selectParameter.Rd

#### 109- maskN #### 
Line 14 :   \item{maskN}{the binary contrast parameter indicating the lesion. \emph{character}. REQUIRED.} # calcTableLesion.Rd

#### 110- max_groups #### 
Line 14 :   \item{max_groups}{the maximum number of groups. \emph{postive integer}.} # calcGroupsCoords.Rd
Line 15 :   \item{max_groups}{the maximum number of groups. \emph{postive integer}.} # calcGroupsCoords_cpp.Rd
Line 13 :   \item{max_groups}{the maximum number of groups. \emph{postive integer}.} # calcGroupsW.Rd
Line 13 :   \item{max_groups}{the maximum number of groups. \emph{postive integer}.} # calcGroupsW_cpp.Rd

#### 111- method #### 
Line 20 :   \item{method}{the distance measure to be used. \emph{character}. This must be one of \code{"euclidean"}, \code{"maximum"}, \code{"minkowski"} or \code{"greatcircle"}.}   # calcW.Rd
Line 23 :   \item{method}{the correlation coefficient which is to be computed. Can be \code{"pearson"}, \code{"kendall"} or \code{"spearman"}.} # heatmapMRIaggr.Rd
Line 21 :   \item{method}{the name of the function that called the initializer. \emph{character}.} # initCol.Rd
Line 12 :   \item{method}{the name of the function that called the initializer. \emph{character}.} # initFilter.Rd
Line 21 :   \item{method}{the name of the function that called the initializer. \emph{character}.} # initGR.Rd
Line 22 :   \item{method}{the name of the function that called the initializer. \emph{character}.} # initIndex.Rd
Line 12 :   \item{method}{the name of the function that called the initializer. \emph{character}.} # initNeighborhood.Rd
Line 20 :   \item{method}{the name of the function that called the initializer. \emph{character}.} # initNum.Rd
Line 23 :   \item{method}{the name of the function that called the initializer. \emph{character}.} # initParameter.Rd
Line 23 :   \item{method}{the name of the function that called the initializer. \emph{character}.} # initWindow.Rd

#### 112- mfrow #### 
Line 19 :   \item{mfrow}{the division of the device in plot region. \emph{numeric vector of size 2}.} # initDisplayWindow.Rd
Line 20 :   \item{mfrow}{the division of the device in plot region. \emph{numeric vector of size 2}.} # initWindow.Rd
Line 69 :   \item{mfrow}{the division of the device in plot region. \emph{numeric vector of size 2} or \code{NULL} leading automatic adjustment.} # multiplot.Rd

#### 113- mgp #### 
Line 33 :   \item{mgp}{the margin line for the axis title, axis labels and axis line. \emph{positive numeric vector of size 3}.} # boxplotMask.Rd
Line 34 :   \item{mgp}{the margin line for the axis title, axis labels and axis line. \emph{positive numeric vector of size 3}.} # calcSigmaGR.Rd
Line 32 :   \item{mgp}{the margin line for the axis title, axis labels and axis line. \emph{positive numeric vector of size 3}.} # heatmapMRIaggr.Rd
Line 23 :   \item{mgp}{the margin line for the axis title, axis labels and axis line. \emph{positive numeric vector of size 3}.} # initDisplayWindow.Rd
Line 71 :   \item{mgp}{the margin line for the axis title, axis labels and axis line. \emph{positive numeric vector of size 3}.} # multiplot.Rd
Line 33 :   \item{mgp}{the margin line for the axis title, axis labels and axis line. \emph{positive numeric vector of size 3}.} # plotDistClass.Rd
Line 27 :   \item{mgp}{the margin line for the axis title, axis labels and axis line. \emph{positive numeric vector of size 3}.} # plotTableLesion.Rd

#### 114- midplane #### 
Line 19 :   \item{midplane}{the position of the mid-sagittal plan. \emph{data.frame}.}   # MRIaggr-class.Rd
Line 45 :   \item{midplane}{should the mid-saggital plan be displayed ? \emph{logical}.} # multiplot.Rd

#### 115- min_dist #### 
Line 14 :   \item{min_dist}{if the distance between the new point and the initial point is inferior to \code{min_dist}, then the definition of the region ends. \emph{numeric}. Only active if \code{sequential} is \code{TRUE}.} # outline.Rd
Line 34 :   \item{min_dist}{if the distance between the new point and the initial point is inferior to \code{min_dist}, then the definition of the region ends. \emph{numeric}. Only active if \code{sequential} is \code{TRUE}.} # outlineMRIaggr.Rd

#### 116- mu #### 
Line 15 :   \item{mu}{should the centering values for the normalization be returned. \emph{logical}.  Active only if \code{type} is \code{"slice"} or \code{"3slices"}.} # selectNormalization.Rd

#### 117- mu_type #### 
Line 16 :   \item{mu_type}{the type of centering. Can be \code{"mean"} or \code{"median"}. } # calcNormalization.Rd

#### 118- n #### 
Line 12 :   \item{n}{maximum number of points to define the outline. \emph{integer}.} # outline.Rd
Line 32 :   \item{n}{maximum number of points to define the outline. \emph{integer}.} # outlineMRIaggr.Rd

#### 119- n.plot #### 
Line 19 :   \item{n.plot}{the number of images to display. \emph{integer}.} # initWindow.Rd

#### 120- n.points #### 
Line 30 :   \item{n.points}{the number of points that represent the mid-saggital plan to computed. \code{positive integer}.} # calcHemisphere.Rd

#### 121- na.rm #### 
Line 16 :   \item{na.rm}{should observations with missing values be removed ? \emph{logical}.} # array2df.Rd
Line 21 :   \item{na.rm}{should observations with missing values in their neighborhood be set to NA ? Otherwise the ponderation is adjusted. \emph{logical}.} # calcFilter.Rd
Line 27 :   \item{na.rm}{should observations with missing values be removed ? \emph{logical}.} # selectContrast.Rd

#### 122- na.value #### 
Line 15 :   \item{na.value}{the value with which NA values are replaced. \emph{numeric} or NA.}     # readMRI.Rd

#### 123- na_rm #### 
Line 15 :   \item{na_rm}{should the observations with missing values in their neighborhood be removed ? Otherwise the ponderation is adjusted.}   # filtrage2D_cpp.Rd
Line 14 :   \item{na_rm}{should the observations with missing values in their neighborhood be removed ? Otherwise the ponderation is adjusted.}   # filtrage2Dmed_cpp.Rd
Line 18 :   \item{na_rm}{should the observations with missing values in their neighborhood be removed ? Otherwise the ponderation is adjusted.}   # filtrage3D_cpp.Rd
Line 16 :   \item{na_rm}{should the observations with missing values in their neighborhood be removed ? Otherwise the ponderation is adjusted.}   # filtrage3Dmed_cpp.Rd

#### 124- name #### 
Line 14 :   \item{name}{the name of the element storing \code{value}. \emph{character}. REQUIRED.} # affectDescStats.Rd
Line 14 :   \item{name}{the name of the element to select. \emph{character} or \code{NULL} leading to select all available elements.} # selectDescStats.Rd

#### 125- name_newparam #### 
Line 14 :   \item{name_newparam}{the name of the contrast parameter to which cooresponds \code{array}. \emph{character}.} # array2df.Rd
Line 16 :   \item{name_newparam}{the name of the new distance parameters. \emph{character vector}.} # calcDistMask.Rd
Line 22 :   \item{name_newparam}{the name of the new parameters. \emph{character vector}.} # calcFilter.Rd
Line 27 :   \item{name_newparam}{the name of the new parameters. \emph{character vector}.} # calcRegionalIntensity.Rd
Line 36 :   \item{name_newparam}{the name of the new parameters. \code{character}.} # calcThresholdMRIaggr.Rd
Line 23 :   \item{name_newparam}{the name of the new paramaters containing the probabilistic segmentation. \emph{character vector of size 3}.}   # calcTissueType.Rd
Line 42 :   \item{name_newparam}{the name of the new parameter. \emph{character}.} # outlineMRIaggr.Rd

#### 126- names_coords #### 
Line 15 :   \item{names_coords}{the name of the coordinates. \emph{character vector}.} # array2df.Rd

#### 127- Neighborhood #### 
Line 19 :   \item{Neighborhood}{the type of neighborhood. \emph{character}.} # calcDistMask.Rd
Line 13 :   \item{Neighborhood}{the type of neighborhood. \emph{character}.} # calcGroupsCoords.Rd
Line 13 :   \item{Neighborhood}{the type of neighborhood. \emph{character}.} # calcGroupsCoords_cpp.Rd
Line 11 :   \item{Neighborhood}{the name of neighborhood configuration. \cr Any of \code{"2D_N4"}, \code{"2D_N8"} \code{"3D_N6"} \code{"3D_N10"} \code{"3D_N18"} \code{"3D_N26"}.} # initNeighborhood.Rd
Line 18 :   \item{Neighborhood}{the type of neighborhood used to defined the edges. \emph{character}.} # plotLesion3D.Rd

#### 128- Neighborhood_2D #### 
Line 21 :   \item{Neighborhood_2D}{the type of 2D neighborhood. \emph{character}.} # calcSmoothMask.Rd

#### 129- Neighborhood_3D #### 
Line 25 :   \item{Neighborhood_3D}{the type of 3D neighborhood. \emph{character}.} # calcSmoothMask.Rd

#### 130- Neighborhood_V #### 
Line 31 :   \item{Neighborhood_V}{the type of neighborhood to use for the spatial regularization. \emph{character}.} # calcSmoothMask.Rd

#### 131- niter #### 
Line 16 :   \item{niter}{the number of iterations used by \code{\link{mritc.bayes}}. \emph{positive integer}.} # calcTissueType.Rd

#### 132- nnei #### 
Line 17 :   \item{nnei}{the number of neighbors. \emph{positive integer}.} # calcTissueType.Rd

#### 133- norm #### 
Line 19 :   \item{norm}{should the filtered correspond to a weighted mean over site ? (or a weighted sum). \emph{logical}.} # calcFilter.Rd

#### 134- norm_mu #### 
Line 23 :   \item{norm_mu}{the type of centering to apply on the parameter values. \emph{character}.} # boxplotMask.Rd
Line 22 :   \item{norm_mu}{the type of centering to apply on the parameter values. \emph{character}.} # calcTableHypoReperf.Rd
Line 49 :   \item{norm_mu}{the type of centering to apply on the parameter values. \emph{character}.} # multiplot.Rd
Line 22 :   \item{norm_mu}{the type of centering to apply on the parameter values. \emph{character}.} # plotDistClass.Rd
Line 25 :   \item{norm_mu}{the type of centering to apply on the parameter values. \emph{character}. See the details section.} # selectContrast.Rd

#### 135- norm_sigma #### 
Line 24 :   \item{norm_sigma}{the type of scaling to apply on the parameter values. \emph{character}.} # boxplotMask.Rd
Line 23 :   \item{norm_sigma}{the type of scaling to apply on the parameter values. \emph{character}.} # calcTableHypoReperf.Rd
Line 50 :   \item{norm_sigma}{the type of scaling to apply on the parameter values. \emph{character}.} # multiplot.Rd
Line 23 :   \item{norm_sigma}{the type of scaling to apply on the parameter values. \emph{character}.} # plotDistClass.Rd
Line 26 :   \item{norm_sigma}{the type of scaling to apply on the parameter values. \emph{character}. See the details section.} # selectContrast.Rd

#### 136- normalization #### 
Line 17 :   \item{normalization}{the normalization values for the contrast parameters. \emph{list}.}   # MRIaggr-class.Rd

#### 137- num #### 
Line 21 :   \item{num}{the slices to consider. \emph{numeric vector} or \code{NULL}. REQUIRED.} # boxplotMask.Rd
Line 15 :   \item{num}{the slices to extract. \emph{numeric vector} or \code{NULL}.} # calcControlateral.Rd
Line 15 :   \item{num}{the slices to use. \emph{numeric vector} or \code{NULL}.} # calcDistTissues.Rd
Line 19 :   \item{num}{the slices to use. \emph{numeric vector} or \code{NULL}.} # calcHemisphere.Rd
Line 17 :   \item{num}{the slices to use. \emph{numeric vector} or \code{NULL}.} # calcRegionalIntensity.Rd
Line 22 :   \item{num}{the slices to use. \emph{numeric vector} or \code{NULL}.} # calcW.Rd
Line 15 :   \item{num}{the slices to extract. \emph{numeric vector} or \code{NULL}.} # Carto3D2MRIaggr.Rd
Line 19 :   \item{num}{the slices to use. \emph{numeric vector} or \code{NULL}.} # heatmapMRIaggr.Rd
Line 14 :   \item{num}{the slices to display. \emph{numeric vector} or \code{NULL}.} # initIndex.Rd
Line 16 :   \item{num}{the slice numbers to check or initialize. \emph{numeric vector} or \code{NULL}. See the details section.} # initNum.Rd
Line 43 :   \item{num}{the slices to display. \emph{numeric vector} or \code{NULL}.} # multiplot.Rd
Line 20 :   \item{num}{the slices on which the outline will be drawn. \emph{numeric vector} or \code{NULL}.} # outlineMRIaggr.Rd
Line 20 :   \item{num}{the slices to use. \emph{numeric vector} or \code{NULL}.} # plotDistClass.Rd
Line 19 :   \item{num}{the slices to display. \emph{numeric vector} or \code{NULL}.} # plotTableLesion.Rd
Line 20 :   \item{num}{the slices to extract. \emph{numeric vector} or \code{NULL}.} # selectContrast.Rd
Line 19 :   \item{num}{the slices to extract. \emph{numeric vector} or \code{NULL}.} # selectCoords.Rd
Line 17 :   \item{num}{the slices to use. \emph{numeric vector} or \code{NULL}.} # selectDescStats.Rd
Line 13 :   \item{num}{the slices to consider. \emph{numeric vector} or \code{NULL}.} # selectN.Rd
Line 17 :   \item{num}{the slices to extract. \emph{numeric vector} or \code{NULL}.} # selectNormalization.Rd

#### 138- num.main #### 
Line 79 :   \item{num.main}{should the slice number be written over each plot. \emph{logical}.} # multiplot.Rd

#### 139- object #### 
Line 12 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # affectClinic.Rd
Line 12 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # affectContrast.Rd
Line 12 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # affectDescStats.Rd
Line 12 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # affectHemisphere.Rd
Line 12 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # affectNormalization.Rd
Line 13 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # affectTable.Rd
Line 17 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # boxplotMask.Rd
Line 17 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # calcBrainMask.Rd
Line 13 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # calcControlateral.Rd
Line 14 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # calcDistMask.Rd
Line 12 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # calcDistTissues.Rd
Line 16 :   \item{object}{an \code{array} or an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # calcFilter.Rd
Line 13 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # calcGroupsMask.Rd
Line 17 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # calcHemisphere.Rd
Line 14 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # calcNormalization.Rd
Line 15 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # calcRegionalIntensity.Rd
Line 16 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # calcROCthreshold.Rd
Line 16 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # calcSmoothMask.Rd
Line 14 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # calcTableHypoReperf.Rd
Line 13 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # calcTableLesion.Rd
Line 21 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # calcThresholdMRIaggr.Rd
Line 14 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # calcTissueType.Rd
Line 18 :   \item{object}{a \emph{data.frame} containing the coordinates of the observations or an \code{object} of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # calcW.Rd
Line 12 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # constCompressMRIaggr.Rd
Line 12 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # constReduceMRIaggr.Rd
Line 17 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # heatmapMRIaggr.Rd
Line 12 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}} or \code{NULL} .} # initIndex.Rd
Line 15 :   \item{object}{an object of class \code{\linkS4class{Carto3D}} or \code{\linkS4class{MRIaggr}}.} # initNum.Rd
Line 14 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}.} # initParameter.Rd
Line 40 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}} or \code{\linkS4class{Carto3D}} or a 3 column \code{data.frame} containing the coordinates of the observations in columns. REQUIRED.} # multiplot.Rd
Line 17 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # outlineMRIaggr.Rd
Line 17 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # plotDistClass.Rd
Line 15 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # plotLesion3D.Rd
Line 17 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # plotTableLesion.Rd
Line 12 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}.} # pointsHemisphere.Rd
Line 12 :   \item{object}{an \code{object} of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # selectClinic.Rd
Line 18 :   \item{object}{an \code{object} of class \code{\linkS4class{Carto3D}} or \code{\linkS4class{MRIaggr}}. REQUIRED.} # selectContrast.Rd
Line 16 :   \item{object}{an \code{object} of class \code{\linkS4class{Carto3D}} or \code{\linkS4class{MRIaggr}}. REQUIRED.} # selectCoords.Rd
Line 15 :   \item{object}{an \code{object} of class \code{\linkS4class{Carto3D}} or \code{\linkS4class{MRIaggr}}. REQUIRED.} # selectDefault_value.Rd
Line 13 :   \item{object}{an \code{object} of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # selectDescStats.Rd
Line 12 :   \item{object}{an \code{object} of class \code{\linkS4class{Carto3D}} or \code{\linkS4class{MRIaggr}}. REQUIRED.} # selectHemispheres.Rd
Line 12 :   \item{object}{an \code{object} of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
} # selectHistory.Rd
Line 15 :   \item{object}{an \code{object} of class \code{\linkS4class{Carto3D}} or \code{\linkS4class{MRIaggr}}. REQUIRED.}
} # selectIdentifier.Rd
Line 12 :   \item{object}{an \code{object} of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
} # selectMidplane.Rd
Line 12 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # selectN.Rd
Line 13 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # selectNormalization.Rd
Line 15 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # selectParameter.Rd
Line 12 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # selectTable.Rd
Line 15 :   \item{object}{an \code{object} of class \code{\linkS4class{Carto3D}} or \code{\linkS4class{MRIaggr}}. REQUIRED.}
} # selectVoxelDim.Rd
Line 12 :   \item{object}{an \code{object} of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # selectVoxelSize.Rd
Line 11 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # summary.MRIaggr.Rd
Line 12 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # supprContrast.Rd
Line 12 :   \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.} # supprDescStats.Rd

#### 140- operator #### 
Line 21 :   \item{operator}{should the median absolute deviation be used to estimte the variability of the group contrast (\code{"mad"}) or the standard deviation (\code{"sd"}).}   # GRalgo.Rd

#### 141- operator_index1 #### 
Line 35 :   \item{operator_index1}{the operator to apply between the index1 observations and the outlined observations. Can be \code{"none"} \code{"intersection"} \code{"difference"} or \code{"union"}.} # outlineMRIaggr.Rd

#### 142- overwrite #### 
Line 15 :   \item{overwrite}{if clinical parameters with the same names are already stored in \code{object@clinic}, can they be overwritten ? \emph{logical}.}   # affectClinic.Rd
Line 16 :   \item{overwrite}{if a contrast parameters with the same names are already stored in \code{object@data}, can they be overwritten ? \emph{logical}.}   # affectContrast.Rd
Line 15 :   \item{overwrite}{if an element with the same name is already stored in \code{object@ls_descStats}, can it be overwritten ? \emph{logical}.}   # affectDescStats.Rd
Line 14 :   \item{overwrite}{if the characteristics of a mid-saggital plan are already stored in \code{object}, can they be overwritten ? \emph{logical}.}   # affectHemisphere.Rd
Line 14 :   \item{overwrite}{if normalization values are already stored in \code{object}, can they be overwritten ? \emph{logical}.}   # affectNormalization.Rd
Line 16 :   \item{overwrite}{if tables are already stored in \code{object@normalization}, can they be overwritten ? \emph{logical}.}   # affectTable.Rd
Line 42 :   \item{overwrite}{if a mask is already stored in \code{object@data}, can it be overwritten ? \emph{logical}.} # calcBrainMask.Rd
Line 23 :   \item{overwrite}{if contrast parameters with the same names are already stored in \code{object@data}, can they be overwritten ? \emph{logical}.} # calcControlateral.Rd
Line 23 :   \item{overwrite}{if contrast parameters with the same names are already stored in \code{object@data}, can they be overwritten ? \emph{logical}.} # calcDistMask.Rd
Line 26 :   \item{overwrite}{if contrast parameters with the same names are already stored in \code{object@data}, can they be overwritten ? \emph{logical}.}   # calcFilter.Rd
Line 22 :   \item{overwrite}{if spatial groups are already stored in \code{object@ls_descStats}, can they be overwritten ? \emph{logical}.}   # calcGroupsMask.Rd
Line 43 :   \item{overwrite}{if a midplane is already stored in \code{object@midplane}, can it be overwritten ? \emph{logical}.}   # calcHemisphere.Rd
Line 24 :   \item{overwrite}{if normalization values are already stored in \code{object@normalization}, can they be overwritten ? \emph{logical}.}   # calcNormalization.Rd
Line 29 :   \item{overwrite}{if contrast parameters with the same names are already stored in \code{object} can they be overwritten ? \emph{logical}.}   # calcRegionalIntensity.Rd
Line 33 :   \item{overwrite}{if a threshold analysis is already stored in \code{object@ls_descStats}, can it be overwritten ? \emph{logical}.}   # calcROCthreshold.Rd
Line 35 :   \item{overwrite}{if a mask is already stored in \code{object@data}, can it be overwritten ? \emph{logical}.}   # calcSmoothMask.Rd
Line 28 :   \item{overwrite}{if reperfusion or hypoperfusion values are already stored in \code{object@table_reperfusion} or \code{object@table_hypofusion}, can they be overwritten ? \emph{logical}.}   # calcTableHypoReperf.Rd
Line 20 :   \item{overwrite}{if a lesion table is already stored in \code{object@table_lesion}, can it be overwritten ? \emph{logical}.}   # calcTableLesion.Rd
Line 39 :   \item{overwrite}{contrast parameters with the same names are already stored in \code{object@data}, can it be overwritten ? \emph{logical}.}   # calcThresholdMRIaggr.Rd
Line 25 :   \item{overwrite}{if tissue types are already stored in \code{object@data}, can they be overwritten ? \emph{logical}.}   # calcTissueType.Rd
Line 31 :   \item{overwrite}{if a neighborhood matrix is already stored in \code{object@ls_descStats}, can it be overwritten ? \emph{logical}.} # calcW.Rd
Line 44 :   \item{overwrite}{if a contrast parameter with the same names is already stored in \code{object@data}, can it be overwritten ? \emph{logical}.}   # outlineMRIaggr.Rd

#### 143- p #### 
Line 14 :   \item{p}{the penalization factor. \emph{positive numeric}.} # calcHemi_cpp.Rd
Line 20 :   \item{p}{the type of distance for the penalization. \emph{positive numeric}.} # calcHemisphere.Rd

#### 144- p_data #### 
Line 13 :   \item{p_data}{spatial dimensions of the data.} # filtrage3D_cpp.Rd
Line 12 :   \item{p_data}{spatial dimensions of the data.} # filtrage3Dmed_cpp.Rd

#### 145- p_operateur #### 
Line 15 :   \item{p_operateur}{spatial dimensions of the filter.} # filtrage3D_cpp.Rd
Line 14 :   \item{p_operateur}{spatial dimensions of the filter.} # filtrage3Dmed_cpp.Rd

#### 146- palette #### 
Line 17 :   \item{palette}{the colors or the palette to use when associating colors to intensities. \emph{character vector} or \emph{character}.} # initCol.Rd
Line 12 :   \item{palette}{the colors or the palette to use when associating colors to observation intensities. \emph{character vector}.} # legendMRI.Rd
Line 55 :   \item{palette}{the colors or the palette to use when associating colors to observation intensities. \emph{character vector} or \emph{character}.} # multiplot.Rd
Line 27 :   \item{palette}{the colors or the palette to use when associating colors to observation intensities. \emph{character vector} or \emph{character}.} # outlineMRIaggr.Rd
Line 15 :   \item{palette}{the colors or the palette to use when associating colors to observation intensities. \emph{character vector}.} # plotMRI.Rd

#### 147- param #### 
Line 14 :   \item{param}{the names of the contrast parameters. \emph{character vector} or \code{NULL} leading to use the name of the \code{value} argument.} # affectContrast.Rd
Line 19 :   \item{param}{the contrast parameter(s) associated with the lesion mask(s). \emph{character vector}. REQUIRED.} # boxplotMask.Rd
Line 18 :   \item{param}{the contrast parameter(s) that should be used to identify the brain observations. \emph{character vector}. REQUIRED} # calcBrainMask.Rd
Line 14 :   \item{param}{the contrast parameters to normalize. \emph{character vector}. REQUIRED.} # calcControlateral.Rd
Line 13 :   \item{param}{the contrast parameters to consider. \emph{character vector}. REQUIRED.} # calcDistTissues.Rd
Line 17 :   \item{param}{the contrast parameter to be filtered. \emph{character vector}. REQUIRED.} # calcFilter.Rd
Line 18 :   \item{param}{the contrast parameter that should be used to distinguish the two hemispheres. \emph{character}. REQUIRED.} # calcHemisphere.Rd
Line 15 :   \item{param}{the contrast parameters to normalize. \emph{character vector}. REQUIRED.} # calcNormalization.Rd
Line 16 :   \item{param}{the contrast parameter(s) from which the regional parameter(s) will be computed. \emph{character vector}. REQUIRED.} # calcRegionalIntensity.Rd
Line 17 :   \item{param}{the contrast parameter(s) that should be used to identify the observations inside the mask. \emph{character vector}. REQUIRED.} # calcROCthreshold.Rd
Line 15 :   \item{param}{the perfusion parameter(s). \emph{character vector}. REQUIRED.} # calcTableHypoReperf.Rd
Line 23 :   \item{param}{the contrast parameters to be thresholded. \code{character vector}. REQUIRED.} # calcThresholdMRIaggr.Rd
Line 15 :   \item{param}{the contrast parameter that should be used to distinguish the WM, the GM and the CSF. \emph{character}. REQUIRED.} # calcTissueType.Rd
Line 14 :   \item{param}{the contrast parameter. \emph{character}. REQUIRED.} # constCarto3D.Rd
Line 14 :   \item{param}{the contrast parameters to load in the new MRIaggr object. \emph{character vector} or \code{NULL}.} # constCompressMRIaggr.Rd
Line 14 :   \item{param}{the contrast parameter(s). \emph{character vector} or \code{NULL}. REQUIRED.} # constMRIaggr.Rd
Line 17 :   \item{param}{the names of directories containing the images. \emph{character vector} or \code{NULL} leading to use all directories.} # constSweave.Rd
Line 18 :   \item{param}{the contrast parameters used to compute the correlations. \emph{character vector}. REQUIRED.} # heatmapMRIaggr.Rd
Line 14 :   \item{param}{the contrast parameter to display. \emph{character}.} # initCol.Rd
Line 15 :   \item{param}{the contrast parameters to check or initialize. \emph{character vector} or \code{NULL}.} # initParameter.Rd
Line 41 :   \item{param}{the contrast parameter to display. \emph{character}. REQUIRED.} # multiplot.Rd
Line 18 :   \item{param}{the contrast parameter map on which the outline will be drawn. \emph{character}. REQUIRED.} # outlineMRIaggr.Rd
Line 18 :   \item{param}{the contrast parameter to display. \emph{character}. REQUIRED.} # plotDistClass.Rd
Line 13 :   \item{param}{the clinical parameters to extract. \emph{character vector} or \code{NULL} leading to extract all the clinical parameters.} # selectClinic.Rd
Line 19 :   \item{param}{the contrast parameters to extract. \emph{character vector} or \code{NULL}.} # selectContrast.Rd
Line 16 :   \item{param}{the contrast parameters for which the reference values should be returned. \emph{character vector} or \code{NULL} leading to extract reference values for all available contrast parameters.} # selectDefault_value.Rd
Line 19 :   \item{param}{the contrast parameters for which the normalization values should be extracted. \emph{character vector} or \code{NULL} indicating all available contrast parameters.} # selectNormalization.Rd
Line 12 :   \item{param}{should detailed information be printed for the contrast parameters ? \emph{logical}.} # summary.MRIaggr.Rd

#### 148- param.ref #### 
Line 17 :   \item{param.ref}{the parameter to use as a reference for the identification of the controlateral voxel. \emph{character} or \code{NULL} if no reference parameter available.} # calcControlateral.Rd

#### 149- param.update #### 
Line 26 :   \item{param.update}{which type of parameter should be stored in the object ? Any of \code{"shift"} \code{"reperf"} \code{"reperf_pc"} \code{"deperf"} \code{"deperf_pc"}.} # calcTableHypoReperf.Rd

#### 150- param_ref #### 
Line 18 :   \item{param_ref}{the parameter to be used as a reference to identify the controlateral voxel. \emph{character}.} # calcContro_cpp.Rd

#### 151- parameter #### 
Line 11 :   \item{parameter}{the name of the contrast parameter. \emph{character}.}    # Carto3D-class.Rd

#### 152- path #### 
Line 41 :   \item{path}{the directory where the plot file will be created. \emph{character}.} # boxplotMask.Rd
Line 36 :   \item{path}{the directory where the plot file will be created. \emph{character}.} # calcBrainMask.Rd
Line 37 :   \item{path}{the directory where the plot file will be created. \emph{character}.} # calcHemisphere.Rd
Line 28 :   \item{path}{the directory where the plot file will be created. \emph{character}.} # calcROCthreshold.Rd
Line 41 :   \item{path}{the directory where the plot file will be created. \emph{character}.} # calcSigmaGR.Rd
Line 40 :   \item{path}{the directory where the plot file will be created. \emph{character}.} # heatmapMRIaggr.Rd
Line 13 :   \item{path}{the directory where the plot file will be created. \emph{character}.} # initDisplayWindow.Rd
Line 14 :   \item{path}{the directory where the plot file will be created. \emph{character}.} # initWindow.Rd
Line 91 :   \item{path}{the directory where the plot file will be created. \emph{character}.} # multiplot.Rd
Line 44 :   \item{path}{the directory where the plot file will be created. \emph{character}.} # plotDistClass.Rd
Line 37 :   \item{path}{the directory where the plot file will be created. \emph{character}.} # plotTableLesion.Rd

#### 153- pch #### 
Line 15 :   \item{pch}{the symbol with which the observations will be displayed. \emph{positive integer}.} # initCol.Rd
Line 57 :   \item{pch}{the symbol with which the observations will be displayed. \emph{positive integer} or \code{NULL} leading to use the \code{image} function instead of \code{plot}.} # multiplot.Rd
Line 16 :   \item{pch}{the symbol with which the observations will be displayed. \emph{positive integer}.} # outline.Rd
Line 38 :   \item{pch}{the symbol with which the observations will be displayed. \emph{positive integer}.} # outlineMRIaggr.Rd
Line 35 :   \item{pch}{the symbol with which the distribution will be displayed. \emph{positive integer}.} # plotDistClass.Rd
Line 21 :   \item{pch}{the symbol with which the observations will be displayed. \emph{positive integer}.} # plotMRI.Rd

#### 154- pch.default #### 
Line 19 :   \item{pch.default}{the default label used to plot the observations. \emph{numeric}.} # initIndex.Rd

#### 155- pch.NA #### 
Line 61 :   \item{pch.NA}{the label to use to plot the NAs. \emph{postive integer}.} # multiplot.Rd
Line 25 :   \item{pch.NA}{the label to use to plot the NAs. \emph{postive integer}.} # plotMRI.Rd

#### 156- penalty #### 
Line 22 :   \item{penalty}{the type of objective function. Can be \code{"symmetry"} or \code{"asymmetry"}.} # calcHemisphere.Rd

#### 157- performance #### 
Line 14 :   \item{performance}{an object of class \code{performance} can be supplied instead of arguments \code{x} and \code{y}.}   # calcAUPRC.Rd

#### 158- plot #### 
Line 25 :   \item{plot}{should the results be plotted ? \emph{logical}.} # calcBrainMask.Rd
Line 32 :   \item{plot}{should the results be plotted ? \emph{logical}.} # calcHemisphere.Rd
Line 21 :   \item{plot}{the type of the graphic to display? \emph{character} or \code{FALSE}. See the details section.} # calcROCthreshold.Rd

#### 159- points.values #### 
Line 24 :   \item{points.values}{should the correlation values be printed on the plot ? \emph{logical}.} # heatmapMRIaggr.Rd

#### 160- pos_default_value #### 
Line 16 :   \item{pos_default_value}{the coordinates of the observation that contains the reference value. \emph{numeric vector} with length the number of dimension of array.} # constCarto3D.Rd
Line 16 :   \item{pos_default_value}{the coordinates of the observations that contains the reference value. \emph{numeric vector}.}    # constMRIaggr.Rd

#### 161- power #### 
Line 13 :   \item{power}{the power of the kernel. \emph{numeric}.} # EDK.Rd

#### 162- power_EDK #### 
Line 23 :   \item{power_EDK}{the power of the kernel. \emph{postive numeric}.} # calcRegionalIntensity.Rd

#### 163- pty #### 
Line 21 :   \item{pty}{the type of plot region to be used. Can be \code{"s"} or \code{"m"}.} # initDisplayWindow.Rd
Line 72 :   \item{pty}{the type of plot region to be used. Can be \code{"s"} or \code{"m"}.}   # multiplot.Rd

#### 164- px_hemiL #### 
Line 11 :   \item{px_hemiL}{the coordinates and the contrast of the observations in the left hemisphere for a given slice. \emph{matrix}.} # calcHemi_cpp.Rd

#### 165- px_hemiR #### 
Line 12 :   \item{px_hemiR}{the coordinates and the contrast of the observations in the right hemisphere for a given slice. \emph{matrix}.} # calcHemi_cpp.Rd

#### 166- px_max #### 
Line 25 :   \item{px_max}{the maximum number of points that can be ploted. \emph{integer}.} # plotLesion3D.Rd

#### 167- quantiles #### 
Line 18 :   \item{quantiles}{the quantiles values to display on the legend. \emph{numeric vector of size 5} or \code{NULL}.} # legendMRI.Rd

#### 168- quantiles.legend #### 
Line 82 :   \item{quantiles.legend}{should the quantiles of the data be displayed on the legend ? \emph{logical}.} # multiplot.Rd

#### 169- radius #### 
Line 26 :   \item{radius}{the radius of spheres. \emph{numeric}. See \code{plot3d} for more details.} # plotLesion3D.Rd

#### 170- range #### 
Line 17 :   \item{range}{the range of acceptable contrast values for the growing region group. \code{numeric vector of size 2}.} # calcGR.Rd
Line 19 :   \item{range}{the range of acceptable contrast values for the growing region group. \code{numeric vector of size 2}.} # calcSigmaGR.Rd
Line 16 :   \item{range}{the range of acceptable contrast values for the growing region group. \code{numeric vector of size 2}.} # GRalgo.Rd
Line 14 :   \item{range}{the range of acceptable contrast values for the growing region group. \code{numeric vector of size 2}.} # initGR.Rd

#### 171- range.coords #### 
Line 15 :   \item{range.coords}{the maximum coordinate in each dimension to be considered. \emph{numeric vector} with length equal to the number of columns of \code{coords}.} # df2array.Rd

#### 172- range.seed #### 
Line 18 :   \item{range.seed}{the range of acceptable contrast values for the seeds. \code{numeric vector of size 2}.} # calcGR.Rd
Line 20 :   \item{range.seed}{the range of acceptable contrast values for the seeds. \code{numeric vector of size 2}.} # calcSigmaGR.Rd
Line 15 :   \item{range.seed}{the range of acceptable contrast values for the seeds. \code{numeric vector of size 2}.} # initGR.Rd

#### 173- reorient #### 
Line 18 :   \item{reorient}{is a logical variable (default = TRUE) that enforces Qform/Sform transformations. Only active if \code{format} equals \code{nifti}.} # readMRI.Rd

#### 174- res #### 
Line 43 :   \item{res}{the nominal resolution in ppi which will be recorded in the bitmap file. \emph{positive integer} or \code{NA}.} # boxplotMask.Rd
Line 38 :   \item{res}{the nominal resolution in ppi which will be recorded in the bitmap file. \emph{positive integer}.} # calcBrainMask.Rd
Line 39 :   \item{res}{the nominal resolution in ppi which will be recorded in the bitmap file. \emph{positive integer}.} # calcHemisphere.Rd
Line 30 :   \item{res}{the nominal resolution in ppi which will be recorded in the bitmap file. \emph{positive integer}.} # calcROCthreshold.Rd
Line 43 :   \item{res}{the nominal resolution in ppi which will be recorded in the bitmap file. \emph{positive integer}.} # calcSigmaGR.Rd
Line 42 :   \item{res}{the nominal resolution in ppi which will be recorded in the bitmap file. \emph{positive integer}.} # heatmapMRIaggr.Rd
Line 17 :   \item{res}{the nominal resolution in ppi which will be recorded in the bitmap file. \emph{positive integer}.} # initDisplayWindow.Rd
Line 18 :   \item{res}{the nominal resolution in ppi which will be recorded in the bitmap file. \emph{positive integer}.} # initWindow.Rd
Line 93 :   \item{res}{the nominal resolution in ppi which will be recorded in the bitmap file. \emph{positive integer}.} # multiplot.Rd
Line 46 :   \item{res}{the nominal resolution in ppi which will be recorded in the bitmap file. \emph{positive integer}.} # plotDistClass.Rd
Line 39 :   \item{res}{the nominal resolution in ppi which will be recorded in the bitmap file. \emph{positive integer}.} # plotTableLesion.Rd

#### 175- rm.2Dhole #### 
Line 22 :   \item{rm.2Dhole}{should the 2D wholes inside the mask be removed ? \emph{logical}.} # calcSmoothMask.Rd

#### 176- rm.3Dhole #### 
Line 26 :   \item{rm.3Dhole}{should the 3D wholes inside the mask be removed ? \emph{logical}.} # calcSmoothMask.Rd

#### 177- rm.array #### 
Line 17 :   \item{rm.array}{should the object on which \emph{array} argument points be removed form the global environment ? \code{logical}.} # constCarto3D.Rd

#### 178- rm.Carto3D #### 
Line 13 :   \item{rm.Carto3D}{should the object on which the \code{ls.Carto3D} argument points be removed form the global environment ? \code{logical}.} # Carto3D2MRIaggr.Rd

#### 179- rm.CSF #### 
Line 18 :   \item{rm.CSF}{should the cerebral spinal fluid observations be excluded ? \emph{logical}.} # calcNormalization.Rd
Line 25 :   \item{rm.CSF}{should the cerebral spinal fluid observations be excluded ? \code{logical} or \emph{character}.} # calcThresholdMRIaggr.Rd

#### 180- rm.GM #### 
Line 19 :   \item{rm.GM}{should the grey matter observations be excluded ? \emph{logical}.} # calcNormalization.Rd

#### 181- rm.ls.array #### 
Line 20 :   \item{rm.ls.array}{should the object on which the \code{ls.array} argument points be removed form the global environment ? \emph{logical}.} # constMRIaggr.Rd

#### 182- rm.WM #### 
Line 20 :   \item{rm.WM}{should the white matter observations be excluded ? \emph{logical}.} # calcNormalization.Rd

#### 183- row.norm #### 
Line 27 :   \item{row.norm}{should the resulting matrix be row-normalized ? \code{TRUE/FALSE}.} # calcW.Rd

#### 184- sample #### 
Line 12 :   \item{sample}{the weight of each voxel in the computation of the barycenter. \emph{positive numeric}.} # calcRadius_cpp.Rd

#### 185- scale #### 
Line 26 :   \item{scale}{should the contrast parameters be scaled ? \emph{logical}.} # boxplotMask.Rd
Line 21 :   \item{scale}{should the contrast be scaled ? \code{logical}.}   # calcGR.Rd
Line 23 :   \item{scale}{should the contrast be scaled ? \code{logical}.}   # calcSigmaGR.Rd
Line 21 :   \item{scale}{should the contrast parameters be scaled ? \emph{logical}.} # heatmapMRIaggr.Rd
Line 16 :   \item{scale}{the scaling factor to convert \code{height} and \code{height} to standard unit. \emph{numeric}.} # initDisplayWindow.Rd
Line 18 :   \item{scale}{should the contrast be scaled ? \code{logical}.}   # initGR.Rd

#### 186- sd.robust #### 
Line 23 :   \item{sd.robust}{should the median absolute deviation be used to estimte the variability of the group contrast, or the standard deviation ? \code{logical}.} # calcGR.Rd
Line 25 :   \item{sd.robust}{should the median absolute deviation be used to estimte the variability of the group contrast (\code{TRUE}), or the standard deviation (\code{FALSE}) ? \code{logical}.} # calcSigmaGR.Rd

#### 187- sd_data #### 
Line 13 :   \item{sd_data}{the standard deviation of the parameter. \emph{numeric}.} # calcHemi_cpp.Rd

#### 188- seed #### 
Line 15 :   \item{seed}{the index of the initial seeds or a binary indicator of the initial seeds. \code{positive integer vector} or \code{logical vector}. REQUIRED.} # calcGR.Rd
Line 17 :   \item{seed}{the index of the initial seeds or a binary indicator of the initial seeds. \code{positive integer vector} or \code{logical vector}. REQUIRED.} # calcSigmaGR.Rd
Line 31 :   \item{seed}{the index of the seeds for the growing region algorithm . \code{positive integer vector}. } # calcThresholdMRIaggr.Rd
Line 14 :   \item{seed}{the index of the initial seeds or a binary indicator of the initial seeds. \code{positive integer vector} or \code{logical vector}.} # GRalgo.Rd
Line 13 :   \item{seed}{the index of the initial seeds or a binary indicator of the initial seeds. \code{positive integer vector} or \code{logical vector}.} # initGR.Rd

#### 189- sep #### 
Line 18 :   \item{sep}{the separator between the parameter names and the time points. \emph{character}.} # calcTableHypoReperf.Rd

#### 190- sequential #### 
Line 13 :   \item{sequential}{should the region edge be updated on the graphical device after each point ? \emph{logical}.} # outline.Rd
Line 33 :   \item{sequential}{should the region edge be updated on the graphical device after each point ? \emph{logical}.} # outlineMRIaggr.Rd

#### 191- sigma #### 
Line 15 :   \item{sigma}{the \code{sigma_max} that have been used in the GR algorithm. \code{positive numeric vector}.} # calcCriteriaGR.Rd
Line 18 :   \item{sigma}{the sequence of maximum admissible values for the group variability \code{positive numeric vector}. REQUIRED.} # calcSigmaGR.Rd
Line 16 :   \item{sigma}{should the scaling values for the normalization be returned. \emph{logical}.  Active only if \code{type} is \code{"slice"} or \code{"3slices"}.} # selectNormalization.Rd

#### 192- sigma_max #### 
Line 16 :   \item{sigma_max}{the maximum admissible value for the variability of the group contrast. \code{positive numeric}. REQUIRED.} # calcGR.Rd
Line 15 :   \item{sigma_max}{the maximum admissible value for the variability of the group contrast. \code{positive numeric}.} # GRalgo.Rd

#### 193- sigma_type #### 
Line 17 :   \item{sigma_type}{the type of scaling. Can be \code{"sd"} or \code{"mad"}.} # calcNormalization.Rd

#### 194- size #### 
Line 14 :   \item{size}{should the values in the table correspond to a number of voxels (\code{FALSE}) or a volume (\code{TRUE}).} # selectTable.Rd

#### 195- size_2Dgroup #### 
Line 20 :   \item{size_2Dgroup}{the minimum size of the 2D groups. \emph{positive integer} or \code{"unique"}.} # calcSmoothMask.Rd

#### 196- size_3Dgroup #### 
Line 24 :   \item{size_3Dgroup}{the minimum size of the 3D groups. \emph{positive integer} or \code{"unique"}.} # calcSmoothMask.Rd

#### 197- skull.n_groups #### 
Line 30 :   \item{skull.n_groups}{the number of groups to use in the kmeans algorithm to obtain the skull.} # calcBrainMask.Rd

#### 198- skull.param #### 
Line 29 :   \item{skull.param}{the parameter used to identify the skull. \emph{character}.} # calcBrainMask.Rd

#### 199- slice_var #### 
Line 19 :   \item{slice_var}{the type of slice to extract. \code{"i"} for sagittal, \code{"j"} for coronal and \code{"k"} for transverse. \emph{character}.}   # initNum.Rd
Line 47 :   \item{slice_var}{the type of view to use. \code{"i"} for sagittal view, \code{"j"} for coronal view and \code{"k"} for transverse view. \emph{character}.}   # multiplot.Rd
Line 22 :   \item{slice_var}{the type of slice to extract. \code{"i"} for sagittal, \code{"j"} for coronal and \code{"k"} for transverse. \emph{character}.}   # selectContrast.Rd
Line 22 :   \item{slice_var}{the type of slice to extract. \code{"i"} for sagittal, \code{"j"} for coronal and \code{"k"} for transverse. \emph{character}.}   # selectCoords.Rd
Line 18 :   \item{slice_var}{the type of slice to extract. \code{"i"} for sagittal, \code{"j"} for coronal and \code{"k"} for transverse. \emph{character}.}   # selectDescStats.Rd

#### 200- spatial_res #### 
Line 17 :   \item{spatial_res}{a dilatation factor for the coordinates. \emph{positive numeric vector of size 3}.} # calcDistMask.Rd
Line 16 :   \item{spatial_res}{a dilatation factor for the coordinates. \emph{positive numeric vector of size 3}.} # calcGroupsMask.Rd
Line 20 :   \item{spatial_res}{a dilatation factor for the coordinates. \emph{positive numeric vector of size 3}.} # calcRegionalIntensity.Rd
Line 21 :   \item{spatial_res}{a dilatation factor for the coordinates. \emph{positive numeric vector of size 3}.} # calcW.Rd
Line 20 :   \item{spatial_res}{a dilatation factor for the coordinates. \emph{positive numeric vector of size 3}.} # plotLesion3D.Rd
Line 18 :   \item{spatial_res}{a dilatation factor for the coordinates. \emph{positive numeric vector of size 3}.} # selectCoords.Rd

#### 201- SPM #### 
Line 17 :   \item{SPM}{is a logical variable (default = FALSE) that forces the voxel data values to be rescaled using the funused1 ANALYZE header field. This is an undocumented convention of ANALYZE files processed using the Statistical Parametric Mapping (SPM) software. Only active if \code{format} equals \code{analyse}.} # readMRI.Rd

#### 202- step #### 
Line 19 :   \item{step}{the step between two consecutive breaks. \emph{numeric}.}   # GRalgo.Rd

#### 203- sub #### 
Line 19 :   \item{sub}{if \code{TRUE}, use the higher resolution model; otherwise, use the whole voxel method. \emph{logical}.} # calcTissueType.Rd

#### 204- subdivisions #### 
Line 13 :   \item{subdivisions}{the maximum number of subintervals used for the integration. \emph{positive integer}.} # calcAUPRC.Rd

#### 205- subsection #### 
Line 21 :   \item{subsection}{the names of subsections for the latex document. \emph{character vector} or \code{NULL} leading to use \code{param} for naming the subsections.} # constSweave.Rd

#### 206- subset #### 
Line 12 :   \item{subset}{the subset of observations to use. \emph{positive integer vector} or \code{NULL} leading to use all observations.} # calcGroupsW.Rd
Line 12 :   \item{subset}{the subset of observations to use. \emph{positive integer vector}.} # calcGroupsW_cpp.Rd
Line 21 :   \item{subset}{the subset of observations to use. \emph{positive integer vector} or \code{NULL} leading to use all observations.} # calcHemisphere.Rd
Line 24 :   \item{subset}{the subset of observations to use. \emph{positive integer vector} or \code{NULL} leading to use all observations.} # calcW.Rd
Line 28 :   \item{subset}{the subset of observations to extract. \emph{positive integer vector} or \code{NULL} leading to use all observations} # selectContrast.Rd
Line 21 :   \item{subset}{the subset of observations to extract. \emph{positive integer vector} or \code{NULL} leading to use all observations} # selectCoords.Rd
Line 15 :   \item{subset}{the subset of observations to consider. \emph{positive integer vector} or \code{NULL} leading to consider all observations.} # selectN.Rd

#### 207- subset_bary #### 
Line 14 :   \item{subset_bary}{an indicator of the observations that should be kept ? \emph{logical vector}.} # calcRadius_cpp.Rd

#### 208- subset_W #### 
Line 15 :   \item{subset_W}{the subset of observations to select. \emph{integer vector} or \code{NULL} leading to select all the observations.} # selectDescStats.Rd

#### 209- subsubsection #### 
Line 24 :   \item{subsubsection}{A list containing the names of the subsubsections for the latex document. \emph{list of character vector} or \code{NULL} leading to no subsubsection.} # constSweave.Rd

#### 210- symetrie #### 
Line 15 :   \item{symetrie}{the type of objective function. \code{TRUE} correspond to \code{"symmetry"} and \code{FALSE} to \code{"asymmetry"}.} # calcHemi_cpp.Rd

#### 211- table #### 
Line 18 :   \item{table}{a list of data.frame to display in the table format. \emph{list of data.frame} or \code{NULL} if there is no table to display.} # constSweave.Rd

#### 212- table_hypoperfusion #### 
Line 23 :   \item{table_hypoperfusion}{the volumic hypoperfusion data. \emph{data.frame}.}   # MRIaggr-class.Rd

#### 213- table_lesion #### 
Line 21 :   \item{table_lesion}{the vertical distribution of the lesion volumes. \emph{data.frame}.}   # MRIaggr-class.Rd

#### 214- table_reperfusion #### 
Line 22 :   \item{table_reperfusion}{the volumic reperfusion data. \emph{data.frame}.}   # MRIaggr-class.Rd

#### 215- test #### 
Line 17 :   \item{test}{should the slice numbers be checked ? \emph{logical}.  } # initNum.Rd
Line 16 :   \item{test}{should the parameters be checked ? \emph{logical}.} # initParameter.Rd

#### 216- th.breaks #### 
Line 21 :   \item{th.breaks}{the number of thresholds to use. \emph{postive integer}.} # calcBrainMask.Rd

#### 217- th.select_optima #### 
Line 23 :   \item{th.select_optima}{the rank of the optimum to retain. \emph{postive integer}.} # calcBrainMask.Rd

#### 218- th.smoothing #### 
Line 22 :   \item{th.smoothing}{should the derivative be smoothed ? \emph{logical}.} # calcBrainMask.Rd

#### 219- th.upper #### 
Line 24 :   \item{th.upper}{should the observations above the selected threshold be retained ? Else the observations bellow will the selected threshold be retained. \emph{logical}.} # calcBrainMask.Rd

#### 220- threshold #### 
Line 13 :   \item{threshold}{observations with a \code{sample} value below the value of \code{threshold} are discarded. \emph{numeric}.}  # calcRadius_cpp.Rd
Line 17 :   \item{threshold}{the value of the hypoperfusion thresholds. \emph{numeric vector}.} # calcTableHypoReperf.Rd
Line 26 :   \item{threshold}{the thresholds to be used for the discretization of the contrast parameter. \code{numeric vector}.} # calcThresholdMRIaggr.Rd
Line 16 :   \item{threshold}{the value above which the local mean of the binary parameters is assigned to 1 (and otherwise to 0). \emph{numeric between 0 and 1}.} # constCompressMRIaggr.Rd

#### 221- time #### 
Line 16 :   \item{time}{two time points. \emph{character vector of size 2}. REQUIRED.} # calcTableHypoReperf.Rd

#### 222- title #### 
Line 35 :   \item{title}{the title of the latex document. \emph{character}.} # constSweave.Rd

#### 223- tol #### 
Line 14 :   \item{tol}{numeric precision for the consistency check. \emph{positive numeric}.} # Carto3D2MRIaggr.Rd
Line 18 :   \item{tol}{numeric precision for the consistency check. \emph{positive numeric}.} # constMRIaggr.Rd

#### 224- trace #### 
Line 16 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # affectClinic.Rd
Line 17 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # affectContrast.Rd
Line 16 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # affectDescStats.Rd
Line 15 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # affectHemisphere.Rd
Line 15 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # affectNormalization.Rd
Line 17 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # affectTable.Rd
Line 40 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # calcBrainMask.Rd
Line 23 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # calcContro_cpp.Rd
Line 21 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # calcControlateral.Rd
Line 21 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.}  # calcDistMask.Rd
Line 24 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # calcFilter.Rd
Line 24 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.}   # calcGR.Rd
Line 15 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # calcGroupsCoords.Rd
Line 16 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # calcGroupsCoords_cpp.Rd
Line 20 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # calcGroupsMask.Rd
Line 41 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # calcHemisphere.Rd
Line 22 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # calcNormalization.Rd
Line 15 :   \item{trace}{should the radius of the spatial group be printed ? \emph{logical}.}  # calcRadius_cpp.Rd
Line 26 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # calcRegionalIntensity.Rd
Line 34 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # calcROCthreshold.Rd
Line 32 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.}   # calcSigmaGR.Rd
Line 33 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # calcSmoothMask.Rd
Line 25 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # calcTableHypoReperf.Rd
Line 18 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # calcTableLesion.Rd
Line 37 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # calcThresholdMRIaggr.Rd
Line 22 :   \item{trace}{indicate the level of output as the algorithm runs. \emph{logical}.} # calcTissueType.Rd
Line 29 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # calcW.Rd
Line 16 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # Carto3D2MRIaggr.Rd
Line 17 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # constCompressMRIaggr.Rd
Line 19 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # constMRIaggr.Rd
Line 28 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # constSweave.Rd
Line 19 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.}   # initGR.Rd
Line 41 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # outlineMRIaggr.Rd
Line 14 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # supprContrast.Rd
Line 14 :   \item{trace}{should the execution of the function be traced ? \emph{logical}.} # supprDescStats.Rd

#### 225- trim #### 
Line 30 :   \item{trim}{the length in mm with which the imported images will be cropped (left, bottom, right top). \emph{list of numeric vector of size 4}.} # constSweave.Rd

#### 226- trim.legend #### 
Line 33 :   \item{trim.legend}{the length in mm with which the legend of the images will be cropped. \emph{numeric vector of size 4}.} # constSweave.Rd

#### 227- type #### 
Line 15 :   \item{type}{the type of volumic information. Can be \code{"lesion"} \code{"reperfusion"} \code{"hypoperfusion"}. REQUIRED.} # affectTable.Rd
Line 19 :   \item{type}{the method to use. Can be \code{"threshold"} or \code{"kmeans"}.} # calcBrainMask.Rd
Line 16 :   \item{type}{the method used to compute the controlateral correspondent of each voxel. Can be \code{"mean"}, \code{"median"} or \code{"1NN_penalised"}.} # calcControlateral.Rd
Line 25 :   \item{type}{the type of plot to display. Any of \code{"image"} or \code{"image.plot"} or \code{FALSE} meaning no plot.} # heatmapMRIaggr.Rd
Line 34 :   \item{type}{the type of plot to display. \emph{character}. See \code{\link{plot.default}} for more details.} # plotDistClass.Rd
Line 27 :   \item{type}{the type of item to plot. \emph{character}. See \code{plot3d} for more details.} # plotLesion3D.Rd
Line 20 :   \item{type}{the type of plot to display. Can be \code{"matplot"} or \code{"evolution"}.} # plotTableLesion.Rd
Line 14 :   \item{type}{the type of normalization. Must be on of \code{"global"}, \code{"slice"}, \code{"3slices"} or \code{NULL} leading to select all types of normalization.} # selectNormalization.Rd
Line 16 :   \item{type}{the type of parameter to select. \emph{character}. See the details section.} # selectParameter.Rd
Line 13 :   \item{type}{the table to extract. Can be \code{"lesion"}, \code{"reperfusion"} or \code{"hypoperfusion"}. REQUIRED.} # selectTable.Rd

#### 228- type.breaks #### 
Line 20 :   \item{type.breaks}{should the break points be equally space according the range of data values (\code{"range"}), centered (\code{"range_center"}) or correspond to the quantile values (\code{"quantile"}) ?} # initCol.Rd
Line 54 :   \item{type.breaks}{should the break points be equally space according the range of data values (\code{"range"}), centered (\code{"range_center"}) or correspond to the quantile values (\code{"quantile"}) ?.} # multiplot.Rd

#### 229- type.plot #### 
Line 24 :   \item{type.plot}{the type of plot to be displayed. Can be \code{"plot3d"} or \code{"shapelist3d"}.} # plotLesion3D.Rd

#### 230- type_med #### 
Line 21 :   \item{type_med}{should the median controlateral value be used ? \emph{logical}.} # calcContro_cpp.Rd

#### 231- type_moy #### 
Line 20 :   \item{type_moy}{should the mean controlateral value be used ? \emph{logical}.} # calcContro_cpp.Rd

#### 232- type_NN #### 
Line 22 :   \item{type_NN}{should the closest controlateral voxel according to the reference parameter be used ? \emph{logical}.} # calcContro_cpp.Rd

#### 233- unit #### 
Line 42 :   \item{unit}{the units in which \code{height} and \code{width} are given. \emph{character}.} # boxplotMask.Rd
Line 37 :   \item{unit}{the units in which \code{height} and \code{width} are given. \emph{character}.} # calcBrainMask.Rd
Line 38 :   \item{unit}{the units in which \code{height} and \code{width} are given. \emph{character}.} # calcHemisphere.Rd
Line 29 :   \item{unit}{the units in which \code{height} and \code{width} are given. \emph{character}.} # calcROCthreshold.Rd
Line 42 :   \item{unit}{the units in which \code{height} and \code{width} are given. \emph{character}.} # calcSigmaGR.Rd
Line 41 :   \item{unit}{the units in which \code{height} and \code{width} are given. \emph{character}.} # heatmapMRIaggr.Rd
Line 17 :   \item{unit}{the units in which \code{height} and \code{width} are given. \emph{character}.} # initWindow.Rd
Line 92 :   \item{unit}{the units in which \code{height} and \code{width} are given. \emph{character}.} # multiplot.Rd
Line 45 :   \item{unit}{the units in which \code{height} and \code{width} are given. \emph{character}.} # plotDistClass.Rd
Line 38 :   \item{unit}{the units in which \code{height} and \code{width} are given. \emph{character}.} # plotTableLesion.Rd
Line 13 :   \item{unit}{should the unit be returned ? \emph{logical}.} # selectVoxelSize.Rd

#### 234- unit_angle #### 
Line 29 :   \item{unit_angle}{the unit in which the angle is given. Can be \code{"radian"} or \code{"degree"}.} # calcHemisphere.Rd

#### 235- update.object #### 
Line 41 :   \item{update.object}{should the resulting mask be stored in \code{object} as a \code{"mask"} parameter ? \emph{logical}.} # calcBrainMask.Rd
Line 22 :   \item{update.object}{should the resulting controlateral parameters be stored in \code{object} ? \emph{logical}.} # calcControlateral.Rd
Line 22 :   \item{update.object}{should the resulting distance parameters be stored in \code{object} ? \emph{logical}.} # calcDistMask.Rd
Line 25 :   \item{update.object}{should the resulting filtered parameters be stored in \code{object} ? \emph{logical}.} # calcFilter.Rd
Line 21 :   \item{update.object}{should the resulting spatial groups be stored in \code{object} ? \emph{logical}.} # calcGroupsMask.Rd
Line 42 :   \item{update.object}{should the resulting midplane be stored in \code{object} ? \emph{logical}.} # calcHemisphere.Rd
Line 23 :   \item{update.object}{should the resulting normalization values be stored in \code{object} ? \emph{logical}.} # calcNormalization.Rd
Line 28 :   \item{update.object}{should the resulting regional parameters be stored in \code{object} ? \emph{logical}.} # calcRegionalIntensity.Rd
Line 32 :   \item{update.object}{should the resulting threshold analysis be stored in \code{object@ls_descStats} ? \emph{logical}.} # calcROCthreshold.Rd
Line 34 :   \item{update.object}{should the resulting regularized mask be stored in \code{object} ? \emph{logical}.} # calcSmoothMask.Rd
Line 27 :   \item{update.object}{should the resulting values be stored in \code{object} ? \emph{logical}.} # calcTableHypoReperf.Rd
Line 19 :   \item{update.object}{should the resulting lesion table be stored in \code{object} ? \emph{logical}.} # calcTableLesion.Rd
Line 38 :   \item{update.object}{should the resulting thresholded parameters be stored in \code{object} ? \emph{logical}.} # calcThresholdMRIaggr.Rd
Line 24 :   \item{update.object}{should the resulting tissue types be stored in \code{object} ? \emph{logical}.} # calcTissueType.Rd
Line 30 :   \item{update.object}{should the resulting neighborhood matrix be stored in \code{object} ? \emph{logical}.} # calcW.Rd
Line 43 :   \item{update.object}{should the new parameter be stored in \code{object} ? \emph{logical}.} # outlineMRIaggr.Rd

#### 236- upper #### 
Line 25 :   \item{upper}{should the entire matrix (\code{NULL}) or only the upper-triagonal (\code{TRUE}) or only the lower-triagonal (\code{FALSE}) values be calculated ?} # calcW.Rd

#### 237- value #### 
Line 13 :   \item{value}{the clinical data. A one row \emph{data.frame}. REQUIRED.} # affectClinic.Rd
Line 13 :   \item{value}{the value of each contrast parameter (in columns) at each voxel (in rows). \emph{data.frame}. REQUIRED.} # affectContrast.Rd
Line 13 :   \item{value}{any R object. REQUIRED.} # affectDescStats.Rd
Line 13 :   \item{value}{a \emph{list} of \emph{data.frame}. Names must be among \code{"midplane"}, \code{"hemispheres"} \code{"data"}. See the Details section. REQUIRED.}   # affectHemisphere.Rd
Line 13 :   \item{value}{the normalisation values. A \emph{list} of \emph{data.frame}. REQUIRED.} # affectNormalization.Rd
Line 14 :   \item{value}{the volumic information. \emph{data.frame}. REQUIRED.} # affectTable.Rd
Line 13 :   \item{value}{the name of the parameter(s) that should be removed. \emph{character vector}. REQUIRED.} # supprContrast.Rd
Line 13 :   \item{value}{the name of the element(s) that should be removed. \emph{character vector}. REQUIRED.} # supprDescStats.Rd

#### 238- var_ref #### 
Line 19 :   \item{var_ref}{the variance of the reference parameter. \emph{numeric}.} # calcContro_cpp.Rd

#### 239- Vbackground_max #### 
Line 30 :   \item{Vbackground_max}{background observations with a proportion of neighbors belonging to the mask higher than \code{Vbackground_max} are affected to the mask. \emph{numeric between 0 and 1}.} # calcSmoothMask.Rd

#### 240- Vec_data #### 
Line 12 :   \item{Vec_data}{vector of data to which the filter will be applied.} # filtrage3D_cpp.Rd
Line 11 :   \item{Vec_data}{vector of data to which the filter will be applied.} # filtrage3Dmed_cpp.Rd

#### 241- Vec_operateur #### 
Line 14 :   \item{Vec_operateur}{vector representing the filter to be applied.} # filtrage3D_cpp.Rd
Line 13 :   \item{Vec_operateur}{vector representing the filter to be applied.} # filtrage3Dmed_cpp.Rd

#### 242- Vmask_min #### 
Line 29 :   \item{Vmask_min}{mask observations with a proportion of neighbors belonging to the mask lower than \code{Vmask_min} are affected to the background. \emph{numeric between 0 and 1}.} # calcSmoothMask.Rd

#### 243- voxelDim #### 
Line 13 :   \item{voxelDim}{the spatial dimensions of the lattice containing the observations. \emph{data.frame} with 1 row and 3 columns named \code{"i"} \code{"j"} \code{"k"}.}   # Carto3D-class.Rd
Line 13 :   \item{voxelDim}{the spatial dimensions of the lattice containing the observations. \emph{data.frame}.}   # MRIaggr-class.Rd

#### 244- W #### 
Line 14 :   \item{W}{the neighborhood matrix. \code{dgCMatrix} or \code{NULL} leading to not compute the d1 criterion.} # calcCriteriaGR.Rd
Line 14 :   \item{W}{the neighborhood matrix. \code{dgCMatrix}. REQUIRED.} # calcGR.Rd
Line 18 :   \item{W}{the neighborhood matrix. \emph{dgCMatrix} or \code{"ifany"} leading to use the neighborhood matrix stored in the \code{object} if any and else computate this matrix.} # calcGroupsMask.Rd
Line 11 :   \item{W}{the neighborhood matrix. \emph{dgCMatrix}. REQUIRED.} # calcGroupsW.Rd
Line 11 :   \item{W}{the neighborhood matrix. \emph{dgCMatrix}.} # calcGroupsW_cpp.Rd
Line 19 :   \item{W}{the neighborhood matrix. \emph{dgCMatrix} or \code{"ifany"} leading to use the neighborhood matrix stored in the \code{object} if possible, else to compute this matrix.} # calcRegionalIntensity.Rd
Line 16 :   \item{W}{the neighborhood matrix. \code{dgCMatrix}. REQUIRED.} # calcSigmaGR.Rd
Line 30 :   \item{W}{the neighborhood matrix. \code{dgCMatrix}.} # calcThresholdMRIaggr.Rd
Line 13 :   \item{W}{the neighborhood matrix. \code{dgCMatrix}.} # GRalgo.Rd
Line 12 :   \item{W}{the neighborhood matrix. \code{dgCMatrix}.} # initGR.Rd

#### 245- W.distband #### 
Line 33 :   \item{W.distband}{only distances smaller than delta are recorded (for the computation of \code{W}). \code{postive numeric}.} # calcThresholdMRIaggr.Rd

#### 246- W.spatial_res #### 
Line 34 :   \item{W.spatial_res}{a dilatation factor for the coordinates (for the computation of \code{W}). \emph{positive numeric vector of size 3}.} # calcThresholdMRIaggr.Rd

#### 247- w_contrast #### 
Line 20 :   \item{w_contrast}{should the influence of each neighbor be ponderated by the difference in signal with the considered observation ? \emph{logical}.} # calcFilter.Rd
Line 14 :   \item{w_contrast}{should the influence of each neighbor be ponderated by the difference in signal with the considered observation ?} # filtrage2D_cpp.Rd
Line 17 :   \item{w_contrast}{should the influence of each neighbor be ponderated by the difference in signal with the considered observation ?} # filtrage3D_cpp.Rd

#### 248- what #### 
Line 14 :   \item{what}{an object whose mode will give the mode of the vector to be read, or a character vector of length one describing the mode: one of "numeric", "double", "integer", "int", "logical", "complex", "character", "raw". Only active if \code{format} equals \code{rawb.gz}.} # readMRI.Rd

#### 249- width #### 
Line 39 :   \item{width}{the width of the device used to export the plot. \emph{postive numeric}.} # boxplotMask.Rd
Line 34 :   \item{width}{the width of the device used to export the plot. \emph{postive numeric}.} # calcBrainMask.Rd
Line 35 :   \item{width}{the width of the device used to export the plot. \emph{postive numeric}.} # calcHemisphere.Rd
Line 26 :   \item{width}{the width of the device used to export the plot. \emph{postive numeric}.} # calcROCthreshold.Rd
Line 39 :   \item{width}{the width of the device used to export the plot. \emph{postive numeric}.} # calcSigmaGR.Rd
Line 29 :   \item{width}{the width of each image relative to the linewidth. \emph{list of positive numeric}.} # constSweave.Rd
Line 38 :   \item{width}{the width of the device used to export the plot. \emph{postive numeric}.} # heatmapMRIaggr.Rd
Line 14 :   \item{width}{the width of the device used to export the plot in inches. \emph{postive numeric}.} # initDisplayWindow.Rd
Line 15 :   \item{width}{the width of the device used to export the plot. \emph{postive numeric}.} # initWindow.Rd
Line 89 :   \item{width}{the width of the device used to export the plot. \emph{postive numeric}.} # multiplot.Rd
Line 42 :   \item{width}{the width of the device used to export the plot. \emph{postive numeric}.} # plotDistClass.Rd
Line 35 :   \item{width}{the width of the device used to export the plot. \emph{postive numeric}.} # plotTableLesion.Rd

#### 250- width.legend #### 
Line 32 :   \item{width.legend}{the width of the legend image relative to the linewidth. \emph{numeric between 0 and 1}.} # constSweave.Rd

#### 251- window #### 
Line 29 :   \item{window}{the type of device on which the plot will be displayed. \emph{logical}, \code{NULL} or \code{character}.} # boxplotMask.Rd
Line 32 :   \item{window}{the type of device on which the plot will be displayed. \emph{logical}, \code{NULL} or \code{character}.} # calcBrainMask.Rd
Line 33 :   \item{window}{the type of device on which the plot will be displayed. \emph{logical}, \code{NULL} or \code{character}.} # calcHemisphere.Rd
Line 24 :   \item{window}{the type of device on which the plot will be displayed. \emph{logical}, \code{NULL} or \code{character}.} # calcROCthreshold.Rd
Line 37 :   \item{window}{the type of device on which the plot will be displayed. \emph{logical}, \code{NULL} or \code{character}.} # calcSigmaGR.Rd
Line 29 :   \item{window}{the type of device on which the plot will be displayed. \emph{logical}, \code{NULL} or \code{character}.} # heatmapMRIaggr.Rd
Line 11 :   \item{window}{the type of device on which the plot will be displayed. \emph{logical}, \code{NULL} or \code{character}.} # initDisplayWindow.Rd
Line 12 :   \item{window}{the type of device on which the plot will be displayed. \emph{logical}, \code{NULL} or \code{character}.} # initWindow.Rd
Line 67 :   \item{window}{the type of device on which the plot will be displayed. \emph{logical}, \code{NULL} or \code{character}.} # multiplot.Rd
Line 30 :   \item{window}{the type of device on which the plot will be displayed. \emph{logical}, \code{NULL} or \code{character}.} # plotDistClass.Rd
Line 22 :   \item{window}{the type of device on which the plot will be displayed. \emph{logical}, \code{NULL} or \code{character}.} # plotTableLesion.Rd

#### 252- x #### 
Line 11 :   \item{x}{the biomarker values. \emph{numeric vector}. REQUIRED.} # calcAUPRC.Rd
Line 11 :   \item{x}{the data on which the kernel will be applied. \emph{numeric} or \emph{numeric vector}.} # EDK.Rd

#### 253- x.legend #### 
Line 34 :   \item{x.legend}{the x coordinates of the legend. \emph{numeric} or \emph{character}.} # boxplotMask.Rd
Line 37 :   \item{x.legend}{the x coordinates of the legend. \emph{numeric} or \emph{character}.} # plotDistClass.Rd

#### 254- xlab #### 
Line 76 :   \item{xlab}{a title for the x axis. \emph{character}.} # multiplot.Rd
Line 26 :   \item{xlab}{a title for the x axis. \emph{character}.} # plotMRI.Rd

#### 255- xlim #### 
Line 21 :   \item{xlim}{the x limits of the plot. \emph{numeric vector of size 2}.} # initWindow.Rd
Line 63 :   \item{xlim}{the x limits of the plot. \emph{numeric vector of size 2} or \code{NULL} leading to automatic setting of the x limits.} # multiplot.Rd
Line 24 :   \item{xlim}{the y limits of the plot. \emph{numeric vector of size 2} or \code{NULL} leading to automatic setting of the x limits.} # outlineMRIaggr.Rd
Line 21 :   \item{xlim}{the x limits of the plot. \emph{numeric vector of size 2} or \code{NULL} leading to automatic setting of the x limits.} # plotLesion3D.Rd
Line 19 :   \item{xlim}{the x limits of the plot. \emph{numeric vector of size 2}.} # plotMRI.Rd

#### 256- y #### 
Line 12 :   \item{y}{the class labels. \emph{numeric vector}, \emph{character vector} or \emph{logical vector}. REQUIRED.} # calcAUPRC.Rd

#### 257- y.legend #### 
Line 35 :   \item{y.legend}{the y coordinates of the legend. \emph{numeric} or \emph{character}.} # boxplotMask.Rd
Line 38 :   \item{y.legend}{the y coordinates of the legend. \emph{numeric} or \emph{character}.} # plotDistClass.Rd

#### 258- ylab #### 
Line 77 :   \item{ylab}{a title for the y axis. \emph{character}.} # multiplot.Rd
Line 27 :   \item{ylab}{a title for the y axis. \emph{character}.} # plotMRI.Rd

#### 259- ylim #### 
Line 30 :   \item{ylim}{the y limits of the plot. \emph{numeric vector of size 2} or \code{NULL} leading to automatic setting of the y limits.} # boxplotMask.Rd
Line 22 :   \item{ylim}{the y limits of the plot. \emph{numeric vector of size 2}.} # initWindow.Rd
Line 64 :   \item{ylim}{the y limits of the plot. \emph{numeric vector of size 2} or \code{NULL} leading to automatic setting of the y limits.} # multiplot.Rd
Line 25 :   \item{ylim}{the y limits of the plot. \emph{numeric vector of size 2} or \code{NULL} leading to automatic setting of the y limits.} # outlineMRIaggr.Rd
Line 28 :   \item{ylim}{the y limits of the plot. \emph{numeric vector of size 2} or \code{NULL} leading to automatic setting of the y limits.} # plotDistClass.Rd
Line 22 :   \item{ylim}{the y limits of the plot. \emph{numeric vector of size 2} or \code{NULL} leading to automatic setting of the y limits.} # plotLesion3D.Rd
Line 20 :   \item{ylim}{the y limits of the plot. \emph{numeric vector of size 2}.} # plotMRI.Rd

#### 260- zlim #### 
Line 23 :   \item{zlim}{the z limits of the plot. \emph{numeric vector of size 2} or \code{NULL} leading to automatic setting of the z limits.} # plotLesion3D.Rd


 
################################# 91 files .Rd #################################
affectClinic.Rd
affectContrast.Rd
affectDescStats.Rd
affectHemisphere.Rd
affectNormalization.Rd
affectTable.Rd
array2df.Rd
boxplotMask.Rd
calcAUPRC.Rd
calcBrainMask.Rd
calcContro_cpp.Rd
calcControlateral.Rd
calcCriteriaGR.Rd
calcDistMask.Rd
calcDistTissues.Rd
calcFilter.Rd
calcGR.Rd
calcGroupsCoords.Rd
calcGroupsCoords_cpp.Rd
calcGroupsMask.Rd
calcGroupsW.Rd
calcGroupsW_cpp.Rd
calcHemi_cpp.Rd
calcHemisphere.Rd
calcNormalization.Rd
calcRadius_cpp.Rd
calcRegionalIntensity.Rd
calcROCthreshold.Rd
calcSigmaGR.Rd
calcSmoothMask.Rd
calcTableHypoReperf.Rd
calcTableLesion.Rd
calcThresholdMRIaggr.Rd
calcTissueType.Rd
calcW.Rd
Carto3D-class.Rd
Carto3D2MRIaggr.Rd
constCarto3D.Rd
constCompressMRIaggr.Rd
constMRIaggr.Rd
constReduceMRIaggr.Rd
constSweave.Rd
df2array.Rd
EDK.Rd
filtrage2D_cpp.Rd
filtrage2Dmed_cpp.Rd
filtrage3D_cpp.Rd
filtrage3Dmed_cpp.Rd
GRalgo.Rd
heatmapMRIaggr.Rd
initCol.Rd
initDisplayWindow.Rd
initFilter.Rd
initGR.Rd
initIndex.Rd
initNeighborhood.Rd
initNum.Rd
initParameter.Rd
initWindow.Rd
legendMRI.Rd
MRIaggr-class.Rd
MRIaggr-package.Rd
MRIaggr.Pat1_red.Rd
multiplot.Rd
outline.Rd
outlineMRIaggr.Rd
plotDistClass.Rd
plotLesion3D.Rd
plotMRI.Rd
plotOutline.Rd
plotTableLesion.Rd
pointsHemisphere.Rd
readMRI.Rd
selectClinic.Rd
selectContrast.Rd
selectCoords.Rd
selectDefault_value.Rd
selectDescStats.Rd
selectHemispheres.Rd
selectHistory.Rd
selectIdentifier.Rd
selectMidplane.Rd
selectN.Rd
selectNormalization.Rd
selectParameter.Rd
selectTable.Rd
selectVoxelDim.Rd
selectVoxelSize.Rd
summary.MRIaggr.Rd
supprContrast.Rd
supprDescStats.Rd
