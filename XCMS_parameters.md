XCMS parameters were optimized by the IPO package, and are listed here.
Peak picking was performed using the centwave algorithm within the
findChromPeaks function, which had the following parameters: ppm = 15,
peakdwidth = 3 – 80, mzdiff = 0.00325, prefilter = 3 – 100, noise =
1000, snthresh = 10. Peaks were grouped using the groupChromPeaks
function with the sample grouping (IEM/control) as input. Retention time
correction was performed using the obiwarp algorithm within the
adjustRtime function using all default parameters except: gapInit =
1.2736 and gapExtend = 3.3336. Peaks were grouped again with the
groupChromPeaks with the same parameters as used previously. Peaks were
filled using the fillChromPeaks with all default parameters. An
intensity matrix was extracted using the featureValues function (method
= “medret”, value = “into”).
