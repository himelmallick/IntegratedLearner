# Plot the summary curves produced by an IntegratedLearner object

Plots outcome-appropriate performance summaries for the training set
and, if available, the validation set produced by an `IntegratedLearner`
fit. Depending on the outcome family, this may include ROC curves,
R-squared bar plots, multiclass one-vs-rest ROC curves, or survival AUC
/ Kaplan-Meier panels. Multiclass plots include all enabled single-layer
and fusion outputs, including stacking, majority voting, and learned
probability averaging. Majority ROC curves use vote fractions, which are
not calibrated probabilities.

## Usage

``` r
# S3 method for class 'learner'
plot(
  x,
  y = NULL,
  label_size = 8,
  label_x = 0.3,
  vjust = 0.1,
  rowwise_plot = TRUE,
  ...
)
```

## Arguments

- x:

  Fitted `IntegratedLearner` object.

- y:

  Unused (required for S3 signature).

- label_size:

  Optional numeric label size for subplot tags. Default is 8.

- label_x:

  Optional single value or vector of x positions for subplot labels,
  relative to each subplot. Defaults to 0.3 for all labels.

- vjust:

  Adjusts the vertical position of each label. More positive values move
  the label further down on the plot canvas. Can be a single value
  (applied to all labels) or a vector of values (one for each label).
  Default is 0.1.

- rowwise_plot:

  If both train and test data are available, should the train and test
  plots be arranged row-wise? Default is `TRUE`. If `FALSE`, plots are
  aligned column-wise.

- ...:

  Additional arguments (currently unused).

## Value

A list whose `$plot` entry is a ggplot2/cowplot composite object, along
with the underlying tabular data used to generate the plot.

## Examples

``` r
PRISM_MAE <- readRDS(system.file(
  "extdata", "PRISM_MAE.rds", package = "IntegratedLearner", mustWork = TRUE
))
mae <- PRISM_MAE
se <- MultiAssayExperiment::experiments(mae)[[1]]
head(SummarizedExperiment::assay(se, SummarizedExperiment::assayNames(se)[1]))
#>                                G35127      G35128      G35152       G36347
#> Granulicella_unclassified -0.05253649 -0.05127158 -0.06133085  0.004887447
#> Actinomyces_graevenitzii   1.04668500 -1.32629194 -1.51654615 -3.247989324
#> Actinomyces_johnsonii     -0.70327678 -0.41575776 -0.29326475 -0.314361595
#> Actinomyces_massiliensis  -0.56808952  0.14722099  0.05660884 -1.077235688
#> Actinomyces_naeslundii    -0.49546119 -0.15921604 -0.03146485 -0.354377267
#> Actinomyces_odontolyticus -2.60071428 -2.10648484 -3.10982881 -3.846828858
#>                                 G36348      G36349      G36350      G36351
#> Granulicella_unclassified -0.006164066 -0.06068187  0.01958055 -0.00291914
#> Actinomyces_graevenitzii  -0.717183019 -1.50427168 -1.68488463 -0.65581069
#> Actinomyces_johnsonii     -0.340485318 -0.30116753 -0.33870436 -0.37999919
#> Actinomyces_massiliensis  -0.159240362  0.06245478 -0.50807364 -0.13001064
#> Actinomyces_naeslundii    -0.139758576 -0.03970686 -0.27103949 -0.18096864
#> Actinomyces_odontolyticus -2.238237450 -0.22407833 -2.10747316  7.10740921
#>                                G36352       G36353      G36354      G36355
#> Granulicella_unclassified -0.05776143 -0.007104484 -0.02798171 -0.04964912
#> Actinomyces_graevenitzii  -1.44903659 -2.227858510 -1.16809841  0.28934691
#> Actinomyces_johnsonii     -0.33673001 -0.631955729 -0.69300749 -0.43551469
#> Actinomyces_massiliensis   4.12202564 -0.871627902 -0.47895003  0.16183585
#> Actinomyces_naeslundii    -0.07679592 -0.480124737 -0.41065892 -0.17982107
#> Actinomyces_odontolyticus -2.75380353 -2.700100958 -2.34537624 -1.94465517
#>                                 G36356      G36357      G36358      G36359
#> Granulicella_unclassified  0.001299263 -0.01103146  0.01730911 -0.07698314
#> Actinomyces_graevenitzii  -0.576026671 -0.80924151 -1.72784526 -1.85085413
#> Actinomyces_johnsonii     -0.431367226 -0.28121451 -0.31104465 -0.72086328
#> Actinomyces_massiliensis  -0.092011992 -0.20308495 -0.52853445 -0.20756435
#> Actinomyces_naeslundii    -0.234541714 -0.07794349 -0.24219244 -0.38066431
#> Actinomyces_odontolyticus -1.493820959 -2.72372647 -2.33403470 -2.60201484
#>                                 G36360     G36361       G36362      G36363
#> Granulicella_unclassified -0.005839574 0.01800834  0.009845776 -0.05970839
#> Actinomyces_graevenitzii  -0.711045786 4.42949600 -1.869001603 -1.48585999
#> Actinomyces_johnsonii     -0.344436706 8.47614049 -0.220162746 -0.31302169
#> Actinomyces_massiliensis  -0.156317389 8.86079698 -0.595762822  0.07122370
#> Actinomyces_naeslundii    -0.143879582 6.12525018 -0.147409306 -0.05206988
#> Actinomyces_odontolyticus -2.205871516 7.27449295  4.629397782  4.99684820
#>                                G36365       G36366       G36367      G36368
#> Granulicella_unclassified -0.06295331 -0.004217111 -0.005806513 -0.06003288
#> Actinomyces_graevenitzii  -1.54723231 -0.680359624  1.168632075 -1.49199722
#> Actinomyces_johnsonii     -0.27350782 -0.364193642 -0.647761278 -0.30907030
#> Actinomyces_massiliensis   0.04199397 -0.141702527 -0.859936012  0.06830073
#> Actinomyces_naeslundii    -0.01085982 -0.164484612 -0.496608761 -0.04794888
#> Actinomyces_odontolyticus -3.27165848 -2.044041844 -2.570637221  2.90391516
#>                                G36369      G36370       G36371       G36372
#> Granulicella_unclassified -0.05224506 -0.01694897 -0.004866096 -0.063277806
#> Actinomyces_graevenitzii  -1.34470364 -0.95943251 -0.692634089 -1.553369545
#> Actinomyces_johnsonii     -0.40390360 -0.82735466 -0.356290867 -0.269556429
#> Actinomyces_massiliensis   0.13845207 -0.37956896 -0.147548472  0.039071002
#> Actinomyces_naeslundii    -0.14685302 -0.55077312 -0.156242600 -0.006738816
#> Actinomyces_odontolyticus -2.20358264 -1.24493447  0.791680370  2.332149272
#>                                G36374       G36375      G36378      G36379
#> Granulicella_unclassified -0.06789735  0.009196791  0.01860708 -0.04142710
#> Actinomyces_graevenitzii  -1.67901162 -1.881276067 -1.70329633 -2.59471794
#> Actinomyces_johnsonii     -0.83150213 -0.212259971 -0.32685020 -0.22035883
#> Actinomyces_massiliensis  -0.12572112 -0.601608767 -0.51684256 -0.34483795
#> Actinomyces_naeslundii    -0.49605247 -0.139167294 -0.25867647 -0.08856765
#> Actinomyces_odontolyticus -1.69576868 -3.143183056  6.35613318 -3.56165133
#>                                G36380      G36381      G36382      G36383
#> Granulicella_unclassified -0.04110260 -0.05224506  0.03593207 -0.05938389
#> Actinomyces_graevenitzii   4.94228519  6.48054316 -1.45023842 -1.47972275
#> Actinomyces_johnsonii     -0.22431021 -0.40390360 -0.38312005 -0.31697308
#> Actinomyces_massiliensis  -0.34191498  0.13845207 -0.92963674  0.07414667
#> Actinomyces_naeslundii    -0.09268866 -0.14685302 -0.35199511 -0.05619089
#> Actinomyces_odontolyticus -3.52928540 -2.20358264 -2.94417212 -2.91563320
#>                                G36384      G36385      G36386      G36387
#> Granulicella_unclassified -0.02150616 -0.01135595 -0.05873491 -0.02904613
#> Actinomyces_graevenitzii  -2.29256218 -0.81537874 -1.46744829 -2.47343934
#> Actinomyces_johnsonii     -0.30823977 -0.27726312 -0.32487585 -0.83462299
#> Actinomyces_massiliensis  -0.73424835 -0.20600792  0.07999262 -0.92534723
#> Actinomyces_naeslundii    -0.21485434 -0.07382248 -0.06443290 -0.66707895
#> Actinomyces_odontolyticus -4.04232502  3.86252714 -2.85090133 -2.72536270
#>                                G36389      G36390      G36392      G36394
#> Granulicella_unclassified  0.02814425  0.01633563 -0.02116737 -0.07503618
#> Actinomyces_graevenitzii  -1.59753200 -1.74625695 -1.03921653 -1.81403073
#> Actinomyces_johnsonii     -0.28828676 -0.29919049 -0.77598663  5.27878701
#> Actinomyces_massiliensis  -0.99978808 -0.53730337 -0.41756761 -0.19002652
#> Actinomyces_naeslundii    -0.25309097 -0.22982943 -0.49720004  2.51945997
#> Actinomyces_odontolyticus -3.72095455 -2.43113250  3.65622523  6.23559551
#>                                G36395      G36396      G48777      G48780
#> Granulicella_unclassified  0.01958055 -0.02580121 -0.06133085 -0.05776143
#> Actinomyces_graevenitzii   3.91094074 -2.41206702 -0.43401430  1.28416581
#> Actinomyces_johnsonii     -0.33870436 -0.87413686 -0.29326475  6.10453778
#> Actinomyces_massiliensis  -0.50807364 -0.89611750  0.05660884  2.67421984
#> Actinomyces_naeslundii    -0.27103949 -0.70828901 -0.03146485 -0.07679592
#> Actinomyces_odontolyticus -2.10747316  4.91222115  1.64870221  5.85895296
#>                                G48781        G48783      G48784      G48786
#> Granulicella_unclassified  0.02087852 -0.0006476922 -0.05675489  0.03171367
#> Actinomyces_graevenitzii  -1.66033570 -0.6128500661  1.96124137  5.51434470
#> Actinomyces_johnsonii     -0.35450991  2.0687611082 -0.65190875 -0.33175202
#> Actinomyces_massiliensis  -0.49638175 -0.1095498280 -0.60608817 -0.96763539
#> Actinomyces_naeslundii    -0.28752351 -0.2098156784  1.24942852 -0.29842203
#> Actinomyces_odontolyticus -1.97800942 -1.6880165658  2.93588947  3.91180455
#>                                 G48787      G48789      G48790      G48791
#> Granulicella_unclassified -0.005515081 -0.01103146 -0.01070696 -0.01070696
#> Actinomyces_graevenitzii  -0.704908554  4.39332806  2.07609172 -0.80310427
#> Actinomyces_johnsonii     -0.348388093 -0.28121451 -0.28516590 -0.28516590
#> Actinomyces_massiliensis  -0.153394417 -0.20308495 -0.20016198 -0.20016198
#> Actinomyces_naeslundii    -0.148000588 -0.07794349 -0.08206449 -0.08206449
#> Actinomyces_odontolyticus -2.173505582  3.96560227  5.71608938  1.57677340
#>                                G48792      G48793       G48794      G48795
#> Granulicella_unclassified -0.01103146 -0.01038247 -0.008760007 -0.02604906
#> Actinomyces_graevenitzii  -0.80924151 -0.79696704 -0.766280879  4.07521830
#> Actinomyces_johnsonii      6.92682253 -0.28911728 -0.308874220 -0.25292035
#> Actinomyces_massiliensis  -0.20308495 -0.19723901 -0.182624143 -0.77516997
#> Actinomyces_naeslundii    -0.07794349 -0.08618550 -0.106790529 -0.15716026
#> Actinomyces_odontolyticus -2.72372647 -2.65899460 -2.497164926  1.93283843
#>                                G48796       G48798      G48799     G48802
#> Granulicella_unclassified  0.01601114 -0.008111022  0.02282548 -0.0090845
#> Actinomyces_graevenitzii  -0.75239419  3.214224126  2.84358630 -0.7724181
#> Actinomyces_johnsonii     -0.29523910 -0.316776995 -0.37821824 -0.3049228
#> Actinomyces_massiliensis  -0.54022634 -0.176778198 -0.47884392 -0.1855471
#> Actinomyces_naeslundii    -0.22570842 -0.115032540 -0.31224955 -0.1026695
#> Actinomyces_odontolyticus -2.46349843  3.057708089  3.01664333  3.4824312
#>                                G48804       G48805      G48806      G48807
#> Granulicella_unclassified  0.01893157  0.006326602 -0.05321854 -0.01297841
#> Actinomyces_graevenitzii  -1.69715909 -2.048447396  8.79860850  2.71750390
#> Actinomyces_johnsonii     -0.33080159 -0.640808932 -0.39204943 -0.25750619
#> Actinomyces_massiliensis  -0.51391959 -1.319497754  0.12968315 -0.22062279
#> Actinomyces_naeslundii    -0.26279747 -0.523991309 -0.13449000 -0.05321745
#> Actinomyces_odontolyticus  1.52083872 -3.828093337 -2.30068045  3.32493041
#>                                G48808      G48810      G48812      G48814
#> Granulicella_unclassified -0.01005798 -0.04110260 -0.02650044  0.02752833
#> Actinomyces_graevenitzii   4.35102026  3.96964561  1.92658333  1.55278289
#> Actinomyces_johnsonii     -0.29306867 -0.22431021 -0.40212264 -0.58370856
#> Actinomyces_massiliensis  -0.19431603 -0.34191498 -0.21038121  1.49469774
#> Actinomyces_naeslundii    -0.09030650 -0.09268866 -0.27813393 -0.59757813
#> Actinomyces_odontolyticus  1.96683686  3.67370141  6.32296792  3.39267094
#>                                G48815       G48817      G48819       G48820
#> Granulicella_unclassified -0.05419201  0.002921726  0.04577656  0.007624573
#> Actinomyces_graevenitzii  -1.38152703 -0.545340509  2.34867916 -2.023898466
#> Actinomyces_johnsonii     -0.38019527 -0.451124162 -0.18772112 -0.656614481
#> Actinomyces_massiliensis   0.12091423 -0.077397130  1.96581859 -1.307805863
#> Actinomyces_naeslundii    -0.12212698 -0.255146744 -0.28134673 -0.540475333
#> Actinomyces_odontolyticus -2.39777825 -1.331991287  2.74634317  0.348781581
#>                                 G48821        G48822      G48823        G48826
#> Granulicella_unclassified -0.008435515 -0.0003231996  1.76377322 -0.0009721848
#> Actinomyces_graevenitzii   3.291973864  6.2483293085 -1.46744829 -0.6189872986
#> Actinomyces_johnsonii      2.229733223 -0.4116102892 -0.32487585 -0.4037075146
#> Actinomyces_massiliensis   0.885925692  0.8933731446  0.07999262 -0.1124728006
#> Actinomyces_naeslundii    -0.110911535 -0.2139366844 -0.06443290 -0.2056946724
#> Actinomyces_odontolyticus  5.007733793  5.4681324015 -1.11390082  4.8903973023
#>                                 G48827      G48831      G48832     G48836
#> Granulicella_unclassified -0.063277806 -0.01103146 -0.05321854  0.0258728
#> Actinomyces_graevenitzii  -1.553369545 -0.80924151 -1.36311534  4.0891743
#> Actinomyces_johnsonii     -0.269556429 -0.28121451  3.35116173 -0.2606270
#> Actinomyces_massiliensis   0.039071002 -0.20308495  0.12968315  1.3085295
#> Actinomyces_naeslundii    -0.006738816 -0.07794349 -0.13449000  0.7757561
#> Actinomyces_odontolyticus  3.712401655 -2.72372647 -1.30068045 -1.8098504
#>                                 G48838     G48839      G48840      G48841
#> Granulicella_unclassified  0.006925342  0.9223679 -0.05707939 -0.03683396
#> Actinomyces_graevenitzii  -1.924236695 -0.8213358  3.92099183  3.02140896
#> Actinomyces_johnsonii     -0.184600260 -0.7129605  6.23945313 -0.73978969
#> Actinomyces_massiliensis  -0.622069575 -0.2134103  3.86113030 -0.99549857
#> Actinomyces_naeslundii    -0.110320252 -0.3724223  3.41760737 -0.56817481
#> Actinomyces_odontolyticus  3.955778233  1.6931378  8.39745649 -3.50214512
#>                                G48842      G72670      G72671      G72672
#> Granulicella_unclassified -0.01232943 -0.06262882 -0.01297841 -0.04272507
#> Actinomyces_graevenitzii   4.05319320 -1.54109508 -0.84606490 -0.42135433
#> Actinomyces_johnsonii      6.20258541 -0.27745920 -0.25750619 -0.20455328
#> Actinomyces_massiliensis   4.93831575  0.04491695 -0.22062279 -0.35652984
#> Actinomyces_naeslundii     2.52663443 -0.01498083 -0.05321745 -0.07208363
#> Actinomyces_odontolyticus  3.61065451  2.69300925 -2.91792207  1.72841300
#>                                G72673      G72674      G72675      G72676
#> Granulicella_unclassified  0.01276621 -0.06165534 -0.01759795 -0.02475108
#> Actinomyces_graevenitzii   1.29841674  5.74354939 -0.97170697 -2.35393450
#> Actinomyces_johnsonii     -0.25572523  0.71068664 -0.81945188 -0.26872590
#> Actinomyces_massiliensis  -0.56945607  0.05368586 -0.38541491 -0.76347808
#> Actinomyces_naeslundii    -0.18449836  3.19559594 -0.54253111 -0.17364429
#> Actinomyces_odontolyticus -2.78715778  2.40625530 -1.30966634  2.26918195
#>                                G72677      G72678      G72679       G72680
#> Granulicella_unclassified -0.04207608 -0.05905940 -0.02960418 -0.008726947
#> Actinomyces_graevenitzii  -2.60699241 -1.47358552  3.40721271  4.650455876
#> Actinomyces_johnsonii     -0.21245605 -0.32092446 -0.67325056 -0.612198792
#> Actinomyces_massiliensis  -0.35068390  0.07706965 -0.49356489  2.287404745
#> Actinomyces_naeslundii    -0.08032564 -0.06031189 -0.39005389 -0.459519707
#> Actinomyces_odontolyticus -3.62638320 -2.88326727 -2.50720591 -2.861930630
#>                                G72681      G72682      G72683      G72684
#> Granulicella_unclassified -0.06295331 -0.04867564  0.02120302 -0.04040337
#> Actinomyces_graevenitzii   2.19187324 -1.27719408  1.12921192 -2.68824248
#> Actinomyces_johnsonii     -0.27350782 -0.44736886 -0.35846130 -0.69632443
#> Actinomyces_massiliensis   0.04199397  0.17060477 -0.49345878 -1.02765127
#> Actinomyces_naeslundii    -0.01085982 -0.19218409  4.29854609 -0.52284374
#> Actinomyces_odontolyticus  2.76693068  3.91372245 -1.94564348  2.19493040
#>                                G72685      G72686       G72687        G72688
#> Granulicella_unclassified -0.07763212 -0.02506128 -0.009375932  1.293047e-06
#> Actinomyces_graevenitzii   5.47669458 -1.11286332 -2.270819137 -6.005756e-01
#> Actinomyces_johnsonii     -0.71296051 -0.72856998 -0.604296018 -4.155617e-01
#> Actinomyces_massiliensis  -0.21341030 -0.45264328 -0.892088710 -1.037039e-01
#> Actinomyces_naeslundii    -0.37242229 -0.44774797 -0.451277695 -2.180577e-01
#> Actinomyces_odontolyticus  3.56077819 -2.05408283 -2.926662499 -1.623285e+00
#>                                G72689      G72690     G72691       G72692
#> Granulicella_unclassified -0.07795662  0.02314997 0.02055403 -0.007104484
#> Actinomyces_graevenitzii  -1.86926582 -1.61737507 4.61147184 -2.227858510
#> Actinomyces_johnsonii     -0.70900912 -0.38216962 8.82237644 -0.631955729
#> Actinomyces_massiliensis  -0.21633327 -0.47592095 5.13125652 -0.871627902
#> Actinomyces_naeslundii    -0.36830129 -0.31637055 1.63645976 -0.480124737
#> Actinomyces_odontolyticus  2.98644355 -1.75144788 9.83577984 -2.700100958
#>                                G72693      G72694       G72695       G72696
#> Granulicella_unclassified -0.02747391 -0.02650044  0.001299263  0.009521283
#> Actinomyces_graevenitzii   4.96992615 -2.31240525 -0.576026671 -1.875138835
#> Actinomyces_johnsonii     -0.39026848 -0.40212264 -0.431367226 -0.216211358
#> Actinomyces_massiliensis  -0.21915013 -0.21038121 -0.092011992 -0.598685794
#> Actinomyces_naeslundii    -0.26577091 -0.27813393 -0.234541714 -0.143288300
#> Actinomyces_odontolyticus -2.16991615 -2.07281835 -1.493820959  1.322187876
#>                                 G72697      G72698      G72699       G72701
#> Granulicella_unclassified  0.009845776 -0.06230433 -0.05094709  0.003570712
#> Actinomyces_graevenitzii  -1.869001603 -1.53495785 -1.32015471 -0.533066044
#> Actinomyces_johnsonii     -0.220162746 -0.28141059 -0.41970914 -0.459026936
#> Actinomyces_massiliensis  -0.595762822  0.04783992  0.15014396 -0.071551184
#> Actinomyces_naeslundii    -0.147409306 -0.01910183 -0.16333704 -0.263388756
#> Actinomyces_odontolyticus  2.214319598 -3.20692661 -2.07411891  7.835787130
#>                                G72702       G72703      G72704      G72705
#> Granulicella_unclassified -0.02700824  0.009845776 -0.02742367  0.02444794
#> Actinomyces_graevenitzii  -1.14968672 -1.869001603 -2.44275318 -1.59282614
#> Actinomyces_johnsonii     -0.70486165 -0.220162746 -0.85437992 -0.39797517
#> Actinomyces_massiliensis  -0.47018111 -0.595762822 -0.91073236 -0.46422906
#> Actinomyces_naeslundii    -0.42302194 -0.147409306 -0.68768398 -0.33285458
#> Actinomyces_odontolyticus -2.24827843 -3.078451187 -2.56353302 -1.62198414
#>                               G72706      G72708      G72709      G72710
#> Granulicella_unclassified  0.0130907 -0.07698314  0.01893157  0.02282548
#> Actinomyces_graevenitzii  -1.8076293  4.61528535  2.40557677 -1.62351230
#> Actinomyces_johnsonii     -0.2596766 -0.72086328 -0.33080159 -0.37821824
#> Actinomyces_massiliensis  -0.5665331 -0.20756435 -0.51391959 -0.47884392
#> Actinomyces_naeslundii    -0.1886194 -0.38066431 -0.26279747 -0.31224955
#> Actinomyces_odontolyticus -2.7547918 -2.60201484 -2.17220502  2.61815782
#>                                G72711     G72712      G72713       G72714
#> Granulicella_unclassified  0.01114375  0.0218520 -0.05094709 -0.001945663
#> Actinomyces_graevenitzii  -0.20402204 -1.6419240 -1.32015471  5.192483408
#> Actinomyces_johnsonii     -0.23596829 -0.3663641 -0.41970914 -0.391853353
#> Actinomyces_massiliensis   5.67000638 -0.4876128  0.15014396 -0.121241718
#> Actinomyces_naeslundii    -0.16389333 -0.2998865 -0.16333704  4.319849522
#> Actinomyces_odontolyticus  5.05774842 -1.8809116 -2.07411891 -1.817480303
#>                                G72715      G72716      G72717      G72718
#> Granulicella_unclassified  0.00584216 -0.01038247 -0.01038247 -0.01330290
#> Actinomyces_graevenitzii  -0.49010542 -0.79696704 -0.79696704 -0.85220213
#> Actinomyces_johnsonii     -0.48668665 -0.28911728 -0.28911728 -0.25355480
#> Actinomyces_massiliensis  -0.05109038 -0.19723901 -0.19723901 -0.22354576
#> Actinomyces_naeslundii    -0.29223580 -0.08618550 -0.08618550 -0.04909644
#> Actinomyces_odontolyticus -1.04069788 -2.65899460 -0.26939485 -2.95028801
#>                                G72719      G72720      G72721      G72722
#> Granulicella_unclassified -0.01265392 -0.01200493 -0.01330290 -0.01330290
#> Actinomyces_graevenitzii  -0.83992767  5.15631919 -0.85220213 -0.85220213
#> Actinomyces_johnsonii     -0.26145757 -0.26936035 -0.25355480 -0.25355480
#> Actinomyces_massiliensis  -0.21769981 -0.21185387 -0.22354576 -0.22354576
#> Actinomyces_naeslundii    -0.05733846 -0.06558047 -0.04909644 -0.04909644
#> Actinomyces_odontolyticus  0.58258159  4.26804118 -2.95028801  2.32194955
#>                                G72723      G72724       G72725      G72726
#> Granulicella_unclassified -0.01232943 -0.01232943 -0.006488559 -0.01103146
#> Actinomyces_graevenitzii  -0.83379044 -0.83379044 -0.723320251 -0.80924151
#> Actinomyces_johnsonii     -0.26540896 -0.26540896 -0.336533931 -0.28121451
#> Actinomyces_massiliensis  -0.21477684 -0.21477684 -0.162163335 -0.20308495
#> Actinomyces_naeslundii    -0.06145946 -0.06145946 -0.135637570 -0.07794349
#> Actinomyces_odontolyticus -2.85319020  2.05198151 -2.270603385 -2.72372647
#>                                G72727      G72728      G72729      G72730
#> Granulicella_unclassified -0.01297841 -0.01330290 -0.01297841 -0.01232943
#> Actinomyces_graevenitzii  -0.84606490 -0.85220213 -0.84606490 -0.83379044
#> Actinomyces_johnsonii     -0.25750619 -0.25355480 -0.25750619 -0.26540896
#> Actinomyces_massiliensis  -0.22062279 -0.22354576 -0.22062279 -0.21477684
#> Actinomyces_naeslundii    -0.05321745 -0.04909644 -0.05321745 -0.06145946
#> Actinomyces_odontolyticus -2.91792207  1.23056949 -2.91792207  3.54596978
#>                                G72731       G72732      G72733      G72734
#> Granulicella_unclassified -0.01265392 -0.006813052 -0.01297841 -0.01297841
#> Actinomyces_graevenitzii  -0.83992767 -0.729457484 -0.84606490 -0.84606490
#> Actinomyces_johnsonii     -0.26145757 -0.332582544 -0.25750619 -0.25750619
#> Actinomyces_massiliensis  -0.21769981 -0.165086307 -0.22062279 -0.22062279
#> Actinomyces_naeslundii    -0.05733846 -0.131516564 -0.05321745 -0.05321745
#> Actinomyces_odontolyticus  2.04674580 -2.302969319 -2.91792207 -2.91792207
#>                                G72735        G72736      G72737
#> Granulicella_unclassified -0.01297841  0.0009747709 -0.01362740
#> Actinomyces_graevenitzii  -0.84606490 -0.5821639036 -0.85833937
#> Actinomyces_johnsonii     -0.25750619 -0.4274158382 -0.24960341
#> Actinomyces_massiliensis  -0.22062279 -0.0949349651 -0.22646873
#> Actinomyces_naeslundii    -0.05321745 -0.2304207083 -0.04497544
#> Actinomyces_odontolyticus  0.19754404 -1.5261868938  0.13542465
fit <- IntegratedLearner(
  MAE_train = mae,
  experiment = names(MultiAssayExperiment::experiments(mae))[1:2],
  assay.type = rep("abundance", 2),
  folds = 2, base_learner = "SL.mean",
  run_stacked = FALSE, run_concat = FALSE,
  family = stats::binomial()
)
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
plt <- plot(fit)
names(plt)
#> [1] "plot"      "ROC_table"
```
