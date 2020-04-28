# using DelimitedFiles
# const UCR_DATASETS = """
# 50words,
# Adiac,
# ArrowHead,
# Beef,
# BeetleFly,
# BirdChicken,
# CBF,
# Car,
# ChlorineConcentration,
# CinC_ECG_torso,
# Coffee,
# Computers,
# Cricket_X,
# Cricket_Y,
# Cricket_Z,
# DiatomSizeReduction,
# DistalPhalanxOutlineAgeGroup,
# DistalPhalanxOutlineCorrect,
# DistalPhalanxTW,
# ECG200,
# ECG5000,
# ECGFiveDays,
# Earthquakes,
# ElectricDevices,
# FISH,
# FaceAll,
# FaceFour,
# FacesUCR,
# FordA,
# FordB,
# Gun_Point,
# Ham,
# HandOutlines,
# Haptics,
# Herring,
# InlineSkate,
# InsectWingbeatSound,
# ItalyPowerDemand,
# LargeKitchenAppliances,
# Lighting2,
# Lighting7,
# MALLAT,
# Meat,
# MedicalImages,
# MiddlePhalanxOutlineAgeGroup,
# MiddlePhalanxOutlineCorrect,
# MiddlePhalanxTW,
# MoteStrain,
# NonInvasiveFatalECG_Thorax1,
# NonInvasiveFatalECG_Thorax2,
# OSULeaf,
# OliveOil,
# PhalangesOutlinesCorrect,
# Phoneme,
# Plane,
# ProximalPhalanxOutlineAgeGroup,
# ProximalPhalanxOutlineCorrect,
# ProximalPhalanxTW,
# RefrigerationDevices,
# ScreenType,
# ShapeletSim,
# ShapesAll,
# SmallKitchenAppliances,
# SonyAIBORobotSurface,
# SonyAIBORobotSurfaceII,
# StarLightCurves,
# Strawberry,
# SwedishLeaf,
# Symbols,
# ToeSegmentation1,
# ToeSegmentation2,
# Trace,
# TwoLeadECG,
# Two_Patterns,
# UWaveGestureLibraryAll,
# Wine,
# WordsSynonyms,
# Worms,
# WormsTwoClass,
# synthetic_control,
# uWaveGestureLibrary_X,
# uWaveGestureLibrary_Y,
# uWaveGestureLibrary_Z,
# wafer,
# yoga
# """
#
# """
#     DynamicAxisWarping.download_data([path="/DynamicAxisWarping/data/"])
#
# Downloads the UC Riverside Time Series Classification Archive to
# the specified path.
# """
# function download_ucr()
#
#     UCR_URL = "http://www.cs.ucr.edu/~eamonn/time_series_data/UCR_TS_Archive_2015.zip"
#
#     # the password is "attempttoclassify"
#     @info(
#         """
#         A simple password is required to extract the UCR Time Series Classification. While the
#         data is downloading, please read the following paper
#
#         Bing Hu, Yanping Chen, Eamonn J. Keogh: Time Series Classification under More Realistic
#         Assumptions. SDM 2013: 578-586.
#         http://dx.doi.org/10.1137/1.9781611972832.64
#
#         Find the following sentence (it's near the end of the first page):
#
#         “Every item that we ******* ## @@@@@@@ belongs to exactly one of our well defined classes”
#
#         Enter the three redacted words (without spaces) to extract the data. The purpose of this
#         exercise is to encourage you to read the associated paper, see:
#
#         http://www.cs.ucr.edu/~eamonn/time_series_data/
#
#         """
#     )
#
#     download(UCR_URL, "$DATAPATH/UCR.zip")
#
#     # should unzip the file correctly in platform-specific manner
#     run(unpack_cmd("UCR", DATAPATH , ".zip", ""))
#
#     # delete the zip file when we're done
#     rm("$DATAPATH/UCR.zip")
#
#     @info("Download and extraction successful!")
#
# end
#
# """
#     data,labels = DynamicAxisWarping.traindata(name)
#
# Loads the training set of the specified dataset. Returns a matrix `data` where each column
# holds a 1-dimensional time series. The class label for each column is held `labels`
# which is a vector of length `size(data,2)`.
#
# Available datasets:
#
# $UCR_DATASETS
# """
# function ucr_traindata(name::AbstractString)
#     try
#         Y = readdlm(DATAPATH*"/UCR_TS_Archive_2015/"*name*"/"*name*"_TRAIN", ',')
#         labels = round(Int,vec(Y[:,1]))
#         data = transpose(Y[:,2:end])
#         return data,labels
#     catch err
#         showerror(stdout, err, backtrace());println()
#         @info("You may have recieved this error because you haven't downloaded the database yet.")
#         @info("Try running download_ucr() first.")
#         @info("You may also have mispelled the name of the dataset.")
#     end
# end
#
# """
#     data,labels = DynamicAxisWarping.traindata(name)
#
# Loads the test set of the specified dataset. Returns a matrix `data` where each column
# holds a 1-dimensional time series. The class label for each column is held `labels`
# which is a vector of length `size(data,2)`.
#
# Available datasets:
#
# $UCR_DATASETS
# """
# function ucr_testdata(name::AbstractString)
#     try
#         Y = return readcsv(DATAPATH*"/UCR_TS_Archive_2015/"*name*"/"*name*"_TEST")
#         labels = round(Int,vec(Y[:,1]))
#         data = transpose(Y[:,2:end])
#         return data,labels
#     catch err
#         showerror(stdout, err, backtrace());println()
#         @info("You may have recieved this error because you haven't downloaded the database yet.")
#         @info("Try running download_ucr() first.")
#         @info("You may also have mispelled the name of the dataset.")
#     end
# end
