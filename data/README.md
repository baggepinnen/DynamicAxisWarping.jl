## Time Warping Datasets

To download the datasets locally, run:

```julia
julia> DynamicAxisWarp.Datasets.download_ucr()
```

You will have to enter a password, follow the printed instructions. The datasets will be saved in `.../DynamicAxisWarp/data/UCR_TS_Archive_2015/`. Then you can load the train/test data as follows (using the `50words` dataset as an example):

```julia
traindata = DynamicAxisWarp.traindata("50words");
testdata = DynamicAxisWarp.testdata("50words")
```

#### References

Yanping Chen, Eamonn Keogh, Bing Hu, Nurjahan Begum, Anthony Bagnall, Abdullah Mueen and Gustavo Batista (2015). The UCR Time Series Classification Archive. URL www.cs.ucr.edu/~eamonn/time_series_data/

#### List of all datasets:

* 50words
* Adiac
* ArrowHead
* Beef
* BeetleFly
* BirdChicken
* CBF
* Car
* ChlorineConcentration
* CinC_ECG_torso
* Coffee
* Computers
* Cricket_X
* Cricket_Y
* Cricket_Z
* DiatomSizeReduction
* DistalPhalanxOutlineAgeGroup
* DistalPhalanxOutlineCorrect
* DistalPhalanxTW
* ECG200
* ECG5000
* ECGFiveDays
* Earthquakes
* ElectricDevices
* FISH
* FaceAll
* FaceFour
* FacesUCR
* FordA
* FordB
* Gun_Point
* Ham
* HandOutlines
* Haptics
* Herring
* InlineSkate
* InsectWingbeatSound
* ItalyPowerDemand
* LargeKitchenAppliances
* Lighting2
* Lighting7
* MALLAT
* Meat
* MedicalImages
* MiddlePhalanxOutlineAgeGroup
* MiddlePhalanxOutlineCorrect
* MiddlePhalanxTW
* MoteStrain
* NonInvasiveFatalECG_Thorax1
* NonInvasiveFatalECG_Thorax2
* OSULeaf
* OliveOil
* PhalangesOutlinesCorrect
* Phoneme
* Plane
* ProximalPhalanxOutlineAgeGroup
* ProximalPhalanxOutlineCorrect
* ProximalPhalanxTW
* RefrigerationDevices
* ScreenType
* ShapeletSim
* ShapesAll
* SmallKitchenAppliances
* SonyAIBORobotSurface
* SonyAIBORobotSurfaceII
* StarLightCurves
* Strawberry
* SwedishLeaf
* Symbols
* ToeSegmentation1
* ToeSegmentation2
* Trace
* TwoLeadECG
* Two_Patterns
* UWaveGestureLibraryAll
* Wine
* WordsSynonyms
* Worms
* WormsTwoClass
* synthetic_control
* uWaveGestureLibrary_X
* uWaveGestureLibrary_Y
* uWaveGestureLibrary_Z
* wafer
* yoga
