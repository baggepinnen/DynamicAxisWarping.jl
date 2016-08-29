import BinDeps: unpack_cmd

function fetch_datasets()

    UCR_URL = "http://www.cs.ucr.edu/~eamonn/time_series_data/UCR_TS_Archive_2015.zip"
    DATAPATH = Pkg.Dir.path()*"/TimeWarp/data"

    # the password is "attempttoclassify"
    info(
        """
        A simple password is required to extract the UCR Time Series Classification. While the
        data is downloading, please read the following paper

        Bing Hu, Yanping Chen, Eamonn J. Keogh: Time Series Classification under More Realistic
        Assumptions. SDM 2013: 578-586.
        http://dx.doi.org/10.1137/1.9781611972832.64

        Find the following sentence (it's near the end of the first page):
        
        “Every item that we ******* ## @@@@@@@ belongs to exactly one of our well defined classes”

        Enter the three redacted words (without spaces) to extract the data. The purpose of this
        exercise is to encourage you to read the associated paper, see:

        http://www.cs.ucr.edu/~eamonn/time_series_data/

        """
    )

    download(UCR_URL, "$DATAPATH/UCR.zip")

    # should unzip the file correctly in platform-specific manner
    run(unpack_cmd("UCR", DATAPATH , ".zip", ""))

    # delete the zip file when we're done
    rm("$DATAPATH/UCR.zip")

    info("Download and extraction successful!")

end
