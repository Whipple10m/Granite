# Granite
Granite data acquisition system for Whipple developped by Glenn Sembrowski at Purdue.

Granite is largely written in FORTAN 90 and running on a VAX. Portions of the code were written as `KUMAC` macros, which allowed certain aspects of the data taking to be modified easily by changing the macros. For example, the GUI was implemented in the FORTAN code, but certain button presses on the GUI resulted in calls to the kumac macros, which would then direct the setting up of the runs etc. The DAQ electronics was based on CAMAC modules, directed by HYTEC list processors (LPs) and read out via HYTEC ethernet crate-contollers (ECCs). Granite builds the LP program for each CAMAC crate on the fly when it starts, uploads it into the LPs and starts them running. It polls the LPs for data, downloads it and packs in the data files in the GRANITE data format (GDF, [see the GDF repository for details](https://github.com/Whipple10m/GDF)). The LP programs built by Granite contain special tags to allow the datastream to be resynchronised if necessary; these were constructed by forcing the LPs to read their own instruction pointers into the data, which could be verified by Granite, and searched for if resynchronisation was necessary.

This is a version I downloaded and used to study the DAQ, and is incomplete, missing (at least) all the scripted kumacs.

## Screenshots of GUI 

![Window selection](https://github.com/Whipple10m/Granite/blob/main/Assets/granite1.png)

![Command entry and log](https://github.com/Whipple10m/Granite/blob/main/Assets/granite2.png)

![Run control](https://github.com/Whipple10m/Granite/blob/main/Assets/granite3.png)

![Display control](https://github.com/Whipple10m/Granite/blob/main/Assets/granite4.png)

![Histogram](https://github.com/Whipple10m/Granite/blob/main/Assets/granite5.png)
