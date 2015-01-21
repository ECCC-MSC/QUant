(You may want to turn Word Wrap on: Format/Word Wrap, in order to better read this file).

You will find 3 MATLAB scripts in this zip folder, that will do LISST-100 data processing for you: tt2mat.m, getscat.m and invert.p

Copy the 3 files to your MATLAB folder, then type help getscat:
> %  Usage:		[scat,tau,zsc,data,cscat] = getscat('datafile','zscfile',X,'ringarefile');
> %
> % Set X=1 for LISST-100X data, all other values for LISST-100 format % 
> % Output:
> % scat is the raw scattering signature in digital counts (i.e. NOT 
> corrected % with the ringarea file!!) % tau is the optical 
> transmission % zsc is the loaded zscat % data is the raw data file, 
> converted from binary, or it can be a .LOG % file % cscat is the 
> corrected scattering, computed if the ringarea file is % supplied as 
> input
>
You will need to locate the ringarea file for your instrument. It is named ringarea_xxxx.asc, where xxxx is the serial number of your instrument.  It is supplied on your ship disk. If you do not ahve it, contact Sequoia at info@sequoiasci.com and ask to have a copy sent to you. You must supply your serial number for us to locate the file!


Run getscat in order to obtain cscat; this is needed by the invert.p script. getscat.m calls tt2mat.m, as this is the function that reads the binary .DAT file from the LISST.


The syntax for invert.p is: (see also this article on our website: http://sequoiasci.com/Articles/ArticlePage.aspx?pageId=128)

[vd dias]=invert(cscat,instrument_type,ST,RANDOM,SHARPEN,GREEN,WAITBARSHOW);
where:
cscat is the fully corrected scattering data in n x 32 format, obtained using getscat.
instrument_type is
	1 for type A, SPHERICAL MATRIX ONLY
	2 for type B
	3 for type C
	4 for FLOC, SPHERICAL MATRIX ONLY
ST = 1 if the data are to be inverted in LISST-ST format (8 size bins), otherwise 0 for -100/-100X format (32 size bins)

RANDOM = 1 if matrices based on scattering from randomly shaped particles are to be used

SHARPEN = 1 causes the routine to check if the size distribution is narrow and, if so, sharpens it. 

GREEN = 1 if inversion is for a green laser unit (only type B as of 8/12/2010)

WAITBARSHOW = 1 if user wants a waitbar to show during processing.

Outputs are:

vd - volume distribution (NOT CALIBRATED WITH VCC)

dias - the midpoint of the size bins for the 8 / 32 size classes for the appropriate instrument, inversion type and laser color


Finally, you need to convert the volume distribution to calibrated units. Locate the InstrumentData.txt file for your instrument and look up the VCC (see this FAQ for details: http://sequoiasci.com/faq/faqquestion.aspx?faqquestion=90). 

The InstrumentData.txt file is on your shipdisk. If you can't find your shipdisk, contact Sequoia via email: info@sequoiasci.com and ask for a copy of your InstrumentData.txt file. You must supply the serial number of your instrument for us to help you.

Then divide the vd with VCC: vd = vd/VCC;