# TCCON .nc checker

## About
A small stand-alone Python script to give an GUI interface to check data points and plot associated spectra  

## Requirements

Numpy, PyQT5, Matplotlib, xarray


## Usage
Run like this:
	python tccon_nc_checker.py path/to/private.nc-file path/to/spectra_rootfolder

In the GUI you can choose the variable of the .nc-file in the top left drop down menu.

Clicking on a point displays addtional information in the textbox and shows the underlying spectrum

Checking the boy day-by-day will go through each day, where there are spectra, consecutively.

You can mark a whole day for later review, saved as a date in a textfile.

'Delete Spectrum' will save the spectrum name to a textfile to be used for flagging later.
