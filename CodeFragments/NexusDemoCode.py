#!/usr/bin/env python
"""
Writes a simple NeXus HDF5 file using h5py with links
according to the example from Figure 2.1 in the Design chapter

https://manual.nexusformat.org/examples/python/simple_example_write2/index.html

"""

from pathlib import Path
import h5py
import numpy

filename = str(Path(__file__).absolute().parent.parent / "simple_example.dat")
buffer = numpy.loadtxt(filename).T
tthData = buffer[0]  # float[]
countsData = numpy.asarray(buffer[1], "int32")  # int[]

with h5py.File("simple_example_write2.hdf5", "w") as f:  # create the HDF5 NeXus file
    f.attrs["default"] = "entry"

    nxentry = f.create_group("entry")
    nxentry.attrs["NX_class"] = "NXentry"
    nxentry.attrs["default"] = "data"

    nxinstrument = nxentry.create_group("instrument")
    nxinstrument.attrs["NX_class"] = "NXinstrument"

    nxdetector = nxinstrument.create_group("detector")
    nxdetector.attrs["NX_class"] = "NXdetector"

    # store the data in the NXdetector group
    ds_tth = nxdetector.create_dataset("two_theta", data=tthData)
    ds_tth.attrs["units"] = "degrees"
    ds_counts = nxdetector.create_dataset("counts", data=countsData)
    ds_counts.attrs["units"] = "counts"

    # create the NXdata group to define the default plot
    nxdata = nxentry.create_group("data")
    nxdata.attrs["NX_class"] = "NXdata"
    nxdata.attrs["signal"] = "counts"
    nxdata.attrs["axes"] = "two_theta"
    nxdata.attrs["two_theta_indices"] = [
        0,
    ]

    source_addr = "/entry/instrument/detector/two_theta"  # existing data
    target_addr = "two_theta"  # new location
    ds_tth.attrs["target"] = source_addr  # a NeXus API convention for links
    nxdata[target_addr] = f[source_addr]  # hard link
    # nxdata._id.link(source_addr, target_addr, h5py.h5g.LINK_HARD)

    source_addr = "/entry/instrument/detector/counts"  # existing data
    target_addr = "counts"  # new location
    ds_counts.attrs["target"] = source_addr  # a NeXus API convention for links
    nxdata[target_addr] = f[source_addr]  # hard link
    # nxdata._id.link(source_addr, target_addr, h5py.h5g.LINK_HARD)



#!/usr/bin/env python
"""Writes a NeXus HDF5 file using h5py and numpy"""

from pathlib import Path
import datetime
import h5py  # HDF5 support
import numpy

print("Write a NeXus HDF5 file")
filename = "simple_example_basic.nexus.hdf5"
timestamp = datetime.datetime.now().astimezone().isoformat()

# load data from two column format
data_filename = str(Path(__file__).absolute().parent.parent / "simple_example.dat")
data = numpy.loadtxt(data_filename).T
mr_arr = data[0]
i00_arr = numpy.asarray(data[1], "int32")

# create the HDF5 NeXus file
with h5py.File(filename, "w") as f:
    # point to the default data to be plotted
    f.attrs["default"] = "entry"
    # give the HDF5 root some more attributes
    f.attrs["file_name"] = filename
    f.attrs["file_time"] = timestamp
    f.attrs["instrument"] = "APS USAXS at 32ID-B"
    f.attrs["creator"] = "simple_example_basic_write.py"
    f.attrs["NeXus_version"] = "4.3.0"
    f.attrs["HDF5_Version"] = h5py.version.hdf5_version
    f.attrs["h5py_version"] = h5py.version.version

    # create the NXentry group
    nxentry = f.create_group("entry")
    nxentry.attrs["NX_class"] = "NXentry"
    nxentry.attrs["default"] = "mr_scan"
    nxentry.create_dataset("title", data="1-D scan of I00 v. mr")

    # create the NXentry group
    nxdata = nxentry.create_group("mr_scan")
    nxdata.attrs["NX_class"] = "NXdata"
    nxdata.attrs["signal"] = "I00"  # Y axis of default plot
    nxdata.attrs["axes"] = "mr"  # X axis of default plot
    nxdata.attrs["mr_indices"] = [
        0,
    ]  # use "mr" as the first dimension of I00

    # X axis data
    ds = nxdata.create_dataset("mr", data=mr_arr)
    ds.attrs["units"] = "degrees"
    ds.attrs["long_name"] = "USAXS mr (degrees)"  # suggested X axis plot label

    # Y axis data
    ds = nxdata.create_dataset("I00", data=i00_arr)
    ds.attrs["units"] = "counts"
    ds.attrs["long_name"] = "USAXS I00 (counts)"  # suggested Y axis plot label

print("wrote file:", filename)



# The script I provided earlier iterates only over the top-level entries in the NeXus file. 
# To handle nested entries up to three layers deep, you'll need to modify the script to recursively search through the groups in the file. Here's an updated version of the script that handles nested entries:

import h5py

def find_nxcanSAS_entries(group, path=''):
    nxcanSAS_entries = []
    
    for name, item in group.items():
        current_path = f"{path}/{name}" if path else name
        
        # Check if the item is a group
        if isinstance(item, h5py.Group):
            # Check if the group has the attribute "NXcanSAS"
            if 'NXcanSAS' in item.attrs:
                nxcanSAS_entries.append(current_path)
            
            # Recursively search within the group
            nxcanSAS_entries.extend(find_nxcanSAS_entries(item, current_path))
    
    return nxcanSAS_entries

def get_nxcanSAS_entries(file_path):
    with h5py.File(file_path, 'r') as nexus_file:
        return find_nxcanSAS_entries(nexus_file)

# Example usage
file_path = 'your_nexus_file.nxs'  # Replace with your NeXus file path
entries_with_nxcanSAS = get_nxcanSAS_entries(file_path)

print("Entries with NXcanSAS attribute:")
for entry in entries_with_nxcanSAS:
    print(entry)
