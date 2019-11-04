from pycovis.postprocess import runtime

with runtime.Runtime() as pp:

    print("Looking for covis_path.json at: %s" % pp.find_input_file("covis_bathy.json") )
