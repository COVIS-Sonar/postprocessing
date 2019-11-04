from pycovis.postprocess import runtime
import json

with runtime.Runtime() as pp:

    print("Looking for covis_path.json at: %s" % pp.find_input_file("covis_bathy.json") )

    with open( pp.find_input_file("covis_bathy.json") ) as f:
        j = json.load( f )

    print(j)
