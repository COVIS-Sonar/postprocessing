from pycovis.postprocess import runtime
import json

with runtime.Runtime() as covis:

    print("Looking for covis_path.json at: %s" % covis.find_input_file("covis_bathy.json") )

    with open( covis.find_input_file("covis_bathy.json") ) as f:
        j = json.load( f )

    print(j)
