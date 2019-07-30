"""
This file contains definition of paths to enable testing which points to default models
"""
import os.path as op

path_to_olga = op.dirname(op.realpath(__file__))
assert op.isdir( path_to_olga )

path_to_olga_default_models = op.join(path_to_olga , "default_models")
assert op.isdir( path_to_olga_default_models  )
