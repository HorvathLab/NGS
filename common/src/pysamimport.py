import pkg_resources
# can't find it in the frozen library.zip...
# pkg_resources.require("pysam>=0.8.1")
import pysam
NEEDED_VERSION = "0.8.1"
assert pkg_resources.parse_version(pysam.version.__version__) >= pkg_resources.parse_version(NEEDED_VERSION), "PySam version at least %s required"%(NEEDED_VERSION,)
