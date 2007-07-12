# --------------------------------------------------------------------
# Specify your configuration below.
# See documentation for hints.
# --------------------------------------------------------------------

PYRTICLE_CONF_TEMPLATE_VERSION = 1

# --------------------------------------------------------------------
# Path options
# --------------------------------------------------------------------

BOOST_INCLUDE_DIRS = []
BOOST_LIBRARY_DIRS = [] 
BPL_LIBRARIES = ["boost_python-gcc41-mt"]

BOOST_BINDINGS_INCLUDE_DIRS = ["/home/andreas/work/boost-sandbox"]

# --------------------------------------------------------------------
# Compiler flags
# --------------------------------------------------------------------
#EXTRA_COMPILE_ARGS = ["-fmessage-length=0", "-Wno-sign-compare"]
EXTRA_COMPILE_ARGS = ["-Wno-sign-compare"]

