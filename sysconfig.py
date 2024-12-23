import os, sys
from distutils import sysconfig
def cygpath(wpath):
    s = os.popen('cygpath -u "%s"' % wpath)
    path = s.read().strip()
    s.close()
    return path
incdir = sysconfig.get_python_inc()
keys = (
    'BLDLIBRARY',
    'LDFLAGS',
    'LDLAST',
    'LDLIBRARY',
    'LIBDIR',
    'LIBP',
    'LIBPL',
    'LIBS',
    'LINKFORSHARED',
    'MODLIBS',
    'SYSLIBS',
    'LA_LDFLAGS',
)
if os.name == "nt":
    # We are running under Python for Windows (the real one...
    # not Cygwin Python, under which 'os.name' is 'posix').
    # We assume that we are still in the Cygwin POSIX environment,
    # however (this is 'configure', after all); so we convert
    # all Windows pathnames to POSIX pathnames using 'cygpath'.
    incdir = cygpath(incdir)
    vars = {}
    libs = os.path.join(sys.prefix, "libs")
    libs = cygpath(libs)
    version = sysconfig.get_python_version()
    version = version.replace('.', '')
    vars['BLDLIBRARY'] = "-L%s -lpython%s" % (libs, version)
else:
    vars = sysconfig.get_config_vars()
    # transform AIX's python.exp
    vars['LINKFORSHARED'] = vars['LINKFORSHARED'].replace('Modules',vars['LIBPL'])
    if vars['LDLIBRARY'] == vars['LIBRARY']:
        # "On systems without shared libraries, LDLIBRARY is the same as LIBRARY"
        vars['BLDLIBRARY'] = "-L%(LIBPL)s -lpython%(VERSION)s" % vars
    elif vars['BLDLIBRARY']:
        #     The 'mpicc' wrapper for LAM/MPI isn't very smart about "-L"
        # options.  Adding "-L/usr/lib" can cause "-lmpi" to be found in /usr/lib
        # instead of LAM's 'lib' directory.  Of course, avoiding "-L/usr/lib"
        # doesn't really fix the problem, but it does make it much less likely;
        # and "-L/usr/lib" is redundant and potentially problematic anyway.
        #     Python 2.4 and later puts a symlink to libpython.so in LIBPL
        # (/usr/lib/python2.x/config), which makes adding "-L"
        # (in addition to "-L") completely redundant.
        #     But we still support Python 2.3, and we prefer shared to static,
        # so we still add "-L" when Python is installed in a non-standard
        # location.  Note that the linker will still prefer shared over static
        # with only "-L/usr/lib/python2.3/config" on the link line.
        libdir = ""
        if vars['LIBDIR'] != "/usr/lib":
            libdir = "-L%(LIBDIR)s "
        # Important: on Cygwin, the import library for libpython.dll is
        # nested inside Python's 'config' directory (see Issue39).  This means
        # that the linker always needs help finding "-lpython2.x" (in the form
        # of "-L"), even for the "system" Python installed under /usr.
        vars['BLDLIBRARY'] = (libdir + "-L%(LIBPL)s -lpython%(VERSION)s") % vars
    else:
        # "On Mac OS X frameworks, BLDLIBRARY is blank"
        # See also Issue39.
        framework = "%(PYTHONFRAMEWORKDIR)s/Versions/%(VERSION)s/%(PYTHONFRAMEWORK)s" % vars
        PYTHONFRAMEWORK = vars.get('PYTHONFRAMEWORK', 'Python')
        vars['LINKFORSHARED'] = vars['LINKFORSHARED'].replace(framework, "-framework " + PYTHONFRAMEWORK)
        vars['LA_LDFLAGS'] = "-Wl,-framework,%s" % PYTHONFRAMEWORK
vars['LDFLAGS'] = '' # only causes trouble (e.g., "-arch i386 -arch ppc" on Mac) -- see issue97
print 'PYTHON_INCDIR="%s"' % incdir
for key in keys:
    print 'PYTHON_%s="%s"' % (key, vars.get(key, ''))

