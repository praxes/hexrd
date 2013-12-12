"""Windows-specific part of the installation"""

import os, sys, shutil

def mkshortcut(target, description, link_file, *args, **kw):
    """make a shortcut if it doesn't exist, and register its creation"""

    create_shortcut(target, description, link_file, *args, **kw)
    file_created(link_file)

def install():
    """Routine to be run by the win32 installer with the -install switch."""

    # Get some system constants
    prefix = sys.prefix
    python = os.path.join(prefix, 'python.exe')
    # Lookup path to common startmenu ...
    start_dir = os.path.join(
        get_special_folder_path('CSIDL_COMMON_PROGRAMS'),
        'HEXRD'
    )
    scripts_dir = os.path.join(prefix, 'Scripts')
    lib_dir = os.path.join(prefix, 'Lib', 'site-packages', 'hexrd')

    # Create entry ...
    if not os.path.isdir(start_dir):
        os.mkdir(start_dir)
        directory_created(start_dir)

    # Create program shortcuts ...
    script = '"%s"' % os.path.join(scripts_dir, 'hexrd_app.py')
    f = os.path.join(start_dir, 'hexrd.lnk')
    mkshortcut(python, 'hexrd', f, script, "%HOMEDRIVE%%HOMEPATH%")

def remove():
    """Routine to be run by the win32 installer with the -remove switch."""
    pass

# main()
if len(sys.argv) > 1:
    if sys.argv[1] == '-install':
        install()
    elif sys.argv[1] == '-remove':
        remove()
    else:
        print "Script was called with option %s" % sys.argv[1]
