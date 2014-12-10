"""Windows-specific part of the installation"""

import os, sys, shutil


def install():
    """Routine to be run by the win32 installer with the -install switch."""

    # Get some system constants
    target = os.path.join(sys.prefix, 'Scripts', 'hexrd_gui.exe')
    # Lookup path to common startmenu ...
    start_dir = os.path.join(
        get_special_folder_path('CSIDL_COMMON_PROGRAMS'),
        'HEXRD'
    )

    # Create entry ...
    if not os.path.isdir(start_dir):
        os.mkdir(start_dir)
        directory_created(start_dir)

    # Create program shortcuts ...
    link_file = os.path.join(start_dir, 'hexrd.lnk')
    create_shortcut(target, 'hexrd', link_file, "gui", "%HOMEDRIVE%%HOMEPATH%")
    file_created(link_file)


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
