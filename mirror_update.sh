#!/bin/bash
# mirror on a separate SVN is needed for some projects
# (SVN changes the inodes)
#
# Jean-Philippe.Lambert@ircam.fr 2008


mirror_dir=rta_mirror # a symbolic link is ok

echo "mirror directory is ${mirror_dir}"

for source in Doxyfile *.{c,h} matlab/*.{c,h,m} doc/README.txt documentation/*.{tex,pdf,graffle} ; do
#    echo "${source}"
#    diff "${source}"  "${mirror_dir}/${source}"
    cat "${source}" > "${mirror_dir}/${source}"
done



# cd ${mirror_dir} && svn commit