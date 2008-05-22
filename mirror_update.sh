#!/bin/bash
# mirror on a separate SVN is needed for some projects
# (SVN changes the inodes)
#
# Jean-Philippe.Lambert@ircam.fr 2008


mirror_dir=rta_mirror # a symbolic link is ok


for source in $( ls *.{c,h} matlab/*.{c,h,m} doc/README.txt documentation/*.{tex,pdf,graffle} ); do
#    diff ${source}  ${mirror_dir}/${source}
    cat ${source} > ${mirror_dir}/${source}
done

# cd ${mirror_dir} && svn commit