/**
@file	miftools.h
@author	Diemo Schwarz
@date	23.02.2010
@brief	Tools for Metric Inverted File Index

Copyright (C) 2008 - 2009 by IRCAM-Centre Georges Pompidou, Paris, France.
*/

/** \page miftools MIF Commandline Tools

\section mifindex mifindex

Index input files, store index in sqlite database.
Data files are in DISCO format.

<b>usage</b>: mifindex nref ki database-name input-file-names...

(-1 for defaults)


\section mifquery mifquery

Query nquery first vectors from query-file against index in index-file, write results to result-file or stdout.
Query-file is in DISCO format, index is an sqlite database generated with mifindex.

<b>usage</b>: mifquery ks mpd  index-file-name query-file-name [nquery] [K to query]  [result-file-name]

(-1 for defaults)

result-file-name will be created as a text file which contains nquery lines of results: 
<pre>
[audioID.segID of the query1] [K] [audioID1.segID1] [dist] ... [audioIDk.segIDk] [dist] 
[audioID.segID of the query2] [K] [audioID1.segID1] [dist] ... [audioIDk.segIDk] [dist] 
...
[audioID.segID of the queryN] [K] [audioID1.segID1] [dist] ... [audioIDk.segIDk] [dist] 
</pre>


\section discotrunc discotrunc

Write nvec first vectors from input-file to output-file or stdout.
All files are in DISCO format.

<b>usage</b>: discotrunc nvec in-file-name [out-file-name]

*/
