#!/usr/bin/env python
"""Fetch GenBank entries for given accessions and report publication title and isolation source.

USAGE:
  python acc2title.epost.py A22237 A22239 A32021 A32022 A33397 > out.tsv
or
  cat ids | python acc2title.epost.py > out.tsv

DEPENDENCIES:
Biopython
"""

import sys
from Bio import Entrez, SeqIO

#define email for entrez login
db           = "nuccore"
Entrez.email = "some_email@somedomain.com"
batchSize    = 100
retmax       = 10**9

#load accessions from arguments
if len(sys.argv[1:]) > 1:
  accs = sys.argv[1:]
else: #load accesions from stdin
  accs = [ l.strip() for l in sys.stdin if l.strip() ]
#first get GI for query accesions
sys.stderr.write( "Fetching %s entries from GenBank: %s\n" % (len(accs), ", ".join(accs[:10])))
query  = " ".join(accs)
handle = Entrez.esearch( db=db,term=query,retmax=retmax )
giList = Entrez.read(handle)['IdList']
sys.stderr.write( "Found %s GI: %s\n" % (len(giList), ", ".join(giList[:10])))
#post NCBI query
search_handle     = Entrez.epost(db=db, id=",".join(giList))
search_results    = Entrez.read(search_handle)
webenv,query_key  = search_results["WebEnv"], search_results["QueryKey"]
#print header
sys.stdout.write("#accesion\ttitle\tisolation_source\n")
#fecth all results in batch of batchSize entries at once
for start in range( 0,len(giList),batchSize ):
  sys.stderr.write( " %9i" % (start+1,))
  #fetch entries in batch
  handle = Entrez.efetch(db=db, rettype="gb", retstart=start, retmax=batchSize, webenv=webenv, query_key=query_key)
  #parse gb entries
  for r in SeqIO.parse(handle, 'gb'):
      title = isolation_source = ""
      #get title
      if 'references' in r.annotations:
          ref   = r.annotations['references'][0]
          title = ref.title
      #get source features
      source = r.features[0]
      #get isolation_source
      if 'isolation_source' in source.qualifiers:
          isolation_source = source.qualifiers['isolation_source'][0]
      #print output to stdout
      #sys.stdout.write("%s\t %s\t %s\n" % r.name, title, isolation_source)
      sys.stdout.write("{}\t {}\t {}\n".format(r.name, title, isolation_source))