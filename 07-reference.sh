#!/bin/sh -e

# From Ava email: https://www.axolotl-omics.org/
# "file" command revealed name within gzip data to use with --output
# Firefox "save link as" comes up with .zip, which is wrong.  Should be .gz.
curl --continue-at - --location --output AmexT_v47_cds.fa.gz 'http://www.axolotl-omics.org/api?method=Assembly.getSequences&assembly=47&type=cds'
curl --continue-at - --location --output AmexT_v47_dna.fa.gz 'http://www.axolotl-omics.org/api?method=Assembly.getSequences&assembly=47&type=dna'

