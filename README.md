# RDF Converter for miRBase

RDF converter for [miRBase](http://www.mirbase.org/)

## Usage

 $ ruby rdf_converter_mirbase.rb -p -f miRNA.dat -o <org code>

 $ ruby rdf_converter_mirbase.rb -p -g -f <gff file>


 example:

 $ for i in aae ame ath bmo bta cbr cel cfa cre dme dps dre ebv fru gga hcmv hsa kshv mdo mghv mml mmu osa ptc ptr rno sme tni vvi xtr zma; do ruby rdf_converter_mirbase.rb -p -f dat/release_22/miRNA.dat -o ${i} > ttl/${i}.ttl ; done

 $ for i in aae ame ath bmo bta cbr cel cfa cre dme dps dre ebv fru gga hcmv hsa kshv mdo mghv mml mmu osa ptc ptr rno sme tni vvi xtr zma; do ruby rdf_converter_mirbase.rb -p -f dat/release_22/gff/${i}.gff3 -g > ttl/${i}_gff.ttl ; done 

## Input data

Input files are available at [here](http://www.mirbase.org/ftp.shtml).
