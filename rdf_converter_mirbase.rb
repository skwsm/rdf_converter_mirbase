#!/usr/bin/env ruby
require 'rubygems'
require 'bio'
require 'optparse'

module MiRBase

  OrganismCodes = {
    "aae" => ["7159",  "Aedes aegypti"],
    "ame" => ["7460",  "Apis mellifera"],
    "ath" => ["3702",  "Arabidopsis thaliana"],
    "bmo" => ["7091",  "Bombyx mori"],
    "bta" => ["9913",  "Bos taurus"],
    "cbr" => ["6238",  "Caenorhabditis briggsae"],
    "cel" => ["6239",  "Caenorhabditis elegans"],
    "cfa" => ["9615",  "Canis familiaris"],
    "cre" => ["3055",  "Chlamydomonas reinhardtii"],
    "dme" => ["7227",  "Drosophila melanogaster"],
    "dps" => ["46245", "Drosophila pseudoobscura pseudoobscura"],
    "dre" => ["7955",  "Danio rerio"],
    "ebv" => ["10376", "Epstein Barr virus"],
    "fru" => ["31033", "Takifugu rubripes"],
    "gga" => ["9031",  "Gallus gallus"],
    "hcmv" => ["10359", "Human cytomegalovirus"],
    "hsa" => ["9606",  "Homo sapiens"],
    "kshv" => ["435895", "Kaposi sarcoma-associated herpesvirus"],
    "mdo" => ["13616", "Monodelphis domestica"],
    "mghv" => ["33708", "Mouse gammaherpesvirus 68"],
    "mml" => ["9544",  "Macaca mulatta"],
    "mmu" => ["10090", "Mus musculus"],
    "osa" => ["39947", "Oryza sativa japonica"],
    "ptc" => ["3694",  "Populus trichocarpa"],
    "ptr" => ["9598",  "Pan troglodytes"],
    "rno" => ["10116", "Rattus norvegicus"],
    "sme" => ["79327", "Schmidtea mediterranea"],
    "tni" => ["99883", "Tetraodon nigroviridis"],
    "vvi" => ["29760", "Vitis vinifera"],
    "xtr" => ["8364",  "Xenopus tropicalis"],
    "zma" => ["4577",  "Zea mays"]
  }

  def get_org(file)
    if /^(\w+)\.gff/ =~ File.basename(file)
      org_code = $1
      org_code if OrganismCodes.key?(org_code)
    end
  end
  module_function :get_org

  def prefixes
    ["faldo: <http://biohackathon.org/resource/faldo#>",
     "rdfs: <http://www.w3.org/2000/01/rdf-schema#>",
     "dcterms: <http://purl.org/dc/terms/>",
     "hco: <http://identifiers.org/hco/>",
     "term: <http://rdf.ebi.ac.uk/terms/ensembl/>",
     "identifiers: <http://identifiers.org/>",
     "pubmed: <http://rdf.ncbi.nlm.nih.gov/pubmed/>",
     "pmid: <http://identifiers.org/pubmed/>",
     "obo: <http://purl.obolibrary.org/obo/>",
     "so: <http://purl.obolibrary.org/obo/so#>",
     "sio: <http://semanticscience.org/resource/>",
     "up: <http://purl.uniprot.org/uniprot/>",
     "upid: <http://identifiers.org/uniprot/>",
     "mirbase: <http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=>",
     "taxonomy: <http://identifiers.org/taxonomy/>",
     "rfam: <http://identifiers.org/rfam/>",
     "hgnc: <http://identifiers.org/hgnc/>",
     "ncbigene: <http://identifiers.org/ncbigene/>",
     "skos: <http://http://www.w3.org/2004/02/skos/core#>",
     ": <http://purl.jp/bio/10/mirbase/ontology/>"
    ].each {|uri| print "@prefix #{uri} .\n"}
    print "\n"
  end
  module_function :prefixes
 
  class GFF

    def initialize(file, org_code = nil)
      org_code == nil ? @org_code = MiRBase.get_org(file) : @org_code = org_code
      @file = open(file)
    end

    def rdf
      @file.each_line do |line|
        unless /^#/ =~ line
          ary = line.chomp.split
          if ary[8] != nil
            if @org_code == "hsa"
              /chr(\d+)/  =~ ary[0]
              chr         = $1
            else
              chr         = ary[0]
            end
            feature     = ary[2]
            mirna_start = ary[3]
            mirna_end   = ary[4]
            info        = ary[8]

            mirna_info = Hash[*info.split(";").map{ |e|
                                                    /(\S+)\=(\S+)/ =~ e;
                                                    [$1, $2]
                                                  }.flatten ]

            /MIMAT/ =~ mirna_info["ID"] ? mature = true : mature = false
            print_ttl(chr, feature, mirna_start, mirna_end, mirna_info, mature)
          end
        end
      end
    end

    def print_ttl(chr, feature, mirna_start, mirna_end, mirna_info, mature)
      print "mirbase:#{mirna_info["ID"]}\n"
      print "  a :MatureMicroRNA ;\n"
      print "  a obo:SO_0000276 ;\n"
      print "  a :MicroRNA ;\n"
      print "  a obo:SO_0001265 ;\n"
      print "  skos:altLabel \"#{mirna_info["Alias"]}\" ;\n"
      print "  dcterms:identifier \"#{mirna_info["ID"]}\" ;\n"
      if /miRNA_primary_transcript/ =~ feature
      elsif /miRNA/ =~ feature
        mimat_id = mirna_info["Derives_from"]
        print "\n"
        print "  :derives_from mirbase:#{mirna_info["Derives_from"]} ;\n"
      end
      print "  faldo:location [\n"
      print "    a faldo:Region ;\n"
      print "    faldo:begin [\n"
      print "      a faldo:ExactPosition ;\n"
      print "      faldo:position #{mirna_start} ; \n"
      if @org_code == 'hsa'
        print "      faldo:reference <http://identifiers.org/hco/#{chr}\/GRCh38>\n"
        #QName notation with URL fragment (like below) is not allowed for
        #Virtuoso due to a bug.\n"
        #print "      faldo:reference hco:#{chr}\\#GRCh38 \n"
      else
        print "      faldo:reference \"#{chr}\"\n"
      end
      print "    ] ;\n"
      print "    faldo:end [\n"
      print "      a faldo:ExactPosition ;\n"
      print "      faldo:position #{mirna_end} ;\n"
      if @org_code == 'hsa'
        print "      faldo:reference <http://identifiers.org/hco/#{chr}\/GRCh38>\n"
        #print "      faldo:reference hco:#{chr}\\#GRCh38 \n"
      else
        print "      faldo:reference \"#{chr}\"\n"
      end
      print "    ] \n"
      print "  ] .\n"
      print "\n"
    end
  end

  class DB

    def initialize(file, org_code = nil)
      org_code == nil ? @org_code = MiRBase.get_org(file) : @org_code = org_code
      @org_name = MiRBase::OrganismCodes[org_code][1]
      @tax_id = MiRBase::OrganismCodes[org_code][0]
      @ref_rn2ref = {}
      @file = open(file)
      @experiment2references = {}
    end

    def rdf
      ff = Bio::FlatFile.new(Bio::EMBL, @file)
      ff.each do |entry|
#        entry.methods.sort{|a, b| a<=>b}.each{|m| print "#{m}\n"}
        entry.features.each do |m|
          m.qualifiers.each do |q|
            case q.qualifier
            when "accession"
            when "product"
            when "evidence"
            when "experiment"
               parse_experiment(q.value)
            else
#              p q.qualifier
            end
          end
        end
        entry.references.each do |ref|
          @ref_rn2ref[ref.embl_gb_record_number] = ref
        end
        if /#{@org_name}/ =~ entry.description
          print_turtle(entry)
        end
      end
    end

    def parse_experiment(value)
      @experiment2references = {}
      experiments = value.split(", ")
      
      experiments.each do |experiment|
        if m = experiment.match(/^([\w-]+)\s+\[(.+)\]$/)
          @experiment2references[m[1]] = parse_references(m[2])
        else
        end
      end
    end

    def parse_references(value)
      references = []
      ary = value.split(",")
      ary.each do |e|
        case e
        when /^(\d+)\-(\d+)$/
          ($1.to_i..$2.to_i).each {|i| references << i}
        when /^(\d+)$/
          references << $1.to_i
        else
        end
      end
      references
    end

    def print_turtle(entry)

      print "mirbase:#{entry.accession}\n"
      print "  dcterms:identifier \"#{entry.accession}\" ;\n"
      print "  dcterms:description \"#{entry.description}\"@en ;\n"
      print "  rdfs:label \"#{entry.entry_id}\"@en ;\n"
      print "  obo:RO_0002162 taxonomy:#{@tax_id} ;\n"
      entry.dblinks.each do |link|
        case link.database
        when "HGNC"
          print "  rdfs:seeAlso hgnc:#{link.id} ;\n"
        when "ENTREZGENE"
          print "  rdfs:seeAlso ncbigene:#{link.id} ;\n"
        when "RFAM"
          print "  rdfs:seeAlso rfam:#{link.id} ;\n"
        when "TARGETS:PICTAR-VERT"
        when "MIR"
        when "WORMBASE"
        else
          STDERR.print link.database
          STDERR.print "Unkown Database\n"
        end
      end
      entry.references.each do |ref|
        @ref_rn2ref[ref.embl_gb_record_number] = ref
#        unless ref.pubmed == ""
#          print "  dcterms:references pubmed:#{ref.embl_gb_record_number} ;\n"
#          print "  dcterms:references pubmed:#{ref.pubmed} ;\n"
#          print "  dcterms:references pmid:#{ref.pubmed} ;\n"
#        end
      end
      @experiment2references.each do |experiment, references|
        print "  :evidence [\n"
          print "  :experiment_type \"#{experiment}\" ;\n"
          pmids = []
          references.each do |embl_gb_record_number|
            pmids << "pubmed:#{@ref_rn2ref[embl_gb_record_number].pubmed}"
            pmids << "pmid:#{@ref_rn2ref[embl_gb_record_number].pubmed}"
          end
          print "      dcterms:references #{pmids.join(', ')} ;\n"
        print "  ] ;\n"
      end
      if entry.accession == "MIMAT" 
        print "  a :MatureMicroRNA .\n"
      else
        print "  a :MicroRNA .\n"
     end
     print "\n"
    end
  end
end


params = ARGV.getopts('pf:o:g', 'file:', 'prefixes', 'org:', 'gff')
if params["p"] || params["prefixes"]
  MiRBase.prefixes
end
if params["g"] || params["gff"]
  file_path = params["f"] || file_path = params["file"]
  mirbase_gff = MiRBase::GFF.new(file_path)
  mirbase_gff.rdf
end
if params["o"] || params["org"]
  org_code = params["o"] || org_code = params["org"]
  file_path = params["f"] || file_path = params["file"]
  mirbase = MiRBase::DB.new(file_path, org_code)
  mirbase.rdf
end

#mirbase = MiRBase::DB.new(ARGV.shift)


#gff = MiRBase::GFF.new(ARGV.shift)
#MiRBase.prefixes



